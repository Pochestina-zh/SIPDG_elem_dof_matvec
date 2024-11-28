/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: ipdg.c,v 1.206 2024/04/11 07:52:35 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

#ifndef INTERFACE_TEST	/* test adaptive refinement on an interface */
# define INTERFACE_TEST 0
#endif

#if INTERFACE_TEST	/* test adaptive refinement on an interface */
static void ls_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{x -= .5; y -= .5; z -= .5; *value = x * x + y * y + z * z - .25 * .25;}

static void ls_grad_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad)
{x -= .5; y -= .5; z -= .5; grad[0] = x + x; grad[1] = y + y; grad[2] = z + z;}
#endif	/* TEST_INTERFACE */

/*---------------------------------------------------------------------------
 * Sample IPDG (interior penalty discontinuous Galerkin) code
 * for solving the following Poisson equation:
 *	-\Delta u = f,		in \Omega
 *		u = g_D,	on \Gamma_D
 *	      u_n = g_N\cdot n,	on \Gamma_N
 * Weak form:
 * 		a(u,v) = f(v), \forall v
 * where:
 * 	a(u,v)	= \int_\Omega \grad u . \grad v
 * 		  +
 *		  \sum_{F\in\Gamma_0\cup\Gamma_D} \int_F (
 *			- {\grad u}.n [v]
 *			- \beta [u]{\grad v}.n
 *			+ \gamma_0 p^2 / h_F [u][v]
 *		  )
 *		  +
 *		  \sum_{F\in\Gamma_0\cup\Gamma_N} \int_F (
 *			\gamma_1 h_F / p^2 [\grad u].n [\grad v].n 
 *		  )
 *	f(v)	= \int_\Omega f v
 *		  +
 *		  \sum_{F\in\Gamma_D} \int_F (
 *			g_D (-\beta(\grad v).n + \gamma_0 p^2 / h_F v)
 *		  )
 *		  +
 *		  \sum_{F\in\Gamma_N} \int_F (
 *			g_N.n (v + \gamma_1 h_F / p^2 (\grad v).n)
 *		  )
 * where F: face, h_F: face diameter,
 * 	 \Gamma_0 is the set of interior faces,
 *	 \Gamma = \Gamma_N \cup \Gamma_D = \partial\Omega,
 *
 * Note:
 * 1. "../test/cube_Neumann_on_top.dat" is a test mesh with Neumann boundary.
 * 2. see "../resistivity/resistivity.c" for an example of Robin BC.
 *---------------------------------------------------------------------------*/

/* The IP parameters */
static FLOAT beta = 1.0, gamma0 = 5.0, gamma1 = 0.0;

static BOOLEAN dump_solver = TRUE, dump_reorder = TRUE;
static INT sol_order = -1;

#define PI _F(3.14159265358979323846264338327950288)

/*static*/ void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    assert(sol_order > -3);
#if Dim == 3
    if (sol_order == -2)
	*value = Cos(2.*PI*x) * Cos(2.*PI*y) * Cos(2.*PI*z);
    else if (sol_order == -1)
	*value = Cos(x * y * z);
    else
	*value = 1.0 + Pow(x, sol_order) - 2.0 * Pow(y, sol_order) +
			3.0 * Pow(z, sol_order);
#else
    if (sol_order == -2)
	*value = Cos(2.*PI*x) * Cos(2.*PI*y);
    else if (sol_order == -1)
	*value = Cos(x * y);
    else
	*value = 1.0 + Pow(x, sol_order) - 2.0 * Pow(y, sol_order);
#endif
}

static void
func_grad_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
#if Dim == 3
    if (sol_order == -2) {
	FLOAT Cx = Cos(2.*PI*x), Cy = Cos(2.*PI*y), Cz = Cos(2.*PI*z);
	values[0] = -2. * PI * Sin(2.*PI*x) * Cy * Cz;
	values[1] = -2. * PI * Cx * Sin(2.*PI*y) * Cz;
	values[2] = -2. * PI * Cx * Cy * Sin(2.*PI*z);
    }
    else if (sol_order == -1) {
	FLOAT u = -Sin(x * y * z);
	values[0] = y * z * u;
	values[1] = x * z * u;
	values[2] = x * y * u;
    }
    else if (sol_order == 0) {
	values[0] = values[1] = values[2] = 0.;
    }
    else {
	values[0] = 1.0 * sol_order * Pow(x, sol_order - 1.0);
	values[1] =-2.0 * sol_order * Pow(y, sol_order - 1.0);
	values[2] = 3.0 * sol_order * Pow(z, sol_order - 1.0);
    }
#else
    if (sol_order == -2) {
	FLOAT Cx = Cos(2.*PI*x), Cy = Cos(2.*PI*y);
	values[0] = -2. * PI * Sin(2.*PI*x) * Cy;
	values[1] = -2. * PI * Cx * Sin(2.*PI*y);
    }
    else if (sol_order == -1) {
	FLOAT u = -Sin(x * y);
	values[0] = y * u;
	values[1] = x * u;
    }
    else if (sol_order == 0) {
	values[0] = values[1] = 0.;
    }
    else {
	values[0] = 1.0 * sol_order * Pow(x, sol_order - 1.0);
	values[1] =-2.0 * sol_order * Pow(y, sol_order - 1.0);
    }
#endif
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (sol_order == -2) {
	func_u(x, y, z, value);
	*value *= Dim * 4. * PI * PI;
    }
    else if (sol_order == -1) {
	func_u(x, y, z, value);
#if Dim == 3
	*value *= (x * x * y * y + x * x * z * z + y * y * z * z);
#else
	*value *= (x * x + y * y);
#endif
    }
    else if (sol_order <= 1) {
	*value = 0.;
    }
    else {
	*value = Pow(x, sol_order - 2.) - 2. * Pow(y, sol_order - 2.);
#if Dim == 3
	*value += 3. * Pow(z, sol_order - 2.);
#endif
	*value *= -sol_order * (sol_order - 1.);
    }
}

#if 0
/* FIXME: linking against the static library libphg.a, or enabling the line
 * below, makes the code noticebly faster (>2x) */
# include "../src/quad-cache.c"
#endif

static void
do_face(SOLVER *solver, DOF *u_h, QCACHE *qc, int Q_gD, int Q_gN,
	ELEMENT *e, int face)
{
    GRID *g = u_h->g;
    ELEMENT *e1;
    BYTE face1;
    int n, n1, p;   /* p = polynomial order */
    INT I, J;
    int i, j;
    FLOAT val, nv[Dim], *rule;
    FLOAT G0, G1, h, a;     /* G0 := gamma0*p^2/h, G1 := gamma1*h/p^2 */
    BTYPE bdry_flag;

    n = DofNBas(u_h, e);
    p = DofTypeOrder(u_h, e);
    h = phgGeomGetFaceDiameter(g, e, face);

#if Dim == 3
    rule = phgQuadGetRule2D(g, e, face, 2 * p);
#else
    rule = phgQuadGetRule1D(g, e, face, 2 * p);
#endif
    phgQCSetRule(qc, rule, -1.);
    phgGeomGetFaceOutNormal(g, e, face, nv);
    phgQCSetConstantNormal(qc, nv);

    if ((e1 = phgGetNeighbour(g, e, face)) == NULL) {
	/* boundary face */
	bdry_flag = GetFaceBTYPE(g, e, face);
	if (bdry_flag != NEUMANN)
	    bdry_flag = DIRICHLET;
	n1 = 0;
	face1 = -1;
	G0 = gamma0 * p * p / h;
	G1 = gamma1 * h / (p * (FLOAT)p);
	/* RHS */
	for (i = 0; i < n; i++) {
	    I = phgSolverMapE2G(solver, 0, e, i);

	    if (bdry_flag == DIRICHLET) {
		/* -\beta\int_\Gamma_D g_D (\grad v).n */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD,   PROJ_NONE, 0,
				qc, e->index, face, Q_GRAD, PROJ_DOT,  i);
		val = -beta * a;
		/* G0 \int_\Gamma_D g_D v */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD,  PROJ_NONE, 0,
				qc, e->index, face, Q_BAS, PROJ_NONE, i);
		val += G0 * a;
	    }
	    else {
		/* \int_\Gamma_N g_N v */
		val = phgQCIntegrateFace(
			    qc, e->index, face, Q_gN,  PROJ_DOT,  0,
			    qc, e->index, face, Q_BAS, PROJ_NONE, i);
		/* \int_\Gamma_N G1 g_N (\grad v).n */
		a = phgQCIntegrateFace(
			    qc, e->index, face, Q_gN,   PROJ_DOT, 0,
			    qc, e->index, face, Q_GRAD, PROJ_DOT, i);
		val += G1 * a;
	    }
	    phgSolverAddGlobalRHSEntry(solver, I, val);
	}
    }
    else {
	/* interior face shared by "e" and "e1" */
	bdry_flag = INTERIOR;
	face1 = phgOppositeFace(g, e, face, e1);
	n1 = DofNBas(u_h, e1);
	i = DofTypeOrder(u_h, e1);
	if (p < i)
	    p = i;
	G0 = gamma0 * p * (FLOAT)p / h;
	G1 = gamma1 * h / (p * (FLOAT)p);
    }

    /* loop on {basis funcs in e} \cup {basis funcs in e1} */
    for (i = 0; i < n + n1; i++) {
	SHORT face_i, bas_i;
	ELEMENT *e_i;

	if (i < n) {
	    e_i = e;
	    face_i = face;
	    bas_i = i;
	}
	else {
	    e_i = e1;
	    face_i = face1;
	    bas_i = i - n;
	}
	I = phgSolverMapE2G(solver, 0, e_i, bas_i);
	/* loop on {basis funcs in e} \cup {basis funcs in e1} */
	for (j = 0; j < n + n1; j++) {
	    SHORT face_j, bas_j;
	    ELEMENT *e_j;

	    if (j < n) {
		e_j = e;
		face_j = face;
		bas_j = j;
	    }
	    else {
		e_j = e1;
		face_j = face1;
		bas_j = j - n;
	    }
	    J = phgSolverMapE2G(solver, 0, e_j, bas_j);

            /* skip jumps on interior face for continuous element */
            if (DofFESpace(u_h) == FE_H1 && e1 != NULL
#ifdef PHG_TO_P4EST
		/* e, e1 are equal size (unless dependent DOFs are removed) */
		&& e->generation == e1->generation
#endif	/* defined(PHG_TO_P4EST) */
		)
                continue;

	    /*-----------------------------------------------------------------
	     * Note: If the normal vector is reversed, the sign of [.] changes,
	     *	      while the sign of {.} is unaffected.
	     *----------------------------------------------------------------*/

	    val = 0.0;

	    if (bdry_flag != NEUMANN) {
		/* -\int {\grad u}.n [v] (i<=>v, j<=>u, n<=>e_j) */
		a = phgQCIntegrateFace(
			qc, e_j->index, face_j, Q_GRAD, PROJ_DOT,  bas_j,
			qc, e_i->index, face_i, Q_BAS,  PROJ_NONE, bas_i);
		if (bdry_flag == INTERIOR)
		    a *= (e_i == e ? 0.5 : -0.5);
		val = -a;

		/* -\int \beta [u]{\grad v}.n (i<=>v, j<=>u, n<=>e_i) */
		a = phgQCIntegrateFace(
			qc, e_j->index, face_j, Q_BAS,  PROJ_NONE, bas_j,
			qc, e_i->index, face_i, Q_GRAD, PROJ_DOT,  bas_i);
		if (bdry_flag == INTERIOR)
		    a *= (e_j == e ? 0.5 : -0.5);
		val += -beta * a;

		/* \int G0 [u][v] (i<=>v, j<=>u, n<=>e) */
		a = phgQCIntegrateFace(
			qc, e_j->index, face_j, Q_BAS, PROJ_NONE, bas_j,
			qc, e_i->index, face_i, Q_BAS, PROJ_NONE, bas_i);
		if (bdry_flag == INTERIOR && e_i != e_j)
		    a = -a;
		val += G0 * a;
	    }

	    if (bdry_flag != DIRICHLET) {
		/* \int G1 [\grad u].n [\grad v].n (i<=>v, j<=>u, n<=>e_i) */
		a = phgQCIntegrateFace(
			    qc, e_j->index, face_j, Q_GRAD, PROJ_DOT, bas_j,
			    qc, e_i->index, face_i, Q_GRAD, PROJ_DOT, bas_i);
		if (bdry_flag == INTERIOR && e_i != e_j)
		    a = -a;
		val += G1 * a;
	    }

	    phgSolverAddGlobalMatrixEntry(solver, I, J, val);
	}
    }

    phgFree(rule);
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    QCACHE *qc;
    int Q_f, Q_gD, Q_gN;

    assert(u_h->dim == 1);

    qc = phgQCNew(QD_DEFAULT, u_h);
    Q_f = phgQCAddFEFunction(qc, f_h);
    Q_gD = phgQCAddXYZFunction(qc, func_u, 1);
    Q_gN = phgQCAddXYZFunction(qc, func_grad_u, 3);

#if USE_OMP
#pragma omp parallel for private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	int N = DofNBas(u_h, e);
	INT I, J, eno = e->index;
	int i, j;
	FLOAT val, *rule;
	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
#if Dim == 3
	rule = phgQuadGetRule3D(g, e, 2 * DofTypeOrder(u_h, e) /*- 2*/);
#else
	rule = phgQuadGetRule2D(g, e, 2 * DofTypeOrder(u_h, e) /*- 2*/);
#endif
	phgQCSetRule(qc, rule, -1.);
	for (i = 0; i < N; i++) {
	    I = phgSolverMapE2G(solver, 0, e, i);
	    for (j = 0; j <= i; j++) {
		J = phgSolverMapE2G(solver, 0, e, j);
		/* \int_T \grad u . \grad v */
		val = phgQCIntegrate(qc, eno, Q_GRAD, j, qc, eno, Q_GRAD, i);
		phgSolverAddGlobalMatrixEntry(solver, I, J, val); 
		if (i != j)
		    phgSolverAddGlobalMatrixEntry(solver, J, I, val); 
	    }
	    /* \int_T fv */
	    val = phgQCIntegrate(qc, eno, Q_BAS, i, qc, eno, Q_f, 0);
	    phgSolverAddGlobalRHSEntry(solver, I, val);
	}
	phgFree(rule);

	/* loop on the faces of the elements */
	for (i = 0; i < NFace; i++) {
	    ELEMENT *e1 = phgGetNeighbour(g, e, i);
	    if (e1 != NULL) {
		/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
		if (e->generation < e1->generation)
		    continue;
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index))
		    continue; /* process each interior face just once */
	    }
	    do_face(solver, u_h, qc, Q_gD, Q_gN, e, i);
	}
    } ForAllElementsEnd

    phgQCFree(&qc);
}

int
main(int argc, char *argv[])
{
#if Dim == 3
    char *fn = "../p4est/cube.mesh";	/* hexahedral mesh */
#else
    char *fn = "../p4est/square.mesh";	/* quadrilateral mesh */
#endif
    INT refine = 7, periodicity = 0 /* {1|2|4} (1=x, 2=y, 4=z) */;
#ifdef PHG_TO_P4EST
    INT refine0 = 2, refine_step = 1;
#else	/* defined(PHG_TO_P4EST) */
    INT refine0 = 0, refine_step = 3;
#endif	/* defined(PHG_TO_P4EST) */
    GRID *g;
    DOF *u_h, *u_old, *f_h, *error, *gerror, *gu_h;
    SOLVER *solver;
    size_t mem_peak;
    double t, L2norm, H1norm, L2err, H1err;
    BOOLEAN vtk_flag = FALSE, rel_err = TRUE, corner_flag = FALSE;
    BOOLEAN debug_pre = FALSE;
    int level, total_level;
#ifndef PHG_TO_P4EST
    BOOLEAN hp_test = FALSE;
    phgOptionsRegisterNoArg("-hp_test", "Test h-p DOF", &hp_test);
#endif	/* PHG_TO_P4EST */

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-periodicity", "Set periodicity", &periodicity);
    phgOptionsRegisterInt("-refine0", "Initial refinement levels", &refine0);
    phgOptionsRegisterInt("-refine", "Repeated refinement levels", &refine);
    phgOptionsRegisterInt("-refine_step", "Refinement step", &refine_step);
    phgOptionsRegisterFloat("-beta", "The parameter beta", &beta);
    phgOptionsRegisterFloat("-gamma0", "The parameter gamma0", &gamma0);
    phgOptionsRegisterFloat("-gamma1", "The parameter gamma1", &gamma1);
    phgOptionsRegisterInt("-sol_order", "Analytic solution's polynomial order "
			  "(<0: non-polynomial, -99: FE order)", &sol_order);

    phgOptionsRegisterNoArg("-rel", "Print relative errors", &rel_err);
    phgOptionsRegisterNoArg("-dump_solver", "Output matrix/RHS as .m files",
				&dump_solver);
    phgOptionsRegisterNoArg("-dump_reorder", "Reorder unknowns in the .m files",
				&dump_reorder);
    phgOptionsRegisterNoArg("-vtk", "Output 'ipdg.vtk'", &vtk_flag);
    phgOptionsRegisterNoArg("-corner", "Refine the corner element during "
			    "initial refinements", &corner_flag);
    phgOptionsRegisterNoArg("-debug_pre", "Set initial x0 to analytic solution",
			    &debug_pre);

    phgOptionsPreset("-dof_type DG1 -solver gmres");

#if USE_HYPRE
    /* BoomerAMG seems to work pretty well */
    phgOptionsPreset("-solver hypre -hypre_solver gmres -hypre_pc boomeramg");
#endif

    phgInit(&argc, &argv);

    if (sol_order == -99)
	sol_order = DOF_DEFAULT->order;

    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    t = phgGetTime(NULL);
    total_level = 0;
    while (refine0 > 0) {
	BOOLEAN marked = FALSE;
	level = refine0 > refine_step ? refine_step : refine0;
	if (corner_flag) {
	    /* only refine the (0,0,0) corner element (for debugging) */
	    ELEMENT *e;
	    ForAllElements(g, e) {
		e->mark = 0;
#ifdef PHG_TO_P4EST
		if (e->corners[0][0] != g->bbox[0][0] ||
		    e->corners[0][1] != g->bbox[0][1] ||
		    e->corners[0][2] != g->bbox[0][2])
		    continue;
#else	/* defined(PHG_TO_P4EST) */
		int i;
		for (i = 0; i < NVert; i++)
		    if (g->verts[e->verts[i]][0] == g->bbox[0][0] &&
			g->verts[e->verts[i]][1] == g->bbox[0][1] &&
			g->verts[e->verts[i]][2] == g->bbox[0][2])
			break;
		if (i >= NVert)
		    continue;
#endif	/* defined(PHG_TO_P4EST) */
		e->mark = level;
		marked = TRUE;
	    }
	}
	if (marked)
	    phgRefineMarkedElements(g);
	else
	    phgRefineAllElements(g, level);	/* uniform refinement */
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	refine0 -= level;
	total_level += level;
    }

    /* The discrete solution */
#ifndef PHG_TO_P4EST
    if (hp_test) {
	ELEMENT *e;
	HP_TYPE *hp;
	hp = phgHPNew(g, DOF_DEFAULT->fe_space == FE_H1 ? HP_HB : HP_DG);
	ForAllElements(g, e)
	    e->hp_order = DOF_DEFAULT->order + GlobalElement(g, e->index) % 3;
	phgHPSetup(hp, FALSE);
	u_h = phgHPDofNew(g, hp, 1, "u_h", DofInterpolation);
	phgHPFree(&hp);
    }
    else 
#endif	/* PHG_TO_P4EST */
    {
	u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    }
    phgDofSetDirichletBoundaryMask(u_h, 0);

    while (TRUE) {
#if INTERFACE_TEST
	XFEM_INFO *xi;
	DOF *ls, *ls_grad;
	ls = phgDofNew(g, DOF_ANALYTIC, 1, "ls", ls_func);
	ls_grad = phgDofNew(g, DOF_ANALYTIC, 3, "ls_grad", ls_grad_func);
	xi = phgXFEMInit(ls, ls_grad, 2, 2 * u_h->type->order + 3);
	phgDofFree(&ls);
	phgDofFree(&ls_grad);
	phgXFEMFree(&xi);
#endif	/* INTERFACE_TEST */
	phgSetupHalo(g, HALO_FACE);
	phgPrintf("********** Level %d, %d proc%s, %"
		  dFMT" elem%s, LIF %0.2lf, refine time: %0.4lg\n",
		  total_level, g->nprocs, g->nprocs > 1 ? "s" : "",
		  g->nleaf_global, g->nleaf_global > 1 ? "s" : "",
		  (double)g->lif, phgGetTime(NULL) - t);

#if 0
# warning ****************** using 0 as initial solution!
	phgDofSetDataByValue(u_h, 0.0);
#endif

#if 0
	/* Test DOF processing in phgSetupHalo and phgRemoveHalo */
	phgRemoveHalo(g);
	phgSetupHalo(g, HALO_FACE);
#endif

	phgPrintf("Building linear equations: ");
	if (beta == 1.0 && !phgOptionsIfUsed("-solver_symmetry"))
	    phgOptionsSetOptions("-solver_symmetry=spd");
	t = phgGetTime(NULL);
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	/* RHS function */
	f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f);
	build_linear_system(solver, u_h, f_h);
	phgDofFree(&f_h);
	phgPrintf("%"dFMT" unknowns, build time: %0.4lg\n", solver->nglobal,
				phgGetTime(NULL) - t);
	u_old = phgDofCopy(u_h, NULL, NULL, "u_old");

	char mat_name[20];
	char rhs_name[20];
	char perm_name[20];
	sprintf(mat_name, "3D_A1_%d", total_level);
	sprintf(rhs_name, "3D_b1_%d", total_level);
	sprintf(perm_name, "3D_p1_%d", total_level);
	if (dump_solver)
	    phgSolverDumpMATLAB_(solver, mat_name, rhs_name, perm_name, dump_reorder);
	phgPrintf("Solving linear equations with %s:\n",
					solver->oem_solver->name);
	if (debug_pre)
	    phgDofSetDataByFunction(u_h, func_u);
	t = phgGetTime(NULL);
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("    nits: %d; residual: %lg; peak mem: %0.2lfMB, "
		  "solve time: %0.4lg\n", solver->nits,
		  (double)solver->residual, mem_peak / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);
	if (solver->cond > 0.0)
	    phgPrintf("    Condition number: %0.2le\n", (double)solver->cond);
	phgSolverDestroy(&solver);

	t = phgGetTime(NULL);
#ifndef PHG_TO_P4EST
	if (u_h->hp != NULL)
	    error = phgHPDofNew(g, u_h->hp, 1, "error", func_u);
	else
#endif  /* PHG_TO_P4EST */
	    error = phgDofNew(g, u_h->type, 1, "error", func_u);
	gerror = phgDofGradient(error, NULL, NULL, NULL);
	gu_h = phgDofGradient(u_h, NULL, NULL, NULL);
	if (rel_err) {
	    L2norm = phgDofNormL2(error);
	    H1norm = phgDofNormL2(gerror);
	    H1norm += L2norm;
#if 1
	    /* get max(|u|,|u_h|) (comment out to revert to before 20211209 */
	    t = phgDofNormL2(u_h);
	    if (L2norm < t)
		L2norm = t;
	    t += phgDofNormL2(gu_h);
	    if (H1norm < t)
		H1norm = t;
#endif
	    L2norm = L2norm == 0. ? 1.0 : L2norm;
	    H1norm = H1norm == 0. ? 1.0 : H1norm;
	}
	else {
	    L2norm = H1norm = 1.0;
	}
	phgDofAXPY(-1.0, u_h, &error);
	phgDofAXPY(-1.0, gu_h, &gerror);
	L2err = phgDofNormL2(error);
	H1err = phgDofNormL2(gerror);
	H1err += L2err;
	phgDofAXPY(-1.0, u_h, &u_old);
	phgPrintf("%s err: L2 = %0.10le; H1 = %0.10le; "
		  "|u_h-u_H| = %0.5le\n", rel_err ? "Rel." : "Abs.",
		  L2err / L2norm, H1err / H1norm,
		  (double)phgDofNormL2(u_old) / L2norm);
	if (vtk_flag) {
	    char name[128];
	    sprintf(name, "ipdg-%02d.vtk", total_level);
	    phgExportVTK(g, name, u_h, error, NULL);
	    phgPrintf("\"%s\" created.\n", name);   /* FIXME: DG to VTK? */
	}

	phgDofFree(&u_old);
	phgDofFree(&error);
	phgDofFree(&gerror);
	phgDofFree(&gu_h);

	if (refine <= 0)
	    break;

	t = phgGetTime(NULL);
	level = refine > refine_step ? refine_step : refine;
	phgRefineAllElements(g, level);
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	refine -= level;
	total_level += level;
    }

    phgDofFree(&u_h);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
