/*Get the coordinates of the Gauss-lobatto points of DGQ1 for 2D/3D on [0,1]^d domain*/

#include "phg.h"
#include <math.h>

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
#undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

int main(int argc, char *argv[])
{
#if Dim == 3
    char *mesh = "./mesh/cube.mesh";
#else
    char *mesh = "./mesh/square.mesh";
#endif

#ifdef PHG_TO_P4EST
    INT finest_level = 7;
#endif /* defined(PHG_TO_P4EST) */
    GRID *g;
    ELEMENT *e;
    DOF *u_h;
    SOLVER *solver;
    INT refine_level = 5;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &mesh);

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, mesh, FALSE))
    {
        phgError(1, "can't read file \"%s\".\n", mesh);
    }

    // phgPrintf("The Dim of problem: %d\n", (int)Dim);

    u_h = phgDofNew(g, DOF_DGQn[1], 1, "u_h", DofInterpolation);
    phgDofSetDirichletBoundaryMask(u_h, 0);

    phgRefineAllElements(g, refine_level);  // intial refine
    while (refine_level < finest_level)
    {
        phgRefineAllElements(g, 1);
        refine_level += 1;

        char filename[20];
        sprintf(filename, "%dD_Coord_dofs_level%d.txt", (int)Dim, (int)refine_level);
        FILE *file = fopen(filename, "w");
        if (file == NULL)
        {
            printf("Can't construct file successfuly!\n");
            return 1;
        }

        // printf("mesh_level:%d\n", refine_level);
        // fprintf(file, "Dimension: %d\n", (int)Dim);
        // fprintf(file, "mesh_refinelevel: %d\n", (int)refine_level);
        // fprintf(file, "Nelem_global: %d ", g->nelem_global);

        solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);

        if (u_h->type->tensor_product == TRUE)
        {
            int n = u_h->type->order;
            // fprintf(file, "Ndof_global: %d\n", (g->nelem_global * (int)pow(n + 1, Dim)));
        }

        ForAllElements(g, e)
        {
            int eno = e->index;
            int N = DofNBas(u_h, e);
            int order = u_h->type->orders[0];

            // Converse the element local dofs
            for (int Lindex = 0; Lindex < N; Lindex++)
            {
                MAP *map = phgMapCreate(u_h, NULL);
                int Gindex = phgMapE2G(map, 0, e, Lindex);
                FLOAT xyz[Dim], lambda[Dim + 1];
                phgDofGetElementBasisInfo_(u_h, e, Lindex, NULL, NULL, NULL, xyz, lambda);
                // Output the eno,Lindex,coordinates to 'filename.txt'
                // fprintf(file, "eno: %d ", eno);
                // fprintf(file, "Gindex:%d ", Gindex);
                // fprintf(file, "Lindex:%d ", Lindex);
                // fprintf(file, "x:%lf ", (double)xyz[0]);
                // fprintf(file, "y:%lf ", (double)xyz[1]);
                fprintf(file, "%lf ", (double)xyz[0]);
                fprintf(file, "%lf ", (double)xyz[1]);
#if Dim == 3
                // fprintf(file, "z:%lf ", (double)xyz[2]);
                fprintf(file, "%lf ", (double)xyz[2]);
#endif
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    phgDofFree(&u_h);
    phgFreeGrid(&g);
    phgFinalize();
}
