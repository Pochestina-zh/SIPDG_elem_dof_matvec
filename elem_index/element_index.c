/*Get the geometric infomation of elements with respect to topology of the specific mesh*/

#include "phg.h"

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
    GRID *g;
    ELEMENT *e;
    BTYPE bdry_flag;

    phgPrintf("%s\n", mesh);

    int refine_level = 7;
    int current_level = 5;
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &mesh);

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, mesh, FALSE))
    {
        phgError(1, "can't read file \"%s\".\n", mesh);
    }

    phgPrintf("The Dim of problem: %d\n", (int)Dim);

    phgRefineAllElements(g,current_level);  // initial refine
    while (current_level < refine_level)
    {
        phgRefineAllElements(g, 1);
        current_level += 1;

        char filename[20];
        sprintf(filename, "%dD_element_index_level%d.txt", (int)Dim, (int)current_level);
        FILE *file = fopen(filename, "w");
        if (file == NULL)
        {
            printf("Can't construct file successfuly!\n");
            return 1;
        }

        ForAllElements(g, e)
        {
            for (size_t face = 0; face < NFace; face++)
            {
                ELEMENT *e1 = phgGetNeighbour(g, e, face);

                if (e1 == NULL) /* boundary face */
                {
                    bdry_flag = GetFaceBTYPE(g, e, face);
                    if (bdry_flag != NEUMANN)
                        bdry_flag = DIRICHLET;
                    fprintf(file, "%d ", -1);
                }
                else // interior face
                {
                    fprintf(file, "%d ", (int)e1->index);
                }

//                 // get the coord of facecenter
//                 int v0, v1, v2,v3;
//                 FLOAT pt[Dim];
// #if Dim == 3
//                 v0 = GetFaceVertex(face, 0);
//                 v1 = GetFaceVertex(face, 1);
//                 v2 = GetFaceVertex(face, 2);
//                 v3 = GetFaceVertex(face, 3);
//                 for (size_t jj = 0; jj < Dim; jj++)
//                 {
//                     pt[jj] = (g->elems[e->index]->corners[(v0 >> jj) & 1][jj] +
//                               g->elems[e->index]->corners[(v1 >> jj) & 1][jj] +
//                               g->elems[e->index]->corners[(v2 >> jj) & 1][jj]) /
//                              3.0;
//                 }
//                 fprintf(file, "FaceCenter:(%2lf, %2lf, %2lf) ", (double)pt[0], (double)pt[1], (double)pt[2]);
// #else
//                 v0 = GetEdgeVertex(face, 0);
//                 v2 = GetEdgeVertex(face, 1);
//                 for (size_t jj = 0; jj < Dim; jj++)
//                 {
//                     pt[jj] = (g->elems[e->index]->corners[(v0 >> jj) & 1][jj] +
//                               g->elems[e->index]->corners[(v1 >> jj) & 1][jj]) *
//                              .5;
//                 }
//                 fprintf(file, "FaceCenter:(%2lf, %2lf) ", (double)pt[0], (double)pt[1]);
// #endif
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    phgFreeGrid(&g);
    return 0;
}