#/bin/bash

make clean
make USER_CFLAGS="-DPHG_TO_P4EST -DDim=2" Coord_dofs
#make USER_CFLAGS="-DDim=2" Coord_dofs
