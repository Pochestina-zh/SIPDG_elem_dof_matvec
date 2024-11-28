main=Coord_dofs
default: $(main)
$(main): $(main).o
$(main).o: $(main).c


clean:
	-/bin/rm -f $(main) 
	-/bin/rm -f *.o
cleanerr:
	-/bin/rm -f *.err
cleanlog:
	-/bin/rm -f *.log

# VERSION 0.9.9-20240820
# PHG_MAKEFILE_INC=~/local/phg-0.9.9-20241026-long_double-long/share/phg/Makefile.inc

# VERSION 0.9.9-20241106
# PHG_MAKEFILE_INC=~/local/phg-0.9.9-20241121-long_double-long/share/phg/Makefile.inc	

# VERSION 0.9.7-20220923
#PHG_MAKEFILE_INC=~/software_package/phg/phg-0.9.7-20220923/Makefile.inc

PHG_MAKEFILE_INC=~/software_package/phg/phg-0.9.9-20241106/Makefile.inc

include ${PHG_MAKEFILE_INC}
