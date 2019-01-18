 SINDO_ROOT = $SINDO_ROOT
 #SINDO_ROOT = $SINDO_ROOT
 SINDO_ROOT = /Users/osode/Desktop/dimer-manuscript-files/step/../tools/sindo
 TARGET = 

 SINDO_BINDIR = ${SINDO_ROOT}/bin/
 SINDO_LIBDIR = ${SINDO_ROOT}/lib/

 F77C = gfortran
 #F77OPT= -w -cm -mp1 -O2 -prefetch  -unroll
 F77OPT= -w -cm -mp1 -O3 -tpp7 -prefetch  -unroll
 F77NOOPT= -w -cm -mp1 -O

 F90C = gfortran 
 #F90OPT= -w -cm -mp1 -O2 -prefetch  -unroll #-m
 F90OPT= -static -std=legacy -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
 #F90OPT= -static -std=legacy -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
 F90NOOPT= -w -cm -mp1 -O -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
 #F90NOOPT= -w -cm -mp1 -O -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

 LINKER = gfortran
 LINKFLAG = 
 LIBS = -L/opt/local/lib/ -L/opt/local/lib/libgcc/ -lblas -llapack -lvib -lmisc -lspr 

 AR = ar
 ARFLAG = cr

 MOD = .mod

 RM = rm -f

.f90.o:
	$(F90C) $(F90OPT) -c $< -o $@
 
.SUFFIXES: .f90 .o

