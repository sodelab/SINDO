
include ../config/make.inc

#
# ------------------------------------------------------------------------------
#             Miscllanious Tools
# ------------------------------------------------------------------------------
#

all: mod obj lib

MODULES= nma.o rg98.o races.o rnwchem.o rpolymer.o rdftb.o spline1d.o spline2d.o spline3d.o
OBJECTS= tools.o intrinsic.o

mod: $(MODULES)
obj: $(OBJECTS)
lib: 
	$(AR) $(ARFLAG) libmisc.a $(OBJECTS) $(MODULES)
	mv libmisc.a $(SINDO_LIBDIR)

clean:
	rm -f $(MODULES) $(OBJECTS) *$(MOD) 

