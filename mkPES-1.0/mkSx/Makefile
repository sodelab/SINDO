
include ../../config/make.inc

all: obj bin

OBJECTS= mkSx.o 

obj: $(OBJECTS)
bin:
	$(LINKER) $(LINKFLAG) -o mkSx *.o ../share/*.o ../convert/getfname.o -L$(SINDO_LIBDIR) $(LIBS)
	mv mkSx $(SINDO_BINDIR)
clean: 
	rm -f $(OBJECTS) *$(MOD)

