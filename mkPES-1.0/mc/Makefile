
include ../../config/make.inc

all: obj bin

OBJECTS= mc.o 

obj: $(OBJECTS)
bin:
	$(LINKER) $(LINKFLAG) -o mc *.o ../share/*.o -L$(SINDO_LIBDIR) $(LIBS)
	mv mc $(SINDO_BINDIR)
clean: 
	rm -f $(OBJECTS) *$(MOD)

