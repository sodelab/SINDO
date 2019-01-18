
include ../../config/make.inc

all: obj bin

OBJECTS= mkrpt.o 

obj: $(OBJECTS)
bin:
	$(LINKER) $(LINKFLAG) -o mkrpt *.o ../share/*.o -L$(SINDO_LIBDIR) $(LIBS)
	mv mkrpt $(SINDO_BINDIR)
clean: 
	rm -f $(OBJECTS) *$(MOD)

