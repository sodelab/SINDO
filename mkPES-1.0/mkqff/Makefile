
include ../../config/make.inc

all: obj bin

OBJECTS= mkqff1.o mkqff2.o mkqff3.o mkqff4.o

obj: $(OBJECTS)
bin:
	$(LINKER) $(LINKFLAG) -o mkqff *.o ../share/*.o -L$(SINDO_LIBDIR) $(LIBS)
	mv mkqff $(SINDO_BINDIR)
clean: 
	rm -f $(OBJECTS) *$(MOD)

