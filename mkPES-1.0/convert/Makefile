
include ../../config/make.inc

all: mod rpt rpt2

BINDIR= $(SINDO_BINDIR)

mod: getfname.o

BINRPT= $(BINDIR)/rpt2sx
rpt: $(BINRPT)
$(BINRPT):  rpt2sx.o getfname.o
	$(LINKER) $(LINKFLAG) -o $(BINRPT) rpt2sx.o getfname.o -L$(SINDO_LIBDIR) $(LIBS)

BINRPT2= $(BINDIR)/rpt22sx
rpt2: $(BINRPT2)
$(BINRPT2):  rpt22sx.o getfname.o
	$(LINKER) $(LINKFLAG) -o $(BINRPT2) rpt22sx.o getfname.o -L$(SINDO_LIBDIR) $(LIBS)

clean: 
	rm -f rpt2sx.o getfname.o rpt22sx.o *$(MOD)

