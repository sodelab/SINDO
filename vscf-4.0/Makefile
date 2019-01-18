
include ../config/make.inc

all: mod obj lib sindo
clean: cllib clsindo

#
#   -------------------------------------------------------------------------------
#                           VSCF module
#   -------------------------------------------------------------------------------
#


MODULES= vib_private.o qff_private.o modal_private.o pegrid_private.o \
         vscf_private.o vcor_private.o vci_private.o vpt_private.o vav_private.o \
         spfn1.o spfn2.o spfn3.o ho.o
OBJECTS= vib_public.o qff_public.o modal_public.o pegrid_public.o \
         vscf_public.o vcor_public.o vci_public.o vpt_public.o vav_public.o

mod: $(MODULES)
obj: $(OBJECTS)
lib:
	$(AR) $(ARFLAG) libvib.a $(OBJECTS) $(MODULES)
	mv libvib.a $(SINDO_LIBDIR)
cllib:
	rm -f $(OBJECTS) $(MODULES) *$(MOD)

spfn1.o          :  vib_private.f90
spfn2.o          :  vib_private.f90
spfn3.o          :  vib_private.f90
modal_private.o  :  vib_private.f90
pegrid_private.o :  vib_private.f90
vscf_private.o   :  vib_private.f90
vav_private.o    :  vib_private.f90
vcor_private.o   :  vib_private.f90
vci_private.o    :  vib_private.f90 vcor_private.f90
vpt_private.o    :  vib_private.f90 vcor_private.f90

vib_public.o     :  vib_private.f90  
modal_public.o   :  vib_private.f90 modal_private.f90  
pegrid_public.o  :  vib_private.f90 pegrid_private.f90
qff_public.o     :  vib_private.f90 qff_private.f90
vscf_public.o    :  vib_private.f90 vscf_private.f90  
vav_public.o     :  vib_private.f90 vav_public.f90
vcor_public.o    :  vib_private.f90 vcor_private.f90
vci_public.o     :  vib_private.f90 vcor_private.f90 vci_private.f90
vpt_public.o     :  vib_private.f90 vcor_private.f90 vpt_private.f90

#
#   -------------------------------------------------------------------------------
#                           PROGRAM Sindo
#   -------------------------------------------------------------------------------
#
PROG= sindo
sindo: main.o bin
bin: 
	$(LINKER) $(LINKFLAG) -o $(SINDO_BINDIR)/$(PROG) main.o -L$(SINDO_LIBDIR) $(LIBS)
clsindo:
	rm -f $(PROG) main.o 

