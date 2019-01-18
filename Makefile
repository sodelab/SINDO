
include ./config/make.inc

all: spr mic vscf pes
clean: sprcl miccl vscfcl pescl
veryclean : clean
	/bin/rm -rf lib/$(TARGET) bin/$(TARGET) 

spr: 
	cd spr-2.2; make ; cd ..
sprcl: 
	cd spr-2.2; make clean ; cd ..

vscf:
	cd vscf-4.0; make ; cd ..
vscfcl:
	cd vscf-4.0; make clean ; cd ..

mic:
	cd misc; make; cd ..
miccl:
	cd misc; make clean; cd ..

pes:
	cd mkPES-1.0; make ; cd ..
pescl:
	cd mkPES-1.0; make clean ; cd ..
