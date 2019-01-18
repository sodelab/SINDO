
include ../config/make.inc

all: sh qff rpt conv sx
clean: clsh clqff clrpt clconv clsx clmc

#   ----------------------------------------------------------------------------------------

sh:
	cd share; make
clsh:
	cd share; make clean

#   ----------------------------------------------------------------------------------------

qff: 
	cd mkqff; make
clqff:
	cd mkqff; make clean 

#   ----------------------------------------------------------------------------------------

rpt:
	cd mkrpt; make
clrpt:
	cd mkrpt; make clean

#   ----------------------------------------------------------------------------------------

conv:
	cd convert; make
clconv:
	cd convert; make clean

#   ----------------------------------------------------------------------------------------

sx:
	cd mkSx; make
clsx:
	cd mkSx; make clean

#   ----------------------------------------------------------------------------------------

mkmc:
	cd mc; make
clmc:
	cd mc; make clean

