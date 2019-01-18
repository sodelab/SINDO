
include ../config/make.inc

#
#   ----------------------------------------------------------------------------------------------------------
#                           Super module provides the basic information
#   ----------------------------------------------------------------------------------------------------------
#

all: mod obj lib
mod: spr_private.o 
obj: spr_public.o
lib: 
	$(AR) $(ARFLAG) libspr.a  spr_private.o spr_public.o
	mv libspr.a $(SINDO_LIBDIR)/

clean:
	$(RM) spr_private.o spr_public.o *$(MOD)

