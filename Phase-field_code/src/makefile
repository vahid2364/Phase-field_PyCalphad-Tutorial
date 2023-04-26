## no implicit rules
.SUFFIXES: 

## definitions
FC=gfortran
FCFLAGS=
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c 
MAKEMOD.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c  

## sources for the program
SOURCES=pfm_CH_2D_AgCu.f90

## main and clean targets
main: $(subst .f90,.o,$(SOURCES))
	$(FC) $(FCFLAGS) -o $@ $+

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod *.anc *.png main *.out fort*

## new compilation rules
%.anc: %.f90
	$(MAKEMOD.f08) $<
	@touch $@

%.o: %.anc
	$(COMPILE.f08) -o $*.o $(<:.anc=.f90)
	@touch $@

## target because -fsyntax-only also needs it)
main.anc: one.anc
one@proc.anc: two.anc
two@proc.anc: one.anc

## .o of a submodule -> .smod of its parent module
one@proc.anc: one.anc
two@proc.anc: two.anc

## list what is inside each file
one.anc: one.mod one.smod
two.anc: two.mod two.smod
one.mod one.smod two.mod two.smod:

## includes
one@proc.anc: one_addone.inc

