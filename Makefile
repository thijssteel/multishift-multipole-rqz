
include make.inc

F90SRC := $(wildcard ./src/*.f90)

LIBOBJS := ./obj/common.o ./obj/d_swaps.o \
			./obj/z_swaps.o ./obj/z_rqz_small.o ./obj/z_rqz_sweep.o ./obj/z_rqz_full.o \
			./obj/z_rqz_aed.o ./obj/d_swaps2.o ./obj/d_rqz_small.o ./obj/d_manipulate_poles.o \
			./obj/d_rqz_aed.o ./obj/d_rqz_sweep.o ./obj/d_rqz_full.o ./obj/d_rqz_small2.o

all: librqz.a examples

librqz.a: $(LIBOBJS)
	ar cr librqz.a $(LIBOBJS)

examples : ./examples/z_rqz_example.out ./examples/d_rqz_example.out

# Compilation rules for the lib itself

./obj/common.o: ./src/common.f make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

# Compile complex rqz

./obj/z_swaps.o: ./src/z_rqz/z_swaps.f ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/z_rqz_small.o: ./src/z_rqz/z_rqz_small.f ./obj/z_swaps.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/z_rqz_aed.o: ./src/z_rqz/z_rqz_aed.f ./obj/z_rqz_small.o ./obj/z_swaps.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/z_rqz_sweep.o: ./src/z_rqz/z_rqz_sweep.f ./obj/z_rqz_small.o ./obj/z_swaps.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/z_rqz_full.o: ./src/z_rqz/z_rqz_full.f ./obj/z_rqz_small.o ./obj/z_swaps.o ./obj/z_rqz_aed.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

# Compile real rqz

./obj/d_swaps.o: ./src/d_rqz/d_swaps.f ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_swaps2.o: ./src/d_rqz/d_swaps2.f ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_manipulate_poles.o: ./src/d_rqz/d_manipulate_poles.f ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_rqz_small.o: ./src/d_rqz/d_rqz_small.f ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/d_manipulate_poles.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_rqz_small2.o: ./src/d_rqz/d_rqz_small2.f ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/d_manipulate_poles.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_rqz_aed.o: ./src/d_rqz/d_rqz_aed.f ./obj/d_rqz_small.o ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/d_manipulate_poles.o ./obj/d_rqz_small2.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_rqz_sweep.o: ./src/d_rqz/d_rqz_sweep.f ./obj/d_rqz_small.o ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/d_manipulate_poles.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

./obj/d_rqz_full.o: ./src/d_rqz/d_rqz_full.f ./obj/d_rqz_small.o ./obj/d_swaps.o ./obj/d_swaps2.o ./obj/d_manipulate_poles.o ./obj/d_rqz_sweep.o ./obj/common.o make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ -I ./include -c $<

# Compilation rules for the examples

./examples/%.out: ./examples/%.f librqz.a make.inc
	$(FC) $(FFLAGS) -J ./include -o $@ $< librqz.a $(LIBS) 

# Some phony commands

.PHONY: clean format

clean:
	rm -f ./*/*.mod ./*/*.o ./*/*.out ./*/*.a