# Project directory
INSTALLDIR := .

# fortran compiler and flags
FC := gfortran
FFLAGS := -Ofast -march=native -fPIC

# libraries required by compiler (LAPACK AND BLAS)
LIBS := -lblas -llapack

# Dependent directories
SRCDIR := $(INSTALLDIR)/src
OBJDIR := $(INSTALLDIR)/obj
INCDIR := $(INSTALLDIR)/include
BINDIR := $(INSTALLDIR)/bin