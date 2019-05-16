# This Makefile is used to compile the main.c code correctly.
#
# SYNTAX:
#target: depencencies
#[tab] 	system command

#==================================
# COMPILATION TIME PARAMETERS
#==================================

EXEC=godunov  # executable name
SRCDIR=../src

RIEMANN=HLLC

RIEMANNFILE=riemann.o
ifeq ($(RIEMANN), EXACT)
	RIEMANNFILE=riemann-exact.o
endif
ifeq ($(RIEMANN), TRRS)
	RIEMANNFILE=riemann-trrs.o
endif
ifeq ($(RIEMANN), TSRS)
	RIEMANNFILE=riemann-tsrs.o
endif
ifeq ($(RIEMANN), HLL)
	RIEMANNFILE=riemann-hll.o
endif
ifeq ($(RIEMANN), HLLC)
	RIEMANNFILE=riemann-hllc.o
endif



DEFINES=-D RIEMANN_$(RIEMANN)

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) -static


# flags for gcc compiler
CFLAGS= -g -Og -Wall -pedantic -std=c11 -Wextra  -fsignaling-nans  -ffinite-math-only -fno-signed-zeros -Warray-bounds $(DEFINES) # debug
# CFLAGS= -O3  		 	   # performance


LDFLAGS= -lm -fopenmp



#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(SRCDIR)


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main.o gas.o godunov.o io.o params.o $(RIEMANNFILE)



#==================================
# RECEPIES
#==================================


default: $(EXEC)


$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o:  %.c %.h
	$(CC) $(CFLAGS) -o $@ -c $<
#-c : compile and assemble, but do not link



#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC) riemann*.o

