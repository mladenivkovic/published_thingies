# This Makefile is used to compile the main.c code correctly.
#
# SYNTAX:
#target: depencencies
#[tab] 	system command

#==================================
# COMPILATION TIME PARAMETERS
#==================================

EXEC=riemann-vacuum  # executable name
SRCDIR=../src

DEFINES=


#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) -fopenmp



# flags for gcc compiler
CFLAGS= -g -Og -Wall -pedantic -std=c11 -Wextra  -fsignaling-nans  -ffinite-math-only -fno-signed-zeros -Warray-bounds $(DEFINES) # debug
# CFLAGS= -O3  		 	   # performance


LDFLAGS= -lm


#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(SRCDIR)


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
OBJECTS = main.o io.o params.o riemann.o gas.o



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
	rm -f $(OBJECTS) $(EXEC)

