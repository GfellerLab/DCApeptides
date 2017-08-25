#
# Makefile: Makefile for the coev program.
#
# Linda Dib & Wim Hordijk   Last modified: 22 July 2014
#
CXX		= g++
CXXFLAGS	= -ansi -O3 -static -g
CC		= gcc
CCFLAGS		= -O3
LIB		= -lgfortran -lblas -llapack -lm -lgsl -lgslcblas
INC_PATH	= -pg -I /usr/include -I /home/ldib/DCA/externalLib
LIB_PATH	= -L /usr/lib -L /usr/lib/lapack/ -L /home/ldib/DCA/externalLib
OBJS		= int_correl.o

int_correl : $(OBJS)
	$(CXX) $(INC_PATH) $(LIB_PATH) -o $@ $(OBJS) $(LIB)

%.o : %.cpp
	$(CXX) $(INC_PATH) -c $(CXXFLAGS) $<

%.o : %.c
	$(CC) $(INC_PATH) -c $(CCFLAGS) $<




clean:
	rm -f *~ *.o
	rm inter_correl
