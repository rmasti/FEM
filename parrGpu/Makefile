ifndef OCCA_DIR
ERROR:
	@echo "Error, envirnoment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/Makefile

# define variables
VPATH   = ./
HDRDIR  = ./include  ../eigen ${OCCA_DIR}/include
LIBDIR	= ./lib  ../eigen ${OCCA_DIR}/lib
# adjust this for your system

# set options for this machine
# specify which compilers to use for c, fortran and linking
#CC	= g++
#LD	= g++
CC = mpic++
LD = mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(addprefix -I, $(HDRDIR)) -O3 -std=c++0x -g3 -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) $(paths)
#-fno-exceptions
# link flags to be used 
LDFLAGS	= $(addprefix -I, $(HDRDIR)) -L. $(addprefix -L, $(LIBDIR)) -O3 -DOCCA_VERSION_1_0

# libraries to be linked in
LIBS = -L ${OCCA_DIR}/lib $(links) -lm  -locca

# types of files we are going to construct rules for
.SUFFIXES: .cpp 

# rule for .cpp files (they are exported to C)
.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled

OBJS = \
       src/occaMpiMain.o\
       src/occaMap2Eigen.o\
       src/occaMeshBlock.o\
       src/occaMpiSetBC.o\
       src/occaMhdRT_f.o
       #src/mhdRT_m.o\

TESTS = \
	test/unitTest.o
	#src/fp_f.o\

MAPTESTS = \
	src/map2Eigen.o\
	test/map2EigenTest.o
	#src/fp_f.o\

MPITEST = \
	test/mpiExampleEigen.o\
	src/mhdRT_f.o\
	src/map2Eigen.o
	#src/fp_f.o\

MPICARTTEST = \
	test/mpiCartTest.o\
	src/mhdRT_f.o
	#src/fp_f.o\

OCCAADDVECTEST = \
	test/occaAddVector.o
	#src/fp_f.o\

OCCAMAP2EIGEN = \
	test/occaMap2Eigen.o\
	src/mhdRT_f.o\
	src/map2Eigen.o
	#src/fp_f.o\



main:$(OBJS) 
	$(LD)  $(LDFLAGS) -o main $(OBJS) $(LIBS)
	rm -r $(OBJS)

tests:$(TESTS)
	$(LD)  $(LDFLAGS) -o main $(TESTS) $(LIBS)
	rm -r $(TESTS)

maptest:$(MAPTESTS)
	$(LD)  $(LDFLAGS) -o main $(MAPTESTS) $(LIBS)
	rm -r $(MAPTESTS)

mpitest:$(MPITEST)
	$(LD)  $(LDFLAGS) -o main $(MPITEST) $(LIBS)
	rm -r $(MPITEST)

mpicarttest:$(MPICARTTEST)
	$(LD)  $(LDFLAGS) -o main $(MPICARTTEST) $(LIBS)
	rm -r $(MPICARTTEST)

occaaddvectest:$(OCCAADDVECTEST)
	$(LD)  $(LDFLAGS) -o main $(OCCAADDVECTEST) $(LIBS)
	#rm -r $(OCCAADDVECTEST)

occamap2eigen:$(OCCAMAP2EIGEN)
	$(LD)  $(LDFLAGS) -o main $(OCCAMAP2EIGEN) $(LIBS)
	#rm -r $(OCCAMAP2EIGEN)


# what to do if user types "make clean"
clean :
	rm -r $(OBJS)

cleanall :
	rm -r $(OBJS)
	rm main unittest
