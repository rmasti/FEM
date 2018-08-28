# define variables
VPATH   = ./
HDRDIR  = ./include 
LIBDIR	= ./lib  
# adjust this for your system

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= g++
LD	= g++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(addprefix -I, $(HDRDIR)) -O3 #-std=c++0x
#-fno-exceptions
# link flags to be used 
LDFLAGS	= $(addprefix -I, $(HDRDIR)) -L. $(addprefix -L, $(LIBDIR)) -O3

# libraries to be linked in
LIBS = -lm 

# types of files we are going to construct rules for
.SUFFIXES: .cpp 

# rule for .cpp files (they are exported to C)
.cpp.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled

OBJS = \
       src/mhdRT_m.o\
       src/mhdRT_f.o

TESTS = \
	test/unitTest.o
	#src/fp_f.o\

main:$(OBJS) 
	$(LD)  $(LDFLAGS) -o main $(OBJS) $(LIBS)
	rm -r $(OBJS)

tests:$(TESTS)
	$(LD)  $(LDFLAGS) -o unittest $(TESTS) $(LIBS)
	rm -r $(TESTS)


# what to do if user types "make clean"
clean :
	rm -r $(OBJS)

cleanall :
	rm -r $(OBJS)
	rm main unittest
