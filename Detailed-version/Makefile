GCC = g++
ICC = icpc

DBGFLAGS = -O0 -Wall -g
GCCFLAGS = -O3 -march=native
ICCFLAGS = -O3 -xS

# set CPP and CPPFLAGS to the desired values
CPP = $(GCC)
CPPFLAGS = $(GCCFLAGS)

LINKER = $(CPP)
LNFLAGS = -lm -lGL -lGLU -lglut
EXE = run

OBJS = global.o streamingvisitor.o connectneighbors.o coeffs_haar.o octtree.o prune.o main.o traverse.o marchingcubes.o reconstruct.o matrix.o readpts.o blur.o eval_haar.o parse.o plane.o eval_daub4.o average.o octnode.o coeffs_daub4.o

all: $(EXE)

clean:
	rm -f $(EXE) $(OBJS)

$(EXE): $(OBJS)
	$(LINKER) $(LNFLAGS) -o $(EXE) $(OBJS)

global.o: global.cpp vect.h global.h parse.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c global.cpp

streamingvisitor.o: streamingvisitor.cpp vect.h common.h
	$(CPP) $(CPPFLAGS) -c streamingvisitor.cpp

connectneighbors.o: connectneighbors.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c connectneighbors.cpp

coeffs_haar.o: coeffs_haar.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c coeffs_haar.cpp

octtree.o: octtree.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c octtree.cpp

prune.o: prune.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c prune.cpp

main.o: main.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c main.cpp

traverse.o: traverse.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c traverse.cpp

marchingcubes.o: marchingcubes.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c marchingcubes.cpp

reconstruct.o: reconstruct.cpp common.h vect.h global.h parse.h array2d.h matrix.h
	$(CPP) $(CPPFLAGS) -c reconstruct.cpp

matrix.o: matrix.cpp vect.h matrix.h
	$(CPP) $(CPPFLAGS) -c matrix.cpp

readpts.o: readpts.cpp vect.h global.h parse.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c readpts.cpp

blur.o: blur.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c blur.cpp

eval_haar.o: eval_haar.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c eval_haar.cpp

parse.o: parse.cpp parse.h
	$(CPP) $(CPPFLAGS) -c parse.cpp

plane.o: plane.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c plane.cpp

eval_daub4.o: eval_daub4.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c eval_daub4.cpp

average.o: average.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c average.cpp

octnode.o: octnode.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c octnode.cpp

coeffs_daub4.o: coeffs_daub4.cpp vect.h global.h common.h array2d.h
	$(CPP) $(CPPFLAGS) -c coeffs_daub4.cpp
