CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall  -pedantic-errors -g -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing
#-O3 -g -Wall -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing
SDIR = $(wildcard ../src)

SRCS = $(SDIR)/*.cpp
OBJS = ${SRCS:%.cpp=${CURDIR}/*.o}
OBJD = ${SRCS:%.cpp=%.o}

MAIN = main

all: ${MAIN}
	echo   Compilation Successful.................

${MAIN}: ${OBJD}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c ${SRCS}

clean:
	${RM} ${PROGS} ${OBJS} *.o *~


	
