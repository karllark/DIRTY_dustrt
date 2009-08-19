#LIBPATH=/home/misselt/Pro/cpp
#LIBPATH2=./lib

NP_SRCS = test_de.cpp \
	  ISRF.cpp \
	  ../DirtyFailure.cpp \
	  ../StringManip.cpp \
          ../ConfigFile.cpp \
	  ../Grain.cpp \
	  ../ComputeDustEmission.cpp \
	  ../StochasticHeating.cpp \
          ../GrainModel.cpp

NP_OBJS = ${NP_SRCS:.cpp=.o}

RM = /bin/rm -f 

#SOURCEDIR=./
LIBDIR=./include
LIBDIR2=../include
INCDIR_FEDORA = /usr/include/cfitsio

#CCFLAGS = -I$(LIBDIR)/lib -I$(LIBDIR) -I$(LIBPATH)/$(LIBDIR) -O3 -pg 
CCFLAGS = -I$(LIBDIR) -I$(LIBDIR2) -I$(INCDIR_FEDORA) -O3 -pg -g
#LDFLAGS = -O3 -L$(LIBDIR) -lConfigFile -lDataFile -lNumUtils -lStringManip -lcfitsio
LDFLAGS = -lcfitsio -pg 
#LDFLAGS = -O3 -L$(LIBDIR) -l90Prime -lcfitsio

CC = g++

all: testeh

testeh: $(NP_OBJS)
	$(CC) -o $(@) $(NP_OBJS) $(CCFLAGS) $(LDFLAGS)

clean: 
	$(RM) $(NP_OBJS) testeh

.cpp.o : $(HDS)
	$(CC) -o ${*}.o $(CCFLAGS) -c $<
