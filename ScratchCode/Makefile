#LIBPATH=/home/misselt/Pro/cpp
#LIBPATH2=./lib

NP_SRCS = ce_idl.cpp \
	  ../StringManip.cpp \
          ../ConfigFile.cpp \
	  ../Grain.cpp \
	  ComputeDustEmission.cpp \
	  ISRF.cpp \
	  StochasticHeating.cpp \
          ../GrainModel.cpp

NP_OBJS = ${NP_SRCS:.cpp=.o}

RM = /bin/rm -f 

#SOURCEDIR=./
LIBDIR=./include
LIBDIR2=../include

#CCFLAGS = -I$(LIBDIR)/lib -I$(LIBDIR) -I$(LIBPATH)/$(LIBDIR) -O3 -pg 
CCFLAGS = -I$(LIBDIR) -I$(LIBDIR2) -O3 -pg -g
#LDFLAGS = -O3 -L$(LIBDIR) -lConfigFile -lDataFile -lNumUtils -lStringManip -lcfitsio
LDFLAGS = -lcfitsio -pg 
#LDFLAGS = -O3 -L$(LIBDIR) -l90Prime -lcfitsio

CC = g++

all: ce_idl

ce_idl: $(NP_OBJS)
	$(CC) -o $(@) $(NP_OBJS) $(CCFLAGS) $(LDFLAGS)

clean: 
	$(RM) $(NP_OBJS) ce_idl

.cpp.o : $(HDS)
	$(CC) -o ${*}.o $(CCFLAGS) -c $<
