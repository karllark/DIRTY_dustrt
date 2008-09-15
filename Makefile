# Makefile for the Dirty code
#
# * 20 Mar 2008 Karl Misselt <misselt@as.arizona.edu> 
# - Made this the default Makefile
# added -g for debugging symbols to track possible leaks. 
# efence library ifn' you want to link against it. 

RM = /bin/rm -f
CPP = g++

INCDIR = include
#LIBS = -lcfitsio -lefence
LIBS = -lcfitsio

PROGRAM = dirty

SOURCES = $(wildcard *.cpp)

OBJECTS = ${SOURCES:.cpp=.o}
#CPPFLAGS = -I${INCDIR} -O2 -Wall -Wextra -g
CPPFLAGS = -I${INCDIR} -O2 -Wall -Wextra 

BUILDOBJ = $(CPP) $(CPPFLAGS)
BUILDLNK = $(CPP) $(LDFLAGS)

all:  ${PROGRAM}

$(PROGRAM) : $(OBJECTS)
	$(BUILDLNK) $^ $(LIBS) -o $(@)

-include $(SOURCES:%.cpp=%.d)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CPP) -MM $(CPPFLAGS) $(LDFLAGS) $< > $(DEPDIR)$@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $(DEPDIR)$@.$$$$ > $(DEPDIR)$@; \
	rm -f $@.$$$$

clean:
	$(RM) $(OBJECTS)
	$(RM) dirty

.cpp.o : 
	$(BUILDOBJ) -c $< -o $@