# Makefile for the Dirty code
#
# * 20 Mar 2008 Karl Misselt <misselt@as.arizona.edu> 
# - Made this the default Makefile

RM = /bin/rm -f
CPP = g++

INCDIR = include
LIBS = -lcfitsio

PROGRAM = dirty

SOURCES = $(wildcard *.cpp)

OBJECTS = ${SOURCES:.cpp=.o}
CPPFLAGS = -I${INCDIR} -O3 

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
