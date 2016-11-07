# Makefile for the Dirty code
#
# * 20 Mar 2008 Karl Misselt <misselt@as.arizona.edu> 
# - Made this the default Makefile
# added -g for debugging symbols to track possible leaks. 
# efence library ifn' you want to link against it. 

RM = /bin/rm -f
CPP = g++

# standard
INCDIR = include
LIBS = -lcfitsio

# STSCI
INCDIR_STSCI = /astro/dust_kg3/kgordon/Dirty_v2/cfitsio
#LIBS = -L/astro/dust_kg3/kgordon/Dirty_v2/cfitsio -lcfitsio

# Fedora with dist cfitsio
INCDIR_FEDORA = /usr/include/cfitsio

PROGRAM = dirty

SOURCES = $(wildcard *.cpp)

OBJECTS = ${SOURCES:.cpp=.o}

#CPPFLAGS = -I${INCDIR} -O2 -Wall -Wextra -g
CPPFLAGS = -I${INCDIR} -I${INCDIR_FEDORA} -O2 -Wall -Wextra 
#CPPFLAGS = -I${INCDIR} -I${INCDIR_STSCI} -O2 -Wall -Wextra 
#CPPFLAGS = -I${INCDIR} -Wall -Wextra 

LDFLAGS = ${LIBS}

BUILDOBJ = $(CPP) $(CPPFLAGS)
BUILDLNK = $(CPP)

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
	$(RM) *.d
	$(RM) dirty

.cpp.o : 
	$(BUILDOBJ) -c $< -o $@
