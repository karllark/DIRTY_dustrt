# Makefile for the Dirty code
#
# * 20 Mar 2008 Karl Misselt <misselt@as.arizona.edu> 
# - Made this the default Makefile
# added -g for debugging symbols to track possible leaks. 
# efence library ifn' you want to link against it. 

RM = /bin/rm -f
CPP = g++

INCDIR = include 
#INCDIR_CFITSIO = /user/kgordon/MIPS/Pipeline/DAT/lib/CFitsIO/
# Fedora with dist cfitsio
#INCDIR_FEDORA = /usr/include/cfitsio
LIBS = -lcfitsio
#LIBS = -L/grp/software/Linux/RH5/x86_64/cfitsio/lib -lcfitsio
#LIBS = -L/user/kgordon/MIPS/Pipeline/DAT/lib/CFitsIO/ -lcfitsio

PROGRAM = write_indiv_grain_data

SOURCES = $(wildcard *.cpp)

OBJECTS = ${SOURCES:.cpp=.o}
CPPFLAGS = -I${INCDIR} -I${INCDIR_CFITSIO} -O2 -Wall -Wextra -g
#CPPFLAGS = -I${INCDIR} -I${INCDIR_FEDORA} -O2 -Wall -Wextra 
#CPPFLAGS = -I${INCDIR} -Wall -Wextra 

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
	$(RM) *.d
	$(RM) dirty

.cpp.o : 
	$(BUILDOBJ) -c $< -o $@
