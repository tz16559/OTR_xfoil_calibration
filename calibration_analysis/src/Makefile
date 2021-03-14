#-----------------------------------------
# This is for Linux
OSFLAGS = -DOS_LINUX -Dextname
LIBS = -lm -lz -lutil -lpthread -lpng

#-------------------------------------------------------------------
# The following lines define direcories. Adjust it if necessary.
#
#INC_DIR	= $(ROOTSYS)/include
INC_DIR	= ../findHoles 
INC_DIR2 =../differencePlotter

####################################################################
# Lines below here should not be edited
####################################################################

CFLAGS = -I$(INC_DIR) -I$(INC_DIR2)
CXXFLAGS = `root-config --cflags --glibs` -Wall
LDLIBS = `root-config --libs` -lHistPainter

OBJS = ../lib/FindHoles.so ../lib/differencePlotter.so

all: holeFinder

holeFinder: holeFinder.cpp
	$(CXX) $(OBJS) $(CFLAGS) $(CXXFLAGS) $(OSFLAGS) -o ../bin/$@ $< $(LDLIBS) $(LIBS)

clean:
	rm holeFinder
