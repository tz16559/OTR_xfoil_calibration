#-----------------------------------------
# This is for Linux
OSFLAGS = -DOS_LINUX -Dextname
LIBS = -lm -lz -lutil -lpthread -lpng

#-------------------------------------------------------------------
# The following lines define direcories. Adjust it if necessary.
#
#INC_DIR	= $(ROOTSYS)/include

####################################################################
# Lines below here should not be edited
####################################################################

CXXFLAGS = `root-config --cflags --glibs` -Wall
LDLIBS = `root-config --libs` -lHistPainter


all: stack_images

holeFinder: stack_images.cc
	$(CXX) $(OBJS) $(CFLAGS) $(CXXFLAGS) $(OSFLAGS) -o $@ $< $(LDLIBS) $(LIBS)

clean:
	rm stack_images
