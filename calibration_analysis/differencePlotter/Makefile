#-----------------------------------------
# This is for Linux
OSFLAGS = -DOS_LINUX -Dextname
CFLAGS = -c
LIBS =  -shared

#-------------------------------------------------------------------
# The following lines define direcories. Adjust it if necessary.
#                 

HDRS = differencePlotter.h 
SRCS = differencePlotter.cpp
CXXFLAGS=`root-config --cflags --glibs`

all: differencePlotter.so

differencePlotter.so: $(HDRS) $(SRCS) 
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o ../lib/$@ $(SRCS) $(LIBS)

clean:
	rm ../lib/differencePlotter.so
