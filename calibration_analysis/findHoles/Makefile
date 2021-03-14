#-----------------------------------------
# This is for Linux
OSFLAGS = -DOS_LINUX -Dextname
CFLAGS = -c
LIBS =  -shared

#-------------------------------------------------------------------
# The following lines define direcories. Adjust it if necessary.
#                 

HDRS = FindHoles.h 
SRCS = FindHoles.cpp

all: FindHoles.so

FindHoles.so: $(HDRS) $(SRCS) 
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o ../lib/$@ $(SRCS) $(LIBS)

clean:
	rm ../lib/FindHoles.so
