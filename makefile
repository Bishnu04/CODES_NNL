CC = gcc 
CXX = g++ 
CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

#TARGET1=     gcc  ## Bishnu's original
#OBJS1=       gcc.o

#TARGET1=     test2
#OBJS1=       test2.o
TARGET1=     test33
OBJS1=       test33.o

#TARGET2=    s0_test # on Jan 29,2019
#OBJS2=     s0_test.o # 0n jan 29,2019

TARGET2=    s0_test_original # on Jan 29,2019 just for  the test purpose
OBJS2=     s0_test_original.o # on Jan 29,2019cjust for  the test purpose

#TARGET2=     sample
#OBJS2=       sample.o

TARGET3=     trk_momentum  
OBJS3=       trk_momentum.o


#TARGET3=     gcc_test  
#OBJS3=       gcc_test.o

#TARGET4=     gcc_test2  
#OBJS4=       gcc_test2.o




all: $(TARGET1)  $(TARGET2)  $(TARGET3)\

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)


$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

#$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)



.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

