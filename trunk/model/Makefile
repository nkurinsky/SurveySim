LIBS = -lCCFits -lcfitsio -lgsl -lgslcblas -lm -lm -lgslcblas -lgsl -lcfitsio -lCCFits

includedir = /usr/local/include
libdir = /usr/local/lib

CXX := clang++
SOURCES = test_model.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXES = $(OBJECTS:.o=)

all: $(OBJECTS) $(EXES)

%: %.o
	$(CXX) -I/usr/local/include -o $@ $< $(LIBS)

clean:
	-rm -f $(EXES) $(OBJECTS) *~