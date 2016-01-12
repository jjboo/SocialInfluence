CC=g++
CFLAGS=-c -std=c++11 -O3 -w
LDFLAGS=
SOURCES=./src/anyoption.cpp ./src/graph.cpp ./src/mc.cpp ./src/infgraph.cpp ./src/timgraph.cpp ./src/comic.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=comic

all: $(SOURCES) $(EXECUTABLE)
	    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *o
