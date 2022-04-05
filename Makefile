CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall -D_GLIBCXX_DEBUG
LDFLAGS=-g

all: raycasting

clean:
	rm -f *.o *.h.gch raycasting raycasting.bmp

test: raycasting
	./raycasting

.PHONY: all clean test

raycasting: raycasting.o
	$(CXX) $(LDFLAGS) -o $(@) $(^)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $(@) $(<)
