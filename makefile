CXXFLAGS = -Wall -g -Iinclude -Llib
LDFLAGS = -lfreeglut -lopengl32 -Wl,--subsystem,windows
OBJS = main.o graph.o sparsematrix.o

GraphMaker: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

main.o: src/main.cpp src/graph.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp

graph.o: src/graph.cpp src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/graph.cpp

sparsematrix.o: src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/sparsematrix.cpp

clean:
	-rm -f *.o core *.core
