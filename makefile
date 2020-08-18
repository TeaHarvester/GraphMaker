CXXFLAGS = -Wall -g -Iinclude -Ifreeglut/include/GL -Iglew/include/GL -Lfreeglut/lib -Lglew/lib
LDFLAGS = -lfreeglut -lglew32
OBJS = main.o graph.o sparsematrix.o

GraphMaker: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $@

main.o: src/main.cpp src/graph.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp

graph.o: src/graph.cpp src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/graph.cpp

sparsematrix.o: src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/sparsematrix.cpp

clean:
	-rm -f *.o core *.core
