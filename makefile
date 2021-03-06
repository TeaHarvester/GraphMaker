CXXFLAGS = -Wall -g -Iinclude -Llib
LDFLAGS = -lfreeglut -lopengl32 -lglew32 #-Wl,--subsystem,windows
OBJS = main.o graph.o graphicobject.o sparsematrix.o #renderer.o

GraphMaker: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

main.o: src/main.cpp src/graph.cpp src/graphicobject.cpp #src/renderer.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp

graph.o: src/graph.cpp src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/graph.cpp

graphicobject.o: src/graphicobject.cpp src/graph.cpp
	$(CXX) $(CXXFLAGS) -c src/graphicobject.cpp

sparsematrix.o: src/sparsematrix.cpp
	$(CXX) $(CXXFLAGS) -c src/sparsematrix.cpp

# renderer.o: src/renderer.cpp src/graphicobject.cpp
# 	$(CXX) $(CXXFLAGS) -c src/renderer.cpp

clean:
	-rm -f *.o core *.core
