CFLAG = -Wall -g
OBJS = main.o graph.o sparsematrix.o

GraphMaker: $(OBJS)
	$(CXX) $(CFLAG) $(OBJS) -o $@

main.o: main.cpp graph.cpp
	$(CXX) $(CFLAG) -c main.cpp

graph.o: graph.cpp sparsematrix.cpp
	$(CXX) $(CFLAG) -c graph.cpp

sparsematrix.o: sparsematrix.cpp
	$(CXX) $(CFLAG) -c sparsematrix.cpp

clean:
	-rm -f *.o core *.core
