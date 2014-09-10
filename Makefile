OBJS = testinstance2.o testset.o individual2.o population.o

FLAGS = -O2

default: 
	@echo "please choose mpi or singlenode"

singlenode: SingleNode.cpp SingleNode.h $(OBJS) config.h
	g++ $(FLAGS) -o singlenode $(OBJS) SingleNode.cpp
testset.o: testset.cpp testset.h testinstance2.o config.h
	g++ $(FLAGS) -c testset.cpp
testinstance2.o: testinstance2.cpp testinstance2.h config.h
	g++ $(FLAGS) -c testinstance2.cpp
population.o: population.cpp population.h individual2.o testset.o config.h
	g++ $(FLAGS) -c population.cpp
individual2.o: individual2.cpp individual2.h testset.o config.h
	g++ $(FLAGS) -c individual2.cpp

mpi mpiisland: MPIisland.cpp MPIisland.h $(OBJS) config.h
	/opt/mpich/gnu/bin/mpicxx $(FLAGS) -o mpiisland $(OBJS) MPIisland.cpp; \
	rm MPIisland.o

clean:
	rm $(OBJS)