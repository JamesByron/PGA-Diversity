OBJS = krktestinstance.o dataset.o individual2.o population.o

FLAGS = -O2

default: 
	@echo "please choose mpi or singlenode"

singlenode: SingleNode.cpp SingleNode.h krktestinstance.o $(OBJS) config.h
	g++ $(FLAGS) -o singlenode $(OBJS) SingleNode.cpp
dataset.o: dataset.cpp dataset.h config.h
	g++ $(FLAGS) -c dataset.cpp
krktestinstance.o: krktestinstance.cpp krktestinstance.h testinstance.o config.h
	g++ $(FLAGS) -c krktestinstance.cpp
testinstance.o: testinstance.cpp dataset.o individual2.o config.h
	g++ $(FLAGS) -c testinstance.cpp
population.o: population.cpp population.h individual2.o dataset.h config.h
	g++ $(FLAGS) -c population.cpp
individual2.o: individual2.cpp individual2.h config.h
	g++ $(FLAGS) -c individual2.cpp

mpi mpiisland: MPIisland.cpp MPIisland.h $(OBJS) config.h
	/opt/mpich/gnu/bin/mpicxx $(FLAGS) -o mpiisland $(OBJS) MPIisland.cpp; \
	rm MPIisland.o

clean:
	rm $(OBJS)