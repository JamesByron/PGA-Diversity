OBJS = krktestinstance.o dataset.o individual2.o population.o

FLAGS = -O2

default: 
	@echo "please choose mpi or singlenode"

singlenode: SingleNode.cpp SingleNode.h $(OBJS) config.h
	g++ $(FLAGS) -o singlenode $(OBJS) SingleNode.cpp
dataset.o: dataset.cpp dataset.h krktestinstance.o config.h
	g++ $(FLAGS) -c dataset.cpp
krktestinstance.o: krktestinstance.cpp krktestinstance.h config.h
	g++ $(FLAGS) -c krktestinstance.cpp
population.o: population.cpp population.h individual2.o dataset.o config.h
	g++ $(FLAGS) -c population.cpp
individual2.o: individual2.cpp individual2.h dataset.o config.h
	g++ $(FLAGS) -c individual2.cpp

mpi mpiisland: MPIisland.cpp MPIisland.h $(OBJS) config.h
	/opt/mpich/gnu/bin/mpicxx $(FLAGS) -o mpiisland $(OBJS) MPIisland.cpp; \
	rm MPIisland.o

clean:
	rm $(OBJS)