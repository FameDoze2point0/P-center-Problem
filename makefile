# compiler
CC = g++ -std=c++11

# module name
NAME = pcenter

# basic directory
DIR = ./

# debug switches
#SW = -Wall -ggdb3
SW = -w -O3

# default target - - - - - - - - - - - - - - - - - - - - - -
default : $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - -
clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define include the necessary libs - - - - - - - - - - - -

CONCERTDIR = /home/antonin/Bureau/IBM_cplex/concert
CPLEXDIR = /home/antonin/Bureau/IBM_cplex/cplex

SYSTEM = x86-64_linux
LIBFORMAT = static_pic

CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINC = -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include

# Flags - - - - - - - - - - - - - - - - - - - - - - - - - -
CPLEXFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl


# makefile
OBJ = main.o data.o structure.o graph.o solver.o perturbation.o search.o
./main.o: ./main.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW)
solver.o: solver.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW)
data.o: data.cpp
	$(CC) -c $< -o $@ $(SW)
structure.o: structure.cpp
	$(CC) -c $< -o $@ $(SW)
graph.o: graph.cpp
	$(CC) -c $< -o $@ $(SW)
perturbation.o: perturbation.cpp
	$(CC) -c $< -o $@ $(CPLEXINC) $(SW)
search.o: search.cpp
	$(CC) -c $< -o $@ $(SW)

$(NAME): $(OBJ)
	$(CC) -o $(NAME) $(OBJ) $(CPLEXFLAGS) $(SW)	