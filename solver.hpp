#ifndef solver_hpp
#define solver_hpp

#include "data.hpp"

#include "ilcplex/ilocplex.h"

#define HISTORIC_SIZE 3

class Solver {
    public:
        unsigned int nnodes, ncenters;
        const Graph *graph;
        const Data *data;

        int * history[HISTORIC_SIZE];
        int * actual_solution;

        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloNumVarArray x,y;
        IloNumVar w;
        IloObjective obj;



        Solver(const Graph & graph_): cplex(env), model(env), y(env){
            graph = &graph_;
            nnodes = graph->nnodes;
            ncenters = graph->current_p;
        }

        Solver(const Graph & graph_, const Data & pcenter_): cplex(env), model(env), y(env), obj(env){
            graph = &graph_;
            nnodes = graph->nnodes;
            ncenters = graph->current_p;
            data = &pcenter_;
            for (int i = 0; i < HISTORIC_SIZE; i++){
                history[i] = new int[nnodes]();
            }
            actual_solution = new int[nnodes]();
        }

        void set_parameters();
        void create_model();
        void solve();
        void freeMemory();
};

#endif