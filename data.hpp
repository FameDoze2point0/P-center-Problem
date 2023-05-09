#ifndef data_hpp
#define data_hpp

#include "structure.hpp"
#include "graph.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>


class Data{
    public:
        // Structs
        unsigned int nnodes;
        unsigned int nedges;
        unsigned int p;
        unsigned int perturbation;

        std::vector<edge> edges;
        std::vector<double> distances;
        std::vector<double> capacity;
        std::vector<double> demand;

        // Constructor
        Data(){
            nnodes = 0;
            nedges = 0;
            p = 0;
        }

        // Methods
        void read_data_pmed(const std::string &f);
        void read_data_tsp(const std::string &f);
        void read_data_lorena(const std::string &f);
        void read_data_galvao(const std::string &f);
        void read_data_beasley(const std::string &f);

        void create_induced_graph(Graph &induced_graph, double max_length);
};

#endif