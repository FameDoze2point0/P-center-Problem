#ifndef graph_hpp
#define graph_hpp

#include "structure.hpp"

class Graph{
    public:
        unsigned int nnodes;
        unsigned int nedges;
        unsigned int current_p;

        std::vector<edge> edges;
        std::vector<node> nodes;

        std::vector<unsigned int> idx_edges;
        std::vector<unsigned int> idx_nodes;

        Graph(){
            nnodes = 0;
            nedges = 0;
            current_p = 0;
        }
        Graph(const Graph &g){
            nnodes = g.nnodes;
            nedges = g.nedges;
            current_p = g.current_p;
            edges = g.edges;
            nodes = g.nodes;
        }
};

#endif