#ifndef structure_hpp
#define structure_hpp

#include <vector>
#include <stddef.h>

struct edge{
    // contains the first and second node of the edge and the cost
    unsigned int i, j;
    double cost;

    edge(){
        i = 0;
        j = 0;
        cost = 0;
    }
    
    edge(unsigned int i, unsigned int j, double cost){
        this->i = i;
        this->j = j;
        this->cost = cost;
    }
    edge(const edge& e){
        this->i = e.i;
        this->j = e.j;
        this->cost = e.cost;
    }
    edge& operator=(const edge& e){
        this->i = e.i;
        this->j = e.j;
        this->cost = e.cost;
        return *this;
    }
};

struct node{
    unsigned id, idcenter;

    node(){
        id = 0;
        idcenter = 0;
    }
    node(unsigned int id){
        this->id = id;
        idcenter = 0;
    }
    node(unsigned int id, unsigned int idcenter){
        this->id = id;
        this->idcenter = idcenter;
    }
    node(const node& n){
        this->id = n.id;
       this->idcenter = n.idcenter;
    }
    node& operator=(const node& n){
        this->id = n.id;
        this->idcenter = n.idcenter;
        return *this;
    }
};

void compute_distance(std::vector<double>& distances, const std::vector<edge>& edges, size_t nnodes, size_t nedges, double unique_cost);
#endif