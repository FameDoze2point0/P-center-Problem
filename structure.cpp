#include "structure.hpp"

void compute_distance(std::vector<double>& distances, const std::vector<edge>& edges, size_t nnodes, size_t nedges, double unique_cost){

    size_t i, j, k;
    double cost;
    distances.resize(nnodes*nnodes, 1e30);
    for (size_t e = 0; e < nedges; e++){
        i = edges[e].i - 1;
        j = edges[e].j - 1;

        if(unique_cost > 0) cost = unique_cost;
        else cost = edges[e].cost;

        distances[i*nnodes + j] = cost;
        distances[j*nnodes + i] = cost;
    }

    for (size_t n = 0; n < nnodes; n++)
        distances[n*nnodes+n] = 0;
    
    for (k = 0; k < nnodes; k++){
        for ( i = 0; i < nnodes; i++){
            for ( j = 0; j < nnodes; j++){
                if (distances[i*nnodes+j] > distances[i*nnodes+k] + distances[k*nnodes+j])
                    distances[i*nnodes+j] = distances[i*nnodes+k] + distances[k*nnodes+j];
            }
        }        
    }
    
    
    
}