#include "data.hpp"
#include <cmath>

void Data::read_data_pmed(const std::string &str){

    std::ifstream f(str);

    if(!f.is_open()){
        std::cout << "Error opening file " << str << std::endl;
        exit(1);
    }

    f >> nnodes >> nedges >> p;

    std::vector <unsigned int> idx_edges(nnodes * nnodes, 0);
    size_t i, j, k, e;
    double c;
    unsigned int count = 0;
    for (size_t e_ = 0; e_ < nedges; e_++){
        
        f >> i >> j >> c;
        e = idx_edges[(i-1)*nnodes+j-1];
        if ( e > 0){
            edges[e-1].cost = c;
        }else{
            edges.push_back(edge(i, j, c));
            count++;
            idx_edges[(i-1)*nnodes+j-1] = count;
            idx_edges[(j-1)*nnodes+i-1] = count;
        }
    }
    f.close();
    nedges = count;

    capacity.resize(nnodes,1e30);
    demand.resize(nnodes,1e30);

    for (size_t i = 0; i < nnodes; i++){
        capacity[i] = nnodes;
        demand[i] = 1;
    }
    
    compute_distance(distances, edges, nnodes, nedges, 0);
    
}

void Data::read_data_tsp(const std::string &str){

    std::ifstream f(str);
    std::string s;

    do{
        f >> s;
    } while (s != "DIMENSION");
    f >> s;
    f >> nnodes;
    do{
        f >> s;
    } while (s != "NODE_COORD_SECTION");
    
    p = ceil(nnodes * 0.2);
    std::vector<std::pair<double, double>> coord(nnodes);
    size_t id;
    for (size_t n = 0; n < nnodes; n++){
        f >> id >> coord[n].first >> coord[n].second;
    }

    edges.reserve(nnodes * nnodes);
    std::vector<unsigned int> idx_edges(nnodes * nnodes, 0);
    unsigned int count = 0;
    for (size_t i = 0; i < nnodes; i++){

        for (size_t j = i+1; j < nnodes; j++){
            unsigned int cost = round(sqrt(pow(coord[j].first - coord[i].first, 2) + pow(coord[j].second - coord[i].second, 2)));
            edges.push_back(edge(i+1, j+1, cost));
            count++;
            idx_edges[i*nnodes+j] = count;
            idx_edges[j*nnodes+i] = count;
        }
    }
    f.close();

    nedges = count;
    compute_distance(distances, edges, nnodes, nedges, 0);
}
void Data::read_data_lorena(const std::string &str){

    std::ifstream f(str);
    f >> nnodes >> p;
    std::vector<std::pair<double, double>> coord(nnodes);
    capacity.resize(nnodes,1e30);
    demand.resize(nnodes,1e30);

    for (size_t i = 0; i < nnodes; i++){
        f >> coord[i].first >> coord[i].second >> capacity[i] >> demand[i];
    }
    f.close();
    distances.resize(nnodes*nnodes, 1e30);

    for (size_t i = 0; i < nnodes; i++)
    {
        for (size_t j = 0; j < nnodes; j++)
        {
            distances[i*nnodes+j] = round(sqrt(pow(coord[i].first - coord[j].first, 2) + pow(coord[i].second - coord[j].second, 2)));
        }
    }
    nedges = nnodes * nnodes;

    
}
void Data::read_data_galvao(const std::string &str){

    std::ifstream f(str);
    f >> nnodes >> p;

    capacity.resize(nnodes,1e30);
    demand.resize(nnodes,1e30);
    distances.resize(nnodes*nnodes, 1e30);

    for (size_t i = 0; i < nnodes; i++){
        f >> capacity[i];
    }
    for (size_t i = 0; i < nnodes; i++){
        f >> demand[i];
    }
    for (size_t i = 0; i < nnodes; i++){
        for (size_t j = 0; j < nnodes; j++){
            f >> distances[i*nnodes+j];
        }
    }
    f.close();
}

void Data::read_data_beasley(const std::string &str){

    std::ifstream f(str);
    int best = 0, c, e;
    f >> best;
    f >> nnodes >> p >> c;
    std::vector <std::pair<double, double>> coord(nnodes);

    capacity.resize(nnodes,1e30);
    demand.resize(nnodes,1e30);
    distances.resize(nnodes*nnodes, 1e30);

    for (size_t i = 0; i < nnodes; i++){
        f >> e >> coord[i].first >> coord[i].second >> demand[i];
        capacity[i] = c;
    }

    for (size_t i = 0; i < nnodes; i++){

        for (size_t j = 0; j < nnodes; j++){
            distances[i*nnodes+j] = round(sqrt(pow(coord[i].first - coord[j].first, 2) + pow(coord[i].second - coord[j].second, 2)));
        }
    }
}


void Data::create_induced_graph(Graph &induced_graph, double max_length){

    induced_graph.nodes.resize(nnodes);
    induced_graph.idx_edges.resize(nnodes*nnodes, 0);
    induced_graph.idx_nodes.resize(nnodes, 0);

    induced_graph.nnodes = nnodes;
    induced_graph.current_p = p;

    for (size_t i = 0; i < nnodes; i++)
    {
        induced_graph.idx_nodes[i] = i;
        induced_graph.nodes[i].id = i+1;
        for (size_t j = 0; j < nnodes; j++)
        {
            if (distances[i*nnodes+j] != distances[j * nnodes+i]){
                std::cerr << "Error: distance matrix is not symmetric" << std::endl;
                std::abort();
            }
            if (distances[i*nnodes+j] <= max_length){
                induced_graph.edges.push_back(edge(i+1, j+1, distances[i*nnodes+j]));
                ++induced_graph.nedges;
                induced_graph.idx_edges[i*nnodes+j] = induced_graph.nedges;
                induced_graph.idx_edges[j*nnodes+i] = induced_graph.nedges;
            }
        }
    }
}
