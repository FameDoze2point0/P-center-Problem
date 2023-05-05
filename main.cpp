#include <iostream>

#include "data.hpp"
#include "solver.hpp"

int main(int argc, char const *argv[]){

    if ( argc < 4){
        std::cout << "Usage : " << argv[0] << " <instance_path> <maximal_distance> <type of instance>" << std::endl;
        return 1;
    }   

    std::string instance = argv[1];
    Data pcenter;

    switch (*argv[3]){
    case 'p':
        pcenter.read_data_pmed(instance);
        std::cout << "Instance type : pmed" << std::endl;
        break;
    case 't':
        pcenter.read_data_tsp(instance);
        std::cout << "Instance type : tsp" << std::endl;
        break;
    case 'l':
        pcenter.read_data_lorena(instance);
        std::cout << "Instance type : Lorena" << std::endl;
        break;
    case 'g':
        pcenter.read_data_galvao(instance);
        std::cout << "Instance type : Galvao" << std::endl;
        break;
    case 'b':
        pcenter.read_data_beasley(instance);
        std::cout << "Instance type : Beasley" << std::endl;
        break;
    default:
        std::cout << "insert a valid instance type\n\t- p for pmed\n\t- t for tsp\n\t- l for Lorena\n\t- g for Galvao\n\t- b for Beasley" << std::endl;
        return 1;
    }

    Graph induced_graph;
    pcenter.create_induced_graph(induced_graph, 1e30);

    Solver cplex_solver(induced_graph, pcenter);
    cplex_solver.solve();

    return 0;
}