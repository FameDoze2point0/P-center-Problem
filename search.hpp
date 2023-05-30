#ifndef search_hpp
#define search_hpp

#include "data.hpp"
#include <vector>
#include <algorithm>
#include <tuple>
#define ALPHA 0.1

#include "search.hpp"


void getK(int * actual, std::vector<int> &K, int nnodes);

void getQ(int * actual, std::vector<int> &Q, int nnodes, int k);

void updateRho(std::vector<int> Q, std::vector<double> &rho, int k, const Data * data);

int selectNextDelete(std::vector<int> Q, std::vector<double> rho);

void algoTri(std::vector<std::pair<int,int>> &W, const Data * data);

double getDistMaxSearch(int nnodes, const Data *pcenter_data);

double getPhi(double distMax, int i, int k, double TG[], const Data *data);

int getSubSet(double distMax,int w, std::vector<int> K,double rest[], const Data * data);


double local_search(int * actual_solution, const Data * data);




#endif