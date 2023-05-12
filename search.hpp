#ifndef search_hpp
#define search_hpp

#include "data.hpp"
#include <vector>
#include <algorithm>
#include <tuple>
#define ALPHA 0.5

#include "search.hpp"

//fonction de récupération de centre 
void getP(int * actual, std::vector<int> &K, int nnodes);

// récupération d'un sous-graphe à partir de la matrice 
int getSubsetNow(int * actual, std::vector<std::pair<int,double>> &subset, int nnodes, int k);

void updateRho(std::vector<std::pair<int,double>> &Q, const Data *data, int k);

int chooseID(std::vector<std::pair<int,double>> Q);

// tri sélection
void selectionSort(std::vector<std::pair<int,double>> &Q);

double getPhi(double distMax, int i, int k, const Data *data, double TG[]);

int getSubSetNext(int id, std::vector<int> P, const Data *data, double distMax, double TG[]);

double getDistMax(const Data* data);


void UpdateCenterSearch(int k, std::vector<int> &P, double TG[], int nnodes, int * actualInt, const Data *data);

double getW(std::vector<int> K, int * actual_solution, const Data * data);

// Recherche locale
double local_search(int * actual_solution, const Data *data);




#endif