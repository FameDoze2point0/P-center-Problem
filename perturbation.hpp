#ifndef perturbation_hpp
#define perturbation_hpp

#include <iostream>
#include "data.hpp"
#include "ilcplex/ilocplex.h"
#include <algorithm>
#include <string.h>

void choosePerturbation(int num, IloNumArray y_, const Data *data,  int * actual_solution, int &k);

void perturbation1(IloNumArray y_, const Data *data, int * actual_solution);

void perturbation2(IloNumArray y_, const Data *data, int * actual_solution, int k);


int chooseCenter(double piArr[]);
double getDistMax(const Data *data);
double getPhi(double distMax, int i, int k, double restArr[], const Data *data);
int getSubset(double distMax, int id, std::vector<int> &P, double restArr[], const Data *data);
void updateCenter(int k, std::vector<int> &P, double restArr[], int * actualSolution, const Data *data);
void updateGamma(bool isCenter[], double gammaArr[], int id, const Data *data);
void updatePi(bool isCenter[], double piArr[], double gammaArr[], const Data *data);
void perturbation3(int * actualInt, int nnodes, int ncenters, const Data *pcenter_data);

#endif