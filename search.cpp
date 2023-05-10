#include "search.hpp"

//fonction de récupération de centre 
void getP(int * actual, std::vector<int> &K, int nnodes){
    for (size_t i = 0; i < nnodes; i++){
        if (actual[i]){
            K.emplace_back(i);
        } 
    }
}

// récupération d'un sous-graphe à partir de la matrice 
int getSubsetNow(int * actual, std::vector<std::pair<int,double>> &subset, int nnodes, int k){
    int nbElement = 1;
    for (size_t i = 0; i < nnodes; i++){
        if (actual[nnodes + nnodes * k + i] && k != i){
            subset.emplace_back(std::make_pair(i,0.0));
            nbElement++;
        } 
    }
    return nbElement;
}

void updateRho(std::vector<std::pair<int,double>> &Q, const Data *data, int k){
    
    int nnodes = data->nnodes;
    // calcul de la somme des distances entre k et chaque sommet du sous-graphe
    double sum = 0.0;
    for (auto s : Q){
        sum += data->distances[s.first*nnodes + k];
    }

    for (size_t i = 0; i < Q.size(); i++){
        Q[i].second = (double)data->distances[Q[i].first*nnodes + k] / sum;
    }
}

int chooseID(std::vector<std::pair<int,double>> Q){
    
    int arg = 0;
    double fin = Q.at(0).second, rng = (double)rand() / (double)RAND_MAX;
    while (rng > fin){
        arg++;
        fin += Q.at(arg).second;
    }
    return arg;
}

// tri sélection
void selectionSort(std::vector<std::pair<int,double>> &Q){
    int max;
    std::pair<int,double> tmp;
    for (size_t i = 0; i < Q.size() - 1; i++){
        max = i;
        for (size_t j = i+1; j < Q.size(); j++){
            if (Q.at(j).second > Q.at(max).second){
                max = j;
            }
        }
        std::swap(Q.at(i), Q.at(max));
    }
}


double getPhi(double distMax, int i, int k, const Data *data, double TG[]){
    if ((double)data->demand[i] <= TG[k]){
        return (double)data->distances[i*data->nnodes+k] / distMax;
    }else{
        return (double)data->demand[i] - TG[k];
    }
}

int getSubSetNext(int id, std::vector<int> P, const Data *data, double distMax, double TG[]){

    int min = P.at(0);
    double phi = __DBL_MAX__, tmp;

    for (auto k : P){
        
        tmp = getPhi(distMax, id, k, data, TG );
        if(tmp < phi){
            phi = tmp;
            min = k;
        }
    }
    return min;
} 

double getDistMax(const Data* data){
    double max = 0.0;
    for (size_t i = 0; i < data->nnodes; i++){
        for (size_t j = i; j < data->nnodes; j++){
            if (data->distances[i*data->nnodes + j] > max){
                max = data->distances[i*data->nnodes + j];
            }
        }
    }
    return max;
}

// Recherche locale
void local_search(int * actual_solution, const Data *data){
    
    std::vector <std::pair<int, double>> Q; // (id, rho)
    std::vector <int> K; // (id)
    std::vector <std::pair<int, int>> W; // (k, c(l(k)) )
    int nbElement = 0, del, arg, nnodes = data->nnodes;

    // récupération des centres
    getP(actual_solution, K, nnodes);

    for ( auto k : K){

        Q.clear();
        if ((nbElement = getSubsetNow(actual_solution, Q, nnodes, k)) == 1) continue;
        updateRho(Q, data, k);
        del = 0;
        while ( del/(nbElement-1) < ALPHA ){
            arg = chooseID(Q);
            auto pair = std::make_pair(Q.at(arg).first, k);
            if ( std::find(W.begin(), W.end(), pair) == W.end() ){
                del++;
                W.emplace_back(pair);
            }
        }
    }

    selectionSort(Q);

    // X <- X \ W
    for (auto w : W){
        actual_solution[nnodes + nnodes * w.second + w.first] = 0;
    }
    
    double TG[nnodes] = {0.0};
    int k;
    double distanceMax = getDistMax(data) + 1;
    for (size_t i = 0; i < nnodes; i++){
        TG[i] = data->capacity[i] - data->demand[i];
    }
    

    for (auto w : W){

        k = getSubSetNext(w.first, K, data, distanceMax, TG);

    }



}











// VND