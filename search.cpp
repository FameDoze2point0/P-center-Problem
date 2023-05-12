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

// void updateRho(std::vector<std::pair<int,double>> &Q, const Data *data, int k){
    
//     int nnodes = data->nnodes;
//     // calcul de la somme des distances entre k et chaque sommet du sous-graphe
//     double sum = 0.0;
//     for (auto s : Q){
//         sum += data->distances[s.first*nnodes + k];
//     }

//     for (size_t i = 0; i < Q.size(); i++){
//         Q[i].second = (double)data->distances[Q[i].first*nnodes + k] / sum;
//     }
// }

int updateRho(int * actual, std::tuple<int, double, bool, int> tab[] , const Data *data, int k, int nnodes){
    
    // calcul de la somme des distances entre k et chaque sommet du sous-graphe
    double sum = 0.0;
    int nb = 0; 
    for (size_t i = 0; i < nnodes; i++){
        if (std::get<3>(tab[i]) == k){
            sum += data->distances[i*nnodes + k];
            nb++;
        }
    }
    if ( nb != 1)
        for (size_t i = 0; i < nnodes; i++){
            if (std::get<3>(tab[i]) == k){
                std::get<1>(tab[i]) = (double) data->distances[i*nnodes + k] / sum;
            }
        }
    return nb;
}


// int chooseID(std::vector<std::pair<int,double>> Q){
    
//     int arg = 0;
//     double fin = Q.at(0).second, rng = (double)rand() / (double)RAND_MAX;
//     while (rng > fin){
//         arg++;
//         fin += Q.at(arg).second;
//     }
//     return arg;
// }

int idToDel(std::tuple<int, double,bool,int> tab[], int k, int nnodes){
    
    std::vector<int> R;
    for (size_t i = 0; i < nnodes; i++){
        if (std::get<3>(tab[i]) == k && i!=k){
            R.emplace_back(i);   
        }
    }
    
    double debut, rng;
    int arg;
    while (1){

        rng = (double)rand() / (double) RAND_MAX;
        debut = std::get<1>(tab[R.at(0)]);
        arg = 1;
        while ( debut < rng && arg < R.size()){
            debut += std::get<1>(tab[R.at(arg)]);
            arg++;
        }
        if( std::get<2>(tab[R.at(arg-1)]) == false) return R.at(arg-1);;
    }
}

// tri sélection
void selectionSort(std::vector<std::pair<int,int>> &W, const Data *data){
    int max;
    double distMax = 0;
    std::pair<int,int> tmp;
    for (size_t i = 0; i < W.size() - 1; i++){
        max = i;
        distMax = 0;
        for (size_t j = i+1; j < W.size(); j++){
            if (distMax < data->distances[ W.at(j).first * data->nnodes + W.at(j).second ]){
                max = j;
                distMax = data->distances[ W.at(j).first * data->nnodes + W.at(j).second ];
            }
        }
        tmp = W.at(i);
        W.at(i) = W.at(max);
        W.at(max) = tmp;
        // std::swap(Q.at(i), Q.at(max));
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


void UpdateCenterSearch(int k, std::vector<int> &P, double TG[], int nnodes, int * actualInt, const Data *data){
    std::tuple<int,double,double> T[nnodes]; // <id, capacited rest, distance max>
    int nb = 0;
    double demands = 0;
    // on récupère les sommets de X_k
    for (size_t i = 0; i < nnodes; i++){
        if (actualInt[nnodes+k*nnodes+i]){
            std::get<0>(T[nb]) = i;
            // get<1>(T[nb]) = data->demand[i]; 
            demands += data->demand[i];
            nb++;
        }
    }
    //on recupère le rayon maximal pour pour chaque sommet
    double max = 0;
    for (size_t i = 0; i < nb; i++){
        max = 0;
        for (size_t j = 0; j < nb; j++){
            if ( max < data->distances[std::get<0>(T[i])*nnodes+std::get<0>(T[j])])
                max = data->distances[std::get<0>(T[i])*nnodes+std::get<0>(T[j])]; 
        }
        std::get<2>(T[i]) = max;
    }

    // on recupère le sommet qui respecte les contraintes de capacité et ayant le plus petit rayon
    int arg = k;
    double dist = __DBL_MAX__;
    for (size_t i = 0; i < nb; i++){
        
        if ( data->capacity[i] - demands > 0 && dist > std::get<2>(T[i])){
            arg = std::get<0>(T[i]);
            dist = std::get<2>(T[i]);
        }
    }

    // if ( arg != k){
        for (size_t i = 0; i < nb; i++){
            actualInt[nnodes+k*nnodes+std::get<0>(T[i])] = 0;
            actualInt[nnodes+arg*nnodes+std::get<0>(T[i])] = 1;
        }
    // }
    
    // on enlève l'ancien centre de notre liste des centres, et on place le nouveau
    if ( arg != k){
        for (std::vector<int>::iterator it = P.begin(); it != P.end(); it++){
            if (*(it) == k){
                std::swap(*it, P.back());
                P.pop_back();
                break;
            }        
        }
        P.emplace_back(arg);
    }
    // on actualise les capacités
    TG[k] = data->capacity[k];
    TG[arg] = data->capacity[arg] - demands;
}

double getW(std::vector<int> K, int * actual_solution, const Data * data){

    int nnodes = data->nnodes;
    double max = 0;
    for ( auto k : K){
        for (size_t i = 0; i < nnodes; i++){
            
            if (actual_solution[nnodes + nnodes * k + i] && data->distances[k * nnodes + i] > max){
                max = data->distances[k* nnodes + i];
            }
        }
    }
    return max;
}




// // Recherche locale
// double local_search(int * actual_solution, const Data *data){
    
//     std::vector <std::pair<int, double>> Q; // (id, rho)
//     std::vector <int> K; // (id)
//     std::vector <std::pair<int, int>> W; // (k, c(l(k)) )
//     int nbElement = 0, del, arg, nnodes = data->nnodes;

//     // récupération des centres
//     getP(actual_solution, K, nnodes);

//     for ( auto k : K){

//         if (!Q.empty()) Q.clear();
//         if ((nbElement = getSubsetNow(actual_solution, Q, nnodes, k)) == 1) continue;
//         updateRho(Q, data, k);
//         del = 0;
//         while ( (double)del/ (double)(nbElement-1) < ALPHA ){
//             arg = chooseID(Q);
//             auto pair = std::make_pair(Q.at(arg).first, k);
//             if ( std::find(W.begin(), W.end(), pair) == W.end() ){
//                 del++;
//                 W.emplace_back(pair);
//             }
//         }
//     }
//     selectionSort(W,data);
//     // X <- X \ W
//     for (auto w : W){
//         actual_solution[nnodes + nnodes * w.second + w.first] = 0;
//     }
    
//     double TG[nnodes] = {0.0};
//     int k;
//     double distanceMax = getDistMax(data) + 1;
//     for (size_t i = 0; i < nnodes; i++){
//         TG[i] = data->capacity[i] - data->demand[i];
//     }
    
//     for (auto w : W){

//         k = getSubSetNext(w.first, K, data, distanceMax, TG);
//         actual_solution[nnodes + nnodes * k + w.first] = 1;
//     }

//     for (auto k : K){
//         UpdateCenterSearch(k,K,TG,nnodes, actual_solution, data);
//     }

//     return getW(K,actual_solution,data);
// }

int sortW1(const void * first, const void * second){
    bool firstBool =  std::get<2>(* (const std::tuple<int, double,bool,int> *) first);
    bool secondBool = std::get<2>(* (const std::tuple<int, double,bool,int> *) second);
    return secondBool - firstBool;
}

void bubleSort(std::tuple<int, double, bool, int> tab[], int w){

    std::tuple <int, double, bool, int> tmp;
    for (size_t i = w-1; i >= 0; i++){
        for (size_t j = 0; j < i-1; j++){
            if (std::get<1>(tab[j+1]) > std::get<1>(tab[j])){
                tmp = tab[j+1];
                tab[j+1] = tab[j];
                tab[j] = tmp;
            }
        }
    }   
}

double local_search_bis(int * actual_solution, const Data * data){

    int nnodes = data->nnodes, nbEl, arg, W = 0;
    std::vector<int> K;
    std::tuple<int, double, bool, int> tab[nnodes]; // rho, isDelete, c(l(k))

    getP(actual_solution,K,nnodes);    

    // initialisation
    for (size_t i = 0; i < nnodes; i++){
        
        std::get<0>(tab[i]) = i;
        std::get<1>(tab[i]) = 0.0;
        std::get<2>(tab[i]) = false;
        std::get<3>(tab[i]) = -1;

        for (size_t j = 0; j < nnodes; j++){
            if (actual_solution[nnodes + nnodes * j + i]){
                std::get<3>(tab[i]) = j;
                break;
            }
        }
    }

    for (auto k : K){
        
        std::cout << "====== " << k << " ======" << std::endl;
        nbEl = updateRho(actual_solution, tab, data, k, nnodes);
        if (nbEl == 1) continue;

        int del = 0;
        while ( (double)del/ (double)(nbEl-1) < ALPHA){
            arg = idToDel(tab, k, nnodes);
            std::get<2>(tab[arg]) = true;
            std::cout << "Je supprime " << arg << std::endl;
            del++;
            W++;
        }
    }
    
    qsort(tab, nnodes, sizeof(std::tuple<int, double,bool,int>), sortW1);

    for (size_t i = 0; i < nnodes; i++){
        std::cout << std::get<0>(tab[i]) << " " << std::get<2>(tab[i]) << " " << std::get<1>(tab[i]) << std::endl;
    }

    std::cout << "============" << W << std::endl;   

    bubleSort(tab,W);
    
    for (size_t i = 0; i < nnodes; i++){
        std::cout << std::get<0>(tab[i]) << " " << std::get<2>(tab[i]) << " " << std::get<1>(tab[i]) << std::endl;
    }
    
    
    return 0;
}







