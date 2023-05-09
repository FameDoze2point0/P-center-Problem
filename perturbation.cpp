#include "perturbation.hpp"

void choosePerturbation(int num, IloNumArray y_, const Data *data, int * actual_solution, int & k){
    switch (num){
    case 1:
        perturbation1(y_, data, actual_solution);
        break;
    case 2:
        perturbation2(y_, data, actual_solution, k);
        break;
    case 3:
        perturbation3(actual_solution, data->nnodes, data->p, data);
        break;
    default:
        std::cout << "This perturbation doesn't exist" << std::endl;
        break;
    }
}
/************************************* Perturbation 1 *******************************************/

/// This perturbation flips a random number of bits in the solution
void perturbation1(IloNumArray y_, const Data *data, int * actual_solution){
    int nnodes = data->nnodes;
    int indexes[nnodes];

    for (size_t i = 0; i < nnodes; i++){
        indexes[i] = i;
    }
    
    std::random_shuffle(indexes, indexes + nnodes);

    int flip = nnodes - rand() % (2 * nnodes / 3), id;
    double rng;
    for (size_t i = 0, j = 0; i < nnodes && j < flip ; i++){
        rng = std::max(0.0, (double)rand() / (double)RAND_MAX - 0.3);
        id = indexes[i];
        if ( std::abs(y_[id] - actual_solution[id]) + rng > 0.5){
            actual_solution[id] ^= 1;
            j++;
        }
    }
}

/************************************* Perturbation 2 *******************************************/

/// This perturbation flips k bits in the solution, and then flips the bits that are at 1 until we have ncentres centers
void perturbation2(IloNumArray y_, const Data *data, int * actual_solution, int k){

    int nnodes = data->nnodes;
    int indexes[nnodes];
    // Indexes va contenir à gauches les indices des noeuds qui sont à 0 dans la solution, et à droite ceux qui sont à 1

    int left = 0, right = nnodes - 1;
    for (size_t i = 0; i < nnodes; i++){
        if (actual_solution[i] == 0){
            indexes[left] = i;
            left++;
        }
        else{
            indexes[right] = i;
            right--;
        }
    }

    // on shuffle chaque côté de manière indépendante
    std::random_shuffle(indexes, indexes + left);
    std::random_shuffle(indexes + left, indexes + nnodes);

    // on va transformer k indices dont la solution est à 0 en 1 
    // si k > left, on va transformer left indices
    int flip = std::min(k, left), id;

    // on flip les k premiers indices dont la solution est à 0 à partir de left
    for (size_t i = left - 1, j = 0; i >= 0 && j < flip; i--){
        id = indexes[i];
        actual_solution[id] = 1;
        j++;
        left--;
        right--;
    }

    // on flip les entiers dont la solution est à 1 tant jusqu'à avoir ncentres centres
    while ( nnodes - right - 1 < data->p){
        id = indexes[right];
        actual_solution[id] = 0;
        right++;
    }

    k += nnodes/10; // k évolutif, il augmente lorsqu'il y a plusieurs perturbations à la suite
    std::cout << "k = " << k << std::endl;
}

/************************************* Perturbation 3 *******************************************/





void updateGamma(bool T[],double TG[],int i, int nnodes, const Data *pcenter_data){
    for (size_t j = 0; j < nnodes; j++)
        if (!T[j])
            TG[j] = pcenter_data->distances[i*nnodes+j] * pcenter_data->capacity[j];
}

void updatePi(bool T[], double TP[], double TG[], int nnodes){
    
    double sum = 0;

    for (size_t j = 0; j < nnodes; j++){
        if (!T[j])
            sum += TG[j];
    }
    
    for (size_t j = 0; j < nnodes; j++){
        if (!T[j])
            TP[j] = TG[j] / sum;
    }
}

int chooseId(double TP[]){
    int arg = 0;
    double fin = TP[0], rng = (double)rand()/RAND_MAX;

    while (fin < rng){
        arg++;
        fin += TP[arg];
    }
    return arg;
}

double getDistMax(int nnodes, const Data *pcenter_data){
    double dist = 0;
    for (size_t i = 0; i < nnodes; i++){
        for (size_t j = i; j < nnodes; j++){
            if (pcenter_data->distances[i*nnodes+j] > dist){
                dist = pcenter_data->distances[i*nnodes+j];
            } 
        }
    }
    return dist;
}

double getPhi(double distMax, int i, int k, double TG[], const Data *pcenter_data, int nnodes){
    if ((double)pcenter_data->demand[i] <= TG[k]){
        return (double)pcenter_data->distances[i*nnodes+k] / distMax;
    }else{
        return (double)pcenter_data->demand[i] - TG[k];
    }
}

int getSubSet(double distMax, int i, std::vector<int> P, double TG[], const Data *pcenter_data, int nnodes){
    int min = P.at(0);
    double phi = __DBL_MAX__, tmp;

    for (auto k : P){
        
        tmp = getPhi(distMax, i, k, TG,pcenter_data, nnodes );
        if(tmp < phi){
            phi = tmp;
            min = k;
        }
    }
    return min;
} 

void UpdateCenter(int k, std::vector<int> &P, double TG[], int nnodes, int * actualInt, const Data *pcenter_data){
    std::tuple<int,double,double> T[nnodes]; // <id, capacited rest, distance max>
    int nb = 0;
    double demands = 0;
    // on récupère les sommets de X_k
    for (size_t i = 0; i < nnodes; i++){
        if (actualInt[nnodes+k*nnodes+i]){
            std::get<0>(T[nb]) = i;
            // get<1>(T[nb]) = pcenter_data->demand[i]; 
            demands += pcenter_data->demand[i];
            nb++;
        }
    }
    //on recupère le rayon maximal pour pour chaque sommet
    double max = 0;
    for (size_t i = 0; i < nb; i++){
        max = 0;
        for (size_t j = 0; j < nb; j++){
            if ( max < pcenter_data->distances[std::get<0>(T[i])*nnodes+std::get<0>(T[j])])
                max = pcenter_data->distances[std::get<0>(T[i])*nnodes+std::get<0>(T[j])]; 
        }
        std::get<2>(T[i]) = max;
    }

    // on recupère le sommet qui respecte les contraintes de capacité et ayant le plus petit rayon
    int arg = k;
    double dist = __DBL_MAX__;
    for (size_t i = 0; i < nb; i++){
        
        if ( pcenter_data->capacity[i] - demands > 0 && dist > std::get<2>(T[i])){
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
    TG[k] = pcenter_data->capacity[k];
    TG[arg] = pcenter_data->capacity[arg] - demands;
}

void perturbation3(int * actualInt, int nnodes, int ncenters, const Data *pcenter_data){
    
    //reset : X <- ensemble vide
    for (size_t i = 0; i < nnodes*nnodes+nnodes; i++)
        actualInt[i] = 0;
    
    std::vector<int> P;
    bool T[nnodes] = {0};
    double TG[nnodes] = {0}, TP[nnodes] = {0};

    // choose i* in V randomly
        int id = rand()%nnodes;
    
    // P <- P U {i*}
        T[id] = actualInt[id] = actualInt[nnodes + id*nnodes+id] = 1; 
        P.emplace_back(id);
    
    // Update gamme and pi
        updateGamma(T,TG,id, nnodes, pcenter_data);
        updatePi(T,TP,TG, nnodes);

    // stage (a)
    for (size_t cpt = 1; cpt < ncenters; cpt++){ // while |P| < p
        // Choose i* in V\P randomly using the probability pi
            id = chooseId(TP);
        // P <- P U {i*}
            P.emplace_back(id);
            T[id] = actualInt[id] = actualInt[nnodes + id*nnodes+id] = 1; 
            TG[id] = TP[id] = 0;
        //Update gamme and pi
            updateGamma(T,TG,id, nnodes, pcenter_data);
            updatePi(T,TP,TG, nnodes);
    }

    // initialisation des capacité (on part du principe qu'ils sont auto suffisant)
    for (size_t i = 0; i < nnodes; i++)
        TG[i] = pcenter_data->capacity[i] - pcenter_data->demand[i];


    // stage (b)
    int k;
    double distanceMax = getDistMax(nnodes, pcenter_data) + 1; // d barre
    for (size_t i = 0; i < nnodes; i++){ // for all  j in V \ P
        if (!T[i]){
            // k* <- arg min_{k from K} phi(j,k)
                k =  getSubSet(distanceMax, i, P, TG, pcenter_data, nnodes);
            // X_k* <- X_k* U {j}
                T[i] = actualInt[nnodes+k*nnodes+i] = 1;
            // Update C(k*)
                UpdateCenter(k,P, TG, nnodes, actualInt, pcenter_data);
        }
    }
}












































































































// void updateGamma(bool isCenter[], double gammaArr[], int id, const Data *data){
//     int nnodes = data->nnodes;

//     for (size_t i = 0; i < nnodes; i++){
//         if (!isCenter[i]){
//             gammaArr[i] = data->distances[id*nnodes + i] * data->capacity[i];
//         }
//     }
// }

// void updatePi(bool isCenter[], double piArr[], double gammaArr[], const Data *data){
//     int nnodes = data->nnodes;

//     double sum = 0;

//     // récupération de la somme des gamma
//     for (size_t i = 0; i < nnodes; i++){
//         if (!isCenter[i]){
//             sum += gammaArr[i];
//         }
//     }

//     // calcul des pi (soit gamma normalisé)
//     for (size_t i = 0; i < nnodes; i++){
//         if (!isCenter[i]){
//             piArr[i] = gammaArr[i] / sum;
//         }
//     }
// }

// int chooseCenter(double piArr[]){
//     int id = 0;
//     double fin = piArr[0], rng = (double)rand() / (double)RAND_MAX;

//     while (fin < rng){
//         id++;
//         fin += piArr[id];
//     }
//     return id; 
// }

// double getDistMax(const Data *data){

//     int nnodes = data->nnodes;
//     double max = 0;
//     for (size_t i = 0; i < nnodes; i++){
//         for (size_t j = 0; j < nnodes; j++){
            
//             if (data->distances[i*nnodes + j] > max){
//                 max = data->distances[i*nnodes + j];
//             }
//         }
//     }
//     return max;
// }

// double getPhi(double distMax, int i, int k, double restArr[], const Data *data){
//     if ((double)data->demand[i] <= restArr[k]){
//         return (double)data->distances[i*data->nnodes + k] / distMax;
//     }else{
//         return (double) data->demand[i] - restArr[k];
//     }
// }

// int getSubset(double distMax, int id, std::vector<int> &P, double restArr[], const Data *data){

//     int min = P.at(0);
//     double phi = __DBL_MAX__, tmp;

//     for (auto k : P){
//         tmp = getPhi(distMax, id, k, restArr, data);
//         if (tmp < phi){
//             phi = tmp;
//             min = k;
//         }
//     }
//     return min;
// }

// void updateCenter(int k, std::vector<int> &P, double restArr[], int * actualSolution, const Data *data){

//     int nnodes = data->nnodes;
//     std::vector<int> nodesX;
//     std::vector<double> radiusMax; // donne le rayon max pour chaque noeud de X
//     double demands = 0;

//     // récupération des noeuds dans l'ensemble et de leur demande
//     for (size_t i = 0; i < nnodes; i++){
//         if ( actualSolution[nnodes + k*nnodes+i]){
//             nodesX.emplace_back(i);
//             demands += data->demand[i];
//         }
//     }

//     // on récupère ensuite le rayon maximal de chaque noeud dans l'ensemble
//     double max = 0;
//     for (auto i : nodesX){
//         max = 0;
//         for (auto j : nodesX){
//             if (data->distances[i*nnodes + j] > max){
//                 max = data->distances[i*nnodes + j];
//             }
//         }
//         radiusMax.push_back(max);
//     }    

//     // on récupère le sommet qui respecte les contraintes de capacité et qui a le rayon maximal le plus petit
//     int arg = k;
//     double dist = __DBL_MAX__;
//     for (size_t i = 0; i < nodesX.size(); i++){
//         if ( data->capacity[i] - demands >= 0 && dist > radiusMax[i]){
//             dist = radiusMax[i];
//             arg = nodesX[i];
//         }
//     }
    
//     // on inverse les liens entre l'ancien k et le nouveau k
//     for (size_t i = 0; i < nodesX.size(); i++){
//         actualSolution[nnodes + k*nnodes + nodesX[i]] = 0;
//         actualSolution[nnodes + arg*nnodes + nodesX[i]] = 1;
//     }

//     if ( arg != k ){
//         for ( std::vector<int>::iterator it = P.begin(); it != P.end(); ++it){
//             if (*it == k){
//                 std::swap(*it, P.back());
//                 P.pop_back();
//                 break;
//             }
//         }
//         P.emplace_back(arg);
//     }

//     // on actualise les capacités
//     restArr[k] = data->capacity[k];
//     restArr[arg] = data->capacity[arg] - demands;
// }

// void perturbation3(const Data *data, int * actual_solution){

//     int nnodes = data->nnodes;
    
//     std::vector <int> P; // P contient les indices des noeuds qui sont centre
//     bool isUsed[nnodes]; // isCenter[i] = true si i est centre, false sinon
//     double gammaArray[nnodes] = {0}, piArray[nnodes] = {0}, restArray[nnodes] = {0}; 

//     // on reset tout
//     for (size_t i = 0; i < nnodes+nnodes*nnodes; i++){
//         actual_solution[i] = 0;
//     }

//     // on choisit un noeud au hasard et on le met à 1 (il est donc centre)
//     int id = rand() % nnodes;
//     isUsed[id] = actual_solution[id] = 1; // id devient un centre
//     P.push_back(id); // on ajoute id à P
//     actual_solution[nnodes + id*nnodes+id] = 1; // on part du principe que id est autosuffisant

//     // on calcule gammaArray et piArray
//     updateGamma(isUsed, gammaArray, id, data);
//     updatePi(isUsed, piArray, gammaArray, data);

//     for (size_t n = 0; n < data->p; n++){
//         // on choisit un noeud au hasard selon piArray
//         id = chooseCenter(piArray);
//         isUsed[id] = actual_solution[id] = 1;
//         P.push_back(id);
//         actual_solution[nnodes + id*nnodes+id] = 1;

//         // on met à jour gammaArray et piArray
//         gammaArray[id] = piArray[id] = 0;

//         updateGamma(isUsed, gammaArray, id, data);
//         updatePi(isUsed, piArray, gammaArray, data);
//     }
    
//     // initialisation de restArray (chaque noeud prend sa capacité - sa demande car autosuffisant)
//     for (size_t i = 0; i < nnodes; i++){
//         restArray[i] = data->capacity[i] - data->demand[i];
//     }

//     // on va maintenant choisir les noeuds qui vont être assignés à un centre
//     int k;
//     double distMax = getDistMax(data) + 1; // d barre
//     for (size_t i = 0; i < nnodes; i++){
        
//         if (!isUsed[i]){
//             // k* <- arg min_{k from K} phi(j,k)
//             k = getSubset(distMax,i,P,restArray,data);

//             isUsed[i] = actual_solution[nnodes + k * nnodes + i] = 1; // on assigne i à k

//             // on met à jour les centres et les restArray
//             updateCenter(k, P, restArray, actual_solution, data);
//         }
//     }
// }
