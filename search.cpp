#include "search.hpp"

void getK(int * actual, std::vector<int> &K, int nnodes){
    for (size_t i = 0; i < nnodes; i++){
        if (actual[i]){
            K.emplace_back(i);
        }
    }
}

void getQ(int * actual, std::vector<int> &Q, int nnodes, int k){

    for (size_t i = 0; i < nnodes; i++){
        if (actual[nnodes + k * nnodes + i] && i != k){
            Q.emplace_back(i);
        }
    }
}

void updateRho(std::vector<int> Q, std::vector<double> &rho, int k, const Data * data){

    double sumD = 0;
    // calcul de la somme des distances dans le sous ensemble
    for (auto q : Q){
        sumD += data->distances[q * data->nnodes + k];
    }
    for (auto q : Q){
        rho.emplace_back(data->distances[q * data->nnodes + k] / sumD);
    }
}

int selectNextDelete(std::vector<int> Q, std::vector<double> rho){
    double debut = rho.at(0), rng = (double)rand() / (double) RAND_MAX;;
    int arg = 1, rep;

    while (debut < rng){
        debut += rho.at(arg);
        arg++;
    }
    rep = Q.at(arg-1);
    return rep;
}

void algoTri(std::vector<std::pair<int,int>> &W, const Data * data){
    int i, j;
    std::vector<std::pair<int,int>> WBis;
    std::pair<int,int> temp;
    double dist;
    int arg = 0;

    while (!W.empty()){  
        dist = 0;
        arg = 0;

        for (size_t i = 0; i < W.size(); i++){
            if ( dist < data->distances[data->nnodes * W.at(i).first + W.at(i).second]){
                dist = data->distances[data->nnodes * W.at(i).first + W.at(i).second];
                arg = i;
            }
        }
        WBis.emplace_back(W.at(arg));
        W.erase(W.begin() + arg );
    }
    for (auto w : WBis){
        W.emplace_back(w);
    }
}

double getDistMaxSearch(int nnodes, const Data *pcenter_data){
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

double getPhi(double distMax, int i, int k, double TG[], const Data *data){
    if ((double)data->demand[i] <= TG[k]){
        return (double)data->distances[i*data->nnodes+k] / distMax;
    }else{
        return (double)data->demand[i] - TG[k];
    }
}

int getSubSet(double distMax,int w, std::vector<int> K,double rest[], const Data * data){
    int min = K.at(0);
    double phi = __DBL_MAX__, tmp;

    for (auto k : K){
        
        tmp = getPhi(distMax, w, k, rest,data);
        if(tmp < phi){
            phi = tmp;
            min = k;
        }
    }
    return min;
} 

bool find(std::vector<std::pair<int,int>> W, int id){

    for (auto w : W){
        if ( w.first == id)
            return true;
    }
    return false;
}

void updateCenterSearch(std::vector<int> &K, int *actual_solution, const Data * data){

    std::vector<std::pair<int,double>> X; // id , distance max
    double demands = 0;
    int nnodes = data->nnodes, nb = 0;
    std::vector<int> KBis;

    for (auto k : K){
        X.clear();
        demands = nb = 0;
        for (size_t i = 0; i < nnodes; i++){
            if (actual_solution[nnodes+k*nnodes+i]){
                demands += data->demand[i];
                nb++;
                X.emplace_back(std::make_pair(i,0));
            }
        }

        double max = 0;
        for (size_t i = 0; i< X.size();i++){
            max = 0;
            for (size_t j = 0; j< X.size();j++){
                if (max < data->distances[X.at(i).first*nnodes + X.at(j).first])
                    max = data->distances[X.at(i).first*nnodes + X.at(j).first];
            }
            X.at(i).second = max;
        }

        int arg = k;
        double dist = __DBL_MAX__;
        for (auto x : X){
            if ( data->capacity[x.first] - demands > 0 && dist > x.second){
                arg = x.first;
                dist = x.second;
            }
        }
    
        for (auto x : X){
            actual_solution[nnodes + k * nnodes + x.first] = 0;
            actual_solution[nnodes + arg * nnodes + x.first] = 1;
        }

        KBis.emplace_back(arg);
    }

    for (auto k : K){
        actual_solution[k] = 0;
    }
    

    K.clear();
    for (auto k : KBis){
        actual_solution[k] = 1;
        K.emplace_back(k);
    }
    
}

double getW(int * actual_solution, const Data* data, std::vector<int> K){
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

double local_search(int * actual_solution, const Data * data){
    std::vector<std::pair<int,int>> W; // ensemble des sommets supprimés
    std::vector<int> K; // ensemble des sommets centres 
    std::vector<int> Q; // ensemble des sommets Q sur lequel nous travaillons actuellement
    std::vector<double> rho;
    int deleted = 0, id = 0;

    getK(actual_solution, K, data->nnodes);

    for (auto k : K){
                
        getQ(actual_solution, Q, data->nnodes, k);
        
        if(Q.empty()){
            continue;
        }
        updateRho(Q,rho,k,data);
            
        deleted = 0;
        while ( ((double)deleted / (double)Q.size()) < ALPHA){
            id = selectNextDelete(Q,rho);
            if ( find(W,id)) continue;
            W.emplace_back(std::make_pair(id, k));
            deleted++;

        }
        rho.clear();
        Q.clear();
    }
    
    algoTri(W,data);


    for (auto w : W){
        actual_solution[data->nnodes + w.first + w.second * data->nnodes] = 0;
    }

    double distMax = getDistMaxSearch(data->nnodes,data) + 1;
    double rest[data->nnodes];

    for (size_t i = 0; i < data->nnodes; i++){
        rest[i] = data->capacity[i] - data->demand[i];
    }
    

    for (auto w : W){

        id = getSubSet(distMax,w.first, K, rest, data);
        actual_solution[data->nnodes + data->nnodes * id + w.first] = 1;
        rest[id] -= data->demand[w.first];
    }    
    updateCenterSearch(K, actual_solution, data);
    
    double dist = getW(actual_solution, data, K);
    return dist;

}