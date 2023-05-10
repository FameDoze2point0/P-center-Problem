#include "solver.hpp"

#define TIME_LIMIT_IN_SECS 7200
#define NUMBER_OF_THREADS 12
#define NUMBER_MAX_ITERATIONS 15000
#define EPSILON 1e-5


ILOSTLBEGIN

using namespace std;

void Solver::set_parameters(){
    cplex.setParam(IloCplex::Param::ClockType, 1);
    cplex.setParam(IloCplex::Param::Threads, 12);
    //cplex.setParam(IloCplex::Param::Emphasis::MIP,5);
    //cplex.setParam(IloCplex::Param::MIP::Strategy::Dive,2);
    /* cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,1.5);
    cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, 3);
    cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur,2); */
    //cplex.setParam(IloCplex::Param::Advance, 2);
    cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_IN_SECS);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-10);
    //cplex.use(InfoMIPCallBack(env, intermediateLB, intermediateUB, requestIntermediate));
    cplex.setOut(env.getNullStream());
}


void Solver::create_model(){
    
    // variable disant à quel centre sont connectés chaque noeud [0,1]
    x = IloNumVarArray(env, nnodes * nnodes, 0.0, 1.0, ILOFLOAT);
    // variable disant si un noeud est centre ou non (1 est centre, 0 n'est pas centre) [0,1]
    y = IloNumVarArray(env, nnodes, 0.0, 1.0, ILOFLOAT);
    // Distance maximale (sera à minimiser)
    w = IloNumVar(env); // correspond au W (attention != w)
        
    // ======== début fonction objective ========
        IloExpr objFct(env);
        objFct += w;
        obj.setExpr(IloMinimize(env, objFct));
        model.add(obj);   
        objFct.end();
    // ======== fin fonction objective ========

    // ======== début contraintes ========

        // Contrainte 1 : On veut que chaque noeud ne soit relié qu'à un seul centre
        // Somme x_{ij} = 1 tel que i et j appartiennent à V  
            IloExprArray constraints_center_linked(env, nnodes); // on crée un tableau de contrainte, afin d'avoir une contrainte par noeud
            for (size_t j = 0; j < nnodes; j++){

                constraints_center_linked[j] = IloExpr(env); // on crée la contrainte pour un seul noeud
                for (size_t i = 0; i < nnodes; i++)
                    constraints_center_linked[j] += x[i*nnodes + j]; 
                
                model.add(constraints_center_linked[j] == 1);
                constraints_center_linked[j].end();
            }
        // fin contrainte 1


        // Contrainte 2 : On veut que le nombre de centre à la fin soit égal au nombre de centre donné
        // Somme y_i = p tel que i appartient à V 
            IloExpr constraints_nb_center(env);
            for (size_t i = 0; i < nnodes; i++)
                constraints_nb_center += y[i];
            
            model.add(constraints_nb_center == ncenters);
            constraints_nb_center.end();
        // fin contrainte 2


        // Contrainte 3 : On veut que la demande totale pour un centre soit inférieure à la capacité qu'il peut produire
            // Somme w_j * x_{ij} - s_i * y_i <= 0 tel que i et j appartiennent à V 
            double pourcentage = 100 / 100; // peut être modifié
            //int capacity = ceil(nnodes*pourcentage); // correspond au s dans l'inéquation (3) (pour un début)
            //int demands = 1; // correspond au w dans l'inéquation (3) (pour un début)
        
            IloExprArray constraints_demands_capacity(env,nnodes);
            for (size_t i = 0; i < nnodes; i++){

                constraints_demands_capacity[i] = IloExpr(env);
                for (size_t j = 0; j < nnodes; j++)
                    constraints_demands_capacity[i] += (data->demand[j]*x[i*nnodes+j]);
                constraints_demands_capacity[i] -= data->capacity[i]*y[i];
                model.add(constraints_demands_capacity[i] <= 0);
                constraints_demands_capacity[i].end();   
            }
        // fin contrainte 3

        // Contrainte 4 : Nous voulons minimiser la distance entre un noeud et son centre associé (et ça pour tous les noeuds)
        // Somme d_{ij} * x{ij} - W <= 0 tel que i et j appartiennent à V 
            IloExprArray constraints_distance_to_center(env,nnodes);
            for (size_t j = 0; j < nnodes; j++){
                
                constraints_distance_to_center[j] = IloExpr(env);
                for (size_t i = 0; i < nnodes; i++)
                    constraints_distance_to_center[j] += (data->distances[i*nnodes+j] * x[i*nnodes+j]);
                constraints_distance_to_center[j] -= w;
                model.add(constraints_distance_to_center[j] <= 0);
                constraints_distance_to_center[j].end();
            }
        // fin contrainte 4

        // Contrainte 5 : 
            IloExprArray contraintes5(env, nnodes*nnodes);
            for (size_t j = 0; j < nnodes; j++){
                for (size_t i = 0; i < nnodes; i++){
                    contraintes5[i*nnodes+j] = IloExpr(env);
                    contraintes5[i*nnodes+j] = x[i*nnodes+j]-y[i];
                    model.add(contraintes5[i*nnodes+j] <= 0);
                    contraintes5[i*nnodes+j].end();
                }
            }
    // ======== fin contraintes ========

}

bool Solver::isRounded(IloNumArray x_, IloNumArray y_){
    for (size_t i = 0; i < nnodes; i++){
        if (y_[i] > 1e-8 && y_[i] < 1 - 1e-8)
            return false;
    }
    for (size_t i = 0; i < nnodes*nnodes; i++){
        if (x_[i] > 1e-8 && x_[i] < 1 - 1e-8)
            return false;
    }
    return true;
}

void Solver::roundXY(IloNumArray x_, IloNumArray y_){
    for (size_t i = 0; i < nnodes; i++){
        actual_solution[i] = round(y_[i]);
    }

    for (size_t i = 0; i < nnodes*nnodes; i++){
        actual_solution[nnodes+i] = round(x_[i]);
    }  
}

bool Solver::isVerified(IloNum w_){
    double cpt = 0;
    // printYX(y_,x_);
    //vérification contrainte 1
    for (size_t j = 0; j < nnodes; j++){
        cpt = 0;
        
        for (size_t i = 0; i < nnodes; i++)
            cpt += actual_solution[nnodes+i*nnodes + j]; 
        // cout << "cptX != 1 : " << (cpt != 1) << "    cpt = " << cpt << endl; 
        if (cpt != 1){
            // save << "La contrainte 1 n'a pas été respectée : un sommet est relié à plusieurs centre\n" << endl;
            return false;
        }    
    }
        
    // vérification contrainte 2
    cpt = 0;
    for (size_t i = 0; i < nnodes; i++)
        cpt += actual_solution[i];
    // cout << "cptY != "<< ncenters << " : " << (cpt != 1) << "    cpt = " << cpt << endl;
    if (cpt != (double)ncenters){
        // save << "La contrainte 2 n'a pas été respectée : le nombre de centre n'est pas égal à "<< ncenters << "\n" << endl;
        return false;
    }  

     // Contrainte 3 : On veut que la demande totale pour un centre soit inférieure à la capacité qu'il peut produire
    // Somme w_j * x_{ij} - s_i * y_i <= 0 tel que i et j appartiennent à V 
    for (size_t i = 0; i < nnodes; i++){

        cpt = 0;
        for (size_t j = 0; j < nnodes; j++)
            cpt += (data->demand[j]*actual_solution[nnodes+i*nnodes+j]);
        cpt -= data->capacity[i]*actual_solution[i];
        if ( cpt > 0)
            return false;
    }
    // fin contrainte 3

    // Contrainte 4 : Nous voulons minimiser la distance entre un noeud et son centre associé (et ça pour tous les noeuds)
    // Somme d_{ij} * x{ij} - W <= 0 tel que i et j appartiennent à V 
    for (size_t j = 0; j < nnodes; j++){
        cpt = 0;
        for (size_t i = 0; i < nnodes; i++)
            cpt += (data->distances[i*nnodes+j] * actual_solution[nnodes+i*nnodes+j]);
        cpt -= w_;
        if (cpt > 0)
            return false;
    }

    // save << "Toutes les contraintes ont été respectées !\n" << endl;
    return true;
}

bool Solver::isInHistory(){

    bool flag = false;
    int cpt = 0;
    for (size_t i = 0; i < 3 && !flag; i++){
        
        cpt = 0;
        // vérification y
        for (size_t j = 0; j < nnodes; j++){
            if (actual_solution[j] == history[i][j]) // si ce sont les mêmes on ajoute 1 au compteur
                cpt++;
        }
        for (size_t j = 0; j < nnodes; j++){

            for (size_t k = 0; k < nnodes; k++){

                if (actual_solution[nnodes+j*nnodes+k] == history[i][nnodes+j*nnodes+k])
                    cpt++;
            }
        }
        if (cpt == (nnodes+nnodes*nnodes)){
            flag = true;
            // save << "Le point entier actuel est déjà présent dans l'historique : historique[" << i << "]\n" << endl;
        }
    }
    return flag;
}

void Solver::addHistory(int index){
    for (size_t i = 0; i < nnodes; i++){
        history[index][i] = actual_solution[i];
    }
    
    for (size_t i = 0; i < nnodes; i++){
        
        for (size_t j = 0; j < nnodes; j++){
            history[index][nnodes+nnodes*i+j] = actual_solution[nnodes+nnodes*i + j];
        }
    }
}

void Solver::updateObjective(){

    IloExpr objFct(env);

    for (size_t i = 0; i < nnodes; i++){
        
        if (actual_solution[i])
            objFct += (1-y[i]);
        else
            objFct += y[i];
    }

    for (size_t i = 0; i < nnodes; i++){
        
        for (size_t j = 0; j < nnodes; j++)
        {
            if (actual_solution[nnodes+i*nnodes+j])
                objFct += 1-x[i*nnodes+j];
            else
                objFct += x[i*nnodes+j];
        }
    }
    IloObjective fo = IloMinimize(env, objFct);
    obj.setExpr(fo);
    fo.end();
    objFct.end();
}

void Solver::freeMemory(){
    for (size_t i = 0; i < HISTORIC_SIZE; i++){
        free(history[i]);
    }
    free(actual_solution);
}

void Solver::solve(){

    double timeSpend = 0, t_initial = cplex.getCplexTime();
    // clock_t begin, end;
    auto start = std::chrono::high_resolution_clock::now();

    cout << "Perturbation Type : " << data->perturbation << endl;
    set_parameters();
    create_model();
    cplex.extract(model);
    
    // begin = clock();
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible){
        cout << "cplex.getStatus() == IloAlgorithm::Infeasible" << endl;
    }else if (cplex.getStatus() == IloAlgorithm::Unbounded){
        cout << "cplex.getStatus() == IloAlgorithm::Unbounded" << endl;
    }else if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal){

        bool isFound = false;
        int round = 0, perturbation = 0, k = 10;
        IloNumArray x_(env);
        IloNumArray y_(env);
        IloNum w_;

        cplex.getValues(x,x_);
        cplex.getValues(y,y_);
        w_ = cplex.getValue(w);
        srand(time(NULL));

        if (!isRounded(x_,y_)){
        
            while (round < NUMBER_MAX_ITERATIONS && cplex.getCplexTime() - t_initial < TIME_LIMIT_IN_SECS){

                roundXY(x_,y_);

                if (isVerified(w_)) { 
                    isFound = true;
                    break;
                }

                if (isInHistory()){
                    choosePerturbation(data->perturbation, y_, data, actual_solution, k);
                    //perturbation
                    perturbation++;
                }else{
                    k = 1;
                }
            
                addHistory(round%HISTORIC_SIZE);

                updateObjective();

                cplex.solve();

                cplex.getValues(x,x_);
                cplex.getValues(y,y_);
                w_ = cplex.getValue(w);

                round++;
            }
        }else {
            isFound = true;
        }
        
        // end = clock();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        if (isFound){
            cout << "Solution found in " << round << " rounds\n";
            cout << "Nombre de perturbations : " << perturbation << "\n";
            // cout << "Time : " << (double)(end-begin)/(double)CLOCKS_PER_SEC/NUMBER_OF_THREADS << "\n";
            cout << "Time : " <<  (double)duration.count()/1000.0 << " seconds\n";
            cout << "Objective value : " << cplex.getValue(w) << "\n\n\n" << endl;
        }else{
            cout << "\n\nSolution not found : " << round << " rounds\n";
            // cout << "Time : " << (double)(end-begin)/(double)CLOCKS_PER_SEC/NUMBER_OF_THREADS << "\n";
            cout << "Nombre de perturbations : " << perturbation << "\n\n\n" << endl;
        }
    }
    
    freeMemory();
}
