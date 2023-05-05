#include "solver.hpp"

#define TIME_LIMIT_IN_SECS 7200
#define NUMBER_OF_THREADS 12
#define NUMBER_MAX_ITERATIONS 15000
#define EPSILON 1e-5

ILOSTLBEGIN

using namespace std;

void Solver::set_parameters(){
    cplex.setParam(IloCplex::Param::ClockType, 1);
    cplex.setParam(IloCplex::Param::Threads, NUMBER_OF_THREADS);
    //cplex.setParam(IloCplex::Param::Emphasis::MIP,5);
    //cplex.setParam(IloCplex::Param::MIP::Strategy::Dive,2);
    /* cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,1.5);
    cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, 3);
    cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur,2); */
    //cplex.setParam(IloCplex::Param::Advance, 2);
    cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_IN_SECS);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-10);
    //cplex.use(InfoMIPCallBack(env, intermediateLB, intermediateUB, requestIntermediate));
    //cplex.setOut(env.getNullStream());
}


void Solver::create_model(){

    x = IloNumVarArray(env, nnodes * nnodes, 0.0, 1.0, ILOFLOAT);
    y = IloNumVarArray(env, nnodes, 0.0, 1.0, ILOFLOAT);
    w = IloNumVar(env);

    // fonction objective
        IloExpr objFct(env);
        objFct += w;
        obj.setExpr(IloMinimize(env, objFct));
        model.add(obj);
        objFct.end();

    // contrainte 1
        IloExprArray constraints_center_linked(env,nnodes);
        for (size_t j = 0; j < nnodes; j++){
            
            constraints_center_linked[j] = IloExpr(env);
            for (size_t i = 0; i < nnodes; i++){
                constraints_center_linked[j] += x[i*nnodes + j];
            }
            model.add(constraints_center_linked[j] == 1);
            constraints_center_linked[j].end();
        }
    
    // contrainte 2
        IloExpr constraint_nb_center(env);
        for (size_t i = 0; i < nnodes; i++){
            constraint_nb_center += y[i];
        }
        model.add(constraint_nb_center == ncenters);
        constraint_nb_center.end();

    // contrainte 3
        IloExprArray constraints_demands_capacity(env,nnodes);
        for (size_t i = 0; i < nnodes; i++){

            constraints_demands_capacity[i] = IloExpr(env);
            for (size_t j = 0; j < nnodes; j++){
                constraints_demands_capacity[i] += x[i*nnodes + j] * data->demand[j];
            }
            model.add(constraints_demands_capacity[i] <= data->capacity[i] * y[i]);
            constraints_demands_capacity[i].end();
        }

    // contrainte 4
        IloExprArray constraints_distance_to_center(env,nnodes);
        for (size_t j = 0; j < nnodes; j++){

            constraints_distance_to_center[j] = IloExpr(env);
            for (size_t i = 0; i < nnodes; i++){
                constraints_distance_to_center[j] += x[i*nnodes + j] * data->distances[i*nnodes+j];
            }
            constraints_distance_to_center[j] -= w;
            model.add(constraints_distance_to_center[j] <= 0);
            constraints_distance_to_center[j].end();
        }

    // contrainte 5
        IloExprArray contrainte5(env,nnodes*nnodes);
        for (size_t i = 0; i < nnodes; i++){
            for (size_t j = 0; j < nnodes; j++){
                contrainte5[i*nnodes + j] = IloExpr(env);
                contrainte5[i*nnodes + j] = x[i*nnodes + j] - y[j];
                model.add(contrainte5[i*nnodes + j] <= 0);
                contrainte5[i*nnodes + j].end();
            }
        }
}

void Solver::freeMemory(){
    for (size_t i = 0; i < HISTORIC_SIZE; i++){
        free(history[i]);
    }
    free(actual_solution);
}

void Solver::solve(){
    set_parameters();
    create_model();
    cplex.extract(model);
    
    cplex.solve();

    freeMemory();
}
