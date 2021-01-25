#include "branchCB.h"

IloCplex::Callback branchCallback(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloNumVar &tVar,
                                  IloExpr &dy, LazyData &lazyData, double *blb, double *bub, bool *isended, Follower &follower, LeaderFollower &leaderFollower) {
  return (IloCplex::Callback(new (env) branchCallbackI(env, xVars, yVars, tVar, dy, lazyData, blb, bub, isended, follower, leaderFollower)));
}
void branchCallbackI::main() {
    try {
        
        // if (isIntegerFeasible())
            // return;

        if (!*isEnded) /* if branching on d^T y is not ended, continue */
        {
            cout << "blb: " << *blb << ", bub: " << *bub << endl;
            
            IloConstraintArray cons_down = IloConstraintArray(getEnv(), 1);
            IloConstraintArray cons_up = IloConstraintArray(getEnv(), 1);

            // cons_down[0] = dy_ <= follower_.getLB() + fUB_range / pow(2,k) * (2 * l + 1);
            // cons_up[0] = dy_ >= follower_.getLB() + fUB_range / pow(2,k) * (2 * l + 1);
            cons_down[0] = dy_ <= (*bub - *blb) / 2;
            cons_up[0] = dy_ >= (*bub - *blb) / 2;

            cout << "halve the range using : " << (*bub - *blb) / 2 << endl;
            
            NodeData *node_down = new NodeData(*blb, *blb + (*bub - *blb) / 2);
            NodeData *node_up = new NodeData(*blb + (*bub - *blb) / 2, *bub);                
            makeBranch(cons_down, lazyData_.current_master_objVal, node_down);
            makeBranch(cons_up, lazyData_.current_master_objVal, node_up);
                    
            cons_down.end();
            cons_up.end();

        }
    } catch (IloException &e) {
        cerr << "CPLEX found the following exception: " << e << " in branchCB.cpp" << endl;
        e.end();
    }
    catch (...) {
        cerr << "The following unknown exception was found in branchCB.cpp: " << endl;
    }

    return;
}

