#include "branchCBfObj.h"
// #define DEBUG_branchCB

IloCplex::Callback branchCallbackfObj(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Follower &follower, bool *nodeType, IloInt64 *nodeId) 
{
  return (IloCplex::Callback(new (env) branchCallbackfObjI(env, xVars, yVars, dy, lazyData, follower, nodeType, nodeId)));
}
void branchCallbackfObjI::main() {
    try {
        
        // if (isIntegerFeasible())
            // return;

        if (*node_type) /* if branching on d^T y */
        {
            /* Get the current solution */
            current_master_objVal = getObjValue();
            getValues(barx, xVars_);
            getValues(bary, yVars_);
            for (i = 0; i < n_l; i++) {
                xVals[i] = barx[i];
            }
            
            follower_.updateProblem(xVals);
            follower_.solve();

            if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {
                fObjVal_ystar = follower_.getObjVal();
                fObjVal_ytilde = follower_.getDyVal(bary);

                if (fObjVal_ystar + 1e-8 <= fObjVal_ytilde) 
                {
                    double branching_dy_val = floor(fObjVal_ystar * 1 / 4 + fObjVal_ytilde * 3 / 4); //floor(fObjVal_ystar + fObjVal_ytilde) / 2; fObjVal_ystar
                    #ifdef DEBUG_branchCB
                    cout << "ystar: " << fObjVal_ystar << ", ytilde: " << fObjVal_ytilde << ", branching: " << branching_dy_val << endl;
                    #endif

                    IloConstraint con_down = dy_ <= branching_dy_val;
                    IloConstraint con_up = dy_ >= branching_dy_val;

                    UserNodeData *node_down = new UserNodeData(-100000, branching_dy_val, -1*(*node_id));
                    UserNodeData *node_up = new UserNodeData(branching_dy_val, 100000, *node_id);      
                    #ifdef DEBUG_branchCB
                    cout << "node_id at branchCB: " << *node_id << endl;          
                    #endif
                    makeBranch(con_down, current_master_objVal, node_down);
                    makeBranch(con_up, current_master_objVal, node_up);
                } 
            }
        }
    } catch (IloException &e) {
        cerr << "CPLEX found the following exception: " << e << " in branchCBfObj.cpp" << endl;
        e.end();
    }
    catch (...) {
        cerr << "The following unknown exception was found in branchCBfObj.cpp: " << endl;
    }

    return;
}

