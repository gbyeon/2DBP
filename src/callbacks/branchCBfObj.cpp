#include "branchCBfObj.h"

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
// cout << "ystar: " << fObjVal_ystar << ", ytilde: " << fObjVal_ytilde << endl;


                if (fObjVal_ystar + 1e-6 <= fObjVal_ytilde) 
                {
                    IloConstraint con_down = dy_ <= floor((fObjVal_ystar + fObjVal_ytilde) / 2); /* (fObjVal_ystar + fObjVal_ytilde) / 2 vs fObjVal_ystar */
                    IloConstraint con_up = dy_ >= floor((fObjVal_ystar + fObjVal_ytilde) / 2); 

                    NodeData *node_down = new NodeData(-1);
                    NodeData *node_up = new NodeData(*node_id);      
                    // cout << "node_id at branchCB: " << *node_id << endl;          
                    makeBranch(con_down, current_master_objVal, node_down);
                    makeBranch(con_up, current_master_objVal, node_up);
                } else 
		{
		//	/* bilevel feasible */
	//		cout << "found a bilevel feasible node. Prune the node..." << endl;
	//	if (isIntegerFeasible())	prune();
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

