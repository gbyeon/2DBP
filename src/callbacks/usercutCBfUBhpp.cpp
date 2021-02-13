#include "usercutCBfUBhpp.h"
// #define CALLBACK_DEBUG

IloCplex::Callback BendersUserCallbackHpp(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, Hpp &hpp){
  return (IloCplex::Callback(new (env) BendersUserCallbackHppI(env, xVars, yVars, dy, lazyData, hpp)));
}
void BendersUserCallbackHppI::main(){

    // Skip the separation if not at the end of the cut loop
    if (!isAfterCutLoop())
        return;
    
    cout << getNodeId() << endl;
    // if (getCurrentNodeDepth() < 30)
        // return;
    
    int i;
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif
    getLBs(lazyData_.xLBs_ilo, xVars_);
    getUBs(lazyData_.xUBs_ilo, xVars_);
    
    getValues(lazyData_.bary, yVars_);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting xBounds: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    hpp_.updateUBProblem (lazyData_.xUBs_ilo, lazyData_.xLBs_ilo);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    hpp_.solve(); 

    double FUB = INFINITY;
    if (hpp_.getStatus() == IloAlgorithm::Status::Optimal || hpp_.getStatus() == IloAlgorithm::Status::Feasible) {
        FUB = hpp_.getObjVal();
    }

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    hpp_.getCplex().getObjective().setSense(IloObjective::Sense::Minimize);

    hpp_.solve(); 

    double FLB = -INFINITY;
    if (hpp_.getStatus() == IloAlgorithm::Status::Optimal || hpp_.getStatus() == IloAlgorithm::Status::Feasible) {
        FLB = hpp_.getObjVal();
    }

    double dyVal = hpp_.getDyVal(lazyData_.bary);
    
        // if (dyVal > FUB) {
            
            // getValues(lazyData_.bary, yVars_);
            // double Dy = follower_.getDy(lazyData_.bary);
            
            // if (FUB >= Dy)
            // if (FUB >= follower_.getObjVal()) 
            // {   
                // cout << "fUB: " << FUB <<", Dy: " << Dy << endl;
                /* do not add cuts when it happens */   
            // } 
            // else 
            // {
                // cout << "add local cut" << endl;
                if (FUB < INFINITY)
                    addLocal(dy_ <= FUB);
                if (FLB > -INFINITY)
                    addLocal(dy_ >= FLB);
                cout << FLB << ", " << dyVal << ", " << FUB << endl;
            // }
        // }
    // } else {
        // cout << "fUBhpp may be infeasible!" << endl;
    // }
    
    // abortCutLoop();

    return;

} 