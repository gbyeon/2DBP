#include "usercutCBfUB.h"
// #define CALLBACK_DEBUG
IloCplex::Callback UserCallbackfUB(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy,
                        LazyData &lazyData, Follower &follower){
  return (IloCplex::Callback(new (env) UserCallbackfUBI(env, xVars, yVars, dy, lazyData, follower)));
}
void UserCallbackfUBI::main(){

    // Skip the separation if not at the end of the cut loop
    if (!isAfterCutLoop())
        return;
    
    // if (getCurrentNodeDepth() < 10)
    //     return;
    /* or in CPLEX 12.7, we can use getNnodes < 300 */

    int i;

#ifdef CALLBACK_DEBUG
   auto start_t = chrono::system_clock::now();
#endif

    getLBs(lazyData_.xLBs_ilo, xVars_);
    getUBs(lazyData_.xUBs_ilo, xVars_);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting xBounds: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    for (i = 0; i < lazyData_.n_l; i++) {
        lazyData_.xLBs[i] = lazyData_.xLBs_ilo[i];
        lazyData_.xUBs[i] = lazyData_.xUBs_ilo[i];
    }

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in constructing xBoundsArrays: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    follower_.updateUBProblem (lazyData_.xUBs, lazyData_.xLBs);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    follower_.solve(); 

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
   cout << "status: " << follower_.getStatus() << endl;
#endif

    if (follower_.getStatus() == IloAlgorithm::Status::Optimal || follower_.getStatus() == IloAlgorithm::Status::Feasible) {

        double FUB = follower_.getObjVal();
        double dyVal = follower_.getDyVal(lazyData_.bary);
        cout << "fUB: " << FUB  << ", dyVal: " << dyVal << endl;

        if (dyVal > FUB) {
        // if (follower_.getbigM() > FUB) {
            
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
                addLocal(dy_ <= FUB);
            // }
        }
    } /* otherwise, dy <= +infinity */ 
    
    abortCutLoop();

    return;

} 