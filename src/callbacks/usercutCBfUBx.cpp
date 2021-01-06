#include "usercutCBfUBx.h"
// #define CALLBACK_DEBUG

IloCplex::Callback BendersUserCallbackX(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy,
                        LazyData &lazyData, FollowerX &followerx){
  return (IloCplex::Callback(new (env) BendersUserCallbackXI(env, xVars, yVars, dy, lazyData, followerx)));
}
void BendersUserCallbackXI::main(){

    // Skip the separation if not at the end of the cut loop
    if (!isAfterCutLoop())
        return;
    
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

    followerx_.updateUBProblem (lazyData_.xUBs_ilo, lazyData_.xLBs_ilo);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif
    followerx_.solve(); 

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving fUB: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (followerx_.getStatus() == IloAlgorithm::Status::Optimal) {

        double FUB = followerx_.getObjVal();
        double dyVal = followerx_.getDyVal(lazyData_.bary);
        cout << "fUB: " << FUB  << ", dyVal: " << dyVal << endl;

        if (dyVal > FUB) {
            
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
    } else {
        cout << "fUBx may be infeasible!" << endl;
    }
    
    abortCutLoop();

    return;

} 