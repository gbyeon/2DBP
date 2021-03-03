#include "lazyCBfUBnu.h"
// #define CALLBACK_DEBUG

IloCplex::Callback LazyCallbackfUBNu(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy,
                        LazyData &lazyData, Follower &follower, Hpp &hpp, Data &data){
  return (IloCplex::Callback(new (env) LazyCallbackfUBNuI(env, xVars, yVars, dy, lazyData, follower, hpp, data)));
}
void LazyCallbackfUBNuI::main(){

    // // Skip the separation if not at the end of the cut loop
    // if (!isAfterCutLoop())
    //     return;
    
    
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif

    getValues(barx_, xVars_);
    for (i = 0; i < n_l; i++) {
        xVals_[i] = barx_[i];
    }

    double fUB = 0;
    /* construct nu values */
    for (i = 0; i < m_f; i++) {
        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            if (fC_lV_coef_[i][j] > 0){
                nu1_[i][j] = xVals_[fC_lV_ind_[i][j]] * fC_lV_coef_[i][j];
                nu4_[i][j] = 0;
                fUB += nu1_[i][j] * dualUb_;// * xVals_[fC_lV_ind_[i][j]];
            } else {
                nu1_[i][j] = 0;
                nu4_[i][j] = 0;
                // nu4_[i][j] = (1-xVals_[fC_lV_ind_[i][j]]) * fabs(fC_lV_coef_[i][j]);
            }
        }
    }

    getLBs(lazyData_.xLBs_ilo, xVars_);
    getUBs(lazyData_.xUBs_ilo, xVars_);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting xBounds: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    /* set objective function: minimized cx + t */
    int nnodes = getNnodes();
    hpp_.updateObjNu(nu1_, nu4_, lazyData_.xLBs_ilo, lazyData_.xUBs_ilo, nnodes);
    hpp_.solve();

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting xBounds: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    follower_.updateProblem (xVals_);

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
#endif

    if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {
// cout << "nu: " << fUB << endl;
cout << "fobj: " << follower_.getObjVal() << " hpp obj: " << hpp_.getObjVal() << endl;
        fUB += follower_.getObjVal() - hpp_.getObjVal();
        
        cout << "fUB: " << fUB << endl;

        if (follower_.getObjVal() > fUB) {
            
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
                addLocal(dy_ <= fUB);
            // }
        }
    } else {
        cout << "fUB may be infeasible!" << endl;
    }
    // abortCutLoop();

    return;

} 