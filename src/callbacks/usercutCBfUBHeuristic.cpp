#include "usercutCBfUBHeuristic.h"
// #define CALLBACK_DEBUG

IloCplex::Callback UserCallbackfUBHeuristic(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy,
                        LazyData &lazyData, Follower &follower, Hpp &hpp, Data &data){
  return (IloCplex::Callback(new (env) UserCallbackfUBHeuristicI(env, xVars, yVars, dy, lazyData, follower, hpp, data)));
}
void UserCallbackfUBHeuristicI::main(){

    // Skip the separation if not at the end of the cut loop
    if (!isAfterCutLoop())
        return;
    
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif

    getLBs(lazyData_.xLBs_ilo, xVars_);
    getUBs(lazyData_.xUBs_ilo, xVars_);

    /* set objective function: minimized ax */
    // hpp_.updateObjax(lazyData_.xLBs_ilo, lazyData_.xUBs_ilo, barx_);
    getValues(barx_, xVars_);
    previous_fobj_ = 1e+100;
    iter_ = 0;
    while (1) {
        
        for (i = 0; i < n_l; i++) {
            xVals_[i] = barx_[i];
        }

        follower_.updateProblem(xVals_);
        follower_.solve(); 

        current_fobj_ = follower_.getObjVal();
        // cout << iter_ << ": " << current_fobj_ << endl;

        // if (fabs(current_fobj_ - previous_fobj_) < 0.01)
        //     break;

        hpp_.updateObjaxpsi(lazyData_.xLBs_ilo, lazyData_.xUBs_ilo, barx_, follower_.getpsiValsPtr());

        bool is_updated = false;
        for (i = 0; i < n_l; i++) {
            if (fabs(xVals_[i] - barx_[i]) > 1e-5){
                is_updated = true;
                break;
            }
        }
        if (!is_updated) break;

        previous_fobj_ = current_fobj_;
        iter_++;
    }

    // getValues(barx_, xVars_);
    // /* round franctional sol */
    // for (i = 0; i < n_l; i++) {
    //     xVals_[i] = round(barx_[i]);
    //     // cout << xVals_[i] << endl;
    // }

//     double fUB = 0;
//     /* construct nu values */
//     for (i = 0; i < m_f; i++) {
//         for (j = 0; j < fC_lV_cnt_[i]; j++) {
//             if (fC_lV_coef_[i][j] > 0){
//                 nu1_[i][j] = xVals_[fC_lV_ind_[i][j]] * fC_lV_coef_[i][j];
//                 nu4_[i][j] = 0;
//                 fUB += nu1_[i][j] * dualUb_;
//             } else {
//                 nu1_[i][j] = 0;
//                 nu4_[i][j] = (1-xVals_[fC_lV_ind_[i][j]]) * fabs(fC_lV_coef_[i][j]);
//             }
//         }
//     }

//     int nnodes = getNnodes();
//     hpp_.updateObjNu(nu1_, nu4_, lazyData_.xLBs_ilo, lazyData_.xUBs_ilo, nnodes);

// #ifdef CALLBACK_DEBUG
//    ticToc_ = (chrono::system_clock::now() - start_t);
//    cout << "time in getting xBounds: " << ticToc_.count() << endl;
//    start_t = chrono::system_clock::now();
// #endif

//     follower_.updateProblem (xVals_);

// #ifdef CALLBACK_DEBUG
//    ticToc_ = (chrono::system_clock::now() - start_t);
//    cout << "time in updating fUB: " << ticToc_.count() << endl;
//    start_t = chrono::system_clock::now();
// #endif
//     follower_.solve(); 

// #ifdef CALLBACK_DEBUG
//    ticToc_ = (chrono::system_clock::now() - start_t);
//    cout << "time in solving fUB: " << ticToc_.count() << endl;
//    start_t = chrono::system_clock::now();
// #endif

//     if (hpp_.solve() == 0){
//         if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {

            getValues(lazyData_.bary, yVars_);
            dyhat_ = follower_.getDyVal(lazyData_.bary);

//             fUB += follower_.getObjVal() - hpp_.getObjVal();
//             // cout << "dy: " << Dy << ", fobj: " << follower_.getObjVal() << " hpp obj: " << hpp_.getObjVal() << ", fUB: " << fUB << endl;
            
            
            // if (Dy > fUB){
                if (dyhat_ > current_fobj_){
                    num_local_cuts_added_++;
                    cout << "fub approx: " << current_fobj_ << ", dyhat: " << dyhat_ << endl;
                addLocal(dy_ <= current_fobj_);
//                 cout << "usercut: fUB: " << fUB <<", Dy: " << Dy << endl;
                // addLocal(dy_ <= fUB);
            }
//         }
//     }

    abortCutLoop();

    return;

} 