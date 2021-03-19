#include "lazyCBBenders.h"
// #define C_SOLVE_DEBUG
// #define addLocalCut

IloCplex::Callback BendersLazyCallback(IloEnv env, Follower &follower, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) 
{
    return (IloCplex::Callback(new (env) BendersLazyCallbackI(env, follower, leaderFollower, xVars, tVar, dy, lazyData)));
}
void BendersLazyCallbackI::main(){
    
    try {
        
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif
        /* Get the current solution */
        getValues(barx_, xVars_);
        tVal_ = getValue(tVar_);
        lazyData_.current_master_objVal = getObjValue();
        for (i = 0; i < n_l; i++) {
            xVals_[i] = barx_[i];
        }

#ifdef CALLBACK_DEBUG
   cout << "current_master_objVal: " << lazyData_.current_master_objVal << endl;
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting master Result: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

        /* solve subproblems: follower and leaderFollower */
        solveSubs();

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving subs: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

        /* generate cut expression */
        addBendersCuts();

        /* get necessary information */
        fcheck_ = follower_.getCheck();
        fobjval_ = follower_.getObjVal();
        lfcheck_ = leaderFollower_.getCheck();
        lfobjval_ = leaderFollower_.getObjVal();
        wVal_ = leaderFollower_.getwVal();
        assert(fabs(fcheck_ - fobjval_) <= tol_);
        actual_tVal_ = lfcheck_ - wVal_ * fcheck_;
        assert(fabs(lfobjval_ - actual_tVal_) <= tol_);
        
#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in generating Benders: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

#ifdef addLocalCut
        /* update follower UB */
        if (getSolutionSource() != CPX_CALLBACK_MIP_INCUMBENT_HEURSOLN) {

            getLBs(lazyData_.xLBs_ilo, xVars_);
            getUBs(lazyData_.xUBs_ilo, xVars_);

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting xBounds: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

            for (i = 0; i < n_l; i++) {
                lazyData_.xLBs[i] = lazyData_.xLBs_ilo[i];
                lazyData_.xLBs[i] = lazyData_.xLBs_ilo[i];
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
#endif

            if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {

                double FUB = follower_.getCplexPtr()->getObjValue();

                if (follower_.getbigM() > FUB) {

                    // if (FUB >= follower_.getObjVal()) 
                    // {
                        /* do not add cuts when it happens */   
                    // } 
                    // else 
                    // {
                        addLocal(dy_ - FUB <= 0);
                    // }

                    // if (cut_type_ == 1) {
                    //     if (tVal_ + tol_ <= actual_tVal_) {
                    //         addLocal(tVar_ >= termLf_ - wVal_ 
                    //                 * (fobjval_ + (FUB - fobjval_) * indicatorTermx_));
                    //     }
                    // } else if (cut_type_ == 3) {
                    //     if (actual_tVal_ >= tol_) {
                    //         addLocal(termLf_ - wVal_ 
                    //                 * (fobjval_ + (FUB - fobjval_) * indicatorTermx_) <= 0);
                    //     }
                    // } else {
                    //     cout << "Incorrect cut type" << endl;
                    // }
                }
            }
        }
#endif
        if (cut_type_ == 1) {
            if (tVal_ - tol_ <= actual_tVal_) {
#ifdef CALLBACK_DEBUG
               cout << "tVal: " << tVal_ << " vs " << actual_tVal_ << endl;
               cout << "wVal_: " << wVal_ << endl;
#endif
                add(tVar_ >= termLf_ - wVal_ * (fobjval_ + (follower_.getbigM() - fobjval_) * indicatorTermx_));
            }
        } else if (cut_type_ == 2) {
            if (fcheck_ >= tol_)
                add(termfP_ <= 0);
        } else if (cut_type_ == 3) {
            if (actual_tVal_ >= tol_)
                add(termLf_ - wVal_ * (fobjval_ + (follower_.getbigM() - fobjval_) * indicatorTermx_) <= 0);
        } else {

            cout << "Incorrect cut type" << endl;
        }

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in adding Benders: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

        if (cut_type_ == 1) {

            IloNum current_obj_val = lazyData_.current_master_objVal - tVal_ + lfobjval_;

            if (current_obj_val < lazyData_.current_best_ub) {
//        if (current_obj_val - tolerance_ <= current_best_ub_) {
                lazyData_.current_best_ub = current_obj_val;
                cout << "*" << current_obj_val << endl;

                lazyData_.feas_x = barx_;
                lazyData_.feas_t = lfobjval_;

                lazyData_.found_new_incumbent = true;
            }
        }
    }
    catch (IloException &e) {
        cerr << "CPLEX found the following exception: " << e << endl;
        e.end();
    }
    catch (...) {
        cerr << "The following unknown exception was found: " << endl;
    }

    return;
}

int BendersLazyCallbackI::solveSubs() {

   auto start_t = chrono::system_clock::now();

    follower_.updateProblem(xVals_);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating follower: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    follower_.solve();

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving follower: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {

        leaderFollower_.updateProblem (xVals_, follower_.getObjVal());

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating leaderFollower: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

        if (leaderFollower_.solve())
            return 1;

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving leaderFollower: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

        if (leaderFollower_.getStatus() == IloAlgorithm::Status::Infeasible) {

            cut_type_ = 3;

        } else
            cut_type_ = 1;
    } else {
        cut_type_ = 2;
    }

    return 0;
}
int BendersLazyCallbackI::addBendersCuts(){
    termfP_.clear();
    termLf_.clear();

   auto start_t = chrono::system_clock::now();

    follower_.getBendersTerms(termfP_, xVars_, xVals_);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termfP: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (cut_type_ != 2)
        leaderFollower_.getBendersTerms(termLf_, xVars_, xVals_);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termLf: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    indicatorTermx_.clear();

    for (i = 0; i < n_l; i++) {
        if (barx_[i] > 0.5) {
            indicatorTermx_ += 1 - xVars_[i];
        } else {
            indicatorTermx_ += xVars_[i];
        }
    }

#ifdef C_SOLVE_DEBUG
    if (lazyData_.use_normalization) {
        cout << "check follower: " << follower_.getCheck() << endl;
        cout << "check leaderFollower: " << leaderFollower_.getCheck() / leaderFollower_.getwVal() << endl;
    } else {
        cout << "check follower: " << follower_.getCheck() << endl;
        cout << "check leaderFollower: " << leaderFollower_.getCheck() - leaderFollower_.getwVal() * follower_.getObjVal() << endl;
    }
#endif

    return 0;
}