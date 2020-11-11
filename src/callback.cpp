/*
 * callback.cpp
 *
 *  Created on: Nov 8, 2020
 *      Author: geunyeongbyeon
 */

#include "callback.h"
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
        int i;
        int n_l = lazyData_.n_l;
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif
        /* Get the current solution */
        /* is it good to define as a member in this class? */
        getValues(lazyData_.barx, xVars_);

        tVal_ = getValue(tVar_);
        lazyData_.current_master_objVal = getObjValue();
        
        for (i = 0; i < n_l; i++) {
            lazyData_.xVals[i] = lazyData_.barx[i];
        }

#ifdef CALLBACK_DEBUG
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

                    if (FUB >= follower_.getObjVal()) 
                    {
                        /* do not add cuts when it happens */   
                    } 
                    else 
                    {
                        addLocal(dy_ - FUB <= 0);
                    }

                    if (cut_type_ == 1) {
                        if (tVal_ + tol_ <= actual_tVal_) {
                            addLocal(tVar_ >= lazyData_.termLf - wVal_ 
                                    * (fobjval_ + (FUB - fobjval_) * lazyData_.indicatorTermx));
                        }
                    } else if (cut_type_ == 3) {
                        if (actual_tVal_ >= tol_) {
                            addLocal(lazyData_.termLf - wVal_ 
                                    * (fobjval_ + (FUB - fobjval_) * lazyData_.indicatorTermx) <= 0);
                        }
                    } else {
                        cout << "Incorrect cut type" << endl;
                    }
                }
            }
        }
#endif
        if (cut_type_ == 1) {
            if (tVal_ - tol_ <= actual_tVal_) {
#ifdef CALLBACK_DEBUG
               cout << "tVal: " << tVal_ << " vs " << actual_tVal_ << endl;
#endif
                add(tVar_ >= lazyData_.termLf - wVal_ * (fobjval_ + (follower_.getbigM() - fobjval_) * lazyData_.indicatorTermx));
            }
        } else if (cut_type_ == 2) {
            if (fcheck_ >= tol_)
                add(lazyData_.termfP <= 0);
        } else if (cut_type_ == 3) {
            if (actual_tVal_ >= tol_)
                add(lazyData_.termLf - wVal_ * (fobjval_ + (follower_.getbigM() - fobjval_) * lazyData_.indicatorTermx) <= 0);
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

                lazyData_.feas_x = lazyData_.barx;
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

    follower_.updateProblem(lazyData_.xVals);

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

        leaderFollower_.updateProblem (lazyData_.xVals, follower_.getObjVal());

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
    lazyData_.termfP.clear();
    lazyData_.termLf.clear();

   auto start_t = chrono::system_clock::now();

    follower_.getBendersTerms(lazyData_.termfP, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termfP: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (cut_type_ != 2)
        leaderFollower_.getBendersTerms(lazyData_.termLf, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termLf: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    lazyData_.indicatorTermx.clear();

    for (int i = 0; i < lazyData_.n_l; i++) {
        if (lazyData_.barx[i] > 0.5) {
            lazyData_.indicatorTermx += 1 - xVars_[i];
        } else {
            lazyData_.indicatorTermx += xVars_[i];
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

IloCplex::Callback incumbentUpdateCallback(IloEnv env, IloNumVarArray& xVars, LazyData& lazyData) {
  return (IloCplex::Callback(new (env) incumbentUpdateCallbackI(env, xVars, lazyData)));
}
void incumbentUpdateCallbackI::main() {
    try {
        if (lazyData_.found_new_incumbent) {
            if (lazyData_.current_best_ub < getIncumbentObjValue()) {
//                cout << "setSolution()" << endl;
                setSolution(xVars_, lazyData_.feas_x, lazyData_.current_best_ub);
            }
            if (lazyData_.current_best_ub >= getIncumbentObjValue()) {
                lazyData_.found_new_incumbent = false;
            }
        }
    }
    catch (...) {
        throw;
    }
}

IloCplex::Callback BendersUserCallback(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy,
                        LazyData &lazyData, Follower &follower, LeaderFollower &leaderFollower){
  return (IloCplex::Callback(new (env) BendersUserCallbackI(env, xVars, yVars, dy, lazyData, follower, leaderFollower)));
}
void BendersUserCallbackI::main(){

    // Skip the separation if not at the end of the cut loop
    if (!isAfterCutLoop())
        return;
    
    int i;

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
#endif

    if (follower_.getStatus() == IloAlgorithm::Status::Optimal) {

        double FUB = follower_.getObjVal();
        // cout << "fUB: " << FUB  << endl;

        if (follower_.getbigM() > FUB) {
            
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
    }
    
    abortCutLoop();

    return;

} 

IloCplex::Callback BendersLazyCallbackMC(IloEnv env, FollowerMC &followerMC, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) 
{
    return (IloCplex::Callback(new (env) BendersLazyCallbackMCI(env, followerMC, leaderFollower, xVars, tVar, dy, lazyData)));
}
void BendersLazyCallbackMCI::main(){
    
    try {
        int i;
        int n_l = lazyData_.n_l;
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif
        /* Get the current solution */
        /* is it good to define as a member in this class? */
        getValues(lazyData_.barx, xVars_);

        tVal_ = getValue(tVar_);
        lazyData_.current_master_objVal = getObjValue();
        
        for (i = 0; i < n_l; i++) {
            lazyData_.xVals[i] = lazyData_.barx[i];
        }

#ifdef CALLBACK_DEBUG
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
        fcheck_ = followerMC_.getCheck();
        fobjval_ = followerMC_.getObjVal();
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
    if (cut_type_ == 1) {
        if (tVal_ - tol_ <= actual_tVal_) {
#ifdef CALLBACK_DEBUG
            cout << "tVal: " << tVal_ << " vs " << actual_tVal_ << endl;
#endif
            add(tVar_ >= lazyData_.termLf - wVal_ * lazyData_.termfPMC);
        }
    } else if (cut_type_ == 2) {
        if (fcheck_ <= -tol_) {
            add(lazyData_.termfPMC >= 0);
        }
    } else if (cut_type_ == 3) {
        if (actual_tVal_ >= tol_)
            add(lazyData_.termLf - wVal_ * lazyData_.termfPMC <= 0);
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

                lazyData_.feas_x = lazyData_.barx;
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

int BendersLazyCallbackMCI::solveSubs() {

   auto start_t = chrono::system_clock::now();

    followerMC_.updateProblem(lazyData_.xVals);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating followerMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    followerMC_.solve();

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving followerMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (followerMC_.getStatus() == IloAlgorithm::Status::Optimal) {

        leaderFollower_.updateProblem (lazyData_.xVals, followerMC_.getObjVal());

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

        if (leaderFollower_.getStatus() == IloAlgorithm::Status::Infeasible) 
        {
            cut_type_ = 3;
        } else
            cut_type_ = 1;
    } else {
        cut_type_ = 2;
    }

    return 0;
}
int BendersLazyCallbackMCI::addBendersCuts()
{
    lazyData_.termfPMC.clear();
    lazyData_.termLf.clear();

   auto start_t = chrono::system_clock::now();

    followerMC_.getBendersTerms(lazyData_.termfPMC, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termfPMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (cut_type_ != 2) {
        leaderFollower_.getBendersTerms(lazyData_.termLf, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termLf: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif
    }

#ifdef C_SOLVE_DEBUG
    if (lazyData_.use_normalization) {
        cout << "check followerMC: " << followerMC_.getCheck() << endl;
        if (cut_type_ != 2)
            cout << "check leaderFollower: " << leaderFollower_.getCheck() / leaderFollower_.getwVal() << endl;
    } else {
        cout << "check followerMC: " << followerMC_.getCheck() << endl;
        if (cut_type_ != 2)
            cout << "check leaderFollower: " << leaderFollower_.getCheck() - leaderFollower_.getwVal() * followerMC_.getObjVal() << endl;
    }
#endif

    return 0;
}


IloCplex::Callback BendersUserCallbackMC(IloEnv env, FollowerMC &followerMC, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) 
{
    return (IloCplex::Callback(new (env) BendersUserCallbackMCI(env, followerMC, leaderFollower, xVars, tVar, dy, lazyData)));
}
void BendersUserCallbackMCI::main(){
    
    try {
        // Skip the separation if not at the end of the cut loop
        if (!isAfterCutLoop())
            return;

        int i;
        int n_l = lazyData_.n_l;
#ifdef CALLBACK_DEBUG
        auto start_t = chrono::system_clock::now();
#endif
        /* Get the current solution */
        /* is it good to define as a member in this class? */
        getValues(lazyData_.barx, xVars_);

        tVal_ = getValue(tVar_);
        lazyData_.current_master_objVal = getObjValue();
        
        for (i = 0; i < n_l; i++) {
            if (lazyData_.barx[i] <= 0.0001) {
                lazyData_.barx[i] = 1e-4;
            } else if (lazyData_.barx[i] >= 0.9999)
                lazyData_.barx[i] = 1-1e-4;
            lazyData_.xVals[i] = lazyData_.barx[i];
            // if (lazyData_.barx[i] > 1e-6) cout << i << ": " << lazyData_.xVals[i] << endl;     
        }

#ifdef CALLBACK_DEBUG
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
        fcheck_ = followerMC_.getCheck();
        fobjval_ = followerMC_.getObjVal();
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
    if (cut_type_ == 1) {
        if (tVal_ + tol_ <= actual_tVal_) {
#ifdef CALLBACK_DEBUG
            cout << "tVal: " << tVal_ << " vs " << actual_tVal_ << endl;
#endif
            add(tVar_ >= lazyData_.termLf - wVal_ * lazyData_.termfPMC);
        }
    } else if (cut_type_ == 2) {
        if (fcheck_ <= -tol_) {
            add(lazyData_.termfPMC >= 0);
        }
    } else if (cut_type_ == 3) {
        if (actual_tVal_ >= tol_)
            add(lazyData_.termLf - wVal_ * lazyData_.termfPMC <= 0);
    } else {
        cout << "Incorrect cut type" << endl;
    }

#ifdef CALLBACK_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in adding Benders: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

//         if (cut_type_ == 1) {

//             IloNum current_obj_val = lazyData_.current_master_objVal - tVal_ + lfobjval_;

//             if (current_obj_val < lazyData_.current_best_ub) {
// //        if (current_obj_val - tolerance_ <= current_best_ub_) {
//                 lazyData_.current_best_ub = current_obj_val;
//                 cout << "*" << current_obj_val << endl;

//                 lazyData_.feas_x = lazyData_.barx;
//                 lazyData_.feas_t = lfobjval_;

//                 lazyData_.found_new_incumbent = true;
//             }
//         }
    abortCutLoop();
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

int BendersUserCallbackMCI::solveSubs() {

   auto start_t = chrono::system_clock::now();

    followerMC_.updateProblem(lazyData_.xVals);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in updating followerMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    followerMC_.solve();

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in solving followerMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (followerMC_.getStatus() == IloAlgorithm::Status::Optimal) {

        leaderFollower_.updateProblem (lazyData_.xVals, followerMC_.getObjVal());

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

        if (leaderFollower_.getStatus() == IloAlgorithm::Status::Infeasible) 
        {
            cut_type_ = 3;
        } else
            cut_type_ = 1;
    } else {
        cut_type_ = 2;
    }

    return 0;
}
int BendersUserCallbackMCI::addBendersCuts()
{
    lazyData_.termfPMC.clear();
    lazyData_.termLf.clear();

   auto start_t = chrono::system_clock::now();

    followerMC_.getBendersTerms(lazyData_.termfPMC, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termfPMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (cut_type_ != 2) {
        leaderFollower_.getBendersTerms(lazyData_.termLf, xVars_, lazyData_.barx);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termLf: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif
    }

#ifdef C_SOLVE_DEBUG
    if (lazyData_.use_normalization) {
        cout << "check followerMC: " << followerMC_.getCheck() << endl;
        if (cut_type_ != 2)
            cout << "check leaderFollower: " << leaderFollower_.getCheck() / leaderFollower_.getwVal() << endl;
    } else {
        cout << "check followerMC: " << followerMC_.getCheck() << endl;
        if (cut_type_ != 2)
            cout << "check leaderFollower: " << leaderFollower_.getCheck() - leaderFollower_.getwVal() * followerMC_.getObjVal() << endl;
    }
#endif

    return 0;
}