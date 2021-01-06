#include "usercutCBBendersMC.h"
// #define CALLBACK_DEBUG

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
        
        /* in out */
        double alpha = 0.5;
        double epsilon = 1e-4;
        for (i = 0; i < n_l; i++) {
            // lazyData_.xVals[i] = lazyData_.barx[i];
            lazyData_.xVals[i] = alpha * lazyData_.feas_x[i] + (1-alpha) * lazyData_.barx[i];
        }
        /* perturbation */
        for (i = 0; i < n_l; i++) {
            if (lazyData_.xVals[i] <= epsilon) {
                lazyData_.xVals[i] = epsilon;
            } else if (lazyData_.xVals[i] >= 1-epsilon)
                lazyData_.xVals[i] = 1-epsilon;
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

    followerMC_.getBendersTerms(lazyData_.termfPMC, xVars_, lazyData_.xVals);

#ifdef C_SOLVE_DEBUG
   ticToc_ = (chrono::system_clock::now() - start_t);
   cout << "time in getting termfPMC: " << ticToc_.count() << endl;
   start_t = chrono::system_clock::now();
#endif

    if (cut_type_ != 2) {
        leaderFollower_.getBendersTerms(lazyData_.termLf, xVars_, lazyData_.xVals);

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