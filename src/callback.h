/*
 * callback.h
 *
 *  Created on: Nov 8, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_CALLBACK_H
#define BILEVEL_CALLBACK_H

#include <ilcplex/ilocplex.h>
#include "follower.h"
#include "followerMC.h"
#include "leaderfollower.h"

ILOSTLBEGIN

/* data or data structure or setting used in lazy callback */
struct LazyData {
    
    IloExpr termfP;
    IloExpr termLf;
    IloExpr termfPMC;
    IloExpr indicatorTermx;
    
    int n_l;
    int n_f;
  
    IloNumArray barx;
    IloNumArray bary;

    IloNumArray feas_x;
    double feas_t;
    bool found_new_incumbent;

    double current_master_objVal;
    double current_best_ub;
    double * xVals;
    double * xLBs;
    double * xUBs;
    IloNumArray xLBs_ilo, xUBs_ilo;

    // double tolerance;
    // bool outputConsole;
    bool use_normalization;

    // int cut_type;
    LazyData() {
        n_l = -1;
    }

    LazyData (int n_l_, int n_f_, IloEnv env){
        n_l = n_l_;
        n_f = n_f_;

        barx = IloNumArray(env, n_l);
        bary = IloNumArray(env, n_f);
        feas_x = IloNumArray(env, n_l);
            
        termfP = IloExpr (env);
        termfPMC = IloExpr (env);
        termLf = IloExpr (env);
        indicatorTermx = IloExpr (env);

        xVals = new double [n_l];
        xLBs_ilo = IloNumArray (env, n_l);
        xUBs_ilo = IloNumArray (env, n_l);
        xLBs = new double [n_l];
        xUBs = new double [n_l]; 

        // tolerance = 1e-4;
        // outputConsole= true;
        use_normalization = false;   
    }

    ~LazyData() {
        // if (n_l > 0) {
        //     barx.end();
        //     feas_x.end();

        //     termfP.end();
        //     termfPMC.end();
        //     termLf.end();
        //     indicatorTermx.end();

        //     xLBs_ilo.end();
        //     xUBs_ilo.end();

        //     delete[] xVals;
        //     delete[] xLBs;
        //     delete[] xUBs;
        //     xVals = nullptr;
        //     xLBs = nullptr;
        //     xUBs = nullptr;
        // }
        // if (n_f > 0) {
        //     bary.end();
        // }
    }
};

class BendersLazyCallbackI : public IloCplex::LazyConstraintCallbackI {
  
    /* ref to subproblems */
    Follower &follower_;
    LeaderFollower &leaderFollower_;

    /* ref to master variables */
    IloNumVarArray &xVars_; 
    IloNumVar &tVar_;

    /* ref to data */
    LazyData &lazyData_;

    int cut_type_;

    IloExpr &dy_;

    chrono::duration<double> ticToc_;

    double fcheck_, lfcheck_;
    double fobjval_, lfobjval_;
    double wVal_;
    double tol_;
    double tVal_;
    double actual_tVal_;

public:
    ILOCOMMONCALLBACKSTUFF(BendersLazyCallback);

    /* constructor */
    BendersLazyCallbackI (IloEnv env, Follower &follower, LeaderFollower &leaderFollower, 
                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) : 
                        IloCplex::LazyConstraintCallbackI(env), follower_(follower), leaderFollower_(leaderFollower), xVars_(xVars), tVar_(tVar), dy_(dy), lazyData_(lazyData) {
                            tol_ = 1e-4;
                        };
    void main();
    int addBendersCuts();
    int solveSubs();
};
IloCplex::Callback BendersLazyCallback(IloEnv env, Follower &follower, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData);

class BendersLazyCallbackMCI : public IloCplex::LazyConstraintCallbackI {
  
    /* ref to subproblems */
    FollowerMC &followerMC_;
    LeaderFollower &leaderFollower_;

    /* ref to master variables */
    IloNumVarArray &xVars_; 
    IloNumVar &tVar_;

    /* ref to data */
    LazyData &lazyData_;

    int cut_type_;

    IloExpr &dy_;

    chrono::duration<double> ticToc_;

    double fcheck_, lfcheck_;
    double fobjval_, lfobjval_;
    double wVal_;
    double tol_;
    double tVal_;
    double actual_tVal_;

public:
    ILOCOMMONCALLBACKSTUFF(BendersLazyCallbackMC);

    /* constructor */
    BendersLazyCallbackMCI (IloEnv env, FollowerMC &followerMC, LeaderFollower &leaderFollower, 
                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) : 
                        IloCplex::LazyConstraintCallbackI(env), followerMC_(followerMC), leaderFollower_(leaderFollower), xVars_(xVars), tVar_(tVar), dy_(dy), lazyData_(lazyData) {
                            tol_ = 1e-4;
                        };
    void main();
    int addBendersCuts();
    int solveSubs();
};
IloCplex::Callback BendersLazyCallbackMC(IloEnv env, FollowerMC &followerMC, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData);

class incumbentUpdateCallbackI : public IloCplex::HeuristicCallbackI {
  IloNumVarArray &xVars_;
  LazyData &lazyData_;
public:
  ILOCOMMONCALLBACKSTUFF(incumbentUpdateCallback);
  incumbentUpdateCallbackI(IloEnv env, IloNumVarArray& xVars, LazyData& lazyData) :
                            IloCplex::HeuristicCallbackI(env), 
                            xVars_(xVars), lazyData_(lazyData) {};
  void main();
};
IloCplex::Callback incumbentUpdateCallback(IloEnv env, IloNumVarArray& xVars, LazyData& lazyData);


class BendersUserCallbackI : public IloCplex::UserCutCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    Follower &follower_;
    LeaderFollower &leaderFollower_;
    IloExpr &dy_;
    
public:
  ILOCOMMONCALLBACKSTUFF(BendersUserCallback) \
  BendersUserCallbackI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Follower &follower, LeaderFollower &leaderFollower)
    : IloCplex::UserCutCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), follower_(follower), leaderFollower_(leaderFollower) {};
  void main();
};
IloCplex::Callback BendersUserCallback(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, Follower &follower, LeaderFollower &leaderFollower);


class BendersUserCallbackMCI : public IloCplex::UserCutCallbackI {
  
    /* ref to subproblems */
    FollowerMC &followerMC_;
    LeaderFollower &leaderFollower_;

    /* ref to master variables */
    IloNumVarArray &xVars_; 
    IloNumVar &tVar_;

    /* ref to data */
    LazyData &lazyData_;

    int cut_type_;

    IloExpr &dy_;

    chrono::duration<double> ticToc_;

    double fcheck_, lfcheck_;
    double fobjval_, lfobjval_;
    double wVal_;
    double tol_;
    double tVal_;
    double actual_tVal_;

public:
    ILOCOMMONCALLBACKSTUFF(BendersUserCallbackMC);

    /* constructor */
    BendersUserCallbackMCI (IloEnv env, FollowerMC &followerMC, LeaderFollower &leaderFollower, 
                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData) : 
                        IloCplex::UserCutCallbackI(env), followerMC_(followerMC), leaderFollower_(leaderFollower), xVars_(xVars), tVar_(tVar), dy_(dy), lazyData_(lazyData) {
                            tol_ = 1e-4;
                        };
    void main();
    int addBendersCuts();
    int solveSubs();
};
IloCplex::Callback BendersUserCallbackMC(IloEnv env, FollowerMC &followerMC, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData);
#endif //BILEVEL_CALLBACK_H
