/*
 * branchCB.h
 *
 *  Created on: Jan 23, 2021
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_BRANCHCB_H
#define BILEVEL_BRANCHCB_H

#include "callback.h"
#include "follower.h"
#include "leaderfollower.h"

class branchCallbackI : public IloCplex::BranchCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    IloNumVar &tVar_;
    LazyData &lazyData_;
    Follower &follower_;
    LeaderFollower &leaderFollower_;
    IloExpr &dy_;
    chrono::duration<double> ticToc_;
    int i, n_l;
    double tVal_;
    double fObjVal_;
    
    /* branch information */
    double *blb, *bub;
    bool *isEnded;
public:
  ILOCOMMONCALLBACKSTUFF(branchCallback);

  branchCallbackI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData, double *bLb, double *bUb, bool *isended, Follower &follower, LeaderFollower &leaderFollower)
    : IloCplex::BranchCallbackI(env), xVars_(xVars), yVars_(yVars), tVar_(tVar), dy_(dy), lazyData_(lazyData), blb(bLb), bub(bUb), follower_(follower), leaderFollower_(leaderFollower) 
    {
        n_l = lazyData.n_l;
        isEnded = isended;
    };
  void main();
};

IloCplex::Callback branchCallback(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloNumVar &tVar,
                                  IloExpr &dy, LazyData &lazyData, double *blb, double *bub, bool *isended, Follower &follower, LeaderFollower &leaderFollower);

#endif //BILEVEL_BRANCHCB_H