/*
 * lazyCBBenders.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_LAZYCBBENDERS_H
#define BILEVEL_LAZYCBBENDERS_H

#include "callback.h"
#include "follower.h"
#include "leaderfollower.h"

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
    void main () override;
    int addBendersCuts();
    int solveSubs();
};
IloCplex::Callback BendersLazyCallback(IloEnv env, Follower &follower, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData);

#endif //BILEVEL_LAZYCBBENDERS_H