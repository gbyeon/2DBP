/*
 * usercutCBBendersMC.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_USERCUTCBBENDERSMC_H
#define BILEVEL_USERCUTCBBENDERSMC_H

#include "callback.h"
#include "followerMC.h"
#include "leaderfollower.h"

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
    void main() override;
    int addBendersCuts();
    int solveSubs();
};
IloCplex::Callback BendersUserCallbackMC(IloEnv env, FollowerMC &followerMC, 
                                        LeaderFollower &leaderFollower, 
                                        IloNumVarArray &xVars, IloNumVar &tVar, IloExpr &dy, LazyData &lazyData);
#endif //BILEVEL_USERCUTCBBENDERSMC_H