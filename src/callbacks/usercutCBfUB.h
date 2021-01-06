/*
 * usercutCBfUB.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_USERCUTCBFUB_H
#define BILEVEL_USERCUTCBFUB_H

#include "callback.h"
#include "follower.h"

class BendersUserCallbackI : public IloCplex::UserCutCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    Follower &follower_;
    IloExpr &dy_;
    
public:
  ILOCOMMONCALLBACKSTUFF(BendersUserCallback) \
  BendersUserCallbackI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Follower &follower)
    : IloCplex::UserCutCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), follower_(follower) {};
  void main();
};
IloCplex::Callback BendersUserCallback(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, Follower &follower);

#endif //BILEVEL_USERCUTCBFUB_H