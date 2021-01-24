/*
 * usercutCBfUBX.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_USERCUTCBFUBX_H
#define BILEVEL_USERCUTCBFUBX_H

#include "callback.h"
#include "followerx.h"

class BendersUserCallbackXI : public IloCplex::UserCutCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    FollowerX &followerx_;
    IloExpr &dy_;
    chrono::duration<double> ticToc_;
    
public:
  ILOCOMMONCALLBACKSTUFF(BendersUserCallbackX); 
  BendersUserCallbackXI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, FollowerX &followerx)
    : IloCplex::UserCutCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), followerx_(followerx) {};
  void main();
};
IloCplex::Callback BendersUserCallbackX(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, FollowerX &followerx);

#endif //BILEVEL_USERCUTCBFUB_H