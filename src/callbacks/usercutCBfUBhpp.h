/*
 * usercutCBfUBX.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_USERCUTCBFUBHPP_H
#define BILEVEL_USERCUTCBFUBHPP_H

#include "callback.h"
#include "hpp.h"

class BendersUserCallbackHppI : public IloCplex::UserCutCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    Hpp &hpp_;
    IloExpr &dy_;
    chrono::duration<double> ticToc_;
    
public:
  ILOCOMMONCALLBACKSTUFF(BendersUserCallbackHpp);
  BendersUserCallbackHppI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Hpp &hpp)
    : IloCplex::UserCutCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), hpp_(hpp) {};
  void main();
};
IloCplex::Callback BendersUserCallbackHpp(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, Hpp &hpp);

#endif //BILEVEL_USERCUTCBFUBHPP_H