/*
 * heuristicCBIncumbentUpdate.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_HEURISTICCBINCUMBENTUPDATE_H
#define BILEVEL_HEURISTICCBINCUMBENTUPDATE_H

#include "callback.h"

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

#endif //BILEVEL_HEURISTICCBINCUMBENTUPDATE_H