/*
 * nodeCB.h
 *
 *  Created on: Jan 23, 2021
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_NODECB_H
#define BILEVEL_NODECB_H

#include "callback.h"

class nodeSelectCallbackI : public IloCplex::NodeCallbackI {
public:
  ILOCOMMONCALLBACKSTUFF(nodeSelectCallback);
  nodeSelectCallbackI(IloEnv env, double *bLb, double *bUb, bool *isended, double fLb, double fRange) : IloCplex::NodeCallbackI(env) {
      call_ = 0;
      k = 1; l = 0; j = 0;
      blb = bLb;
      bub = bUb;
      flb = fLb;
      frange = fRange;
      isEnded = isended;
  };
  void main();
  int call_;
  int k, l, j;
  double *blb, *bub;
  double flb, frange;
  bool *isEnded;
};
IloCplex::Callback nodeSelectCallback(IloEnv env, double *bLb, double *bUb, bool *isended, double fLb, double fRange);

#endif //BILEVEL_NODECB_H