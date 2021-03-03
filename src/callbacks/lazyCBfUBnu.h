/*
 * usercutCBfUBnu.h
 *
 *  Created on: March 2, 2021
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_LAZYCBBENDERSNU_H
#define BILEVEL_LAZYCBBENDERSNU_H

#include "callback.h"
#include "follower.h"
#include "hpp.h"

class LazyCallbackfUBNuI : public IloCplex::LazyConstraintCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    Follower &follower_;
    Hpp &hpp_;
    IloExpr &dy_;
    Data &data_;
    chrono::duration<double> ticToc_;

    /* leader variable coeffcients */
    int * fC_lV_cnt_;
    double ** fC_lV_coef_;
    int ** fC_lV_ind_;

    IloNumArray barx_;
    double * xVals_;

    int i, j;
    int n_l, m_f;
    double **nu1_, **nu4_;

    /* dual ub */
    double dualUb_;

  ControlCallbackI::IntegerFeasibilityArray stat;
    
public:
  ILOCOMMONCALLBACKSTUFF(LazyCallbackfUBNu); 
  LazyCallbackfUBNuI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Follower &follower, Hpp &hpp, Data &data)
    : IloCplex::LazyConstraintCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), follower_(follower), hpp_(hpp), data_(data) {
      
      n_l = data.n_l_;
      m_f = data.m_f_;
      dualUb_ = 100000;

      barx_ = IloNumArray(env, n_l);
      xVals_ = new double [n_l];

      // stat = ControlCallbackI::IntegerFeasibilityArray(env);

      fC_lV_cnt_ = data.fC_lV_cnt_;
      fC_lV_coef_ = data.fC_lV_coef_;
      fC_lV_ind_ = data.fC_lV_ind_;

      nu1_ = new double * [m_f];
      nu4_ = new double * [m_f];

      for (int i = 0; i < m_f; i++) {
        if (fC_lV_cnt_[i] > 0) {
          nu1_[i] = new double[fC_lV_cnt_[i]]; 
          nu4_[i] = new double[fC_lV_cnt_[i]]; 
        }
      }
    };
  void main() override;
};
IloCplex::Callback LazyCallbackfUBNu(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                        IloExpr &dy, 
                        LazyData &lazyData, Follower &follower, Hpp &hpp, Data &data);

#endif //BILEVEL_LAZYCBBENDERSNU_H