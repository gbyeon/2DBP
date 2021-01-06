/*
 * callback.h
 *
 *  Created on: Nov 8, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_CALLBACK_H
#define BILEVEL_CALLBACK_H

#include <ilcplex/ilocplex.h>

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

#endif //BILEVEL_CALLBACK_H
