/*
 * branchCBfObj.h
 *
 *  Created on: Jan 24, 2021
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_BRANCHCBFOBJ_H
#define BILEVEL_BRANCHCBFOBJ_H

#include "callback.h"
#include "follower.h"

class branchCallbackfObjI : public IloCplex::BranchCallbackI {
    IloNumVarArray& xVars_;
    IloNumVarArray& yVars_;
    LazyData &lazyData_;
    Follower &follower_;
    IloExpr &dy_;
    chrono::duration<double> ticToc_;
    int i, n_l, n_f;
    double fObjVal_ystar; // ystar: opt follower sol for the follower
    double fObjVal_ytilde; // ytilde: opt follower sol for the leader
    double current_master_objVal;
    
    /* branch information */
    IloInt64 *node_id;
    bool *node_type;

    IloNumArray barx;
    IloNumArray bary;
    double * xVals;
    
public:
    ILOCOMMONCALLBACKSTUFF(branchCallbackfObj);

    branchCallbackfObjI(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars, IloExpr &dy, LazyData &lazyData, Follower &follower, bool *nodeType, IloInt64 *nodeId) : IloCplex::BranchCallbackI(env), xVars_(xVars), yVars_(yVars), dy_(dy), lazyData_(lazyData), follower_(follower), node_type(nodeType), node_id(nodeId) 
    {
        n_l = lazyData.n_l;
        n_f = lazyData.n_f;
        barx = IloNumArray(env, n_l);
        bary = IloNumArray(env, n_f);
        xVals = new double [n_l];  
    };
    // ~branchCallbackfObjI() 
    // {
    //     delete xVals;
        
    //     barx.end();
    //     bary.end();
    // };
    void main() override;
};

IloCplex::Callback branchCallbackfObj(IloEnv env, IloNumVarArray& xVars, IloNumVarArray& yVars,
                                  IloExpr &dy, LazyData &lazyData, Follower &follower, bool *nodeType, IloInt64 *nodeId);

#endif //BILEVEL_BRANCHCBFOBJ_H