/*
 * followerx.h
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_FOLLOWERX_H
#define BILEVEL_FOLLOWERX_H

#include "data.h"

class FollowerX {
public:
	/* default constructor */
	FollowerX();

    /* constructor */
    FollowerX(int n_l,
            int n_f,
            int m_f,
            int ylb_cnt,
            int * ylb_ind,
            double * ylb_coef,
            int yub_cnt,
            int * yub_ind,
            double * yub_coef,
            double * fObj,
            double* fC_rhs,
            int * fC_fV_cnt,
            double ** fC_fV_coef,
            int ** fC_fV_ind,
            int * fC_lV_cnt,
            double ** fC_lV_coef,
            int ** fC_lV_ind);

	/* copy constructor */
	FollowerX(const FollowerX & rhs);

	/* default destructor */
	~FollowerX();

public:
    void loadProblem(Data &data);
    void createProblem();
    
    void updateUBProblem(IloNumArray &xUBs, IloNumArray &xLBs);
    
    int solve();

    void getResults();
    IloCplex * getCplexPtr() {return &cplex_;};
    double getObjVal() {return objVal_;};
    double getStatus() {return status_;};
    IloExpr * getdyExprPtr() {return &dy_expr_;};
    double getDyVal(IloNumArray &yVals);

private:

    struct fXVars {

        IloNumVarArray y;
        IloNumVarArray x;

        fXVars (IloEnv * env) {
            y = IloNumVarArray (*env);
            x = IloNumVarArray (*env);
        };

        fXVars() {};
    };

    struct fXConstrs {

    IloRangeArray f; 
    IloRangeArray yLBs, yUBs;
    IloRangeArray xBds;

    fXConstrs (IloEnv * env){

        f = IloRangeArray(*env);
        yLBs = IloRangeArray(*env);
        yUBs = IloRangeArray(*env);
        xBds = IloRangeArray(*env);
    };

    fXConstrs() {};
    };  

    /* data */
    int n_l_;
    int n_f_;
    int m_f_;

    /* y bounds data */
    int ylb_cnt_;
    int * ylb_ind_;
    double * ylb_coef_;
    int yub_cnt_;
    int * yub_ind_;
    double * yub_coef_;

    /* x bounds data */
    vector<vector<int>> ind_col_;
    vector<vector<double>> scale_col_;
    vector<double> xlb_;
    vector<double> xub_;
    map<int,int> map_varind_to_lvarind_;

    /* integer var indicator for leader var */
    int * is_integer_;
    
    /* objective function data */
    double * fObj_;     /* d for dy */

    /* follower constraint data */
    double* fC_rhs_;
    
    /* follower variable coeffcients */
    int * fC_fV_cnt_;
    double ** fC_fV_coef_;
    int ** fC_fV_ind_;
    
    /* leader variable coeffcients */
    int * fC_lV_cnt_;
    double ** fC_lV_coef_;
    int ** fC_lV_ind_;

    /* opt members */
    fXVars vars_;
    fXConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray dy_;
    IloExpr dy_expr_;

    IloAlgorithm::Status status_;
    IloNum objVal_;
    IloNumArray yVals_;

    chrono::duration<double> ticToc_;
};

#endif //BILEVEL_FOLLOWERX_H
