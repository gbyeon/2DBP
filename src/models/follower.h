/*
 * follower.h
 *
 *  Created on: Nov 7, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_FOLLOWER_H
#define BILEVEL_FOLLOWER_H

#include "data.h"

class Follower {
public:
	/* default constructor */
	Follower();

    /* constructor */
    Follower(int n_l,
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
	Follower(const Follower & rhs);

	/* default destructor */
	~Follower();

public:
    void loadProblem(Data &data);
    void createProblem();
    
    void updateProblem(double* &xVals);
    void updateUBProblem(double* &xUBs, double* &xLBs);
    
    int solve();
    
    void getBendersTerms(IloExpr &termsfP, IloNumVarArray &xVars, double * xVals);

    void getResults();
    double getDyVal(IloNumArray &yVals);
    IloCplex * getCplexPtr() {return &cplex_;};
    double getObjVal() {return objVal_;};
    IloAlgorithm::Status getStatus() {return status_;};
    double getbigM() {return M_;};
    double getLB() {return fLB_;};
    double getCheck() {return check_;};
    IloNumArray * getpsiValsPtr() {return &psiVals_;};
    IloNumArray * getpsiValsPtr_yLBs() {return &psiVal_yLBs_;};
    IloNumArray * getpsiValsPtr_yUBs() {return &psiVal_yUBs_;};
    IloExpr * getdyExprPtr() {return &dy_expr_;};

    void setfUb(double M) {M_ = M;};
    void setLB(double fLB) {fLB_ = fLB;};

private:

    struct fVars {

        IloNumVarArray y;

        fVars (IloEnv * env) {
            y = IloNumVarArray (*env);
        };

        fVars() {};
    };

    struct fConstrs {

    IloRangeArray fF; 
    IloRangeArray yLBs, yUBs;

    fConstrs (IloEnv * env){

        fF = IloRangeArray(*env);
        yLBs = IloRangeArray(*env);
        yUBs = IloRangeArray(*env);
    };

    fConstrs() {};
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
    fVars vars_;
    fConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray dy_;
    IloExpr dy_expr_;

    IloAlgorithm::Status status_;
    IloNum objVal_;
    IloNumArray yVals_;
    IloNumArray psiVals_;
    IloNumArray psiVal_yLBs_, psiVal_yUBs_;

    double M_, fLB_;

    double * xbar_coef_;
    IloNum check_, constant_;

    chrono::duration<double> ticToc_;
};

#endif //BILEVEL_FOLLOWER_H
