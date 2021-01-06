/*
 * followerMC.h
 *
 *  Created on: Nov 9, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_FOLLOWERMC_H
#define BILEVEL_FOLLOWERMC_H

#include "data.h"

class FollowerMC {
public:
	/* default constructor */
	FollowerMC();

    /* constructor */
    FollowerMC(int n_l,
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
	FollowerMC(const FollowerMC & rhs);

	/* default destructor */
	~FollowerMC();

public:
    void loadProblem(Data &data);
    void createProblem();
    
    void updateProblem(double* xVals);
    
    int solve();
    
    void getBendersTerms(IloExpr &termsfP, IloNumVarArray &xVars, double * xVals);

    void getResults();
    IloCplex * getCplexPtr() {return &cplex_;};
    double getObjVal() {return objVal_;};
    double getStatus() {return status_;};
    double getCheck() {return check_;};
    IloNumArray * getyValsPtr() {return &yVals_;};
    vector<IloNumArray> * getnuValsVector() {return &nuVals_;};
    IloExpr * getdyExprPtr() {return &dy_expr_;};
    double getpsiUB() {return psiUB_;};

    void setpsiUB(double ub) {psiUB_ = ub;};

private:

    struct fMCVars {

        IloNumVarArray y; 
        vector<IloNumVarArray> nu;
        IloNumVar r; 

        fMCVars (IloEnv * env) {
            y = IloNumVarArray (*env);
            nu = vector<IloNumVarArray> (4);
            for (int i = 0; i < 4; i++)
                nu[i] = IloNumVarArray (*env);
        };

        fMCVars() {};
    };

    struct fMCConstrs {

    IloRangeArray fFMC; 
    IloRangeArray yLBs, yUBs;
    IloRangeArray nu;
    IloRange r, normRay;

    fMCConstrs (IloEnv * env){

        fFMC = IloRangeArray(*env);
        yLBs = IloRangeArray(*env);
        yUBs = IloRangeArray(*env);
        nu = IloRangeArray(*env);
    };

    fMCConstrs() {};
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

    /* ub on psiVals */
    double psiUB_;

    /* opt members */
    fMCVars vars_;
    fMCConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray dy_;
    IloExpr dy_expr_;

    IloAlgorithm::Status status_;
    IloNum objVal_;
    IloNumArray yVals_;
    vector<IloNumArray> nuVals_;

    double * xbar_coef_;
    IloNum check_, constant_;

    chrono::duration<double> ticToc_;
};

#endif //BILEVEL_FOLLOWERMC_H
