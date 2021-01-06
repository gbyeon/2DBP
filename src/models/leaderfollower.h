/*
 * leaderfollower.h
 *
 *  Created on: Nov 5, 2020
 *      Author: geunyeongbyeon
 */
#ifndef BILEVEL_LEADERFOLLOWER_H
#define BILEVEL_LEADERFOLLOWER_H

#include "data.h"

class LeaderFollower {
public:
	/** default constructor */
	LeaderFollower();

    /* constructor */
    LeaderFollower(int n_l,
                int n_f,
                int m_l,
                int m_f,
                int ylb_cnt,
                int * ylb_ind,
                double * ylb_coef,
                int yub_cnt,
                int * yub_ind,
                double * yub_coef,
                double * fObj,
                double * lfObj,
                double * fC_rhs,
                int * fC_fV_cnt,
                double ** fC_fV_coef,
                int ** fC_fV_ind,
                int * fC_lV_cnt,
                double ** fC_lV_coef,
                int ** fC_lV_ind,
                double * lC_rhs,
                int * lC_fV_cnt,
                double ** lC_fV_coef,
                int ** lC_fV_ind,
                int * lC_lV_cnt,
                double ** lC_lV_coef,
                int ** lC_lV_ind);

	/** copy constructor */
	LeaderFollower(const LeaderFollower & rhs);

	/** default destructor */
	~LeaderFollower();

public:
    void loadProblem(Data &data);
    void createProblem ();
    void updateProblem (double* xVals, double O);
    int solve ();

    void getBendersTerms (IloExpr &lFterms, IloNumVarArray &xVars, double * xVals);
    void getResults ();

    IloCplex * getCplexPtr() {return &cplex_;};
    double getStatus() {return status_;};
    double getCheck() {return check_;};
    double getwVal() {return wVal_;};
    double getObjVal() {return objVal_;};

private:

    struct lfVars {

        IloNumVarArray y, psi, s;
        IloNumVar slack;

        lfVars (IloEnv * env) {
            y = IloNumVarArray (*env);
            psi = IloNumVarArray (*env);
            s = IloNumVarArray (*env);
        };

        lfVars() {};
    };

    struct lfConstrs {

        IloRangeArray lF, fF;
        IloRange fObj;     

        IloRangeArray yLBs, yUBs;

        lfConstrs (IloEnv * env) {

            lF = IloRangeArray(*env);
            fF = IloRangeArray(*env);

            yLBs = IloRangeArray(*env);
            yUBs = IloRangeArray(*env);
        };

        lfConstrs (){};
    };

    /* data */
    int n_l_;
    int n_f_;
    int m_l_;
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
    double * lfObj_;    /* c for cy */

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

    /* leader constraint data */
    double* lC_rhs_;
    
    /* follower variable coeffcients */
    int * lC_fV_cnt_;
    double ** lC_fV_coef_;
    int ** lC_fV_ind_;
    
    /* leader variable coeffcients */
    int * lC_lV_cnt_;
    double ** lC_lV_coef_;
    int ** lC_lV_ind_;

    /* opt members */
    lfVars vars_;
    lfConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray cy_;
    IloExpr cy_expr_, obj_slack_expr_;

    // results
    IloAlgorithm::Status status_;
    IloNum objVal_;
    IloNumArray psiVal_l_, psiVal_f_;
    IloNum wVal_;
    IloNumArray psiVal_yLBs_, psiVal_yUBs_;

    double * xbar_coef_;
    IloNum check_, constant_;

    chrono::duration<double> ticToc_;
};

#endif //BILEVEL_LEADERFOLLOWER_H
