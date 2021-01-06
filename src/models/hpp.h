/*
 * hpp.h
 *
 *  Created on: Nov 6, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_HPP_H
#define BILEVEL_HPP_H

#include "data.h"

class Hpp {
public:

	/** default constructor */
	Hpp();

	/** copy constructor */
	Hpp(const Hpp & rhs);

	/** default destructor */
	~Hpp();

public:
    /* get problem data */
    void loadProblem (Data &data);

    /* create hpp problem */
    void createProblem (); 

    int solve();

    /* get functions */
    /* get environment */
    IloEnv * getEnvPtr() {return env_;}
    /* get current objective value */
    double getObjVal() {return objVal_;};
    /* get time duration */
    chrono::duration<double> getTicToc(){return ticToc_;};

    bool solvefUb();
    double getfUb(){return fUB_;};

    /* set functions */
    void setTimeLimit(double timelimit) {
        timelimit_ = timelimit;
        cplex_.setParam(IloCplex::Param::TimeLimit, timelimit_);
    };

private:

    struct hppVars {

        IloNumVarArray x, y;

        hppVars (IloEnv * env) {
            x = IloNumVarArray (*env);
            y = IloNumVarArray (*env);
        };

        hppVars() {};
    };

    struct hppConstrs {

        IloRangeArray l, f;
        IloRangeArray yLBs, yUBs; // for bounds of y
        IloRangeArray xBds;

        hppConstrs (IloEnv * env) {
            l = IloRangeArray(*env);
            f = IloRangeArray(*env);
            xBds = IloRangeArray(*env);
            yLBs = IloRangeArray(*env);
            yUBs = IloRangeArray(*env);
        };

        hppConstrs (){};
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
    double * llObj_;    /* c for cx */
    double * lfObj_;    /* c for cy */
    double lObj_constant_;

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
    hppVars vars_;
    hppConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray xVals_;

    IloNumArray cx_;

    IloNumArray dy_;
    IloExpr dy_expr_;
    IloNumArray cy_;
    IloExpr cy_expr_;

    double objVal_;
    double gap_;

    chrono::duration<double> ticToc_;
    double timelimit_;

    double fUB_;
};

#endif //BILEVEL_HPP_H
