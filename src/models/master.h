/*
 * master.h
 *
 *  Created on: Oct 29, 2020
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_MASTER_H
#define BILEVEL_MASTER_H

#include "data.h"
#include "callback.h"
#include "lazyCBBenders.h"
#include "lazyCBfUBnu.h"
#include "lazyCBBendersMC.h"
#include "heuristicCBIncumbentUpdate.h"
#include "usercutCBfUB.h"
#include "usercutCBfUBHeuristic.h"
#include "usercutCBBendersMC.h"
#include "branchCB.h"
#include "branchCBfObj.h"
#include "nodeCB.h"
#include "nodeCBfObj.h"

class Master {
public:

	/** default constructor */
	Master();

	/** copy constructor */
	Master(const Master & rhs);

	/** default destructor */
	~Master();

public:
    /* get problem data */
    void loadProblem (Data &data);

    /* create master problem */
    void createProblem (); 

    int solve();
    void solveCallback(Follower &follower, FollowerMC &followerMC, LeaderFollower &leaderFollower, Data &data);

    /* get functions */
    /* get environment */
    IloEnv * getEnvPtr() {return env_;}
    /* get current objective value */
    double getObjVal() {return objVal_;};
    /* get current best obj bound */
    double getBestObjVal() {return bestObjVal_;};
    /* get optimality gap */
    double getGap() {return gap_;};
    int getNumLocalCutsAdded() {return num_local_cuts_added_;};
    int getNumUserBranches() {return num_user_branches_;};
    double getFUbGap() {return fub_gap_;};
    /* get current x solution */
    IloNumArray getxVals() {return xVals_;};
    /* get time duration */
    chrono::duration<double> getTicToc(){return ticToc_;};

    /* set functions */
    void setTimeLimit(double timelimit) {
        timelimit_ = timelimit;
        cplex_.setParam(IloCplex::Param::TimeLimit, timelimit_);
    };
    void setParams(bool s, bool h, double ht){
        use_numerically_stable_cut_cb_ = s;
        use_heuristic_cb_ = h;
        heuristic_time_limit_ = ht;
    }

private:

    struct masterVars {

        IloNumVarArray x, y;
        IloNumVar t;

        masterVars (IloEnv * env) {
            x = IloNumVarArray (*env);
            y = IloNumVarArray (*env);
        };

        masterVars() {};
    };

    struct masterConstrs {

        IloRangeArray l, f;
        IloRangeArray yLBs, yUBs; 
        IloRangeArray xBds;
        // IloRange fObj;

        // For Benders & cut class = d
        IloRangeArray optCuts, feasCuts; // for master

        masterConstrs (IloEnv * env) {

            l = IloRangeArray(*env);
            f = IloRangeArray(*env);
            xBds = IloRangeArray(*env);
            yLBs = IloRangeArray(*env);
            yUBs = IloRangeArray(*env);

            optCuts = IloRangeArray(*env);
            feasCuts = IloRangeArray(*env);
        };

        masterConstrs (){};
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
    masterVars vars_;
    masterConstrs constrs_;

    IloEnv * env_;
    IloModel m_;
    IloCplex cplex_;

    IloNumArray xVals_;
    IloNum tVal_;

    IloNumArray cx_;

    IloNumArray dy_;
    IloExpr dy_expr_;
    IloNumArray cy_;
    IloExpr cy_expr_;

    /* record */
    double objVal_;
    double bestObjVal_;
    double gap_;
    int num_local_cuts_added_;
    int num_user_branches_;
    double fub_gap_;

    chrono::duration<double> ticToc_;
    double timelimit_;

    /* data for lazy callback */
    LazyData lazyData_;

    bool use_numerically_stable_cut_cb_;
    bool use_heuristic_cb_;
    double heuristic_time_limit_;
    double elapsed_time_;
};

#endif //BILEVEL_MASTER_H
