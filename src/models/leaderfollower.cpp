/*
 * leaderfollower.cpp
 *
 *  Created on: Nov 5, 2020
 *      Author: geunyeongbyeon
 */
#undef NDEBUG
// #define LF_SOLVE_DEBUG

#include "leaderfollower.h"
#include "macro.h"

LeaderFollower::LeaderFollower()
{
    env_ = new IloEnv,
    m_ = IloModel(*env_);
    vars_ = lfVars(env_);
    constrs_ = lfConstrs(env_);
    psiVal_l_ = IloNumArray (*env_);
    psiVal_f_ = IloNumArray (*env_);
    psiVal_yLBs_ = IloNumArray (*env_);
    psiVal_yUBs_ = IloNumArray (*env_);
    cy_ = IloNumArray (*env_);
    cy_expr_ = IloExpr (*env_);
    obj_slack_expr_ = IloExpr (*env_);
}

LeaderFollower::LeaderFollower(int n_l,
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
                int ** lC_lV_ind) : 
n_l_(n_l), n_f_(n_f), m_l_(m_l), m_f_(m_f), ylb_cnt_(ylb_cnt), yub_cnt_(yub_cnt)
{
    if (ylb_cnt_ > 0) {
        ylb_ind_ = new int [ylb_cnt];
        ylb_coef_ = new double [ylb_cnt];

        ylb_ind_ = ylb_ind;
        ylb_coef_ = ylb_coef;
    }

    if (yub_cnt_ > 0) {
        yub_ind_ = new int [yub_cnt];
        yub_coef_ = new double [yub_cnt];

        yub_ind_ = yub_ind;
        yub_coef_ = yub_coef;
    }

    if (n_f_ > 0) {
        fObj_ = new double [n_f];
        fObj_ = fObj;

        lfObj_ = new double [n_f];
        lfObj_ = lfObj;
    }

    if (m_f > 0) {
        
        fC_rhs_ = new double [m_f];
        fC_rhs_ = fC_rhs;

        fC_fV_cnt_ = new int [m_f];
        fC_lV_cnt_ = new int [m_f];

        fC_fV_cnt_ = fC_fV_cnt;
        fC_lV_cnt_ = fC_lV_cnt;

        fC_fV_coef_ = new double * [m_f];
        fC_fV_ind_ = new int * [m_f];
        
        fC_lV_coef_ = new double * [m_f];
        fC_lV_ind_ = new int * [m_f];

        for (int i = 0; i < m_f; i++)
        {
            if (fC_fV_cnt[i] > 0) {
                fC_fV_coef_[i] = new double [fC_fV_cnt[i]];
                fC_fV_ind_[i] = new int [fC_fV_cnt[i]];

                fC_fV_coef_[i] = fC_fV_coef[i];
                fC_fV_ind_[i] = fC_fV_ind[i];
            }

            if (fC_lV_cnt[i] > 0) {
                fC_lV_coef_[i] = new double [fC_lV_cnt[i]];
                fC_lV_ind_[i] = new int [fC_lV_cnt[i]];

                fC_lV_coef_[i] = fC_lV_coef[i];
                fC_lV_ind_[i] = fC_lV_ind[i];
            }
        }
    }

    if (m_l > 0) {
        
        lC_rhs_ = new double [m_l];
        lC_rhs_ = lC_rhs;

        lC_fV_cnt_ = new int [m_l];
        lC_lV_cnt_ = new int [m_l];

        lC_fV_cnt_ = lC_fV_cnt;
        lC_lV_cnt_ = lC_lV_cnt;

        lC_fV_coef_ = new double * [m_l];
        lC_fV_ind_ = new int * [m_l];
        
        lC_lV_coef_ = new double * [m_l];
        lC_lV_ind_ = new int * [m_l];

        for (int i = 0; i < m_l; i++)
        {
            if (lC_fV_cnt[i] > 0) {
                lC_fV_coef_[i] = new double [lC_fV_cnt[i]];
                lC_fV_ind_[i] = new int [lC_fV_cnt[i]];

                lC_fV_coef_[i] = lC_fV_coef[i];
                lC_fV_ind_[i] = lC_fV_ind[i];
            }

            if (lC_lV_cnt[i] > 0) {
                lC_lV_coef_[i] = new double [lC_lV_cnt[i]];
                lC_lV_ind_[i] = new int [lC_lV_cnt[i]];

                lC_lV_coef_[i] = lC_lV_coef[i];
                lC_lV_ind_[i] = lC_lV_ind[i];
            }
        }
    }

    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = lfVars(env_);
    constrs_ = lfConstrs(env_);
    psiVal_l_ = IloNumArray (*env_);
    psiVal_f_ = IloNumArray (*env_);
    psiVal_yLBs_ = IloNumArray (*env_);
    psiVal_yUBs_ = IloNumArray (*env_);
    cy_ = IloNumArray (*env_);
    cy_expr_ = IloExpr (*env_);
    obj_slack_expr_ = IloExpr (*env_);

    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
}

/* copy constructor */
LeaderFollower::LeaderFollower(const LeaderFollower & rhs) :
env_(rhs.env_),
m_(rhs.m_),
vars_(rhs.vars_),
constrs_(rhs.constrs_),
psiVal_l_(rhs.psiVal_l_),
psiVal_f_(rhs.psiVal_f_),
psiVal_yLBs_(rhs.psiVal_yLBs_),
psiVal_yUBs_(rhs.psiVal_yUBs_),
cy_(rhs.cy_),
cy_expr_(rhs.cy_expr_),
obj_slack_expr_(rhs.obj_slack_expr_),
xbar_coef_(rhs.xbar_coef_)
{
    /* nothing to do */
}

LeaderFollower::~LeaderFollower()
{
    psiVal_l_.end();
    psiVal_f_.end();
    psiVal_yLBs_.end();
    psiVal_yUBs_.end();
    cy_.end();
    cy_expr_.end();
    obj_slack_expr_.end();
    
    vars_.y.end();
    vars_.psi.end();
    vars_.s.end();
      
    constrs_.lF.end();
    constrs_.fF.end();
    constrs_.yLBs.end();
    constrs_.yUBs.end();
    
    cplex_.end();
    m_.end();
    env_->end();

    delete[] xbar_coef_;
}

void LeaderFollower::loadProblem (Data &data) {

    n_l_ = data.n_l_;
    n_f_ = data.n_f_;
    m_l_ = data.m_l_;
    m_f_ = data.m_f_;

    ylb_cnt_ = data.ylb_cnt_;
    yub_cnt_ = data.yub_cnt_;

    if (ylb_cnt_ > 0) {
        ylb_ind_ = new int [ylb_cnt_];
        ylb_coef_ = new double [ylb_cnt_];

        ylb_ind_ = data.ylb_ind_;
        ylb_coef_ = data.ylb_coef_;
    }

    if (yub_cnt_ > 0) {
        yub_ind_ = new int [yub_cnt_];
        yub_coef_ = new double [yub_cnt_];

        yub_ind_ = data.yub_ind_;
        yub_coef_ = data.yub_coef_;
    }

    if (n_f_ > 0) {
        fObj_ = new double [n_f_];
        fObj_ = data.fObj_;

        lfObj_ = new double [n_f_];
        lfObj_ = data.lfObj_;
    }

     if (m_f_ > 0) {
        
        fC_rhs_ = new double [m_f_];
        fC_rhs_ = data.fC_rhs_;

        fC_fV_cnt_ = new int [m_f_];
        fC_lV_cnt_ = new int [m_f_];

        fC_fV_cnt_ = data.fC_fV_cnt_;
        fC_lV_cnt_ = data.fC_lV_cnt_;

        fC_fV_coef_ = new double * [m_f_];
        fC_fV_ind_ = new int * [m_f_];
        
        fC_lV_coef_ = new double * [m_f_];
        fC_lV_ind_ = new int * [m_f_];

        for (int i = 0; i < m_f_; i++)
        {
            if (fC_fV_cnt_[i] > 0) {
                fC_fV_coef_[i] = new double [fC_fV_cnt_[i]];
                fC_fV_ind_[i] = new int [fC_fV_cnt_[i]];

                fC_fV_coef_[i] = data.fC_fV_coef_[i];
                fC_fV_ind_[i] = data.fC_fV_ind_[i];
            }

            if (fC_lV_cnt_[i] > 0) {
                fC_lV_coef_[i] = new double [fC_lV_cnt_[i]];
                fC_lV_ind_[i] = new int [fC_lV_cnt_[i]];

                fC_lV_coef_[i] = data.fC_lV_coef_[i];
                fC_lV_ind_[i] = data.fC_lV_ind_[i];
            }
        }
    }

    if (m_l_ > 0) {
        
        lC_rhs_ = new double [m_l_];
        lC_rhs_ = data.lC_rhs_;

        lC_fV_cnt_ = new int [m_l_];
        lC_lV_cnt_ = new int [m_l_];

        lC_fV_cnt_ = data.lC_fV_cnt_;
        lC_lV_cnt_ = data.lC_lV_cnt_;

        lC_fV_coef_ = new double * [m_l_];
        lC_fV_ind_ = new int * [m_l_];
        
        lC_lV_coef_ = new double * [m_l_];
        lC_lV_ind_ = new int * [m_l_];

        for (int i = 0; i < m_l_; i++)
        {
            if (lC_fV_cnt_[i] > 0) {
                lC_fV_coef_[i] = new double [lC_fV_cnt_[i]];
                lC_fV_ind_[i] = new int [lC_fV_cnt_[i]];

                lC_fV_coef_[i] = data.lC_fV_coef_[i];
                lC_fV_ind_[i] = data.lC_fV_ind_[i];
            }

            if (lC_lV_cnt_[i] > 0) {
                lC_lV_coef_[i] = new double [lC_lV_cnt_[i]];
                lC_lV_ind_[i] = new int [lC_lV_cnt_[i]];

                lC_lV_coef_[i] = data.lC_lV_coef_[i];
                lC_lV_ind_[i] = data.lC_lV_ind_[i];
            }
        }
    }

    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
}

void LeaderFollower::createProblem() {

    /* tic */
    auto start_t = chrono::system_clock::now();

    addyVars();

#ifdef LF_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in adding vars to LF: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* add follower constraints without x: By >= b */
    addfFConstr();

#ifdef LF_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addfFConstr to LF: " << ticToc_.count() << endl;
    /* tic */    
    start_t = chrono::system_clock::now();
#endif

    /* add leader constraints without x: Gy >= h */
    addlFConstr();
    
#ifdef LF_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addlFConstr to LF: " << ticToc_.count() << endl;
    /* tic */    
    start_t = chrono::system_clock::now();
#endif

    /* add follower objective constraint: dy <= 0 */
    addfObjConstr();

#ifdef LF_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addfObjConstr to LF: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* add lb and ub constraints on y */
    addyBdsConstr();

#ifdef LF_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addyBdsConstr to LF: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    // slack for infeasible case ===
    vars_.slack = IloNumVar (*env_, 0, IloInfinity, ILOFLOAT, "slack");
    obj_slack_expr_.setLinearCoef(vars_.slack, 1);
    
    /* define cplex */
    cplex_ = IloCplex(m_);

    /* set objective function: minimized cy */
    m_.add(IloMinimize(*env_));
    for (int i = 0; i < n_f_; i++) {
        cy_.add(lfObj_[i]);
    }
    cy_expr_.setLinearCoefs(vars_.y, cy_);
    cplex_.getObjective().setExpr(cy_expr_);

    /* turn off display */
    cplex_.setOut(env_->getNullStream());
    // cplex_.setParam(IloCplex::Param::TimeLimit, sett.timeLimit);

    // cplex_.exportModel("leaderFollower.lp");

}

void LeaderFollower::updateProblem (double* xVals, double O) {

    double rhsChg;

    int i, j;

    for (i = 0; i < m_l_; i++) {

        rhsChg = lC_rhs_[i];

        for (j = 0; j < lC_lV_cnt_[i]; j++) {

            rhsChg -= lC_lV_coef_[i][j] * xVals[lC_lV_ind_[i][j]];
        }

        constrs_.lF[i].setLB(rhsChg);
    }

    for (i = 0; i < m_f_; i++) {

        rhsChg = fC_rhs_[i];

        for (j = 0; j < fC_lV_cnt_[i]; j++) {

            rhsChg -= fC_lV_coef_[i][j] * xVals[fC_lV_ind_[i][j]];

        }

        constrs_.fF[i].setLB(rhsChg);
    }

    constrs_.fObj.setLB(-O);
}

int LeaderFollower::solve () {

    /* tic */
    auto start_t = chrono::system_clock::now();

    if (!cplex_.solve()) {
//        env_->error() << "Failed to optimize s2." << endl;
//        throw (-1);
    }

    status_ = cplex_.getStatus();

    if (status_ == IloAlgorithm::Status::Infeasible) {

        cplex_.getObjective().setExpr(obj_slack_expr_);

        int i;
        for (i = 0; i < m_f_; i++)
            constrs_.fF[i].setLinearCoef(vars_.slack, 1);
        for (i = 0; i < m_l_; i++)
            constrs_.lF[i].setLinearCoef(vars_.slack, 1);
        constrs_.fObj.setLinearCoef(vars_.slack, 1);

        if (!cplex_.solve()) {
            env_->error() << "Failed to optimize leaderFollower." << endl;
            throw (-1);
        }

        getResults();

        for (i = 0; i < m_f_; i++)
            constrs_.fF[i].setLinearCoef(vars_.slack, 0);
        for (i = 0; i < m_l_; i++)
            constrs_.lF[i].setLinearCoef(vars_.slack, 0);
        constrs_.fObj.setLinearCoef(vars_.slack, 0);

        cplex_.getObjective().setExpr(cy_expr_);

        #ifdef LF_SOLVE_DEBUG
            ticToc_ = (chrono::system_clock::now() - start_t);
            cout << "LeaderFollower is infeasible, LeaderFollower bd opt: " << objVal_ << ", time: " << ticToc_.count() << endl;
        #endif

    } else if (status_ == IloAlgorithm::Status::Optimal) {

        getResults();

        #ifdef LF_SOLVE_DEBUG
            ticToc_ = (chrono::system_clock::now() - start_t);
            cout << "LeaderFollower has optimum: " << objVal_ << ", time: " << ticToc_.count() << endl;
        #endif


    } else {
        cout << "LeaderFollower terminated with status " << cplex_.getStatus() << endl;
        return 1;
    }

    return 0;
}

void LeaderFollower::getResults () {

    objVal_ = cplex_.getObjValue();

    wVal_ = cplex_.getDual(constrs_.fObj);

    cplex_.getDuals(psiVal_l_, constrs_.lF);
    cplex_.getDuals(psiVal_f_, constrs_.fF);

    cplex_.getDuals(psiVal_yLBs_, constrs_.yLBs);
    cplex_.getDuals(psiVal_yUBs_, constrs_.yUBs);
}

void LeaderFollower::getBendersTerms (IloExpr &termLf, IloNumVarArray &xVars, double *xVals) {

    int num;
    int i, j; 

    double psiVal_l [m_l_];
    double psiVal_f [m_f_];

    double psiVal_yLBs [ylb_cnt_];
    double psiVal_yUBs [yub_cnt_];

    for (i = 0; i < m_l_; i++) {
        psiVal_l[i] = psiVal_l_[i];
    }
    for (i = 0; i < m_f_; i++) {
        psiVal_f[i] = psiVal_f_[i];
    }
    for (i = 0; i < ylb_cnt_; i++) {
        psiVal_yLBs[i] = psiVal_yLBs_[i];
    }
    for (i = 0; i < yub_cnt_; i++) {
        psiVal_yUBs[i] = psiVal_yUBs_[i];
    }

    for (i = 0; i < n_l_; i++)
        xbar_coef_[i] = 0;
    constant_ = 0;

    for (i = 0; i < m_l_; i++) {

        constant_ += lC_rhs_[i] * psiVal_l[i];

        for (j = 0; j < lC_lV_cnt_[i]; j++) {
            xbar_coef_[lC_lV_ind_[i][j]] -= lC_lV_coef_[i][j] * psiVal_l[i];
        }
    }

    for (i = 0; i < m_f_; i++) {
        
        constant_ += fC_rhs_[i] * psiVal_f[i];

        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            xbar_coef_[fC_lV_ind_[i][j]] -= fC_lV_coef_[i][j] * psiVal_f[i];
        }
    }

    for (i = 0; i < ylb_cnt_; i++) {
        constant_ += ylb_coef_[i] * psiVal_yLBs[i];
    }
    for (i = 0; i < yub_cnt_; i++) {
        constant_ += -yub_coef_[i] * psiVal_yUBs[i];
    }

    for (j = 0; j < n_l_; j++) {
        termLf += xbar_coef_[j] * xVars[j];
    }
    termLf += constant_;

//    if (sett.outputConsole) {
        check_ = 0;
        for (j = 0; j < n_l_; j++) {
            check_ += xbar_coef_[j] * xVals[j];
        }
        check_ += constant_;
//    }

}