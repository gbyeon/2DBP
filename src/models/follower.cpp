/*
 * follower.cpp
 *
 *  Created on: Oct 29, 2020
 *      Author: geunyeongbyeon
 */
#undef NDEBUG
// #define F_SOLVE_DEBUG

#include "follower.h"
#include "macro.h"

Follower::Follower()
{
    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = fVars(env_);
    constrs_ = fConstrs(env_);
    yVals_ = IloNumArray (*env_);
    psiVals_ = IloNumArray (*env_);
    psiVal_yLBs_ = IloNumArray (*env_);
    psiVal_yUBs_ = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);
}

Follower::Follower(int n_l,
                int n_f,
                int m_f,
                int ylb_cnt,
                int * ylb_ind,
                double * ylb_coef,
                int yub_cnt,
                int * yub_ind,
                double * yub_coef,
                double * fObj,
                double * fC_rhs,
                int * fC_fV_cnt,
                double ** fC_fV_coef,
                int ** fC_fV_ind,
                int * fC_lV_cnt,
                double ** fC_lV_coef,
                int ** fC_lV_ind) : 
n_l_(n_l), n_f_(n_f), m_f_(m_f), ylb_cnt_(ylb_cnt), yub_cnt_(yub_cnt)
{
    ylb_ind_ = ylb_ind;
    ylb_coef_ = ylb_coef;
    yub_ind_ = yub_ind;
    yub_coef_ = yub_coef;
    
    fObj_ = fObj;
    fC_rhs_ = fC_rhs;

    fC_fV_cnt_ = fC_fV_cnt;
    fC_lV_cnt_ = fC_lV_cnt;

    fC_fV_coef_ = fC_fV_coef;
    fC_fV_ind_ = fC_fV_ind;
    fC_lV_coef_ = fC_lV_coef;
    fC_lV_ind_ = fC_lV_ind;

    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = fVars(env_);
    constrs_ = fConstrs(env_);
    yVals_ = IloNumArray (*env_);
    psiVals_ = IloNumArray (*env_);
    psiVal_yLBs_ = IloNumArray (*env_);
    psiVal_yUBs_ = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);

    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
    else xbar_coef_ = NULL;
}

/* copy constructor */
Follower::Follower(const Follower & rhs) :
env_(rhs.env_),
m_(rhs.m_),
vars_(rhs.vars_),
constrs_(rhs.constrs_),
yVals_(rhs.yVals_),
psiVals_(rhs.psiVals_),
psiVal_yLBs_(rhs.psiVal_yLBs_),
psiVal_yUBs_(rhs.psiVal_yUBs_),
dy_(rhs.dy_),
dy_expr_(rhs.dy_expr_),
xbar_coef_(rhs.xbar_coef_)
{
    /* nothing to do */
}

Follower::~Follower()
{
    yVals_.end();
    psiVals_.end();
    psiVal_yLBs_.end();
    psiVal_yUBs_.end();
    dy_.end();
    dy_expr_.end();
    
    vars_.y.end();
    constrs_.fF.end();
    constrs_.yLBs.end();
    constrs_.yUBs.end();

    cplex_.end();
    m_.end();
    env_->end();

    freeArrayPtr(xbar_coef_);

    delete env_;
}

void Follower::loadProblem (Data &data) {

    n_l_ = data.n_l_;
    n_f_ = data.n_f_;
    m_f_ = data.m_f_;
    
    ylb_cnt_ = data.ylb_cnt_;
    yub_cnt_ = data.yub_cnt_;
    
   ylb_ind_ = data.ylb_ind_;
    ylb_coef_ = data.ylb_coef_;
    yub_ind_ = data.yub_ind_;
    yub_coef_ = data.yub_coef_;
    
    fObj_ = data.fObj_;
    fC_rhs_ = data.fC_rhs_;

    fC_fV_cnt_ = data.fC_fV_cnt_;
    fC_lV_cnt_ = data.fC_lV_cnt_;

    fC_fV_coef_ = data.fC_fV_coef_;
    fC_fV_ind_ = data.fC_fV_ind_;
    fC_lV_coef_ = data.fC_lV_coef_;
    fC_lV_ind_ = data.fC_lV_ind_;
    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
    else xbar_coef_ = NULL;
}

void Follower::createProblem () {
#ifdef DATA_DEBUG
    for (int i = 0; i < m_f_; i++) {
        for (int j = 0; j < fC_fV_cnt_[i]; j++) {
            cout << fC_fV_ind_[i][j] << ", " << fC_fV_coef_[i][j] << ", ";
        }
        for (int j = 0; j < fC_lV_cnt_[i]; j++) {
            cout << fC_lV_ind_[i][j] << ", " << fC_lV_coef_[i][j] << ", ";
        }
        cout << endl;
    }
    cout << "rhs: ";
    for (int i = 0; i < m_f_; i++) {
        cout << fC_rhs_[i] << ", ";
    }
    cout << endl;
#endif
    
    /* add follower variable y */
    addyVars();
    
    /* add follower constraint without x: By >= b */
    addfFConstr();
    
    /* add bound constraints on y */
    addyBdsConstr();
    
    /* initialize cplex */
    cplex_ = IloCplex(m_);

    /* set objective function: minimize dy */
    int i;
    for (i = 0; i < n_f_; i++)
        dy_.add(fObj_[i]);
    dy_expr_.setLinearCoefs(vars_.y, dy_);
    m_.add(IloMinimize(*env_, dy_expr_));

    /* turn off display */
    cplex_.setOut(env_->getNullStream());
    // cplex_.setParam(IloCplex::Param::Sifting::Display,2);

//    cplex_.setParam(IloCplex::Param::TimeLimit, 10);
    // cplex_.exportModel("follower.lp");

//    // dummy r
//    vars_.r = IloNumVar (env_->, 0, 0, ILOFLOAT, "r");
//    m_.add(vars_.r);
//    IloInt i;
//    for (i = 0; i < constrs_.fF.getSize(); i++)
//        fF[i].setLinearCoef(vars_.r, 1);

}

void Follower::updateProblem (double* &xVals) {

    double rhsChg;

    int i, j;

    for (i = 0; i < m_f_; i++) {

        rhsChg = fC_rhs_[i];

        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            rhsChg -= fC_lV_coef_[i][j] * xVals[fC_lV_ind_[i][j]];
        }

        constrs_.fF[i].setLB(rhsChg);
    }
}

void Follower::updateUBProblem (double* &xUBs, double* &xLBs) {
#ifdef F_SOLVE_DEBUG
cout << "solve follower UB: ";
#endif 
    double rhsChg;

    int i, j;

    // for (i=0; i < n_l_; i++) {
    //     if (xUBs[i] > 1)
    //         xUBs[i] = 1;
    //     if (xLBs[i] < 0)
    //         xLBs[i] = 0;
    // }

    for (i = 0; i < m_f_; i++) {

        rhsChg = fC_rhs_[i];

        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            if (fC_lV_coef_[i][j] * xLBs[fC_lV_ind_[i][j]] < fC_lV_coef_[i][j] * xUBs[fC_lV_ind_[i][j]])
                rhsChg -= fC_lV_coef_[i][j] * xLBs[fC_lV_ind_[i][j]];
            else
                rhsChg -= fC_lV_coef_[i][j] * xUBs[fC_lV_ind_[i][j]];
        }
        constrs_.fF[i].setLB(rhsChg);
    }
//    cplex_.exportModel("follower_fUB.lp");
}

int Follower::solve () {
#ifdef F_SOLVE_DEBUG
    /* tic */
    auto start_t = chrono::system_clock::now();
#endif
    /* solve */
    if (!cplex_.solve()) {
        // env_->error() << "Failed to optimize follower problem." << endl;
        // cout << "status: " << status_ << endl;
        // return 1;
    }
    
    status_ = cplex_.getStatus();
    getResults();

#ifdef F_SOLVE_DEBUG
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "Follower has optimum: " << objVal_ << ", time: " << ticToc_.count() << endl;
#endif

    return 0;
}

void Follower::getResults () {

    objVal_ = cplex_.getObjValue();

    cplex_.getValues(yVals_, vars_.y);
    cplex_.getDuals(psiVals_, constrs_.fF);
    cplex_.getDuals(psiVal_yLBs_, constrs_.yLBs);
    cplex_.getDuals(psiVal_yUBs_, constrs_.yUBs);
}

double Follower::getDyVal(IloNumArray &yVals) {
    int i; 
    double Dy = 0; 
    for (i = 0; i < n_f_; i++)
    {
        Dy += fObj_[i] * yVals[i];
    }

    return Dy;
}

void Follower::getBendersTerms(IloExpr &termsfP, IloNumVarArray &xVars, double * xVals) {
    
    int num;
    int i, j; 
    
    double psiVals[m_f_];
    double psiVal_yLBs [ylb_cnt_];
    double psiVal_yUBs [yub_cnt_];

    for (i = 0; i < m_f_; i++) {
        psiVals[i] = psiVals_[i];
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

    for (i = 0; i < m_f_; i++) {

        constant_ += fC_rhs_[i] * psiVals[i];

        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            xbar_coef_[fC_lV_ind_[i][j]] -= fC_lV_coef_[i][j] * psiVals[i];
        }
    }

    for (i = 0; i < ylb_cnt_; i++) {
        constant_ += ylb_coef_[i] * psiVal_yLBs[i];
    }
    for (i = 0; i < yub_cnt_; i++) {
        constant_ += -yub_coef_[i] * psiVal_yUBs[i];
    }

    for (j = 0; j < n_l_; j++) {
        termsfP += xbar_coef_[j] * xVars[j];
    }
    termsfP += constant_;

    check_ = 0;
    for (j = 0; j < n_l_; j++) {
        check_ += xbar_coef_[j] * xVals[j];
    }
    check_ += constant_;

}