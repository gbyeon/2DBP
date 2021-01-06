/*
 * followerMC.cpp
 *
 *  Created on: Nov 9, 2020
 *      Author: geunyeongbyeon
 */
#undef NDEBUG
// #define FMC_SOLVE_DEBUG

#include "followerMC.h"
#include "macro.h"

FollowerMC::FollowerMC()
{
    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = fMCVars(env_);
    constrs_ = fMCConstrs(env_);
    yVals_ = IloNumArray (*env_);
    nuVals_ = vector<IloNumArray> (4);
    for (int i = 0; i < 4; i++)
        nuVals_[i] = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);
}

FollowerMC::FollowerMC(int n_l,
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

    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = fMCVars(env_);
    constrs_ = fMCConstrs(env_);
    yVals_ = IloNumArray (*env_);
    nuVals_ = vector<IloNumArray> (4);
    for (int i = 0; i < 4; i++)
        nuVals_[i] = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);

    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
}

/* copy constructor */
FollowerMC::FollowerMC(const FollowerMC & rhs) :
env_(rhs.env_),
m_(rhs.m_),
vars_(rhs.vars_),
constrs_(rhs.constrs_),
yVals_(rhs.yVals_),
nuVals_(rhs.nuVals_),
dy_(rhs.dy_),
dy_expr_(rhs.dy_expr_),
xbar_coef_(rhs.xbar_coef_)
{
    /* nothing to do */
}

FollowerMC::~FollowerMC()
{
    yVals_.end();
    for (int i = 0; i < 4; i++)
        nuVals_[i].end();
    dy_.end();
    dy_expr_.end();
    
    vars_.y.end();
    constrs_.fFMC.end();
    constrs_.yLBs.end();
    constrs_.yUBs.end();

    cplex_.end();
    m_.end();
    env_->end();

    delete[] xbar_coef_;
}

void FollowerMC::loadProblem (Data &data) {

    n_l_ = data.n_l_;
    n_f_ = data.n_f_;
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
    if (n_l_ > 0) xbar_coef_ = new double [n_l_];
}

void FollowerMC::createProblem () {
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

    /* add dual variables nu for McComick constraints */
    addnuVars();
    
    /* add follower constraint without x and with nu: By >= b + K^T nu */
    addfFMCConstr();

    /* add nu constraints: nu_0 + nu_1 - nu_2 - nu_3 == 1 */
    addnuConstr();
    
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

    /* to distinguish infeasible vs unboundedness */
    /* no unboundedness case as master has the follower constraint */
    // cplex_.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
    /* setting feas tol */
    cplex_.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);

//    cplex_.setParam(IloCplex::Param::TimeLimit, 10);
    // cplex_.exportModel("followerMC.lp");

   // dummy r for unbounded case
    vars_.r = IloNumVar (*env_, -IloInfinity, IloInfinity, ILOFLOAT, "r");
    m_.add(vars_.r);
    for (i = 0; i < constrs_.fFMC.getSize(); i++) {
       constrs_.fFMC[i].setLinearCoef(vars_.r, -constrs_.fFMC[i].getLB());
       constrs_.fFMC[i].setLB(0);
    } for (i = 0; i < constrs_.nu.getSize(); i++) {
       constrs_.nu[i].setLinearCoef(vars_.r, -constrs_.nu[i].getLB());
       constrs_.nu[i].setLB(0);
       constrs_.nu[i].setUB(0);
    } for (i = 0; i < constrs_.yUBs.getSize(); i++) {
       constrs_.yUBs[i].setLinearCoef(vars_.r, -constrs_.yUBs[i].getLB());
       constrs_.yUBs[i].setLB(0);
    } for (i = 0; i < constrs_.yLBs.getSize(); i++) {
       constrs_.yLBs[i].setLinearCoef(vars_.r, -constrs_.yLBs[i].getLB());
       constrs_.yLBs[i].setLB(0);
    }
    constrs_.r = IloRange(*env_, 1, vars_.r, 1, "r");
    m_.add(constrs_.r);
    
    IloExpr lhs(*env_);
    for (i = 0; i < n_f_; i++) {
        lhs += vars_.y[i];
    }
    for (i = 0; i < vars_.nu[0].getSize(); i++) {
        lhs += vars_.nu[0][i] + vars_.nu[1][i] + vars_.nu[2][i] + vars_.nu[3][i];
    }
    constrs_.normRay = IloRange(*env_, -1, lhs, 1, "normRay");
    lhs.end();

    // cplex_.exportModel("followerMC2.lp");
}

void FollowerMC::updateProblem (double* xVals) {

    IloNumArray objChg_nu1(*env_);
    IloNumArray objChg_nu2(*env_);

    int i, j;

    /* TODO: do not update all. only for those updated */
    for (i = 0; i < m_f_; i++) {
        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            if (fC_lV_coef_[i][j] > 0) {
                objChg_nu1.add(-fC_lV_coef_[i][j] * psiUB_ * (xVals[fC_lV_ind_[i][j]] - 1));
                objChg_nu2.add(fC_lV_coef_[i][j] * psiUB_ * xVals[fC_lV_ind_[i][j]]);
            } else {
                objChg_nu1.add(fC_lV_coef_[i][j] * psiUB_ * (xVals[fC_lV_ind_[i][j]] - 1));
                objChg_nu2.add(-fC_lV_coef_[i][j] * psiUB_ * xVals[fC_lV_ind_[i][j]]);
            }
        }
    }
    cplex_.getObjective().setLinearCoefs(vars_.nu[1], objChg_nu1);
    cplex_.getObjective().setLinearCoefs(vars_.nu[2], objChg_nu2);

    objChg_nu1.end();
    objChg_nu2.end();

    // cplex_.exportModel("followerMC_updated.lp");
}

int FollowerMC::solve () {
#ifdef FMC_SOLVE_DEBUG
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
    if (status_ == IloAlgorithm::Status::Unbounded) {

        constrs_.r.setLB(0);
        constrs_.r.setUB(0);
        m_.add(constrs_.normRay);

        // cplex_.exportModel("followerMC_ubd.lp");

        if (!cplex_.solve()) {
            env_->error() << "Failed to optimize FollowerMC." << endl;
            throw (-1);
        }

        getResults();

        cplex_.clear();
        constrs_.r.setLB(1);
        constrs_.r.setUB(1);

        m_.remove(constrs_.normRay);

        #ifdef FMC_SOLVE_DEBUG
            ticToc_ = (chrono::system_clock::now() - start_t);
            cout << "FollowerMC is unbounded, FollowerMC bd opt: " << objVal_ << ", time: " << ticToc_.count() << endl;
        #endif

    } else if (status_ == IloAlgorithm::Status::Optimal) {

        getResults();

        #ifdef FMC_SOLVE_DEBUG
            ticToc_ = (chrono::system_clock::now() - start_t);
            cout << "FollowerMC has optimum: " << objVal_ << ", time: " << ticToc_.count() << endl;
        #endif


    } else {
        cout << "FollowerMC terminated with status " << cplex_.getStatus() << endl;
        return 1;
    }
    return 0;
}

void FollowerMC::getResults () {

    objVal_ = cplex_.getObjValue();

    cplex_.getValues(yVals_, vars_.y);
    // for (int i = 0; i < 4; i++)
        // cplex_.getValues(nuVals_[i], vars_.nu[i]);
    cplex_.getValues(nuVals_[1], vars_.nu[1]);
    cplex_.getValues(nuVals_[2], vars_.nu[2]);
}

void FollowerMC::getBendersTerms(IloExpr &termsfPMC, IloNumVarArray &xVars, double * xVals) {
    
    int i, j; 

    for (i = 0; i < n_l_; i++)
        xbar_coef_[i] = 0;
    constant_ = 0;

    int pos = 0;
    for (i = 0; i < m_f_; i++) {
        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            if (fC_lV_coef_[i][j] > 0) {
                xbar_coef_[fC_lV_ind_[i][j]] += fC_lV_coef_[i][j] * psiUB_ * (nuVals_[2][pos] - nuVals_[1][pos]);
                constant_ += nuVals_[1][pos] * fC_lV_coef_[i][j] * psiUB_;
            } else {
                xbar_coef_[fC_lV_ind_[i][j]] += fC_lV_coef_[i][j] * psiUB_ * (nuVals_[1][pos] - nuVals_[2][pos]);
                constant_ -= nuVals_[1][pos] * fC_lV_coef_[i][j] * psiUB_;
            }
            pos++;
        }
    }

    for (i = 0; i < n_f_; i++) {
        constant_ += fObj_[i] * yVals_[i];
    }
    
    for (j = 0; j < n_l_; j++) {
        termsfPMC += xbar_coef_[j] * xVars[j];
    }
    termsfPMC += constant_;

    check_ = 0;
    for (j = 0; j < n_l_; j++) {
        check_ += xbar_coef_[j] * xVals[j];
    }
    check_ += constant_;

}