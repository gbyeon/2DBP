/*
 * hpp.cpp
 *
 *  Created on: Nov 6, 2020
 *      Author: geunyeongbyeon
 */

#include "hpp.h"
#include "macro.h"

Hpp::Hpp()
{
    env_ = new IloEnv,
    m_ = IloModel(*env_);
    vars_ = hppVars(env_);
    constrs_ = hppConstrs(env_);
    xVals_ = IloNumArray (*env_);
    cx_ = IloNumArray (*env_);
    cy_ = IloNumArray (*env_);
    cy_expr_ = IloExpr (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);
    ax_expr_ = IloExpr (*env_);
}

/* copy constructor */
Hpp::Hpp(const Hpp & rhs) :
env_(rhs.env_),
m_(rhs.m_),
vars_(rhs.vars_),
constrs_(rhs.constrs_),
xVals_(rhs.xVals_),
cx_(rhs.cx_),
cy_(rhs.cy_),
cy_expr_(rhs.cy_expr_),
dy_(rhs.dy_),
dy_expr_(rhs.dy_expr_),
ax_expr_(rhs.ax_expr_)
{
    /* nothing to do */
}

Hpp::~Hpp()
{
    xVals_.end();
    cx_.end();
    dy_.end();
    dy_expr_.end();
    cy_.end();
    cy_expr_.end();
    ax_expr_.end();

    vars_.x.end();
    vars_.y.end();
    
    constrs_.l.end();
    constrs_.f.end();
    constrs_.yLBs.end();
    constrs_.yUBs.end();
    constrs_.xBds.end();

    obj_.end();
    obj_expr_.end();

    cplex_.end();
    m_.end();
    env_->end();
    delete env_;
}

void Hpp::loadProblem (Data &data) {

    n_l_ = data.n_l_;
    n_f_ = data.n_f_;
    m_l_ = data.m_l_;
    m_f_ = data.m_f_;

    xlb_ = data.xlb_;
    xub_ = data.xub_;
    ind_col_ = data.ind_col_;
    scale_col_ = data.scale_col_;
    map_varind_to_lvarind_ = data.map_varind_to_lvarind_;
    lObj_constant_ = data.lObj_constant_;

    is_integer_ = data.is_integer_;
    
    llObj_ = data.llObj_;
    fObj_ = data.fObj_;
    lfObj_ = data.lfObj_;

    ylb_cnt_ = data.ylb_cnt_;
    ylb_ind_ = data.ylb_ind_;
    ylb_coef_ = data.ylb_coef_;
    yub_cnt_ = data.yub_cnt_;
    yub_ind_ = data.yub_ind_;
    yub_coef_ = data.yub_coef_;
    
    fC_rhs_ = data.fC_rhs_;
    fC_fV_cnt_ = data.fC_fV_cnt_;
    fC_lV_cnt_ = data.fC_lV_cnt_;
    fC_fV_coef_ = data.fC_fV_coef_;
    fC_fV_ind_ = data.fC_fV_ind_;
    fC_lV_coef_ = data.fC_lV_coef_;
    fC_lV_ind_ = data.fC_lV_ind_;

    lC_rhs_ = data.lC_rhs_;
    lC_fV_cnt_ = data.lC_fV_cnt_;
    lC_lV_cnt_ = data.lC_lV_cnt_;
    lC_fV_coef_ = data.lC_fV_coef_;
    lC_fV_ind_ = data.lC_fV_ind_;
    lC_lV_coef_ = data.lC_lV_coef_;
    lC_lV_ind_ = data.lC_lV_ind_;

    obj_ = IloNumArray(*env_, n_l_);
    obj_expr_ = IloExpr(*env_);
}

void Hpp::createProblem(bool relaxed) {

    int i, j;

#ifdef HPP_BUILD_DEBUG    
    /* tic */
    auto start_t = chrono::system_clock::now();
#endif

    if (!relaxed) {
        addxVars();
    } else {
        addxVars(true);
    }
    addyVars();

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in adding vars to hpp: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* add leader constraints: Gx + Gy >= h */
    addlConstr();

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addlConstr to hpp: " << ticToc_.count() << endl;    
    /* tic*/
    start_t = chrono::system_clock::now();
#endif

    /* add follower constraints with x: Ax + By >= b */
    addfConstr();

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addfConstr to hpp: " << ticToc_.count() << endl;
    /* tic*/
    start_t = chrono::system_clock::now();
#endif

    /* add bound constraints on x */
    addxBdsConstr();

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in addxBdsConstr to hpp: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* add bound constraints on y */
    addyBdsConstr();

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in adding yBdconstrs to hpp: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* define cplex */
    cplex_ = IloCplex(m_);

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in defining hpp cplex: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* set objective function: minimized cx + t */
    m_.add(IloMinimize(*env_));
    for (int i = 0; i < n_l_; i++) {
        cx_.add(llObj_[i]);
    }
    for (i = 0; i < n_f_; i++) {
        cy_.add(lfObj_[i]);
    }
    cy_expr_.setLinearCoefs(vars_.y, cy_);
    cplex_.getObjective().setLinearCoefs(vars_.x, cx_);
    cplex_.getObjective().setLinearCoefs(vars_.y, cy_);
    cplex_.getObjective().setConstant(lObj_constant_);

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in setting hpp obj: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    /* set cplex parameters */
    cplex_.setOut(env_->getNullStream());

    cplex_.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 1e-9);
    cplex_.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
    cplex_.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, CPX_VARSEL_STRONG);

    setTimeLimit(3600);
    // /* set branching priority */
    // /* if x appears in follower problem, increase priority */
    // vector<int> priority(n_l_, 0);
    
    // for (i = 0; i < m_f_; i++) {
    //     for (j = 0; j < fC_lV_cnt_[i]; j++) {
    //         priority[fC_lV_ind_[i][j]] += 1;
    //         cplex_.setPriority(vars_.x[fC_lV_ind_[i][j]], priority[fC_lV_ind_[i][j]]);
    //     }
    // }

#ifdef HPP_BUILD_DEBUG
    /* toc */
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "time in setting branching priority of hpp: " << ticToc_.count() << endl;
    /* tic */
    start_t = chrono::system_clock::now();
#endif

    for (i = 0; i < n_f_; i++)
        dy_.add(fObj_[i]);
    dy_expr_.setLinearCoefs(vars_.y, dy_);

//    cplex_.exportModel("hpp.lp");

    IloNumArray ax = IloNumArray(*env_, n_l_);
    double ax_constant = 0;
    for (int i = 0; i < n_l_; i++) 
        ax[i] = 0;

    for (int i = 0; i < m_f_; i++) {
        for (int j = 0; j < fC_lV_cnt_[i]; j++) {
            if (fC_lV_coef_[i][j] > 0) {
                ax[fC_lV_ind_[i][j]] += fC_lV_coef_[i][j];
            } else {
                ax[fC_lV_ind_[i][j]] -= fabs(fC_lV_coef_[i][j]);
                ax_constant += fabs(fC_lV_coef_[i][j]);
            }
        }
    }
    ax_expr_.setLinearCoefs(vars_.x, ax);
    ax_expr_.setConstant(ax_constant);
}

void Hpp::updateUBProblem (IloNumArray &xUBs, IloNumArray &xLBs) {
    
    for (int i=0; i < n_l_; i++) {
        if (xUBs[i] > 1)
            xUBs[i] = 1;
        if (xLBs[i] < 0)
            xLBs[i] = 0;
    }

    cplex_.getObjective().setSense(IloObjective::Sense::Maximize);
    cplex_.getObjective().setExpr(dy_expr_);

    vars_.x.setBounds(xLBs, xUBs);
    // cplex_.exportModel("hpp_ub_updated.lp");
}

int Hpp::solve() {
    auto start_t = chrono::system_clock::now();

    if (!cplex_.solve()) {
        env_->error() << "Failed to optimize hpp." << endl;
        // throw (-1);
        return 1;
    }

    status_ = cplex_.getStatus();

    if (status_ != IloAlgorithm::Status::Optimal
        && status_ != IloAlgorithm::Status::Unbounded ) {
        cout << "hpp is terminated w/ status " << status_ << endl;
        return 1;
    }

    if (status_ == IloAlgorithm::Status::Unbounded)
        objVal_ = -3e+6;
    else
        objVal_ = cplex_.getObjValue();

    cplex_.getValues(xVals_, vars_.x);
    
    ticToc_ = (chrono::system_clock::now() - start_t);
    // cout << "Hpp obj: " << objVal_ << endl;
    // cout << "time in solving Hpp: " << ticToc_.count() << endl;

    return 0;
}

bool Hpp::solvefUb () {

    // set timelimit and use best upper bounds

    cplex_.getObjective().setSense(IloObjective::Sense::Maximize);
    cplex_.getObjective().setExpr(dy_expr_);
    
    double original_time_limit = timelimit_;
    setTimeLimit(100);
    
//    cplex_.exportModel("fub.lp");

    if (!cplex_.solve()) {
        env_->error() << "Failed to optimize fUB." << endl;
        cout << cplex_.getStatus() << endl;
        fUB_ = 1e+9;
        return true;
        // return false;
    }

    status_ = cplex_.getStatus();
    if (status_ == IloAlgorithm::Status::Optimal || status_ == IloAlgorithm::Status::Feasible) {

        fUB_ = cplex_.getBestObjValue();//cplex_.getObjValue();
        // cout << "f.M: " << fUB_ << endl;

    } else {

        cout << "Follower does not have finite UB" << endl;
        return false;
    }

    setTimeLimit(original_time_limit);
    return true;
}

bool Hpp::solvefLb () {

    // set timelimit and use best upper bounds

    cplex_.getObjective().setSense(IloObjective::Sense::Minimize);
    cplex_.getObjective().setExpr(dy_expr_);
    
    double original_time_limit = timelimit_;
    setTimeLimit(100);
    
//    cplex_.exportModel("fub.lp");

    if (!cplex_.solve()) {
        env_->error() << "Failed to optimize fLB." << endl;
        return false;
    }

    status_ = cplex_.getStatus();
    if (status_ == IloAlgorithm::Status::Optimal || status_ == IloAlgorithm::Status::Feasible) {

        fLB_ = cplex_.getBestObjValue();//cplex_.getObjValue();
        cout << "f.LB: " << fLB_ << endl;

    } else {

        cout << "Follower does not have finite LB" << endl;
        return false;
    }

    setTimeLimit(original_time_limit);
    return true;
}

double Hpp::getDyVal(IloNumArray &yVals) {
    int i; 
    double Dy = 0; 
    for (i = 0; i < n_f_; i++)
    {
        Dy += fObj_[i] * yVals[i];
    }

    return Dy;
}

bool Hpp::updateObjNu(double **nu1, double **nu4, IloNumArray &xlbs, IloNumArray &xubs, int &nnodes) {

    if (nnodes > 0) 
        vars_.x.setBounds(xlbs, xubs); 


    cplex_.getObjective().setSense(IloObjective::Sense::Minimize);

    dualUb_ = 100000;
    
    for (int i = 0; i < n_l_; i++) 
        obj_[i] = 0;

    for (int i = 0; i < m_f_; i++) {
        for (int j = 0; j < fC_lV_cnt_[i]; j++) {
            obj_[fC_lV_ind_[i][j]] += (nu1[i][j] - nu4[i][j]) * dualUb_;
        }
    }

    obj_expr_.setLinearCoefs(vars_.x, obj_);
    cplex_.getObjective().setExpr(obj_expr_);
    
    // double original_time_limit = timelimit_;
    // setTimeLimit(100);
    
//    cplex_.exportModel("hppnu.lp");

    // if (!cplex_.solve()) {
    //     env_->error() << "Failed to optimize fLB." << endl;
    //     return false;
    // }

    // status_ = cplex_.getStatus();
    // if (status_ == IloAlgorithm::Status::Optimal || status_ == IloAlgorithm::Status::Feasible) {

    //     fLB_ = cplex_.getBestObjValue();//cplex_.getObjValue();
    //     cout << "f.LB: " << fLB_ << endl;

    // } else {

    //     cout << "Follower does not have finite LB" << endl;
    //     return false;
    // }

    // setTimeLimit(original_time_limit);

    return true;
}

bool Hpp::updateObjax(IloNumArray &xlbs, IloNumArray &xubs, IloNumArray &barx) {

    vars_.x.setBounds(xlbs, xubs); 

    cplex_.getObjective().setSense(IloObjective::Sense::Minimize);

    cplex_.getObjective().setExpr(ax_expr_);
    
    double original_time_limit = timelimit_;
    setTimeLimit(100);
    
//    cplex_.exportModel("hppnu.lp");

    if (!cplex_.solve()) {
        env_->error() << "Failed to optimize fLB." << endl;
        return false;
    }

    cplex_.getValues(barx, vars_.x);

    // status_ = cplex_.getStatus();
    // if (status_ == IloAlgorithm::Status::Optimal || status_ == IloAlgorithm::Status::Feasible) {

    //     fLB_ = cplex_.getBestObjValue();//cplex_.getObjValue();
    //     cout << "f.LB: " << fLB_ << endl;

    // } else {

    //     cout << "Follower does not have finite LB" << endl;
    //     return false;
    // }

    setTimeLimit(original_time_limit);

    return true;
}

bool Hpp::updateObjaxpsi(IloNumArray &xlbs, IloNumArray &xubs, IloNumArray &barx, IloNumArray *psiVal) {

    vars_.x.setBounds(xlbs, xubs); 

    cplex_.getObjective().setSense(IloObjective::Sense::Minimize);

    IloNumArray psiax = IloNumArray(*env_, n_l_);
    IloExpr psiax_expr = IloExpr(*env_);
    for (int i = 0; i < n_l_; i++) 
        psiax[i] = 0;

    for (int i = 0; i < m_f_; i++) {
        for (int j = 0; j < fC_lV_cnt_[i]; j++) {
            // if (fC_lV_coef_[i][j] > 0) {
                psiax[fC_lV_ind_[i][j]] += fC_lV_coef_[i][j] * (*psiVal)[i];
            // } else {
                // psiax[fC_lV_ind_[i][j]] -= fabs(fC_lV_coef_[i][j]) * (*psiVal)[i];
            // }
        }
    }
    psiax_expr.setLinearCoefs(vars_.x, psiax);
    cplex_.getObjective().setExpr(psiax_expr);
    
    double original_time_limit = timelimit_;
    setTimeLimit(100);
    
//    cplex_.exportModel("hppnu.lp");

    if (!cplex_.solve()) {
        env_->error() << "Failed to optimize fLB." << endl;
        return false;
    }

    cplex_.getValues(barx, vars_.x);

    // status_ = cplex_.getStatus();
    // if (status_ == IloAlgorithm::Status::Optimal || status_ == IloAlgorithm::Status::Feasible) {

    //     fLB_ = cplex_.getBestObjValue();//cplex_.getObjValue();
    //     cout << "f.LB: " << fLB_ << endl;

    // } else {

    //     cout << "Follower does not have finite LB" << endl;
    //     return false;
    // }

    setTimeLimit(original_time_limit);

    return true;
}
