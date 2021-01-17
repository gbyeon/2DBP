/*
 * followerX.cpp
 *
 *  Created on: Nov 13, 2020
 *      Author: geunyeongbyeon
 */
#undef NDEBUG
// #define F_SOLVE_DEBUG

#include "followerx.h"
#include "macro.h"

FollowerX::FollowerX()
{
    env_ = new IloEnv;
    m_ = IloModel(*env_);
    vars_ = fXVars(env_);
    constrs_ = fXConstrs(env_);
    yVals_ = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);
}

FollowerX::FollowerX(int n_l,
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
    vars_ = fXVars(env_);
    constrs_ = fXConstrs(env_);
    yVals_ = IloNumArray (*env_);
    dy_ = IloNumArray (*env_);
    dy_expr_ = IloExpr (*env_);
}

/* copy constructor */
FollowerX::FollowerX(const FollowerX & rhs) :
env_(rhs.env_),
m_(rhs.m_),
vars_(rhs.vars_),
constrs_(rhs.constrs_),
yVals_(rhs.yVals_),
dy_(rhs.dy_),
dy_expr_(rhs.dy_expr_)
{
    /* nothing to do */
}

FollowerX::~FollowerX()
{
    yVals_.end();
    dy_.end();
    dy_expr_.end();
    
    vars_.y.end();
    vars_.x.end();
    constrs_.f.end();
    constrs_.yLBs.end();
    constrs_.yUBs.end();
    constrs_.xBds.end();

    cplex_.end();
    m_.end();
    env_->end();

    delete env_;
}

void FollowerX::loadProblem (Data &data) {

    n_l_ = data.n_l_;
    n_f_ = data.n_f_;
    m_f_ = data.m_f_;

    xlb_ = data.xlb_;
    xub_ = data.xub_;
    ind_col_ = data.ind_col_;
    scale_col_ = data.scale_col_;
    map_varind_to_lvarind_ = data.map_varind_to_lvarind_;

    is_integer_ = data.is_integer_;
    
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

   
}

void FollowerX::createProblem () {
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
    
    /* add leader variable y */
    addxVars();

    /* add follower variable y */
    addyVars();
    
    /* add follower constraint without x: By >= b */
    addfConstr();
    
    /* add bound constraints on y */
    addyBdsConstr();

    /* add bound constraints on x */
    addxBdsConstr();
    
    /* initialize cplex */
    cplex_ = IloCplex(m_);

    /* set objective function: minimize dy */
    int i;
    for (i = 0; i < n_f_; i++)
        dy_.add(fObj_[i]);
    dy_expr_.setLinearCoefs(vars_.y, dy_);
    m_.add(IloMaximize(*env_, dy_expr_));

    /* turn off display */
    cplex_.setOut(env_->getNullStream());

//    cplex_.setParam(IloCplex::Param::TimeLimit, 10);
    cplex_.exportModel("followerx.lp");

//    // dummy r
//    vars_.r = IloNumVar (env_->, 0, 0, ILOFLOAT, "r");
//    m_.add(vars_.r);
//    IloInt i;
//    for (i = 0; i < constrs_.fF.getSize(); i++)
//        fF[i].setLinearCoef(vars_.r, 1);

}

void FollowerX::updateUBProblem (IloNumArray &xUBs, IloNumArray &xLBs) {
#ifdef F_SOLVE_DEBUG
cout << "solve followerX UB: ";
#endif 
    vars_.x.setBounds(xLBs, xUBs);
    cplex_.exportModel("followerx_updated.lp");
}

int FollowerX::solve () {
#ifdef F_SOLVE_DEBUG
    /* tic */
    auto start_t = chrono::system_clock::now();
#endif
    /* solve */
    if (!cplex_.solve()) {
        // env_->error() << "Failed to optimize followerX problem." << endl;
        // cout << "status: " << status_ << endl;
        // return 1;
    }
    
    status_ = cplex_.getStatus();
    getResults();

#ifdef F_SOLVE_DEBUG
    ticToc_ = (chrono::system_clock::now() - start_t);
    cout << "FollowerX has optimum: " << objVal_ << ", time: " << ticToc_.count() << endl;
#endif

    return 0;
}

void FollowerX::getResults () {

    objVal_ = cplex_.getObjValue();
}


double FollowerX::getDyVal(IloNumArray &yVals) {
    int i; 
    double Dy = 0; 
    for (i = 0; i < n_f_; i++)
    {
        Dy += fObj_[i] * yVals[i];
    }

    return Dy;
}