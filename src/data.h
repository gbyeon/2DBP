/* data.h
 * 
 * Created on: Nov 7, 2020
 *     Author: geunyeongbyeon
 */

#ifndef BILEVEL_DATA_H
#define BILEVEL_DATA_H

#include <map>
#include <vector>

#include <math.h>
#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <iostream>

#include <fstream>
#include <chrono>
#include <algorithm> // std::find

#include <macro.h>

class Data {
public:
    int read(string &filename);

public:
    string problemName_;

    /* number of HPP var and constr */
    int n_, m_; 

    /* number of follower var and constr */
    int n_f_, m_f_;

    /* number of leader var and constr */
    int n_l_, m_l_;

    /* objective coefficients */
    double *fObj_; /* d for dy */
    double *lObj_; /* c for cx + cy */
    double *llObj_; /* c for cx */
    double *lfObj_; /* c for cy */
    int lObj_constant_;

    /* integer var indicator for leader problem: 1 if integer */
    int * is_integer_;

    /* follower constraint data */
    double* fC_rhs_;
    int * fC_fV_cnt_;
    double ** fC_fV_coef_;
    int ** fC_fV_ind_;
    int * fC_lV_cnt_;
    double ** fC_lV_coef_;
    int ** fC_lV_ind_;

    /* leader constraint data */
    double* lC_rhs_;
    int * lC_fV_cnt_;
    double ** lC_fV_coef_;
    int ** lC_fV_ind_;
    int * lC_lV_cnt_;
    double ** lC_lV_coef_;
    int ** lC_lV_ind_;

    /* y bounds data */
    int ylb_cnt_;
    int yub_cnt_;
    double * ylb_coef_;
    int * ylb_ind_;
    double * yub_coef_;
    int * yub_ind_;

    /* x bounds data */
    vector<vector<int>> ind_col_;
    vector<vector<double>> scale_col_;
    vector<double> xlb_;
    vector<double> xub_;
    map<int,int> map_varind_to_lvarind_;

    Data() : n_(0), n_l_(0), n_f_(0), m_(0), m_l_(0), m_f_(0), ylb_cnt_(0), yub_cnt_(0), ylb_ind_(NULL), ylb_coef_(NULL), yub_ind_(NULL), yub_coef_(NULL), fObj_(NULL), lObj_(NULL), llObj_(NULL), lfObj_(NULL), is_integer_(NULL), fC_rhs_(NULL), fC_fV_cnt_(NULL), fC_fV_coef_(NULL), fC_fV_ind_(NULL), fC_lV_cnt_(NULL), fC_lV_coef_(NULL), fC_lV_ind_(NULL), lC_rhs_(NULL), lC_fV_cnt_(NULL), lC_fV_coef_(NULL), lC_fV_ind_(NULL), lC_lV_cnt_(NULL), lC_lV_coef_(NULL), lC_lV_ind_(NULL) {};

    ~Data(){     
        freeArrayPtr(fC_rhs_);
        free2DArrayPtr(m_f_, fC_fV_coef_);
        free2DArrayPtr(m_f_, fC_fV_ind_);
        free2DArrayPtr(m_f_, fC_lV_coef_);
        free2DArrayPtr(m_f_, fC_lV_ind_);
        freeArrayPtr(fC_fV_cnt_);
        freeArrayPtr(fC_lV_cnt_);
           
        freeArrayPtr(lC_rhs_);
        free2DArrayPtr(m_l_, lC_fV_coef_);
        free2DArrayPtr(m_l_, lC_fV_ind_);
        free2DArrayPtr(m_l_, lC_lV_coef_);
        free2DArrayPtr(m_l_, lC_lV_ind_);
        freeArrayPtr(lC_fV_cnt_);
        freeArrayPtr(lC_lV_cnt_);

        freeArrayPtr(ylb_coef_);
        freeArrayPtr(yub_coef_);
        freeArrayPtr(ylb_ind_);
        freeArrayPtr(yub_ind_);
        freeArrayPtr(fObj_);
        freeArrayPtr(lObj_);
        freeArrayPtr(lfObj_);
        freeArrayPtr(llObj_);
        freeArrayPtr(is_integer_);

    };
};



    
#endif //BILEVEL_DATA_H
