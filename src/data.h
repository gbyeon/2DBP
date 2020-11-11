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

    Data(){
        ylb_cnt_ = 0;
        yub_cnt_ = 0;
    };

    ~Data(){     
        if (m_f_ > 0) 
        {
            delete[] fC_rhs_;

            for (int i = 0; i < m_f_; i++) 
            {    
                if (fC_fV_cnt_[i] > 0) {
                    delete[] fC_fV_coef_[i];
                    delete[] fC_fV_ind_[i];
                }
                if (fC_lV_cnt_[i] > 0) {
                    delete[] fC_lV_coef_[i];
                    delete[] fC_lV_ind_[i];
                }
            }

            delete[] fC_fV_coef_;
            delete[] fC_fV_ind_;
            delete[] fC_lV_coef_;
            delete[] fC_lV_ind_;

            delete[] fC_fV_cnt_;
            delete[] fC_lV_cnt_;
        }

        if (m_l_ > 0) {
            delete[] lC_rhs_;

            for (int i = 0; i < m_l_; i++) 
            {    
                if (lC_fV_cnt_[i] > 0) {
                    delete[] lC_fV_coef_[i];
                    delete[] lC_fV_ind_[i];
                }
                if (lC_lV_cnt_[i] > 0) {
                    delete[] lC_lV_coef_[i];
                    delete[] lC_lV_ind_[i];
                }
            }
            delete[] lC_fV_coef_;
            delete[] lC_fV_ind_;
            delete[] lC_lV_coef_;
            delete[] lC_lV_ind_;

            delete[] lC_fV_cnt_;
            delete[] lC_lV_cnt_;
        }
    
        if (ylb_cnt_ > 0) {
            delete[] ylb_coef_;
            delete[] yub_coef_;
        }
        if (yub_cnt_ > 0) {
            delete[] ylb_ind_;
            delete[] yub_ind_;
        }
        if (n_f_ > 0) {
            delete[] fObj_;
            delete[] lfObj_;
        }
        if (n_ > 0)
            delete[] lObj_;
        if (n_l_ > 0) {
            delete[] is_integer_;
            delete[] llObj_;
        }
    };
};



    
#endif //BILEVEL_DATA_H
