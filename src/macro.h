/*
 * macro.h
 *
 * Created on: Nov 8, 2020
 *     Author: geunyeongbyeon
 */

#ifndef BILEVEL_MACRO_H_
#define BILEVEL_MACRO_H_

/* ##################################################
 * # functions for adding variables
 * ################################################## */

/* add leader variables */
#define addxVars() {    \
    int i;  \
    for (i = 0; i < (this)->n_l_; i++) {    \
        char cname[80]; \
        sprintf(cname, "x_%d", i); \
        if ((this)->is_integer_[i]) \
            (this)->vars_.x.add(IloNumVar(*((this)->env_), 0, 1, ILOBOOL, cname));  \
        else {  \
            /* if it's continuous, the original indices are researved */    \
            (this)->vars_.x.add(IloNumVar(*((this)->env_), (this)->xlb_[i], (this)->xub_[i], ILOFLOAT, cname));    \
        }   \
    }   \
}

/* add follower variables */
#define addyVars() {    \
    int i;  \
    for (i = 0; i < (this)->n_f_; i++) {    \
        char cname[80]; \
        sprintf(cname, "y_%d", i); \
        (this)->vars_.y.add(IloNumVar(*((this)->env_), -IloInfinity, IloInfinity, ILOFLOAT, cname));  \
    }   \
}

/* add unrestricted auxiliary var t (used in master) */
#define addtVar() {    \
    (this)->vars_.t = IloNumVar(*((this)->env_), -1e+19, IloInfinity, ILOFLOAT, "t");\
}

/* add nu variables */
#define addnuVars() {    \
    int i, j, k;  \
    for (i = 0; i < (this)->m_f_; i++) {    \
        for (j = 0; j < (this)->fC_lV_cnt_[i]; j++) {   \
            for (k = 0; k < 4; k++) {   \
                char cname[80]; \
                sprintf(cname, "nu_%d_%d_%d", i, (this)->fC_lV_ind_[i][j], k); \
                (this)->vars_.nu[k].add(IloNumVar(*((this)->env_), 0, IloInfinity, ILOFLOAT, cname));  \
            }   \
        }   \
    }   \
}

/* ##################################################
 * # functions for adding constraints
 * ################################################## */

/* add leader constraint: G_x x + G_y y >= h */
#define addlConstr() {  \
    IloExpr lhs(*((this)->env_)); \
    int i, j;   \
    for (i = 0; i < (this)->m_l_; i++) {    \
        char rname[80]; \
        sprintf(rname, "l_%d", i); \
        for (j = 0; j < (this)->lC_lV_cnt_[i]; j++) {   \
            lhs += (this)->lC_lV_coef_[i][j] * (this)->vars_.x[(this)->lC_lV_ind_[i][j]];   \
        }   \
        for (j = 0; j < (this)->lC_fV_cnt_[i]; j++) {   \
            lhs += (this)->lC_fV_coef_[i][j] * (this)->vars_.y[(this)->lC_fV_ind_[i][j]];   \
        }   \
        (this)->constrs_.l.add(IloRange(*((this)->env_), (this)->lC_rhs_[i], lhs, IloInfinity, rname));    \
        lhs.clear();    \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.l); \
}

/* add follower constraint: A x + B y >= b */
#define addfConstr() {  \
    IloExpr lhs(*((this)->env_));   \
    int i, j;   \
    for (i = 0; i < (this)->m_f_; i++) {    \
        char rname[80]; \
        sprintf(rname, "f_%d", i); \
        for (j = 0; j < (this)->fC_lV_cnt_[i]; j++) {   \
            lhs += (this)->fC_lV_coef_[i][j] * (this)->vars_.x[(this)->fC_lV_ind_[i][j]];   \
        }   \
        for (int j = 0; j < (this)->fC_fV_cnt_[i]; j++) {   \
            lhs += (this)->fC_fV_coef_[i][j] * (this)->vars_.y[(this)->fC_fV_ind_[i][j]];   \
        }   \
        (this)->constrs_.f.add(IloRange(*((this)->env_), (this)->fC_rhs_[i], lhs, IloInfinity, rname));    \
        lhs.clear();    \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.f); \
}

/* add x bounds constraint */
#define addxBdsConstr() {   \
    int i, j;   \
    IloExpr lhs(*((this)->env_));   \
    for (i = 0; i < (this)->ind_col_.size(); i++) { \
        if ((this)->ind_col_[i].size() > 1) {   \
            char rname[80]; \
            sprintf(rname, "xBds_%d", i); \
            for (j = 0; j < (this)->ind_col_[i].size(); j++) {  \
                lhs += (this)->scale_col_[i][j] * (this)->vars_.x[(this)->map_varind_to_lvarind_[(this)->ind_col_[i][j]]];  \
            }   \
            (this)->constrs_.xBds.add(IloRange(*((this)->env_), 0, lhs, (this)->xub_[i] - (this)->xlb_[i], rname));    \
            lhs.clear();    \
        }   \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.xBds);  \
}

/* add y bounds constraint */
#define addyBdsConstr() {   \
    int i;  \
    for (i = 0; i < (this)->ylb_cnt_; i++){ \
        char rname[80]; \
        sprintf(rname, "yLBs_%d", ylb_ind_[i]); \
        (this)->constrs_.yLBs.add(IloRange(*((this)->env_), (this)->ylb_coef_[i], (this)->vars_.y[(this)->ylb_ind_[i]], IloInfinity, rname));  \
    }   \
    for (i = 0; i < (this)->yub_cnt_; i++){ \
        char rname[80]; \
        sprintf(rname, "yUBs_%d", yub_ind_[i]); \
        (this)->constrs_.yUBs.add(IloRange(*((this)->env_), -(this)->yub_coef_[i], -(this)->vars_.y[(this)->yub_ind_[i]], IloInfinity, rname));    \
    }   \
    (this)->m_.add((this)->constrs_.yLBs);  \
    (this)->m_.add((this)->constrs_.yUBs);  \
}

/* add follower constraint with only y vars: B y >= b */
#define addfFConstr() { \
    IloExpr lhs(*((this)->env_));   \
    int i, j;   \
    for (i = 0; i < (this)->m_f_; i++) {    \
        char rname[80]; \
        sprintf(rname, "fF_%d", i); \
        for (j = 0; j < (this)->fC_fV_cnt_[i]; j++) {   \
            lhs += (this)->fC_fV_coef_[i][j] * (this)->vars_.y[(this)->fC_fV_ind_[i][j]];   \
        }   \
        (this)->constrs_.fF.add(IloRange(*((this)->env_), (this)->fC_rhs_[i], lhs, IloInfinity, rname));   \
        lhs.clear();    \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.fF);    \
}

/* add leader constraint with only y vars: G_y y >= h */
#define addlFConstr() {    \
    IloExpr lhs(*((this)->env_)); \
    int i, j;   \
    for (i = 0; i < (this)->m_l_; i++) {    \
        char rname[80]; \
        sprintf(rname, "lF_%d", i); \
        for (j = 0; j < (this)->lC_fV_cnt_[i]; j++) {   \
            lhs += (this)->lC_fV_coef_[i][j] * (this)->vars_.y[(this)->lC_fV_ind_[i][j]];   \
        }   \
        (this)->constrs_.lF.add(IloRange(*((this)->env_), (this)->lC_rhs_[i], lhs, IloInfinity, rname));   \
        lhs.clear();    \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.lF);    \
}

/* add constraint dy <= O */
#define addfObjConstr() {   \
    IloExpr lhs(*((this)->env_));   \
    int j;  \
    for (j = 0; j < (this)->n_f_; j++) {    \
        lhs += (this)->fObj_[j] * (this)->vars_.y[j];   \
    }   \
    (this)->constrs_.fObj = IloRange(*((this)->env_), 0, -lhs, IloInfinity, "fObj");    \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.fObj);  \
}

/* add follower constraint without x and with nu: By >= b + K^T nu */
#define addfFMCConstr() { \
    IloExpr lhs(*((this)->env_));   \
    int i, j;   \
    int pos = 0;    \
    for (i = 0; i < (this)->m_f_; i++) {    \
        char rname[80]; \
        sprintf(rname, "fFMC_%d", i); \
        for (j = 0; j < (this)->fC_lV_cnt_[i]; j++) {   \
            if ((this)->fC_lV_coef_[i][j] > 0){ \
                lhs += (this)->fC_lV_coef_[i][j] * ((this)->vars_.nu[1][pos] - (this)->vars_.nu[3][pos]);   \
            } else {    \
                lhs += -(this)->fC_lV_coef_[i][j] * ((this)->vars_.nu[1][pos] - (this)->vars_.nu[3][pos]);   \
            }   \
            pos++;  \
        }   \
        for (int j = 0; j < (this)->fC_fV_cnt_[i]; j++) {   \
            lhs += (this)->fC_fV_coef_[i][j] * (this)->vars_.y[(this)->fC_fV_ind_[i][j]];   \
        }   \
        (this)->constrs_.fFMC.add(IloRange(*((this)->env_), (this)->fC_rhs_[i], lhs, IloInfinity, rname));    \
        lhs.clear();    \
    }   \
    lhs.end();  \
    (this)->m_.add((this)->constrs_.fFMC); \
}

/* add nu constraints */
#define addnuConstr() { \
    int i, j;   \
    int pos = 0;    \
    for (i = 0; i < (this)->m_f_; i++) {    \
        for (j = 0; j < (this)->fC_lV_cnt_[i]; j++) {   \
            char rname[80]; \
            sprintf(rname, "nu_%d_%d", i, (this)->fC_lV_ind_[i][j]); \
            if ((this)->fC_lV_coef_[i][j] > 0){ \
                (this)->constrs_.nu.add(IloRange(*((this)->env_), 1, (this)->vars_.nu[0][pos] + (this)->vars_.nu[1][pos] - (this)->vars_.nu[2][pos] - (this)->vars_.nu[3][pos], 1, rname));    \
            } else {    \
                (this)->constrs_.nu.add(IloRange(*((this)->env_), 1, -(this)->vars_.nu[0][pos] -(this)->vars_.nu[1][pos] + (this)->vars_.nu[2][pos] + (this)->vars_.nu[3][pos], 1, rname));    \
            }   \
            pos++;  \
        }   \
    }   \
    (this)->m_.add((this)->constrs_.nu); \
}

/* add nu constraints for fUB:
 * nu_{ij,1} - nu_{ij,2} + nu_{ij,3} - nu_{ij,4} = a_{ij} */
#define addnuConstrfUB() { \
    int i, j;   \
    int pos = 0;    \
    for (i = 0; i < (this)->m_f_; i++) {    \
        for (j = 0; j < (this)->fC_lV_cnt_[i]; j++) {   \
            char rname[80]; \
            sprintf(rname, "nu_%d_%d", i, (this)->fC_lV_ind_[i][j]); \
            (this)->constrs_.nu.add(IloRange(*((this)->env_), (this)->fC_lV_coef_[i][j], (this)->vars_.nu[0][pos] - (this)->vars_.nu[1][pos] + (this)->vars_.nu[2][pos] - (this)->vars_.nu[3][pos], (this)->fC_lV_coef_[i][j], rname));    \
            pos++;  \
        }   \
    }   \
    (this)->m_.add((this)->constrs_.nu); \
}

#define freeArrayPtr(ptr)  \
	if (ptr) {              \
		delete [] ptr;      \
		ptr = NULL;         \
    }

#define free2DArrayPtr(n, ptr)         \
	if (ptr) {                             \
		for (int _i = 0; _i < n; ++_i) { \
			freeArrayPtr(ptr[_i])        \
		}                                  \
		delete [] ptr;                     \
		ptr = NULL;                        \
	}

#endif //BILEVEL_MACRO_H_