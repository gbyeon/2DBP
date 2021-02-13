/* cpp
 * 
 * Created on: Nov 7, 2020
 *     Author: geunyeongbyeon
 */

#include "data.h"
#include "CoinMpsIO.hpp"
// #define DATA_DEBUG

void resize(int * &arr, int size, int newSize) {
    int * newArr = new int [newSize];

    if (size > 0){
        memcpy(newArr, arr, size * sizeof(int));
    }

    // size = newSize;
    freeArrayPtr(arr);

    arr = newArr; 
}
void resize(double * &arr, int size, int newSize) {
    double * newArr = new double [newSize];
    
    if (size > 0) {
        memcpy(newArr, arr, size * sizeof(double));
    }

    // size = newSize;
    freeArrayPtr(arr);
    arr = newArr; 
}
int Data::read(string &filename) {
    
    int i, j, k;

    CoinMpsIO mps;
    mps.readMps(filename.c_str());
    /* number of leader and follower variables */
    n_ = mps.getNumCols(); 
    m_ = mps.getNumRows(); 
    cout << "n: " << n_ << endl;
    cout << "m: " << m_ << endl;

    /* readauxfile */
    vector<bool> is_follower_var (n_, false);
    vector<bool> is_follower_constr (m_, false);
    
    int obj_cnt = 0;
    int fObjSense;
    int objSense = -1;
    
    /* leader obj sense is assumed to be the same as the fObjSense */
    int scale_lobj;

    filename.append(".aux");
    ifstream auxInput (filename.c_str(), ios::in | ios:: binary);
    // ifstream auxInput ("../instances/toy/toy.aux", ios::in | ios:: binary);

    string item;

    while (auxInput >> item) {
        if (item == "N")
        {
            auxInput >> n_f_;
            if (n_f_ > 0)
                fObj_ = new double [n_f_];
        }
        else if (item == "M")
            auxInput >> m_f_;
        else if (item == "LC") 
        {
            auxInput >> i;
            is_follower_var[i] = true;

        } else if (item == "LR") 
        {
            auxInput >> i;
            is_follower_constr[i] = true;

        } else if (item == "LO") 
        {
            auxInput >> fObj_[obj_cnt];
            obj_cnt++;
        } 
        else if (item == "OS") {
            auxInput >> fObjSense;
        } else {
            cout << "Incorrect row in aux file: " << item << endl;
            return 1;
        }
    }
    auxInput.close();

    assert(obj_cnt == n_f_);

    cout << "filename: " << mps.getFileName() << endl;
    problemName_.assign(mps.getProblemName());
    cout << "problem name: " << problemName_.c_str() << endl;
    // cout << "getObjectiveName: " << mps.getObjectiveName() << endl;

    m_l_ = m_ - m_f_;
    n_l_ = n_ - n_f_;
    cout << "m_l: " << m_l_ << ", m_f: " << m_f_ << endl;
    cout << "n_l: " << n_l_ << ", n_f: " << n_f_ << endl;

    /* if there is some leader variable that is integer and not binary, calculate the number of resultant binary vars
    * if there is a contraint that is equal to sense, calculate the number of resultant less than or equal to constrs 
    */ 
    vector<vector<int>> ind_col(n_);
    vector<vector<double>> scale_col(n_);
    xlb_.resize(n_);
    xub_.resize(n_);
    for (j = 0; j < n_; j++) {    

        xlb_[j] = mps.getColLower()[j];
        xub_[j] = mps.getColUpper()[j];
        // cout << "xlb_:" << xlb_[j] << ", xub_:" << xub_[j] << endl;

        if (!is_follower_var[j] && mps.isInteger(j)) {
            if (mps.getColLower()[j] != 0 || mps.getColUpper()[j] != 1) 
            {
                int r = floor(log2(mps.getColUpper()[j] - mps.getColLower()[j])) + 1;
                ind_col[j].resize(r);
                scale_col[j].resize(r);
                
                for (k = 0; k < r; k++) {
                    if (k == 0) {
                        ind_col[j][k] = j;
                        scale_col[j][k] = 1;
                    } else {
                        
                        ind_col[j][k] = is_follower_var.size(); // last index of current variables
                        scale_col[j][k] = pow(2,k);
                        
                        is_follower_var.push_back(false);
                        n_l_++;
                    }
                }
            } else {
                ind_col[j].push_back(j);
                scale_col[j].push_back(1);
            }
        } else {
            ind_col[j].push_back(j);
            scale_col[j].push_back(1);
        }
    }

    n_ = n_l_ + n_f_;
    cout << "new n: " << n_ << endl;
    assert(n_ == is_follower_var.size());

#ifdef DATA_DEBUG
    for (j = 0; j < ind_col.size(); j++) {    
        cout << j << "th var: " << endl; 
        for (k = 0; k < ind_col[j].size(); k++) {
            cout << "        " << k << "th subvar index " << ind_col[j][k];
            cout << ", scale " << scale_col[j][k] << endl;
        }
    }
#endif
    vector<vector<int>> ind_row(m_);
    vector<vector<double>> scale_row(m_);
    for (i = 0; i < m_; i++) {
        
        if (mps.getRowSense()[i] == 'G') {
            ind_row[i].push_back(i);
            scale_row[i].push_back(1);
        }
        else if (mps.getRowSense()[i] == 'E') {
            
            ind_row[i].resize(2);
            scale_row[i].resize(2);

            ind_row[i][0] = i;
            scale_row[i][0] = 1;

            ind_row[i][1] = is_follower_constr.size(); // last index of current constraints
            scale_row[i][1] = -1;
            
            if (is_follower_constr[i]) {
                m_f_++;
                is_follower_constr.push_back(true);
            } else {
                m_l_++;
                is_follower_constr.push_back(false);
            }
        } else if (mps.getRowSense()[i] == 'L') {
            ind_row[i].push_back(i);
            scale_row[i].push_back(-1);
        } else {
            throw "invalid row sense";//, "main", "bilevel.cpp";
        }
    }
#ifdef DATA_DEBUG
    for (j = 0; j < ind_row.size(); j++) {    
        cout << j << "th row. sense " << mps.getRowSense()[j] << endl; 
        for (k = 0; k < ind_row[j].size(); k++) {
            cout << "        " << k << "th subrow index " << ind_row[j][k];
            cout << ", scale " << scale_row[j][k] << endl;
        }
    }
#endif
    cout << "new m: " << m_ << endl;

    m_ = m_l_ + m_f_;
    assert(m_ == is_follower_constr.size());
    
    map<int,int> map_varind_to_fvarind;
    map<int,int> map_varind_to_lvarind;
    int fpos = 0;
    int lpos = 0;
    for (j = 0; j < n_; j++) {
        if (is_follower_var[j])
            map_varind_to_fvarind.insert(pair<int,int>(j,fpos++));
        else 
            map_varind_to_lvarind.insert(pair<int,int>(j,lpos++));    
    }
#ifdef DATA_DEBUG
    /* test map_varind_to_f or lvarind */
    for (j = 0; j < n_; j++) {
        if (is_follower_var[j])
            cout << j  << "th var is " << map_varind_to_fvarind[j] << "th fvar" << endl;
        else 
        cout << j  << "th var is " << map_varind_to_lvarind[j] << "th lvar" << endl;
    }
#endif
    
    assert(n_f_ == map_varind_to_fvarind.size());
    assert(n_l_ == map_varind_to_lvarind.size());

    map<int,int> map_cind_to_fcind;
    map<int,int> map_cind_to_lcind;
    fpos = 0;
    lpos = 0;
    for (j = 0; j < m_; j++) {
        if (is_follower_constr[j])
            map_cind_to_fcind.insert(pair<int,int>(j,fpos++));
        else 
            map_cind_to_lcind.insert(pair<int,int>(j,lpos++));    
    }

#ifdef DATA_DEBUG
    /* test map_cind_to_f or lcind */
    for (j = 0; j < m_; j++) {
        if (is_follower_constr[j])
            cout << j  << "th constr is " << map_cind_to_fcind[j] << "th fcon" << endl;
        else 
        cout << j  << "th constr is " << map_cind_to_lcind[j] << "th lcon" << endl;
    }
#endif
    /* scale follower objective */
    /* change to minimization problem */
    if (fObjSense == -1) {
        for (i = 0; i < n_f_; i++) {
            fObj_[i] = -fObj_[i];
        }
        fObjSense = 1;
        // scale_lobj = 1;   
        scale_lobj = -1;   
    } else scale_lobj = 1; 
#ifdef DATA_DEBUG
    for (i = 0; i < n_f_; i++) {
        cout << "fobj of " << i << " th var: " << fObj_[i] << endl;
    }
#endif

    assert(m_f_ == map_cind_to_fcind.size());
    assert(m_l_ == map_cind_to_lcind.size());

    /* get leader objective */
    /* change to minimization problem */
    if (n_ > 0)
        lObj_ = new double [n_];
    if (n_l_ > 0)
        llObj_ = new double [n_l_];
    if (n_f_ > 0)
        lfObj_ = new double [n_f_];
    lObj_constant_ = 0;
    for (i = 0; i < ind_col.size(); i++) {
        for (k = 0; k < ind_col[i].size(); k++) {
            lObj_[ind_col[i][k]] = pow(2,k) * mps.getObjCoefficients()[i] * scale_lobj;
            if (!is_follower_var[i]) {
                llObj_[map_varind_to_lvarind[ind_col[i][k]]] = lObj_[ind_col[i][k]];
            }
        }
        if (ind_col[i].size() > 1) {
            lObj_constant_ += mps.getColLower()[i] * mps.getObjCoefficients()[i] * scale_lobj;
        }
        if (is_follower_var[i]) {
            lfObj_[map_varind_to_fvarind[i]] = lObj_[i];
        }
    }
#ifdef DATA_DEBUG
    for (i = 0; i < n_; i++) {
        cout << "lobj of " << i << " th var: " << lObj_[i] << endl;
    }
    cout << "lobj const: " << lObj_constant_ << endl;
    for (i = 0; i < n_f_; i++) {
        cout << "lfobj of " << i << " th fvar: " << lfObj_[i] << endl;
    }
#endif

    /* get integer var indicator of leader */
    if (n_l_ > 0)
        is_integer_ = new int [n_l_];
    int pos = 0;
    for (i = 0; i < ind_col.size(); i++) {
        if (!is_follower_var[i]) {
            for (k = 0; k < ind_col[i].size(); k++) {
                is_integer_[pos] = mps.isInteger(i);
                pos++;
            }
        }
    }
#ifdef DATA_DEBUG
    for (i = 0; i < n_l_; i++) {
        cout << "is_integer? " << i << " th lvar: " << is_integer_[i] << endl;
    }
#endif
    /* follower constraint data */
    if (m_f_ > 0) {
        fC_rhs_ = new double [m_f_];
        fC_fV_cnt_ = new int [m_f_];
        fC_fV_coef_ = new double * [m_f_];
        fC_fV_ind_ = new int * [m_f_];
        fC_lV_cnt_ = new int [m_f_];
        fC_lV_coef_ = new double * [m_f_];
        fC_lV_ind_ = new int * [m_f_];
    }
    /* leader constraint data */
    if (m_l_ > 0) {
        lC_rhs_ = new double [m_l_];
        lC_fV_cnt_ = new int [m_l_];
        lC_fV_coef_ = new double * [m_l_];
        lC_fV_ind_ = new int * [m_l_];
        lC_lV_cnt_ = new int [m_l_];
        lC_lV_coef_ = new double * [m_l_];
        lC_lV_ind_ = new int * [m_l_];
    }
    
    /* get default rhs */
    for (j = 0; j < ind_row.size(); j++) {
        
        if (!is_follower_constr[j]) {
            int row_index = map_cind_to_lcind[j];
            if (mps.getRowSense()[j] == 'L') {
                lC_rhs_[row_index] = -mps.getRowUpper()[j];
            } else if (mps.getRowSense()[j] == 'E') {
                lC_rhs_[row_index] = mps.getRowUpper()[j];
                lC_rhs_[map_cind_to_lcind[ind_row[j][1]]] = -mps.getRowUpper()[j];
            } else if (mps.getRowSense()[j] == 'G') {
                lC_rhs_[row_index] = mps.getRowLower()[j];
            }
        } else {
            int row_index = map_cind_to_fcind[j];
            if (mps.getRowSense()[j] == 'L') {
                fC_rhs_[row_index] = -mps.getRowUpper()[j];
            } else if (mps.getRowSense()[j] == 'E') {
                fC_rhs_[row_index] = mps.getRowUpper()[j];
                fC_rhs_[map_cind_to_fcind[ind_row[j][1]]] = -mps.getRowUpper()[j];
            } else if (mps.getRowSense()[j] == 'G') {
                fC_rhs_[row_index] = mps.getRowLower()[j];
            }
        }
    }

    /* get coeff matrix */
    // cout << "ind_row size: " << ind_row.size() << endl;
    assert(mps.getMatrixByRow()->getMajorDim() == ind_row.size());

    for (j = 0; j < ind_row.size(); j++) {
        
        // cout << j << " th row: " << endl;
        const CoinShallowPackedVector row = mps.getMatrixByRow()->getVector(j);

        for (k = 0; k < ind_row[j].size(); k++) {
            
            if (!is_follower_constr[ind_row[j][k]]) {
                // cout << "leader constr. " << endl;
                int row_index = map_cind_to_lcind[ind_row[j][k]];
                
                lC_fV_cnt_[row_index] = 0;
                lC_lV_cnt_[row_index] = 0;

                // cout << "numElement: " << row.getNumElements() << endl;
                for (i = 0; i < row.getNumElements(); i++) {
                    int var_org_index = row.getIndices()[i];

                    for (int l = 0; l < ind_col[var_org_index].size(); l++) {
                        int var_index = ind_col[var_org_index][l];
                        if (is_follower_var[var_index]) {
                            lC_fV_cnt_[row_index]++;
                            // cout << "fv cnt: " << lC_fV_cnt[row_index] << endl;
                        } else {
                            lC_lV_cnt_[row_index]++;
                            // cout << "lv cnt: " << lC_lV_cnt[row_index] << endl;
                        }
                    }
                }  
                if (lC_fV_cnt_[row_index] > 0) {
                    lC_fV_ind_[row_index] = new int [lC_fV_cnt_[row_index]];
                    lC_fV_coef_[row_index] = new double [lC_fV_cnt_[row_index]]; 
                }
                
                if (lC_lV_cnt_[row_index] > 0) {
                    lC_lV_ind_[row_index] = new int [lC_lV_cnt_[row_index]];
                    lC_lV_coef_[row_index] = new double [lC_lV_cnt_[row_index]];
                }
                fpos = 0;
                lpos = 0;
                for (i = 0; i < row.getNumElements(); i++) {
                    int var_org_index = row.getIndices()[i];

                    for (int l = 0; l < ind_col[var_org_index].size(); l++) {
                        
                        int var_index = ind_col[var_org_index][l];
                        
                        if (is_follower_var[var_index]) {
                            lC_fV_ind_[row_index][fpos] = map_varind_to_fvarind[var_index];
                            lC_fV_coef_[row_index][fpos] = row.getElements()[i] * scale_col[var_org_index][l] * scale_row[j][k];
                            fpos++;
                        } else {
                            lC_lV_ind_[row_index][lpos] = map_varind_to_lvarind[var_index];
                            lC_lV_coef_[row_index][lpos] = row.getElements()[i] * scale_col[var_org_index][l] * scale_row[j][k];
                            lpos++;
                        }
                    }

                    if (ind_col[var_org_index].size() > 1) {
                        lC_rhs_[row_index] -= mps.getColLower()[var_org_index] * row.getElements()[i] * scale_row[j][k];
                    }
                }
            } 
            else 
            {
                // cout << "follower constr. " << endl;
                // cout << "ind_row: " << ind_row[j][k] << endl;
                int row_index = map_cind_to_fcind[ind_row[j][k]];
                // cout << "row_index: " << row_index << endl;
                
                fC_fV_cnt_[row_index] = 0;
                fC_lV_cnt_[row_index] = 0;

                // cout << "numElement: " << row.getNumElements() << endl;
                for (i = 0; i < row.getNumElements(); i++) {
                    
                    int var_org_index = row.getIndices()[i];
                    // cout << "var_org_index: " << var_org_index << endl;

                    for (int l = 0; l < ind_col[var_org_index].size(); l++) {
                        
                        int var_index = ind_col[var_org_index][l];
                        // cout << "var_index: " << var_index << endl;
                        
                        if (is_follower_var[var_index]) {
                            fC_fV_cnt_[row_index]++;
                            // cout << "fv cnt: " << fC_fV_cnt_[row_index] << endl;
                        } else {
                            fC_lV_cnt_[row_index]++;
                            // cout << "lv cnt: " << fC_lV_cnt_[row_index] << endl;
                        }
                    }
                }
                if (fC_fV_cnt_[row_index] > 0)
                {
                    fC_fV_ind_[row_index] = new int [fC_fV_cnt_[row_index]];
                    fC_fV_coef_[row_index] = new double [fC_fV_cnt_[row_index]];
                }
                
                if (fC_lV_cnt_[row_index] > 0)
                {
                    fC_lV_ind_[row_index] = new int [fC_lV_cnt_[row_index]];
                    fC_lV_coef_[row_index] = new double [fC_lV_cnt_[row_index]];
                }  
                fpos = 0;
                lpos = 0;
                for (i = 0; i < row.getNumElements(); i++) {
                    int var_org_index = row.getIndices()[i];
                    // cout << "var_org_index: " << var_org_index << endl;

                    for (int l = 0; l < ind_col[var_org_index].size(); l++) {
                        
                        int var_index = ind_col[var_org_index][l];
                        // cout << "var_index: " << var_index << endl;
                        
                        if (is_follower_var[var_index]) {
                            fC_fV_ind_[row_index][fpos] = map_varind_to_fvarind[var_index];
                            // cout << "map_varind_to_fvarind: " << map_varind_to_fvarind[var_index] << endl;
                            fC_fV_coef_[row_index][fpos] = row.getElements()[i] * scale_col[var_org_index][l] * scale_row[j][k];
                            // cout << row.getElements()[i] * scale_col[var_org_index][l] * scale_row[j][k] << endl;
                            // cout << "fc_fv_coef of row " << row_index << " and var " << fC_fV_ind[row_index][fpos] << ": " << fC_fV_coef[row_index][fpos] << endl;
                            fpos++;
                        } else {
                            fC_lV_ind_[row_index][lpos] = map_varind_to_lvarind[var_index];
                            // cout << "map_varind_to_lvarind: " << map_varind_to_lvarind[var_index] << endl;
                            fC_lV_coef_[row_index][lpos] = row.getElements()[i] * scale_col[var_org_index][l] * scale_row[j][k];
                            // cout << "org coeff: " << row.getElements()[i] << ", scale_col: " << scale_col[var_org_index][l] << ", scale_row: " <<  scale_row[j][k] << endl;
                            // cout << "fc_lv_coef of row " << row_index << " and var " << fC_lV_ind[row_index][lpos] << ": " << fC_lV_coef[row_index][lpos] << endl;
                        
                            lpos++;
                        }
                    }

                    if (ind_col[var_org_index].size() > 1) {
                        // cout << row_index << "th fC_rhc: " << fC_rhs[row_index] << endl;
                        fC_rhs_[row_index] -= mps.getColLower()[var_org_index] * row.getElements()[i] * scale_row[j][k];
                        // cout << row_index << "th fC_rhc after : " << fC_rhs[row_index] << endl;
                    }
                }
            }
        }
    }

    for (i = 0; i < n_; i++) {
        if (is_follower_var[i]){
            if (mps.getColLower()[i] != -mps.getInfinity()) 
                ylb_cnt_++;
            if (mps.getColUpper()[i] != mps.getInfinity())
                yub_cnt_++;
        }
    }
    
    // cout << "ylb_cnt: " << ylb_cnt << ", yub_cnt: " << yub_cnt << endl;
    if (ylb_cnt_ > 0) {
        ylb_coef_ = new double [ylb_cnt_];
        ylb_ind_ = new int [ylb_cnt_];
    }
    if (yub_cnt_ > 0) {
        yub_coef_ = new double [yub_cnt_];
        yub_ind_ = new int [yub_cnt_];
    }
    int lbpos = 0;
    int ubpos = 0;
    for (i = 0; i < n_; i++) {
        if (is_follower_var[i]){
            if (mps.getColLower()[i] > -mps.getInfinity()) {
                ylb_ind_[lbpos] = map_varind_to_fvarind[i];
                ylb_coef_[lbpos] = mps.getColLower()[i];
                // cout << map_varind_to_fvarind[i] << "th fvar has lb of " << ylb_coef[lbpos] << endl;
                lbpos++;
            } if (mps.getColUpper()[i] < mps.getInfinity()) {
                yub_ind_[ubpos] = map_varind_to_fvarind[i];
                yub_coef_[ubpos] = mps.getColUpper()[i];
                ubpos++;
            }
        }
    }
#ifdef DATA_DEBUG
    cout << "is_follower_var" << endl;
    for (i = 0; i < n_; i++) {
        cout << i << ", " << is_follower_var[i] << endl;
    }
    cout << "is_follower_constr" << endl;
    for (i = 0; i < m_; i++) {
        cout << i << ", " << is_follower_constr[i] << endl;
    }
    for (i = 0; i < m_f_; i++) {
        for (j = 0; j < fC_fV_cnt_[i]; j++) {
            cout << fC_fV_ind_[i][j] << ", " << fC_fV_coef_[i][j] << ", ";
        }
        for (j = 0; j < fC_lV_cnt_[i]; j++) {
            cout << fC_lV_ind_[i][j] << ", " << fC_lV_coef_[i][j] << ", ";
        }
        cout << endl;
    }
    cout << "rhs: ";
    for (i = 0; i < m_f_; i++) {
        cout << fC_rhs_[i] << ", ";
    }
    cout << endl;
#endif

    ind_col_ = ind_col;
    scale_col_ = scale_col;
    map_varind_to_lvarind_ = map_varind_to_lvarind;

    vector<vector<double>> y_coeffs(n_f_);
    for (i = 0; i < m_f_; i++) {
        for (j = 0; j < fC_fV_cnt_[i]; j++) {
            y_coeffs[fC_fV_ind_[i][j]].push_back(fC_fV_coef_[i][j]);
        }
    }

    /* preprocess follower variable fixing */
    for (i = 0; i < n_; i++) {
        if (is_follower_var[i]){
            int f_ind = map_varind_to_fvarind[i];
            if (fObj_[f_ind] * fObjSense < -1e-4) {
                /* the larger the better: max with positive obj coef */
                bool can_be_fixed = true;

                for (auto &c : y_coeffs[f_ind]) {
                    if (c < 0) {
                        can_be_fixed = false;
                        break;
                    }
                }
            
                if (can_be_fixed) {
                    bool found_index = false;
                    if (mps.getColUpper()[i] >= mps.getInfinity()) {
                        cout << "unbounded problem " << endl;
                        return 1;
                    }
                    for (j = 0; j < ylb_cnt_; j++) {
                        if (ylb_ind_[j] == f_ind) {
                            ylb_coef_[j] = mps.getColUpper()[i];
                            found_index = true;
                            break;
                        }
                    }
                    if (!found_index){
                        resize(ylb_ind_, ylb_cnt_, ylb_cnt_+1);
                        resize(ylb_coef_, ylb_cnt_, ylb_cnt_+1);
                        // ylb_ind_ = (int *) realloc(ylb_ind_, sizeof(int) * (ylb_cnt_+1));
                        // ylb_coef_ = (double *) realloc(ylb_coef_, sizeof(double) * (ylb_cnt_+1));
                        
                       ylb_ind_[ylb_cnt_] = f_ind;
                       ylb_coef_[ylb_cnt_] = mps.getColUpper()[i];
                       ylb_cnt_++;
                    }
                }
                
            } else if (fObj_[f_ind] * fObjSense > 1e-4) {
                /* the smaller the better: max with negative obj coef */
                bool can_be_fixed = true;

                for (auto &c : y_coeffs[f_ind]) {
                    if (c > 0) {
                        can_be_fixed = false;
                        break;
                    }
                }
                
                if (can_be_fixed) {
                    bool found_index = false;
                    if (mps.getColLower()[i] <= -mps.getInfinity()) {
                        cout << "unbounded problem " << endl;
                        return 1;
                    }
                    for (j = 0; j < yub_cnt_; j++) {
                        if (yub_ind_[j] == f_ind) {
                            yub_coef_[j] = mps.getColLower()[i];
                            found_index = true;
                            break;
                        }
                    }
                    if (!found_index){
                        resize(yub_ind_, yub_cnt_, yub_cnt_+1);
                        resize(yub_coef_, yub_cnt_, yub_cnt_+1);

                        // yub_ind_ = (int *) realloc(yub_ind_, sizeof(int) * (yub_cnt_+1));
                        // yub_coef_ = (double *) realloc(yub_coef_, sizeof(double) * (yub_cnt_+1));

                       yub_ind_[yub_cnt_] = f_ind;
                       yub_coef_[yub_cnt_] = mps.getColLower()[i];
                       yub_cnt_++;
                    }
                }
            }
        }
    }

    return 0;
}
