/* bilevel.cpp
 * 
 * Created on: Nov 7, 2020
 *     Author: geunyeongbyeon
 */
#include <iostream>
#undef NDEBUG
#include <assert.h>

#include "master.h"

int getArgValues (int argc, char *argv[], string &filepath, double * &testValue, bool &s, double &t, bool &h, double &ht, double &u){

    for (int i = 0; i < floor((double) argc / (double) 2); i++) {

        int j = 2 * i + 1;
        int v = 2 * (i + 1);

        if (strcmp(argv[j], "-f") == 0) {
            /* instance path/name */
            filepath = string(argv[v]);
        } 
        else if (strcmp(argv[j], "-test") == 0) {
            /* test value */
            testValue = new double;
            *testValue = stod(argv[v]);
            cout << "testing value: " << *testValue << endl;
        }
        else if (strcmp(argv[j], "-s") == 0) {
            /* use_numerically_stable_cut_cb_ */
            s = stoi(argv[v]);
            cout << "use_numerically_stable_cut_cb_: " << s << endl;
        } 
        else if (strcmp(argv[j], "-t") == 0) {
            /* time_limit_ of master */
            t = stod(argv[v]);
            cout << "time limit: " << t << endl;
        } 
        else if (strcmp(argv[j], "-h") == 0) {
            /* use_heuristic_cb_ */
            h = stoi(argv[v]);
            cout << "use_heuristic_cb_: " << h << endl;
        }
        else if (strcmp(argv[j], "-ht") == 0) {
            /* heuristic_time_limit_ */
            ht = stod(argv[v]);
            cout << "heuristic_time_limit_: " << ht << endl;
        }
        else if (strcmp(argv[j], "-u") == 0) {
            /* test value */
            u = stod(argv[v]);
            cout << "dual variable UB: " << u << endl;
        }
    } 

   // check
    cout << "This program was called with " << argv[0] << endl;
    cout << "filepath: " << filepath << endl;
 
    return 0;
}

int main (int argc, char *argv[]) {

    string filepath;
    double * testValue = NULL;
    bool s = 1;
    double t = 3600;
    bool h = 0;
    double ht = 150;
    double u = 100000;
    if (getArgValues(argc, argv, filepath, testValue, s, t, h, ht, u))
        return 1;

    Data data;
    if (data.read(filepath)) {
        cout << "fail to read data" << endl;
        return 1;
    }

    /* create follower problem:
     * min dy : By >= b (-Ax)   (fF)
     *          y bounds        (yLbs, yUbs)
     */
    Follower follower;
    follower.loadProblem(data);
    follower.createProblem();
    
    /* create follower problem with x:
     * min dy : Ax + By >= b      (f)
     *          x, y bounds       (xBds, yLbs, yUbs)
     */
    // FollowerX followerx;
    // followerx.loadProblem(data);
    // followerx.createProblem();
    
    /* create leaderFollower problem:
     * min c_y y : G_y y >= h (-G_x x)  (lF)
     *             B y >= b (-Ax)       (fF)
     *             dy <= O 
     *             y bounds             (yLbs, yUbs)
     */
    LeaderFollower leaderFollower;
    leaderFollower.loadProblem(data);
    leaderFollower.createProblem();

    /* create follower McCormick problem:
     * min dy - nu (k + K^T x): By >= b + K_psi^T nu   (fF)
     *                          y bounds        (yLbs, yUbs)
     */
    FollowerMC followerMC;
    if (true)
    {
        followerMC.loadProblem(data);
        followerMC.createProblem();
        // followerMC.setpsiUB(5000);
        followerMC.setpsiUB(u);
    }

    Master master;
    master.loadProblem(data);
    master.createProblem();
    master.setTimeLimit(t);
    master.setParams(s,h,ht);

    /* solve using Benders implemented with callback */
    master.solveCallback(follower, followerMC, leaderFollower, data);
    
    /* get results */
    cout << "objVal: " << master.getObjVal() << endl;
    cout << "num local cuts added: " << master.getNumLocalCutsAdded() << endl;
    cout << "num user branches: " << master.getNumUserBranches() << endl;
    ofstream bendersResultFile;
    bendersResultFile.open ("Benders_result.csv", ios::out | ios::app);
    bendersResultFile << data.problemName_ << " ," << master.getObjVal() << ", " << master.getBestObjVal() << ", " << master.getGap() * 100 << ", " << master.getTicToc().count() << ", " << master.getNumLocalCutsAdded() << ", " << master.getNumUserBranches()
                  << endl;
    bendersResultFile.close();

    // ofstream result_iter;
    // if (!sett.is_callback) {
    //     if (outputIterResultFile(result_iter, aInfo, instance, sett))
    //         return 1;
    // }

    /* for testing */
    if (testValue != NULL) {
        if (fabs(master.getObjVal() - *testValue) / (1e-10 + fabs(*testValue)) < 1e-4) {
            delete testValue;
            return 0;
        } else {
            delete testValue;
            return 1;
        }
    } else 
        return 0;
}

