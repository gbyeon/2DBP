/* bilevel.cpp
 * 
 * Created on: Nov 7, 2020
 *     Author: geunyeongbyeon
 */
#include <iostream>
#undef NDEBUG
#include <assert.h>

#include "master.h"

int getArgValues (int argc, char *argv[], string &filename, double * &testValue){

    for (int i = 0; i < floor((double) argc / (double) 2); i++) {

        int j = 2 * i + 1;
        int v = 2 * (i + 1);

        if (strcmp(argv[j], "-f") == 0) {
            /* instance path/name */
            filename = string(argv[v]);
        } 
        else if (strcmp(argv[j], "-test") == 0) {
            /* test value */
            testValue = new double;
            *testValue = stod(argv[v]);
            cout << "testing value: " << *testValue << endl;
        }
    } 

   // check
    cout << "This program was called with " << argv[0] << endl;
    cout << "filename: " << filename << endl;
 
    return 0;
}

int main (int argc, char *argv[]) {

    string filename;
    double * testValue = NULL;
    if (getArgValues(argc, argv, filename, testValue))
        return 1;

    Data data;
    if (data.read(filename)) {
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
    // followerMC.loadProblem(data);
    // followerMC.createProblem();
    // followerMC.setpsiUB(5000);

    Master master;
    master.loadProblem(data);
    master.createProblem();
    master.setTimeLimit(3600);

    /* solve using Benders implemented with callback */
    master.solveCallback(follower, followerMC, leaderFollower, data);
    
    /* get results */
    cout << "objVal: " << master.getObjVal() << endl;
    ofstream bendersResultFile;
    bendersResultFile.open ("Benders_result.csv", ios::out | ios::app);
    bendersResultFile << data.problemName_ << " ," << master.getObjVal() << ", " << master.getBestObjVal() << ", " << master.getGap() * 100 << ", " << master.getTicToc().count()
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

