/* bilevel.cpp
 * 
 * Created on: Nov 7, 2020
 *     Author: geunyeongbyeon
 */
#include <iostream>
#undef NDEBUG
#include <assert.h>

#include "master.h"
#include "hpp.h"

int getArgValues (int argc, char *argv[], string &filename, double * &testValue){

    cout << "This program was called with " << argv[0] << endl;

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
        
        // else if (strcmp(argv[j], "-m") == 0) {

        //     // solution method
        //     if (strcmp(argv[v], "MIP") == 0)
        //         sett.is_MIP = true;
        //     else if (strcmp(argv[v], "Benders") == 0)
        //         sett.is_MIP = false;
        //     else if (strcmp(argv[v], "Callback") == 0) {
        //         sett.is_callback = true;
        //         sett.is_MIP = false;
        //     } else {
        //         cout << "-m should be chosen among 'MIP', 'Benders', 'Callback'" << endl;
        //         return 1;
        //     }
        // } 
        // else if (strcmp(argv[j], "-c") == 0) {

        //     // Benders cut type
        //     if (strcmp(argv[v], "o") == 0)      // original
        //         sett.cutStyle = "o";
        //     else if (strcmp(argv[v], "d") == 0) // decomposed
        //         sett.cutStyle = "d";
        //     else if (strcmp(argv[v], "m") == 0) // modified
        //         sett.cutStyle = "m";
        //     else {
        //         cout << "-c should be chosen among 'o', 'd', 'm'" << endl;
        //         return 1;
        //     }
        // } else if (strcmp(argv[j], "-nc") == 0) {

        //     sett.use_nogoodcut = stoi(argv[v]);

        // } else if (strcmp(argv[j], "-n") == 0) {
        //     /* normalized ray */
        //     sett.use_normalization = stoi(argv[v]);

        // } else if (strcmp(argv[j], "-i") == 0) {
        //     /* in-out method */
        //     sett.use_in_out = stoi(argv[v]);

        // } else if (strcmp(argv[j], "-o") == 0) {
        //     sett.outputConsole = stoi(argv[v]);
        // } else if (strcmp(argv[j], "-LBCut") == 0) {
        //     sett.useLBCut = stoi(argv[v]);
        // } else if (strcmp(argv[j], "-p") == 0) {
        //     sett.preprosessing = stoi(argv[v]);
        // }
    }

//    // check
//    cout << "File: " << instance.name << endl;
//    cout << "is_MIP: " <<  sett.is_MIP << endl;
//    cout << "norm: " << sett.use_normalization << endl;
//    cout << "inout: " << sett.use_in_out << endl;

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

    /* create hpp problem:
     * min c_x x + c_y y : G_x x + G_y y >= h   (l)
     *                     Ax + By >= b         (f)
     *                     x bounds, y bounds   (xBds, yLbs, yUbs)
     */
    Hpp hpp;
    hpp.loadProblem(data);
    hpp.createProblem();

    /* create follower problem:
     * min dy : By >= b (-Ax)   (fF)
     *          y bounds        (yLbs, yUbs)
     */
    Follower follower;
    follower.loadProblem(data);
    follower.createProblem();
    
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

    /* calculate upper bound of follower obj val */
    if (hpp.solvefUb())
        follower.setfUb(hpp.getfUb());
    // follower.setfUb(100000);

    Master master;
    master.loadProblem(data);
    master.createProblem();
    master.setTimeLimit(3600);

    /* solve using Benders implemented with callback */
    master.solveCallback(follower, followerMC, leaderFollower);
    
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

