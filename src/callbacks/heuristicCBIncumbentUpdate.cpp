#include "heuristicCBIncumbentUpdate.h"


IloCplex::Callback incumbentUpdateCallback(IloEnv env, IloNumVarArray& xVars, LazyData& lazyData) {
  return (IloCplex::Callback(new (env) incumbentUpdateCallbackI(env, xVars, lazyData)));
}
void incumbentUpdateCallbackI::main() {
    try {
        if (lazyData_.found_new_incumbent) {
            if (lazyData_.current_best_ub < getIncumbentObjValue()) {
//                cout << "setSolution()" << endl;
                setSolution(xVars_, lazyData_.feas_x, lazyData_.current_best_ub);
            }
            if (lazyData_.current_best_ub >= getIncumbentObjValue()) {
                lazyData_.found_new_incumbent = false;
            }
        }
    }
    catch (...) {
        throw;
    }
}
