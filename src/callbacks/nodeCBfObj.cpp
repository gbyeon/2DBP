#include "nodeCBfObj.h"
// #define DEBUG
IloCplex::Callback nodeSelectCallbackfObj(IloEnv env, bool* nodeType, IloInt64 *nodeId) {
  return (IloCplex::Callback(new (env) nodeSelectCallbackfObjI(env, nodeType, nodeId)));
}
void nodeSelectCallbackfObjI::main() {
   
   try
   { 
      if (fabs(current_best_objVal - getBestObjValue()) < 1e-3)
         count++;
      else count = 0;
      current_best_objVal = getBestObjValue();

      IloInt remainingNodes = getNremainingNodes();
      IloInt nNodes = getNnodes();

      if ((count >= 3 && nNodes <= 25000) || (count >= 2 && nNodes > 25000))
      {
         num_user_branches_++;
         #ifdef DEBUG 
         cout << "number of remaining nodes: " << remainingNodes << endl;
         #endif
         
         if (!*node_type) 
         {
            
            bestnode = -1;
            bestbound = INFINITY;
            for (i = 0; i < remainingNodes; i++) 
            {
               #ifdef DEBUG
               cout << i << "th objval: " << getObjValue(i) << " at Depth " << getDepth(i) << endl;
               #endif
               if (getObjValue(i) < bestbound)
               {
                  bestnode = i;
                  bestbound = getObjValue(i);
               }
            }
            #ifdef DEBUG
            cout << "node " << bestnode << " with objVal " << bestbound << " selected" << endl;
            #endif

            if (bestnode > -1)// && !getNodeData(bestnode))
            {   
               selectNode(bestnode);
               *node_type = true;
               *node_id = num_selected_nodes;
               num_selected_nodes++;
            }
               
         } else 
         {
            *node_type = false;
            #ifdef DEBUG
            cout << "node_id at nodeCB: " << *node_id << endl;
            #endif
            for (i = 0; i < remainingNodes; i++) 
            {
               if (getNodeData(i))
               {
                  UserNodeData * nd = (UserNodeData *) getNodeData(i);
                  if (*node_id == nd->id)
                  {
                     selectNode(i);
                     break;
                  }
               }
            }
            count = 0;
         }
      }
   } catch (IloException &e) {
        cerr << "CPLEX found the following exception: " << e << " in nodeCBfObj.cpp" << endl;
        e.end();
    }
    catch (...) {
        cerr << "The following unknown exception was found in nodeCBfObj.cpp: " << endl;
    }

    return;
}