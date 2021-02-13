/*
 * nodeCBfObj.h
 *
 *  Created on: Jan 24, 2021
 *      Author: geunyeongbyeon
 */

#ifndef BILEVEL_NODECBFOBJ_H
#define BILEVEL_NODECBFOBJ_H

#include "callback.h"

class nodeSelectCallbackfObjI : public IloCplex::NodeCallbackI {
public:
  ILOCOMMONCALLBACKSTUFF(nodeSelectCallbackfObj);
  nodeSelectCallbackfObjI(IloEnv env, bool* nodeType, IloInt64 *nodeId) : IloCplex::NodeCallbackI(env) 
  {
      node_type = nodeType;
      node_id = nodeId;
      num_selected_nodes = 0;
      count = 0;
      current_best_objVal = -INFINITY;
  };
  
  void main() override;
  
  IloInt i;
  int count;
  double current_best_objVal;
  IloInt bestnode;
  double bestbound;
  int num_selected_nodes;
  bool *node_type; /* true: node with dy <= fObj is selected */
  IloInt64 *node_id;
};
IloCplex::Callback nodeSelectCallbackfObj(IloEnv env, bool* nodeType, IloInt64 *nodeId);

#endif //BILEVEL_NODECBFOBJ_H