#include "nodeCB.h"

IloCplex::Callback nodeSelectCallback(IloEnv env,  double *bLb, double *bUb, bool *isEnded, double fLb, double fRange) {
  return (IloCplex::Callback(new (env) nodeSelectCallbackI(env, bLb, bUb, isEnded, fLb, fRange)));
}
void nodeSelectCallbackI::main() {
   IloInt remainingNodes = getNremainingNodes();
   
   if (!*isEnded) {
   for (IloInt i = 0; i < remainingNodes; i++) {
      
      if (getNodeData(i))
      {
         UserNodeData * nd = (UserNodeData *) getNodeData(i);
         cout << i << ": " << nd->dy_lb <<  ", " << nd->dy_ub << endl;
      }
   }
   
   int discretization = 5;
   if (k <= discretization) 
   {
      if (l < pow(2,k-1)) 
      {
         /* If k = 3, it will iterate over (k, l) = (1,0), (2,0), (2,1), (3,0), (3,1), (3,2)
               * for each branching node it will branch having (2*l+1) / 2^k at the center
               * i.e., 1/2, [1/4, 3/4], [1/8, 3/8, 5/8, 7/8]
               *      Branching node at level 0: [flb, fub]
               *      Branching node 1 at level 1: [flb, flb + frange * 1/2]
               *      Branching node 2 at level 1: [flb + frange * 1/2, fub]
               *      Branching node 1 at level 2: [flb, flb + frange * 1/4]
               *      Branching node 2 at level 2: [flb + frange * 1/4, flb + frange * 1/2]
               *      Branching node 3 at level 2: [flb + frange * 1/2, flb + frange * 3/4]
               *      Branching node 4 at level 2: [flb + frange * 1/4, fub]
               *      Branching node 1 at level 3: [flb, flb + frange * 1/8]
               *      Branching node 2 at level 3: [flb + frange * 1/8, flb + frange * 2/8]
               *      Branching node 3 at level 3: [flb + frange * 2/8, flb + frange * 3/8]
               *      Branching node 4 at level 3: [flb + frange * 3/8, flb + frange * 4/8]
               *      Branching node 5 at level 3: [flb + frange * 4/8, flb + frange * 5/8]
               *      Branching node 6 at level 3: [flb + frange * 5/8, flb + frange * 6/8]
               *      Branching node 7 at level 3: [flb + frange * 6/8, flb + frange * 7/8]
               *      Branching node 8 at level 3: [flb + frange * 7/8, fub]
               */
         
         *blb = flb + frange / pow(2,k) * (2 * l + j);
         *bub = flb + frange / pow(2,k) * (2 * l + 1 + j);
         cout << "choose with lb: " << *blb  <<" and ub: " << *bub << endl;
         
         j++;

         if (j == 2) {
            l++;
            if (l == pow(2,k-1)) {
               k++;
               if (k > discretization)
                  *isEnded = true;
               l = 0;
            }
            j = 0;
         }
         
         
         for (IloInt i = 0; i < remainingNodes; i++) {
            if (getNodeData(i))
            {
               UserNodeData * nd = (UserNodeData *) getNodeData(i);
               if (*blb == nd->dy_lb && *bub == nd->dy_ub){
                  selectNode(i);
                  cout << "chosen i: " << i << endl;
                  break;
               }
            }
         }
      }
   }
   call_++;
   }
}