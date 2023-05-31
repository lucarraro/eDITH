#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector evalConc2_cpp(S4 OCN, IntegerVector ss, NumericVector AS, double tau, NumericVector p, String str = "RN", bool normalize = false){
  //
  List L = OCN.slot(str);
  double nNodes = L["nNodes"];
 NumericVector leng = L["leng"];
 NumericVector width = L["width"];
 NumericVector velocity = L["velocity"];
 NumericVector depth = L["depth"];
 IntegerVector downNode = L["downNode"];
 // NumericVector AS = leng * width;
 NumericVector  Q = depth * velocity * width;
 
 NumericVector conc (nNodes);

 for (int i{ 0 }; i<nNodes; ++i)
{
 double node = ss[i];
 
 double localP = AS[node-1]*p[node-1];
 double decayFact = exp(-leng[node-1]/velocity[node-1]/tau)/Q[node-1];
 conc[node-1] = (conc[node-1] + localP)*decayFact;
 double flux = conc[node-1]*Q[node-1];
 int dN = downNode[node-1];
 conc[dN-1] = conc[dN-1] + flux;

}
  return(conc);
}