#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector DPLOOP(NumericVector p, NumericVector id)
{
  int n = p.size() / 2;
  NumericVector res(n);
  int maxp = max(p);
  for(int i = 0; i < res.size(); ++i)
    {
      res[i] = maxp;
    }
  for(int i = 0; i < p.size(); ++i)
    {
      int is = i-2;
      if(is < 0){
	is = 0;
      }
      int ie = i+2;
      if(ie > p.size()){
	ie = p.size();
      }
      for(int j = is; j < ie; ++j)
	{
	  if(id[j] != id[i] & res[id[i]] > abs(p[j]-p[i]))
	    {
	      res[id[i]] = abs(p[j]-p[i]);
	    }
	}
    }
  return wrap(res);
}
