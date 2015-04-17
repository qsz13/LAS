//#include <Rcpp.h>
//using namespace Rcpp;


// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//// [[Rcpp::export]]
//double lascoreold(const NumericVector& x, const NumericVector& y, const NumericVector& z) {
//  if(x.size()!=y.size()||x.size()!=z.size())
//  {
//    stop("x, y, z sizes are not the same!");
//  }
//  int size = x.size();
//  NumericVector::const_iterator itx, ity, itz;
//  double sum = 0;
//  for(itx = x.begin(),ity = y.begin(),itz = z.begin();
//      itx!=x.end() && ity!=y.end() && itz!=z.end(); 
//      ++itx,++ity,++itz)
//      {
//          sum+=(*itx)*(*ity)*(*itz);
//      }
//  return sum/size;
//}

//// [[Rcpp::export]]
//NumericVector lascoreXY(NumericVector x, NumericVector y, NumericMatrix m)
//{
//  int nrow = m.nrow();
// NumericVector score(nrow);
//  for (int i = 0; i < nrow; i++) {
//    score[i] = lascore(x, y, m.row(i));
//  }
//  return score;
//}

