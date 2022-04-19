#include <Rcpp.h>
#include <cmath>
#include <fstream>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector range(int first, int last){
  if (first > last) return IntegerVector(0);
  IntegerVector ans(last - first + 1);
  for (int i = 0; i < ans.size(); ++i)
    ans[i] = first + i;
  return ans;
}


//[[Rcpp::export]]
int HammingScore(const IntegerVector& left, const IntegerVector& right){
  if (left.size() != right.size()) return -1;
  int score = 0;
  for (int i = 0; i < left.size(); ++i)
    score += (left[i] == right[i]);
  return score;
}


//[[Rcpp::export]]
IntegerMatrix MHStrSampling(int N, const IntegerVector& targetStr, int alphabetLen, double gamma,
                            const IntegerVector& startStr, const IntegerVector& posChange,
                            const IntegerVector& shift, const NumericVector& unif01) {
  size_t strLen = targetStr.size();
  IntegerMatrix x(N, strLen);
  x(0, _) = startStr;
  IntegerVector curr_x, proposed_x;
  for (int i = 1; i < N; ++i){
    proposed_x = x(i - 1, _);
    proposed_x[posChange[i] - 1] = (proposed_x[posChange[i] - 1] + shift[i]) % alphabetLen;
    int delta = HammingScore(targetStr, proposed_x) - HammingScore(targetStr, x(i - 1, _));
    double A = exp(gamma * delta);
    x(i, _) = (A < unif01[i]) ? x(i - 1, _) : proposed_x;
  }
  return x;
}

//[[Rcpp::export]]
NumericVector MHScoreSampling(int N, const IntegerVector& targetStr, int alphabetLen, double gamma,
                              const IntegerVector& startStr, const IntegerVector& posChange,
                              const IntegerVector& shift, const NumericVector& unif01) {
  //int strLen = targetStr.size();
  NumericVector scores(N);
  IntegerVector currStr = clone(startStr);
  scores[0] = HammingScore(targetStr, currStr);
  
  IntegerVector proposedStr;
  for (int i = 1; i < N; ++i){
    proposedStr = clone(currStr);
    proposedStr[posChange[i] - 1] = (proposedStr[posChange[i] - 1] + shift[i]) % alphabetLen;
    int proposedScore = HammingScore(proposedStr, targetStr);
    int delta = proposedScore - scores[i - 1];
    double A = exp(gamma * delta);
    
    if (A < unif01[i]) {
      scores[i] = scores[i - 1];
    } 
    else {
      scores[i] = proposedScore;
      currStr = clone(proposedStr);
    }
  }
  return scores;
}


//[[Rcpp::export]]
IntegerMatrix Vector2Matrix(const IntegerVector& v, int nrow, int ncol){
  IntegerMatrix matr(nrow, ncol);
  for (int i = 0; i < nrow; ++i)
    for (int j = 0; j < ncol; ++j)
      matr(i, j) = v[i * ncol + j];
  return matr;
}

//[[Rcpp::export]]
void swapRows(IntegerMatrix& matr, int row1, int row2){
  int nrow = matr.nrow();
    for (int j = 0; j < matr.ncol(); ++j)
      std::swap(matr[j * nrow + row1], matr[nrow * j + row2]);
}


//[[Rcpp::export]]
IntegerVector getRow(const IntegerMatrix& matr, int ind){
  IntegerVector row(matr.ncol());
  int nrow = matr.nrow();
  for (int i = 0; i < matr.ncol(); ++i)
    row[i] = matr[i * nrow + ind];
  return row;
}

//[[Rcpp::export]]
double getElement(const IntegerMatrix& matr, int row, int col){
  return matr[col * matr.nrow() + row];
}

//[[Rcpp::export]]
void assign(IntegerMatrix& matr, int row, int col, int val){
  matr[matr.nrow() * col + row] = val;
}

//[[Rcpp::export]]
void assignRow(IntegerMatrix& matr, int row, const IntegerVector& v){
  for (int j = 0; j < matr.ncol(); ++j)
    assign(matr, row, j, v[j]);
}

//[[Rcpp::export]]
void swapElements(IntegerMatrix& matr, int row1, int col1, int row2, int col2){
  std::swap(matr[matr.nrow() * col1 + row1], matr[matr.nrow() * col2+ row2]);
}

//[[Rcpp::export]]
IntegerVector MC3ScoreSample(int N, const IntegerVector& targetStr, int alphabetLen, 
                             const NumericVector& gammas, int step)
{
  
  int strLen = targetStr.size(), nchains = gammas.size();
  IntegerMatrix scores(nchains, N);
  IntegerMatrix currStr = Vector2Matrix(sample(range(0, alphabetLen - 1), nchains * strLen, true), 
                                       nchains, strLen);
  
  for (int i = 0; i < nchains; ++i)
    //scores[i, 0] = HammingScore(targetStr, currStr[i]);
    assign(scores, i, 0, HammingScore(targetStr, currStr(i, _)));
  
  IntegerVector posChange = sample(range(0, strLen - 1), N * nchains, true);
  IntegerVector shift = sample(range(0, alphabetLen - 1), N * nchains, true);
  NumericVector unif01 = runif(N * nchains + 0.5 * (nchains + 1) * (N / step + 1));
  
  
  int unifCounter = 0;
  for (int i = 1; i < N; ++i){
    for (int j = 0; j < nchains; ++j){
      IntegerVector proposedStr(getRow(currStr, j));
      proposedStr[posChange[(i - 1) * nchains + j]] = 
        (proposedStr[posChange[(i - 1) * nchains + j]] + shift[i]) % alphabetLen;
      int proposedScore = HammingScore(targetStr, proposedStr);
      //int delta = proposedScore - scores[j, i - 1];
      int delta = proposedScore - getElement(scores, j, i - 1);
      //if (j > 2) return(IntegerVector(1));
      double A = exp(gammas[j] * delta);
      
      if (A < unif01[unifCounter++]) {
        //scores[j, i] = scores[j, i - 1];
        assign(scores, j, i, getElement(scores, j, i - 1));
      }
      else {
        //scores[j, i] = proposedScore;
        assign(scores, j, i, proposedScore);
        //for (int k = 0; k < strLen; ++k) currStr[j, k] = proposedStr[k];
        assignRow(currStr, j, proposedStr);
        
      }
    }
    
    if (i % step == 0){
      IntegerVector perm = sample(range(0, nchains - 1), nchains, false);
      for (int j = 0; j < nchains / 2; ++j){
        int delta = getElement(scores, perm[j], i) - getElement(scores, perm[nchains - j - 1], i);
        double A = exp((gammas[perm[j]] - gammas[perm[nchains - j - 1]]) * delta);
        if (A > unif01[unifCounter++]) {
          //std::swap(scores[perm[j]], scores[perm[nchains - j - 1]]);
          swapElements(scores, perm[j], i, perm[nchains - j - 1], i);
          //std::swap(currStr[perm[j]], currStr[perm[nchains - j - 1]]);
          swapRows(currStr, perm[j], perm[nchains - j - 1]);
        }
      }
    }
  }
  return getRow(scores, 0);
}
  
//[[Rcpp::export]]
double get(const NumericMatrix& matr, int ind){
  return matr[ind - 1];
}
