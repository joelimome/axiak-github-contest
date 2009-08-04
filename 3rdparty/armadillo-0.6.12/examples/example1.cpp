#include <iostream>

#include "armadillo"

using namespace arma;
using namespace std;


int main(int argc, char** argv)
  {
  // matrix construction from a string representation
  mat A = \
   "\
   0.555950  0.274690  0.540605  0.798938;\
   0.108929  0.830123  0.891726  0.895283;\
   0.948014  0.973234  0.216504  0.883152;\
   0.023787  0.675382  0.231751  0.450332;\
   ";

  // print to the cout stream
  // with an optional string before the contents of the matrix
  A.print("A =");

  // save to disk
  A.save("A.txt", raw_ascii);

  mat B(12,34);
  // .n_rows = number of rows (read only)
  // .n_cols = number of columns (read only)
  cout << "B.n_rows = " << B.n_rows << endl;
  cout << "B.n_cols = " << B.n_cols << endl;
  
  // scalars are treated as a 1x1 matrix,
  // hence the code below will set B to have a size of 1x1
  B = 5.0;

  // load from disk
  B.load("A.txt");


  B += 2.0 * A;

  // print the contents of B to an arbitrary stream
  // (cout in this case) 
  cout << "B =" << endl << B << endl;

  // submatrix types:
  //
  // .submat(first_row, first_column, last_row, last_column)
  // .row(row_number)
  // .col(column_number)
  // .cols(first_column, last_column)
  // .rows(first_row, last_row)

  cout << "B.submat(0,0,3,1) =" << endl;
  cout << B.submat(0,0,3,1) << endl;

  // generate the identity matrix
  mat C = eye(4,4);
  
  C.submat(0,0,3,1) = B.cols(1,2);
  C.print("C =");

  // transpose
  cout << "trans(A) =" << endl;
  cout << trans(A) << endl;

  // maximum from each row
  cout << "max(A) =" << endl;
  cout << max(A) << endl;

  // maximum from each column
  cout << "max(A,1) =" << endl;
  cout << max(A,1) << endl;

  // maximum value in A
  cout << "max(max(A)) = " << max(max(A)) << endl;

  // sum along rows
  cout << "sum(A) =" << endl;
  cout << sum(A) << endl;

  // sum along columns
  cout << "sum(A,1) =" << endl;
  cout << sum(A,1) << endl;

  // sum of all elements
  cout << "sum(sum(A)) = " << sum(sum(A)) << endl;
  cout << "accu(A)     = " << accu(A) << endl;

  // trace = sum along diagonal
  cout << "trace(A)    = " << trace(A) << endl;

  // random matrix;
  // values are uniformly distributed in the [0,1] interval
  mat D = rand(4,4);
  D.print("D =");

  cout << endl;

  // row vectors are treated like a matrix with one row
  rowvec r = "0.59499  0.88807  0.88532  0.19968";

  // column vectors are treated like a matrix with one column
  colvec q = "0.81114  0.06256  0.95989  0.73628";

  r.print("r =");
  q.print("q =");

  // dot or inner product
  cout << "r*q = " << r*q << endl;
  

  // outer product
  cout << "q*r =" << endl;
  cout << q*r << endl;

  // multiply-and-accumulate operation
  // (no temporary matrices are created)
  cout << "accu(A % B) = " << accu(A % B) << endl;

  // sum of four matrices;
  // (no temporary matrices are created)
  mat F = A + B + C + D;
  F.print("F =");


  // two dimensional field of arbitrary length row vectors
  field<rowvec> xyz(3,2);

  xyz(0,0) = rand(1,2);
  xyz(1,0) = rand(1,3);
  xyz(2,0) = rand(1,4);
  xyz(0,1) = rand(1,5);
  xyz(1,1) = rand(1,6);
  xyz(2,1) = rand(1,7);

  cout << "xyz = " << endl;
  cout << xyz << endl;


  // comparison of matrices (element-wise)
  imat AA = "1 2 3; 4 5 6; 7 8 9;";
  imat BB = "3 2 1; 6 5 4; 9 8 7;";
  
  umat ZZ = (AA >= BB);
  ZZ.print("ZZ =");


  return 0;
  }

