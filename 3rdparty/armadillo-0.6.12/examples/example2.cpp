#include <iostream>

#include "armadillo"

using namespace arma;
using namespace std;


int main(int argc, char** argv)
  {
  // matrix construction from a string representation
  mat A = \
   "\
   0.165300  0.454037  0.995795  0.124098  0.047084;\
   0.688782  0.036549  0.552848  0.937664  0.866401;\
   0.348740  0.479388  0.506228  0.145673  0.491547;\
   0.148678  0.682258  0.571154  0.874724  0.444632;\
   0.245726  0.595218  0.409327  0.367827  0.385736;\
   ";

  // print to the cout stream
  // with an optional string before the contents of the matrix
  A.print("A =");

  // determinant
  cout << "det(A) = " << det(A) << endl;

  // inverse
  cout << "inv(A) = " << endl << inv(A) << endl;


  return 0;
  }

