#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};


  Matrix<double> A(3, 3);

  Matrix<double> B(A);

  cout << B;
  
  A.setElem(elem_A);

  fprintf(stdout, "-----     M_ori   -----\n");

  cout << A;

  fprintf(stdout, "----- SVD Example -----\n");

  cout << svd( A )[0];
  cout << svd( A )[1];
  cout << svd( A )[2];

  fprintf(stdout, "----- INV Example -----\n");

  cout << inverse( A );

  fprintf(stdout, "-----     End    -----\n");

  return 0;
}
