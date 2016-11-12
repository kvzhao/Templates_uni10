#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  vector<int> CC(3, 2);

  Matrix<double> A(3, 3);

  Matrix< std::complex<double> > C(3, 3);

  A.setElem(elem_A);

  Matrix<double> B(A);

  cout << A;

  A += A;

  cout << A;

  A -= B;

  cout << A;

  cout << B;
  

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