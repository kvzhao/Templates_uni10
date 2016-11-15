#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  std::complex<double> elem_C[9] = {std::complex<double>(8.1,2.1),std::complex<double>(6.2,12.1),std::complex<double>(6.2,12.1),
                      std::complex<double>(2.1,3.9),std::complex<double>(40.1,3.90),std::complex<double>(101.1,3.9),
                      std::complex<double>(2.1,4.1),std::complex<double>(1.1, 21.0),std::complex<double>(8.1, 1.1)};

  vector<int> CC(3, 2);

  Matrix<double> A(3, 3);

  Matrix<double> B(A);

  Matrix< std::complex<double> > C(3, 3);

  C.setElem(elem_C);

  A.setElem(elem_A);

  cout << det( A ) << endl;

  cout << det( C ) << endl;

  //A += A;

  //cout << A;

  //A -= B;

  //cout << A;

  //cout << B;
  //

  //fprintf(stdout, "-----     M_ori   -----\n");

  //cout << A;

  //fprintf(stdout, "----- SVD Example -----\n");

  //cout << svd( A )[0];
  //cout << svd( A )[1];
  //cout << svd( A )[2];

  //fprintf(stdout, "----- INV Example -----\n");

  //cout << inverse( A );

  //fprintf(stdout, "-----     End    -----\n");

  return 0;
}
