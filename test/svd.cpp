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

  fprintf(stdout, "-----     M_ori   -----\n");

  cout << A;

  fprintf(stdout, "----- SVD Example -----\n");

  vector< Matrix<double> > SVD = svd(A);
  cout << SVD[0];
  cout << SVD[1];
  cout << SVD[2];

  fprintf(stdout, "----- Multi Example -----\n");

  Matrix<double> A_mul;
  dots(A_mul, SVD[0], SVD[1], SVD[2]);

  cout << A_mul;

  fprintf(stdout, "----- SVD Inplace Example -----\n");

  Matrix<double> U, S, VT;
  svd(A, U, S, VT, INPLACE);
  cout << U;
  cout << S;
  cout << VT;

  Matrix<double> A_mul1;
  dots(A_mul1, SVD[0], SVD[1], SVD[2]);

  cout << A_mul1;

  return 0;
}
