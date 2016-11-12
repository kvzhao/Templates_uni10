#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  Matrix<double> A(3, 3);

  Matrix< std::complex<double> > C(3, 3);

  A.setElem(elem_A);

  Matrix<double> B(A);

  fprintf(stdout, "-----     Reszie [3,3]->[3,5]  -----\n");

  cout << A;

  resize( A, 3, 5 );

  cout << A;

  fprintf(stdout, "-----     Reszie [3,5]->[5,5]  -----\n");

  resize( A, 5, 5 );

  cout << A;

  fprintf(stdout, "-----     Reszie [5,5]->[2,2]  -----\n");
  
  resize( A, 2, 2 );

  cout << A;

  fprintf(stdout, "------------------------\n");

  return 0;
}
