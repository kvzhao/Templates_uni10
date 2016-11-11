#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  double elem_T[9] = {0. ,2.1,3.2,
                      2.1,0. ,4.1,
                      2.1,4.1,0. };

  Matrix<double> A(3, 3);
  A.setElem(elem_A);

  fprintf(stdout, "Copy Constructor with the same type.\n");

  Matrix<double> B(A);

  cout << A;
  cout << B;

  A.setElem(elem_T);

  cout << A;
  cout << B;

  fprintf(stdout, "assignment operator with the same type.\n");
  //
  A.setElem(elem_A);
  Matrix<double> C = A;

  cout << A;
  cout << C;

  A.setElem(elem_T);

  cout << A;
  cout << C;

  return 0;

}
