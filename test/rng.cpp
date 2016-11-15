#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  Matrix<double> A(3, 3);

  A.setElem(elem_A);

  cout << A;

  A[1] = 32.1;
  A[2] = 200.1;

  cout << A;

  uni10_rng(A, uni10_mt19937, uni10_uniform_real, -1, 1, -1);

  cout << A;

  
  return 0;
}
