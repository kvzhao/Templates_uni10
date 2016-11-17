#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Matrix<double> A(5, 4);
  uni10_orthoRand(A, uni10_mt19937, uni10_normal, 0, 1, 1);
  cout << A;
  cout << dot(transpose(A), A);

  identity(A);
  cout << A;


  return 0;
}
