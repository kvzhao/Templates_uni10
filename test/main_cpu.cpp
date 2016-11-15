#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  uni10_create();
  uni10_print_env_info();

  uni10_double64 elemA[] = {1., 2., 3, 4.};
  uni10_complex128 elemB[] = {uni10_complex128(1., 4.), uni10_complex128(3, 4.32), uni10_complex128(2.12, 4.12), uni10_complex128(33.12, 4.12)  };

  Matrix<uni10_double64> A(2, 2);
  Matrix<uni10_complex128> B(4, 4, true);

  cout << A;
  cout << B;

  printf("\n------ Matrix A ------\n");
  A.setElem(elemA);
  cout << A;

  printf("\n------ Matrix B ------\n");
  B.setElem(NULL);
  cout << B;

   return 0;

}
