#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  uni10_create();
  uni10_print_env_info();

  double elemA[] = {1., 2., 3, 4.};
  double elemB[] = {1., 4., 3, 4.};
  //cout << "==========" << endl;
  uni10_elem_cpu<double> A(elemA, 2, 2);
  //cout << "==========" << endl;
  A.setElem(elemB);
  //A.print_elem_cpu(2, 2);

  uni10_elem_cpu<double> B(A);
  //B.print_elem_cpu(4, 1);
  B.setElem(NULL);

   return 0;

}
