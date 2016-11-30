#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  Matrix<double> A(3, 3);

  A.setElem(elem_A);

  cout << exph( 2.0, A );

}
