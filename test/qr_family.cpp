#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[16] = {8.1,2.1,3.2,1.2,
                      2.1,3.9,4.1,20,
                      92.1,4.1,21.1,0.2,
                      2.1,41.1,1.1,32
                      };


  Matrix<double> A(4, 4, elem_A);

  cout << A;

  vector< Matrix<double> > QDR = qdr_cpivot(A);

  cout << QDR[0];
  cout << QDR[1];
  cout << QDR[2];

  dots(QDR[0], QDR[1],QDR[2]);

  cout << QDR[0];

  return 0;
}
