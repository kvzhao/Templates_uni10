#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[9] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1};

  Matrix<double> A(3, 3);
  A.setElem(elem_A);

  fprintf(stdout, "-----    M_ori   -----\n");

  cout << A;

  fprintf(stdout, "----- QR Example -----\n");

  cout << qr( A )[0];
  cout << qr( A )[1];

  //fprintf(stdout, "----- RQ Example -----\n");

  //cout << rq( A )[0];
  //cout << rq( A )[1];

  fprintf(stdout, "----- QL Example -----\n");

  cout << ql( A )[0];
  cout << ql( A )[1];

  fprintf(stdout, "----- LQ Example -----\n");

  cout << lq( A )[0];
  cout << lq( A )[1];

  fprintf(stdout, "-----     End    -----\n");

  return 0;
}
