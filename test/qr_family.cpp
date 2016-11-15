#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_A[12] = {8.1,2.1,3.2,
                      2.1,3.9,4.1,
                      2.1,4.1,1.1,
                      2.1,4.1,1.1,
                      };

  Matrix<double> A(4, 3);
  A.setElem(elem_A);

  fprintf(stdout, "-----    M_ori   -----\n");
  cout << A;

  vector<Matrix<double> > QR = qr( A );
  vector<Matrix<double> > QL = ql( A );
  
  fprintf(stdout, "-----    QR   -----\n");
  cout << A;
  cout << QR[0];
  cout << QR[1];
  cout << dot( QR[0], QR[1] );
  cout << transpose( QR[0] );
  cout << dot( transpose( QR[0]), QR[0] );

  fprintf(stdout, "-----    QL   -----\n");
  cout << A;
  cout << QL[0];
  cout << QL[1];
  cout << dot( QL[0], QL[1] );
  cout << dot(transpose( QL[0] ), QL[0]);

  transpose(A, INPLACE);

  vector<Matrix<double> > RQ = rq( A );
  vector<Matrix<double> > LQ = lq( A );

  fprintf(stdout, "-----    RQ   -----\n");

  cout << A;
  cout << RQ[0];
  cout << RQ[1];

  cout << dot( RQ[0], RQ[1] );
  cout << dot( RQ[1],  transpose(RQ[1]) );

  fprintf(stdout, "-----    LQ   -----\n");

  cout << A;
  cout << LQ[0];
  cout << LQ[1];

  cout << dot( LQ[0], LQ[1] );
  cout << dot( LQ[1], transpose( LQ[1] ) );
  exit(0);

  //cout << QR[0];
  //cout << QR[1];

  //fprintf(stdout, "----- RQ Example -----\n");

  //cout << rq( A )[0];
  //cout << rq( A )[1];

  //fprintf(stdout, "----- QL Example -----\n");

  //cout << ql( A )[0];
  //cout << ql( A )[1];

  //fprintf(stdout, "----- LQ Example -----\n");

  //cout << lq( A )[0];
  //cout << lq( A )[1];

  //fprintf(stdout, "-----     End    -----\n");

  return 0;
}
