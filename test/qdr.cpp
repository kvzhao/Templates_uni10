#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Matrix<double> M;

  Matrix<double> A(4, 4);
  uni10_rng(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);

  //vector< Matrix<double> > QDR = qdr( A );
  //vector< Matrix<double> > QR = qr( A );

  vector< Matrix<double> > LDQ = qdr_cpivot( A );
  vector< Matrix<double> > LQ  = qr( A );
  cout << "==== Ori =====\n";

  cout << "   ldq: \n";
  cout <<LDQ[0];
  cout <<LDQ[1];
  cout <<LDQ[2];

  cout << "   lq: \n";
  cout <<LQ[0];
  cout <<LQ[1];

  //cout << "==== qr =====\n";
  cout << dot( LQ[0], LQ[1] );
  dots( M, LDQ[0], LDQ[1], LDQ[2] );
  cout << M;
  cout << A;

  return 0;
}
