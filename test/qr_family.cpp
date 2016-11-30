#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){


  vector< Matrix< double > > A(2, Matrix<double>(4, 4, true));

  for(int i = 0; i < (int)A.size(); i++){
    cout << "i:" << i << endl;
    uni10_rand(A[i], uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  }

  cout << A[0];
  cout << A[1];
  exit(0);

  //vector< Matrix<complex<double> > >QDR = qdr_cpivot(A);

  //cout << QDR[0];
  //cout << QDR[1];
  //cout << QDR[2];

  //dots(QDR[0], QDR[1],QDR[2]);

  //cout << QDR[0];

  return 0;
}
