#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  vector<Bond> bonds1(3, Bond(BD_IN, 3));
  bonds1[2] = Bond(BD_OUT, 2);

  vector<Bond> bonds2(3, Bond(BD_IN, 3));
  bonds2[2] = Bond(BD_OUT, 4);

  Matrix<double> BA(9, 2);
  uni10_rand(BA, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);

  Matrix<double> BA2(9, 4);
  uni10_rand(BA2, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);

  UniTensor<double> A(bonds1);
  A.putBlock(BA);

  UniTensor<double> A2(bonds2);
  A2.putBlock(BA2);

  cout << A;
  cout << A2;

  Network<double> N1("N1.net");
  N1.putTensor(0, A);
  N1.putTensor(1, A2);

  UniTensor<double> C = N1.launch();

  cout << C;
  return 0;
}
