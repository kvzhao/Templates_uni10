#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Block<double> A;
  Matrix<double> buf(3, 9);

  uni10_rand(buf, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

  vector<Bond> bonds(3, Bond(BD_OUT, 3));
  bonds[0] = Bond(BD_IN, 3);
  //UniTensor<double> B(2);
  int label[] = {4, 5, 10};
  UniTensor<double> B(bonds, label);

  cout << B;

  cout << buf;

  double cc[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    ,1, 2, 3, 4, 5, 6, 7, 8, 9, 10
      ,1, 2, 3, 4, 5, 6, 7
  };

  B.setElem(cc);
  cout << B;
  B.putBlock(buf);
  cout << B;
  //cout << B.getBlock();

  int label2[] = {1, 3, 2};
  B.setLabel(label2);
  cout << B;

  return 0;
}
