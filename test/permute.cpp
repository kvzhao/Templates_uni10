#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Block<double> A;
  Matrix<double> buf(3, 9);
  Matrix<double> buf2(3, 9);

  uni10_rand(buf, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);
  uni10_rand(buf2, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

  cout << buf;
  cout << buf2;

  std::cout << "-=====================- \n";

  vector<Bond> bonds(2, Bond(BD_OUT, 9));
  bonds[0] = Bond(BD_IN, 3);
  //UniTensor<double> B(2);
  int label[] = {4, 5};
  UniTensor<double> B(bonds, label);

  B.putBlock(buf);

  UniTensor<double> C(B);

  std::cout << C;

  C.putBlock(buf2);

  std::cout << B;
  std::cout << C;

  UniTensor<double> D(C);
  set_zeros(D);

  std::cout << C;
  std::cout << D;

  int newlabel[] = {5, 4};
 
  std::cout << "-=====================-\n";
  std::cout << B;
  cout << permute(B, newlabel, 1);
  
  cout << 2.0 * B;

  return 0;
}
