#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Matrix<double> Res;

  Matrix<double> A(4, 4);
  uni10_rng(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<double> B(4, 4);
  uni10_rng(B, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> C(4, 4);
  uni10_rng(C, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> D(4, 4);
  uni10_rng(D, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> E(4, 4);
  uni10_rng(E, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> F(4, 4);
  uni10_rng(F, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> G(4, 4);
  uni10_rng(G, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> H(4, 4);
  uni10_rng(H, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );
  Matrix<double> I(4, 4);
  uni10_rng(I, uni10_mt19937, uni10_normal, 0, 1, uni10_clock );

  cout << "\n==========          Dirty dots          ================\n";
  cout << dot(A, dot(B, dot(C, dot(D, dot(E, dot(F, dot(G, dot(H, I))))))));

  cout << "\n==========  assign in empty Matrix Res  ================\n";
  dots(Res, A, B, C, D, E, F, G, H, I);
  cout << Res;

  cout << "\n==========       a = a * b * ,,,        ================\n";
  dots(A, B, C, D, E, F, G, H, I);
  cout << A;

}
