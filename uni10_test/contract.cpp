#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  int dim1 = 7, dim2 = 10, dim3 = 5;
  ///A is matrix to be decomposed
  Matrix<double> A(dim1, dim2);
  Matrix<double> B(dim2, dim3);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  uni10_rand(B, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  
  Bond bond1( BD_IN, dim1);
  Bond bond2i( BD_IN, dim2);
  Bond bond2o( BD_OUT, dim2);
  Bond bond3( BD_OUT, dim3);
  vector<Bond> bds_At( 1, bond1);
  bds_At.push_back( bond2o );
  vector<Bond> bds_Bt( 1, bond2i);
  bds_Bt.push_back( bond3 );

  UniTensor<double> A_t( bds_At);
  UniTensor<double> B_t( bds_Bt);
  vector<int> A_label {1, 2};
  vector<int> B_label {2, 3};
  A_t.setLabel(A_label);
  B_t.setLabel(B_label);
  A_t.putBlock(A);
  B_t.putBlock(B);
  
  UniTensor<double> C_t = contract( A_t, B_t, false);
  Matrix<double> AB;
  dots( AB, A, B );

  cout<<norm((AB - C_t.getBlock()))<<endl;

  return 0;
}
