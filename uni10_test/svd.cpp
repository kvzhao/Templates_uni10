#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  int dim1 = 7, dim2 = 10;
  ///A is matrix to be decomposed
  Matrix<double> A(dim1, dim2);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  ///u, s, vt are answers
  Matrix<double> u, s, vt;
  svd( A, u, s, vt, INPLACE);
  vector<Matrix<double>> svds = svd(A);
  ///u_s_vt multiplication of answers
  Matrix<double> u_s_vt_inplace;
  dots( u_s_vt_inplace, u, s, vt);
  Matrix<double> u_s_vt_noninplace;
  dots( u_s_vt_noninplace, svds.at(0), svds.at(1), svds.at(2));
  ///output error
  cout<<u<<endl;
  cout<<s<<endl;
  cout<<vt<<endl;
  cout<<"error = "<<norm( A - u_s_vt_noninplace)<<endl;
  cout<<"error = "<<norm( A - u_s_vt_inplace)<<endl;

  return 0;
}
