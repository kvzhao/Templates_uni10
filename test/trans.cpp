#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_R[12] = {8.1,  2.1,  3.2,
                       2.1,  3.9,  4.1,
                       2.1,  4.1,  1.1,
                      20.1, 42.1, 10.1};

  std::complex<double> elem_C[12] = 
  {std::complex<double>(8.1, 10.1), std::complex<double>( 2.1, 4.2), std::complex<double>( 3.2, 6.2),
   std::complex<double>( 2.1 , 10.2), std::complex<double>( 3.9, 23.1), std::complex<double>( 4.1, 1.1),
   std::complex<double>( 2.1 , 13.2), std::complex<double> ( 4.1, 12.2), std::complex<double>( 1.1, 2.1),
   std::complex<double>( 20.1, 12.2),std::complex<double>( 42.1, 3.2),std::complex<double>( 3.2, 10.1)};

  Matrix<double> MR(4, 3);
  MR.setElem(elem_R);

  Matrix< std::complex<double> > MC(4, 3);
  MC.setElem(elem_C);
  
  fprintf(stdout, "-----    Real    -----\n");

  fprintf(stdout, "-----    MR_ori   -----\n");

  cout << MR;

  fprintf(stdout, "-----  transpose -----\n");

  cout << transpose( MR );

  fprintf(stdout, "-----    dagger  -----\n");

  cout << dagger( MR );

  fprintf(stdout, "-----    conj    -----\n");

  cout << conj( MR );

  fprintf(stdout, "-----   Complex  -----\n");

  fprintf(stdout, "-----    MC_ori   -----\n");

  cout << MC;

  fprintf(stdout, "-----  transpose -----\n");

  cout << transpose( MC );

  fprintf(stdout, "-----    dagger  -----\n");

  cout << dagger( MC );

  fprintf(stdout, "-----    conj    -----\n");

  cout << conj( MC );

  fprintf(stdout, "-----    end    -----\n");

  cout << norm( MR ) << endl;

  cout << norm( MC ) << endl;

  cout << uni10::sum( MR ) << endl;

  cout << uni10::sum( MC ) << endl;

  return 0;
}
