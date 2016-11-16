#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_R[12] = {8.1,  2.1,  3.2,
                       2.1,  3.9,  4.1,
                       2.1,  4.1,  1.1,
                      20.1, 42.1, 10.1};

  double elem_diag1[4] = {2.12, 3.12, 4.12, 32.1};

  double elem_diag2[3] = {6.5, 8.2, 94.1};

  std::complex<double> elem_C[9] = 
  {std::complex<double>(8.1, 1.1), std::complex<double>( 2.1, 4.2), std::complex<double>( 3.2, 6.2),
   std::complex<double>( 2.1 , 1.2), std::complex<double>( 3.9, 2), std::complex<double>( 4.1, 1.1),
   std::complex<double>( 2.1 , 1.2), std::complex<double> ( 4.1, 1.2), std::complex<double>( 1.1, 2.1)};//,
  // std::complex<double>( 2.1, 1.2),std::complex<double>( 2.1, 3.2),std::complex<double>( 3.2, 1.1)};

  std::complex<double> celem_diag1[4] = {std::complex<double>( 2.12, 4.21 ),std::complex<double>( 3.12, 21.23),std::complex<double>( 4.12, 11.2) ,std::complex<double>( 32.1, 1.23)};

  std::complex<double> celem_diag2[3] = {std::complex<double>( 6.5, 2.11 ), std::complex<double>( 8.2 , 2.13), std::complex<double>( 94.1, 21.2)};

  Matrix<double> MR1(4, 3);
  MR1.setElem(elem_R);

  Matrix<double> MR2(3, 4);
  MR2.setElem(elem_R);

  Matrix<double> diag1(4, 4, true);
  diag1.setElem(elem_diag1);

  Matrix<double> diag2(3, 3, true);
  diag2.setElem(elem_diag2);

  cout << "======= Real =======\n";

  cout << MR1;

  cout << MR2;

  cout << diag1;

  cout << diag2;

  cout << "======= MR1 MR2 =======\n";

  cout << dot(MR1, diag1);

  exit(0);

  cout << "======= diag1 MR1 =======\n";

  cout << dot(diag1, MR1);

  cout << "======= diag2 MR1 =======\n";

  cout << dot(diag2, MR2);

  cout << "======= diag1 diag1 =======\n";

  cout << dot(diag1, diag1);

  cout << "======= Complex =======\n";

  Matrix< std::complex<double> > MC1(3, 3);
  MC1.setElem(elem_C);

  Matrix< std::complex<double> > MC2(3, 4);
  MC2.setElem(elem_C);

  Matrix< std::complex<double> > cdiag1(4, 4, true);
  cdiag1.setElem(celem_diag1);

  Matrix< std::complex<double> > cdiag2(3, 3, true);
  cdiag2.setElem(celem_diag2);

  cout << "======= Real =======\n";

  cout << MC1;

  cout << MC2;

  cout << cdiag1;

  cout << cdiag2;

  cout << "======= MR1 MR2 =======\n";

  cout << dot(MC1, MC2);

  cout << "======= diag1 MR1 =======\n";

  cout << dot(cdiag1, MC1);

  cout << "======= diag2 MR1 =======\n";

  cout << dot(cdiag2, MC2);

  cout << "======= diag1 diag1 =======\n";

  cout << dot(cdiag1, cdiag1);

  cout << "======= Complex =======\n";

  return 0;
}
