#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  ///real, complex
  ///diagonal, nondiagonal

  ///test of different rol and col number for real non diagonal case
    ///declare real nondiagonal matrix
    double elem_R32[6] = {8.1, 2.1,
                         -0.3, 3.9,
                          0.1, 4.1};
    double elem_R23[6] = {8.1, -2.1, -3.1,
                         -1.9,  0.1,  4.1};
    double elem_R22[4] = {4.1, 2.1,
                          0.8, 1.1};
    double elem_R33[9] = {1.1, 2.1, 2.1,
                          3.9, 1.1, 4.1,
                          0.9, 0.1, 3.1};
    Matrix<double> R32_0( 3, 2, false);
    R32_0.setElem( elem_R32);
    Matrix<double> R32_1 = R32_0;
    Matrix<double> R23( 2, 3, false); R23.setElem( elem_R23);
    Matrix<double> R22( 2, 2, false); R22.setElem( elem_R22);
    Matrix<double> R33( 3, 3, false); R33.setElem( elem_R33);
    ///check if matrix declared successfully
    ///correct
    cout<<norm(R32_0)<<endl;
exit(0);
    cout<<"test of different rol and col number for real non diagonal case"<<endl;
    cout<<"R32_0"<<R32_0;
    cout<<"R32_1"<<R32_1;
    cout<<"R23"<<R23;
    cout<<"R22"<<R22;
    cout<<"R33"<<R33;
    ///multuply by a scalar
    cout<<"2.0*R32_0"<<2.0*R32_0;
    cout<<"-1.0*R32_1"<<(-1.0)*R32_1;
    ///check addition
    ///correct, lack of error msg
    cout<<"adding same matrix"<<R32_0+R32_0;
    cout<<"adding two matrix"<<R32_0+R32_1;
    //cout<<"adding matrix with different row num"<<R32_0+R22;
    //cout<<"adding matrix with different col num"<<R32_0+R33;
    //cout<<"adding matrix with different col and row num"<<R32_0+R23;

  ///test of addition of diagonal and nondiagonal cae
    ///declare real diagonal matrix
    double elem_diag2[2] = {2.12, 3.12};
    double elem_diag3[3] = {6.5, 8.2, 94.1};
    Matrix<double> RD2( 2, 2, true); RD2.setElem( elem_diag2);
    Matrix<double> RD3( 3, 3, true); RD3.setElem( elem_diag3);
    cout<<"test of addition of real diagonal and nondiagonal cae"<<endl;
    cout<<"RD2"<<RD2;
    cout<<"RD3"<<RD3;
    cout<<"RD2+R22"<<RD2+R22;
    cout<<"RD3+R33"<<RD3+R33;

  ///test of addition of addition of complex and real case;
    ///declare tensors
    std::complex<double> elem_C33[9] = 
    {std::complex<double>(8.1, 1.1), std::complex<double>( 2.1, 4.2), std::complex<double>( 3.2, 6.2),
     std::complex<double>( 2.1 , 1.2), std::complex<double>( 3.9, 2), std::complex<double>( 4.1, 1.1),
     std::complex<double>( 2.1 , 1.2), std::complex<double> ( 4.1, 1.2), std::complex<double>( 1.1, 2.1)};
    std::complex<double> elem_C22[4] = 
    {std::complex<double>(0.1,  1.1), std::complex<double>( 2.1, -1.2), 
     std::complex<double>(3.2, -0.2), std::complex<double>( 0.1 , 0.8)};
    std::complex<double> elem_CD2[2] = {std::complex<double>( 2.12, 4.21 ),std::complex<double>( 3.12, 21.23)};
    std::complex<double> elem_CD3[3] = {std::complex<double>( 6.5, 2.11 ), std::complex<double>( 8.2 , 2.13), std::complex<double>( 94.1, 21.2)};
    Matrix<std::complex<double>> C33(3, 3), C22( 2, 2), CD2( 2, 2, true), CD3( 3, 3, true);
    C33.setElem(elem_C33); C22.setElem(elem_C22); CD2.setElem(elem_CD2); CD3.setElem(elem_CD3);

    cout<<"C22"<<C22;
    cout<<"C33"<<C33;
    cout<<"CD2"<<CD2;
    cout<<"CD3"<<CD3;
    ///perform addition
    /*
    cout<<"C22+R22"<<C22+R22;
    cout<<"C33+R33"<<C33+R33;
    cout<<"CD2+RD2"<<CD2+RD2;
    cout<<"CD3+RD3"<<CD3+RD3;
    */

  ///can not add complex and real type
  return 0;
}
