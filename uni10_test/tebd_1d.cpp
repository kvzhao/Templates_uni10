#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

Matrix<double> sigma_x(){
  double mat_elem[4] = {\
      0,  1.0,\
    1.0,    0};
  Matrix<double> output(2, 2);
  output.setElem( mat_elem);
  return output;
}

Matrix<double> sigma_z(){
  double mat_elem[4] = {\
    1.0,    0,\
      0, -1.0};
  Matrix<double> output(2, 2);
  output.setElem( mat_elem);
  return output;
}

UniTensor<double> oneD_Ising( double h){
  Matrix<double> sx = sigma_x();
  Matrix<double> sz = sigma_z();
  Matrix<double> I( sx.row(), sx.col(), true);
  identity(I);

  Matrix<double> ham = otimes( sz, sz) + (0.5)*h*( otimes( sx, I) + otimes( I, sx));
  Bond bdi = Bond(BD_IN, 2);
  Bond bdo = Bond(BD_OUT, 2);
  vector<Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);

  UniTensor<double> H( bonds);
  H.putBlock(ham);
  return H;
}

void initialize_gamma_lambda( vector<UniTensor<double> > &gamma_tensors, vector<Matrix<double> > &lambda_matrix, const int Dbond, const int dbond){
  gamma_tensors.clear(); lambda_matrix.clear();

  vector<Bond> gamma_bds(3);
               gamma_bds.at(0) = Bond( BD_IN, Dbond);
               gamma_bds.at(1) = Bond( BD_IN, Dbond);
               gamma_bds.at(2) = Bond( BD_OUT, dbond);
  gamma_tensors = vector<UniTensor<double> > ( 2, UniTensor<double>( gamma_bds));
  lambda_matrix = vector<Matrix<double> > ( 2, Matrix<double>( Dbond, Dbond, true));

  ///create random matrix
  Matrix<double> temp = gamma_tensors.at(0).getBlock();
  uni10_rand( temp, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<double> temp1(Dbond, Dbond, true);
  uni10_rand( temp1, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);

  for (int i=0; i!=2; i++){
    gamma_tensors.at(i).putBlock(temp);
    lambda_matrix[i] = temp1;
  }
}

void bondcat(UniTensor<double>& T, const Matrix<double>& L, int bidx){
  int inBondNum = T.inBondNum();
	vector<int> labels = T.label();
  vector<int> per_labels = labels;
  int l = labels[bidx];
  per_labels.erase(per_labels.begin() + bidx);
	per_labels.insert(per_labels.begin(), l);

  UniTensor<double> T_c = permute(T, per_labels, 1);
  T_c.putBlock((dot(L, T_c.getBlock())));
  T = permute( T_c, labels, inBondNum);
}

void bondrm(UniTensor<double>& T, const Matrix<double>& L, int bidx){
	Matrix<double> invL = L;
  for(unsigned long i=0; i!=L.col(); i++){
    invL[i] = invL[i] == 0 ? 0 : (1 / invL[i]);
	}
	bondcat(T, invL, bidx);
}

/*
int tebd_two_site( UniTensor &gate, vector<UniTensor> &gamma_tensors, vector<Matrix> &lambda_matrix, const int i_start, const bool if_nonunitary){
  int Dbond = lambda_matrix.at(0).col();
  int dbond = gamma_tensors.at(0).bond(2).dim();
  int i_l = i_start, i_r = (i_start+1)%2;
  
  bondcat( gamma_tensors.at(i_l), lambda_matrix.at(i_r), 0);
  bondcat( gamma_tensors.at(i_r), lambda_matrix.at(i_r), 1);
  UniTensor theta; contract_theta( gamma_tensors, lambda_matrix, i_l, theta); 
  theta.setLabel( {1, 2, 3, 4}); gate.setLabel( {3, 4, -3, -4});
  theta = contract( theta, gate);///now theta label { 1, 2; -3, -4}
  
  if (if_nonunitary){
    theta.combineBond( {-3, -4});
    orthogonalization_one_site( theta, lambda_matrix.at(i_r));
    seperBond( theta, 2, dbond);
    theta.setLabel( { 1, 2, -3, -4});
  }

  theta.permute( {1, -3, 2, -4}, 2);
  svd_update( theta, i_l, gamma_tensors, lambda_matrix);
  bondrm( gamma_tensors.at(i_l), lambda_matrix.at(i_r), 0);
  bondrm( gamma_tensors.at(i_r), lambda_matrix.at(i_r), 1);
}
*/

int main(){
  const int Dbond = 5;
  const int dbond = 2;
  const int N_evolution = 1000000;
  const double tau = 1.0e-4;
  const double accuracy = 1.0e-14;
  UniTensor<double> hamiltonian = oneD_Ising( 1.05 );
  UniTensor<double> gate = UniTensor<double>( hamiltonian.bond());
  gate.putBlock( exph( -1.0*tau, hamiltonian.getBlock()));
  hamiltonian.setLabel({4, 5, 6, 7});
  gate.setLabel({4, 5, 6, 7});

  ///declare gamma and lambda
  vector<UniTensor<double> > gamma_tensors;
  vector<Matrix<double> > lambda_matrix;
  initialize_gamma_lambda( gamma_tensors, lambda_matrix, Dbond, dbond);

  for (int i=0; i!=N_evolution; i++){
    int i_l = i%2, i_r = (i+1)%2;

    ///contract theta
    bondcat( gamma_tensors.at(i_l), lambda_matrix.at(i_r), 0);
    bondcat( gamma_tensors.at(i_l), lambda_matrix.at(i_l), 1);
    bondcat( gamma_tensors.at(i_r), lambda_matrix.at(i_r), 1);
    gamma_tensors.at(i_l).setLabel({1, 2, 4});
    gamma_tensors.at(i_r).setLabel({2, 3, 5});
    UniTensor<double> theta; 
    theta = contract( gamma_tensors.at(i_l), gamma_tensors.at(i_r), false );
    theta = contract( theta, gate, false ); ///now theta label 1, 3; 6, 7
    theta = permute( theta, {1, 6, 3, 7}, 2 );

    ///svd and update
    vector<Matrix<double> > usv = svd(theta.getBlock());
    resize( usv.at(0),  usv.at(0).row(),           Dbond);
    resize( usv.at(1),            Dbond,           Dbond);
    resize( usv.at(2),            Dbond, usv.at(2).col());
    lambda_matrix.at(i_l) = usv.at(1);
    lambda_matrix.at(i_l) *= 1.0/norm( lambda_matrix.at(i_l));

    gamma_tensors.at(i_l) = permute(gamma_tensors.at(i_l), {1, 4, 2}, 2); 
    gamma_tensors.at(i_r) = permute(gamma_tensors.at(i_r), {2, 3, 5}, 1);
    gamma_tensors.at(i_l).putBlock( usv.at(0)); 
    gamma_tensors.at(i_r).putBlock( usv.at(2));
    gamma_tensors.at(i_l) = permute( gamma_tensors.at(i_l), {1, 2, 4}, 2); 
    gamma_tensors.at(i_r) = permute( gamma_tensors.at(i_r), {2, 3, 5}, 2);

    bondrm( gamma_tensors.at(i_l), lambda_matrix.at(i_r), 0);
    bondrm( gamma_tensors.at(i_r), lambda_matrix.at(i_r), 1);

    ///measure
    if (i%100==0){
      vector<UniTensor<double> > gamma_now = gamma_tensors;
      bondcat( gamma_now.at(i_l), lambda_matrix.at(i_r), 0);
      bondcat( gamma_now.at(i_l), lambda_matrix.at(i_l), 1);
      bondcat( gamma_now.at(i_r), lambda_matrix.at(i_r), 1);
      theta = contract( gamma_now.at(0), gamma_now.at(1), false); ///now theta label 1, 4; 3, 5
      UniTensor<double> theta_T =  theta;
      UniTensor<double> theta_H = contract( theta, hamiltonian, false ); ///now theta_H label 1, 3; 6, 7
      theta_H.setLabel( {1, 3, 4, 5} );
      //UniTensor<double> norm = contract( theta, theta_T );
      //UniTensor<double> expec = contract( theta, theta_H );
      double norm  = contract( theta, theta_T, false )[0];
      double expec = contract( theta, theta_H, false )[0];
      //cout<<contract( theta, theta_T, false );
      //cout<<contract( theta, theta_H, false );
      printf( "%22.14f\n",  expec/norm );
    }
    else{}
  }
  return 0;
}
