#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  Block<double> A;
  vector<Bond> bonds(3, Bond(BD_IN, 3));
  cout << bonds[0];
  cout << bonds[1];
  UniTensor_Nsy<double> B(bonds);
  cout << B;

  return 0;
}
