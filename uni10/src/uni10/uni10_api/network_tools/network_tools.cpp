#include "uni10/uni10_api/network_tools/network_tools.h"

namespace uni10{

  std::vector<_Swap> recSwap(std::vector<int>& _ord) { //Given the reshape order out to in.
    //int ordF[n];
    int n = _ord.size();
    std::vector<int> ordF(n);
    for(int i = 0; i < n; i++)
      ordF[i] = i;
    return recSwap(_ord, ordF);
  }
  std::vector<_Swap> recSwap(std::vector<int>& _ord, std::vector<int>& ordF) { //Given the reshape order out to in.
    int n = _ord.size();
    std::vector<int> ord = _ord;
    std::vector<_Swap> swaps;
    _Swap sg;
    int tmp;
    for(int i = 0; i < n - 1; i++)
      for(int j = 0; j < n - i - 1; j++)
        if(ord[j] > ord[j + 1]) {
          sg.b1 = ordF[ord[j + 1]];
          sg.b2 = ordF[ord[j]];
          tmp = ord[j];
          ord[j] = ord[j + 1];
          ord[j + 1] = tmp;
          swaps.push_back(sg);
        }
    return swaps;
  }

};
