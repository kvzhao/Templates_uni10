#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type> 
    void set_zeros( UniTensor<uni10_type>& A ){

      tensor_tools::set_zero(A.paras, A.style);

    }

  template<typename uni10_type> 
    void randomize( UniTensor<uni10_type>& A ){

      tensor_tools::randomize(A.paras, A.style);

    }

  template<typename uni10_type>
    void permute( const UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UniTensor<uni10_type>& Tout, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg((*T.status) & UniTensor<uni10_type>::GET_HAVEBOND == 0, "%s", "There is no bond in the tensor(scalar) to permute.");
      uni10_error_msg((T.labels->size() == newLabels.size()) == 0, "%s", "The size of the input new labels does not match for the number of bonds.");

      uni10_int32 bondNum = T.bonds->size();
      std::vector<uni10_int32> rsp_outin(bondNum);
      uni10_int32 cnt = 0;
      for(uni10_int32 i = 0; i < bondNum; i++)
        for(uni10_int32 j = 0; j < bondNum; j++)
          if((*T.labels)[i] == newLabels[j]){
            rsp_outin[j] = i;
            cnt++;
          }
      uni10_error_msg((cnt == newLabels.size()) == 0, "%s", "The input new labels do not 1-1 correspond to the labels of the tensor.");

      uni10_bool inorder = true;

      for(uni10_int32 i = 1; i < bondNum; i++)
        if(rsp_outin[i] != i){
          inorder = false;
          break;
        }
      if(inorder && (*T.RBondNum) == rowBondNum) {
        Tout = T;
        return ;
      }   //do nothing
      else{
        std::vector<Bond> outBonds;
        for(uni10_int32 b = 0; b < T.bonds->size(); b++){
          outBonds.push_back((*T.bonds)[rsp_outin[b]]);
        }
        for(uni10_uint64 b = 0; b < T.bonds->size(); b++){
          if(b < (uni10_uint64)rowBondNum)
            outBonds[b].change(BD_IN);
          else
            outBonds[b].change(BD_OUT);
        }
        Tout.assign(outBonds);
        Tout.setName(*T.name);
        // ON CPU 
        // ON GPU Developping
        if((*T.status) & UniTensor<uni10_type>::GET_HAVEELEM())
          tensor_tools::permute(T.paras, T.style, rsp_outin, Tout.paras, inorder);

        *Tout.status |= UniTensor<uni10_type>::GET_HAVEELEM();
        Tout.setLabel(newLabels);

      }
    }

  template<typename uni10_type>
    void permute( const UniTensor<uni10_type>& T, int* newLabels, int rowBondNum, UniTensor<uni10_type>& Tout, UNI10_INPLACE on){

      std::vector<uni10_int32> _labels(newLabels, newLabels + T.bond().size());
      permute(T, _labels, rowBondNum, Tout, on);

    }

  template<typename uni10_type>
    void permute( const UniTensor<uni10_type>& T, int rowBondNum, UniTensor<uni10_type>& Tout, UNI10_INPLACE on){

      std::vector<uni10_int32> ori_labels = T.label();
      permute(T, ori_labels, rowBondNum, Tout, on);

    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on){
      UniTensor<uni10_type> Tout;
      permute(T, newLabels, rowBondNum, Tout, on);
      T = Tout;
    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, int* newLabels, int rowBondNum, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      permute(T, newLabels, rowBondNum, Tout, on);
      T = Tout;

    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, int rowBondNum, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      permute(T, rowBondNum, Tout, on);
      T = Tout;

    }

  template<typename uni10_type>
    void contract(UniTensor<uni10_type>& Ta, UniTensor<uni10_type>& Tb, UniTensor<uni10_type>& Tc, uni10_uint64 fast, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(!((*Ta.status) & (*Tb.status) & Ta.HAVEELEM), "%s" ,"Cannot perform contraction of two tensors before setting their elements.");

      if(&Ta == &Tb){
        UniTensor<uni10_type> Ttmp(Tb);
        contract(Ta, Ttmp, Tc, fast, on);
      }

      if((*Ta.status) & Ta.HAVEBOND && (*Tb.status) & Ta.HAVEBOND){
        int AbondNum = Ta.bonds->size();
        int BbondNum = Tb.bonds->size();
        std::vector<int> oldLabelA = (*Ta.labels);
        std::vector<int> oldLabelB = (*Tb.labels);
        int oldRnumA = (*Ta.RBondNum);
        int oldRnumB = (*Tb.RBondNum);
        std::vector<int> newLabelA;
        std::vector<int> interLabel;
        std::vector<int> newLabelB;
        std::vector<int> markB(BbondNum, 0);
        std::vector<int> newLabelC;
        bool match;
        for(int a = 0; a < AbondNum; a++){
          match = false;
          for(int b = 0; b < BbondNum; b++)
            if((*Ta.labels)[a] == (*Tb.labels)[b]){
              markB[b] = 1;
              interLabel.push_back((*Ta.labels)[a]);
              newLabelB.push_back((*Tb.labels)[b]);

              uni10_error_msg(!( (*Ta.bonds)[a].dim() == (*Tb.bonds)[b].dim() ), "%s", 
                  "Cannot contract two bonds having different dimensions");

              match = true;
              break;
            }
          if(!match){
            newLabelA.push_back((*Ta.labels)[a]);
            newLabelC.push_back((*Ta.labels)[a]);
          }
        }
        for(size_t a = 0; a < interLabel.size(); a++)
          newLabelA.push_back(interLabel[a]);
        for(int b = 0; b < BbondNum; b++)
          if(markB[b] == 0){
            newLabelB.push_back((*Tb.labels)[b]);
            newLabelC.push_back((*Tb.labels)[b]);
          }
        int conBond = interLabel.size();

        Ta = permute(Ta, newLabelA, AbondNum - conBond);
        Tb = permute(Tb, newLabelB, conBond);

        std::vector<Bond> cBonds;
        for(int i = 0; i < AbondNum - conBond; i++)
          cBonds.push_back((*Ta.bonds)[i]);
        for(int i = conBond; i < BbondNum; i++)
          cBonds.push_back((*Tb.bonds)[i]);

        Tc.assign(cBonds);
        if(cBonds.size())
          Tc.setLabel(newLabelC);
        Block<uni10_type> blockA, blockB, blockC;
        typename std::map<Qnum, Block<uni10_type> >::iterator it;
        typename std::map<Qnum, Block<uni10_type> >::iterator it2;
        for(it = Ta.blocks->begin() ; it != Ta.blocks->end(); it++){
          if((it2 = Tb.blocks->find(it->first)) != Tb.blocks->end()){
            blockA = it->second;
            blockB = it2->second;
            blockC = (*Tc.blocks)[it->first];
            uni10_error_msg(!(blockA.row() == blockC.row() && blockB.col() == blockC.col() && blockA.col() == blockB.row()), 
                "%s", "The dimensions the bonds to be contracted out are different.");
                       
            matrixDot(&blockA.elem_enforce(), &blockA.diag_enforce(), &blockB.elem_enforce(), &blockB.diag_enforce(), 
                &blockA.row_enforce(), &blockB.col_enforce(), &blockA.col_enforce(), &blockC.elem_enforce());
          }
        }
        (*Tc.status) |= Tc.HAVEELEM;

        if(conBond == 0){                     //Outer product
          int idx = 0;
          for(int i = 0; i < oldRnumA; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = 0; i < oldRnumB; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          for(int i = oldRnumA; i < AbondNum; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(int i = oldRnumB; i < BbondNum; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          Tc = permute(Tc, newLabelC, oldRnumA + oldRnumB);
        }

        if(!fast){
          Ta = permute(Ta, oldLabelA, oldRnumA);
          Tb = permute(Tb, oldLabelB, oldRnumB);
        }
        return ;
      }
      else if((*Ta.status) & Ta.HAVEBOND)
        Tc = Ta * Tb[0];
      else if((*Tb.status) & Tb.HAVEBOND)
        Tc = Ta[0] * Tb;
      else
        Tc = UniTensor<uni10_type>(Ta[0] * Tb[0]);

    }

  template<typename uni10_type>
    void otimes(const UniTensor<uni10_type>& Ta, const UniTensor<uni10_type>& Tb, UniTensor<uni10_type>& Tc, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      
      UniTensor<uni10_type> T1 = Ta;
      UniTensor<uni10_type> T2 = Tb;
      std::vector<int> label1(T1.bondNum());
      std::vector<int> label2(T2.bondNum());
      for(size_t i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(size_t i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);

      contract(T1, T2, Tc, true, on);

    }

  template<typename uni10_type>
    void otimes(const Block<uni10_type>& Ta, const Block<uni10_type>& Tb, Matrix<uni10_type>& Mc, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      UniTensor<uni10_type> T1(Ta);
      UniTensor<uni10_type> T2(Tb);
      std::vector<int> label1(T1.bondNum());
      std::vector<int> label2(T2.bondNum());
      for(size_t i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(size_t i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);

      Mc = contract(T1, T2, true).getBlock();

    }

};

#endif
