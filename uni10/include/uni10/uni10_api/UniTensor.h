#ifndef __UNI10_UNITENSOR_H__
#define __UNI10_UNITENSOR_H__

#include <string>

#include "uni10/uni10_api/UniTensor_para.h"
#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  enum contain_type{
    no_sym    = 0,
    blk_sym   = 1,
    spar_sym  = 2 
  };

  struct nsy_paras;

  template<typename uni10_type>
    class UniTensor;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const UniTensor<uni10_type>& _b);

  template <typename uni10_type>
    class UniTensor {

      private:

        contain_type style;
        struct U_para<uni10_type>* paras;

        void init_para();
        void meta_link();
        void init();
        void free_para();
        //void initBlocks();
        //uni10_uint64 grouping();

        // General variables.
        std::string* name;
        std::vector<Bond>* bonds;
        std::vector<uni10_int32>* labels;
        uni10_int32*  RBondNum;       //Row bond number
        uni10_uint64* Rdim;
        uni10_uint64* Cdim;
        uni10_uint64* U_elemNum;
        std::map< Qnum, Block<uni10_type> >* blocks;
        UELEM(uni10_elem, _package, _type)<uni10_type>* U_elem;     // pointer to a real matrix
        uni10_int32*  status;     //Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements

        static uni10_int32 COUNTER;
        static uni10_uint64 ELEMNUM;
        static uni10_uint64 MAXELEMNUM;
        static uni10_uint64 MAXELEMTEN;               //Max number of element of a tensor

        static const uni10_int32 HAVEBOND = 1;        /**< A flag for initialization */
        static const uni10_int32 HAVEELEM = 2;        /**< A flag for having element assigned */

      public:

        //static uni10_int32  GET_COUNTER(){return COUNTER;}
        //static uni10_uint64 GET_ELEMNUM(){return ELEMNUM;}
        //static uni10_uint64 GET_MAXELEMNUM(){return MAXELEMNUM;};
        //static uni10_uint64 GET_MAXELEMTEN(){return MAXELEMTEN;};

        static uni10_int32 GET_HAVEBOND(){return HAVEBOND;}; 
        static uni10_int32 GET_HAVEELEM(){return HAVEELEM;};      /**< A flag for having element assigned */

        explicit UniTensor();
        explicit UniTensor(uni10_type val, contain_type _style = no_sym);
        explicit UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "", contain_type _style = no_sym);
        explicit UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "", contain_type _style = no_sym);
        explicit UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "", contain_type _style = no_sym);
        explicit UniTensor(const UniTensor& UniT);
        //explicit UniTensor(const std::string& fname);
        explicit UniTensor(const Block<uni10_type>& UniT);
        ~UniTensor();

        void putBlock(const Block<uni10_type>& mat);

        void putBlock(const Qnum& qnum, const Block<uni10_type>& mat);
        
        std::vector<int> label()const{return *labels;};

        uni10_int32 label(size_t idx)const{return (*labels)[idx];};

        std::string getName() const{return *name;};

        void setName(const std::string& _name);

        const std::map<Qnum, Block<uni10_type> >& const_getBlocks()const{return *blocks;};

        uni10_uint64 bondNum()const{return bonds->size();};

        uni10_uint64 inBondNum()const{return *RBondNum;};

        std::vector<Bond> bond()const{return *bonds;};

        Bond bond(uni10_uint64 idx)const{return (*bonds)[idx];};

        uni10_uint64 elemNum()const{return (*Rdim)*(*Cdim); };

        uni10_uint64 blockNum()const{return blocks->size();};

        const Block<uni10_type>& const_getBlock()const;

        const Block<uni10_type>& const_getBlock(const Qnum& qnum)const;

        void setLabel(const uni10_int32 newLabel, const uni10_uint64 idx);

        void setLabel(const std::vector<uni10_int32>& newLabels);

        void setLabel(uni10_int32* newLabels);

        std::vector<Qnum> blockQnum()const;

        Qnum blockQnum(uni10_uint64 idx)const;

        std::map< Qnum, Matrix<uni10_type> > getBlocks()const;

        Matrix<uni10_type> getBlock(bool diag = false)const;

        Matrix<uni10_type> getBlock(const Qnum& qnum, bool diag = false)const;

        void setElem(const uni10_type* _elem);

        void setElem(const std::vector<uni10_type>& _elem);

        void printDiagram()const;

        friend std::ostream& operator<< <>(std::ostream& os, const UniTensor& _b);  // --> uni10_elem().print_elem()

        //void set_zero();
        //void set_zero(const Qnum& qnum);
        //void identity();
        //void identity(const Qnum& qnum);
        //void randomize();
        //void orthoRand();
        //void orthoRand(const Qnum& qnum);
        //void save(const std::string& fname) const;
        //UniTensor& transpose();
        //UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
        //UniTensor& permute(int* newLabels, int inBondNum);
        //UniTensor& permute(int inBondNum);
        //friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
        //friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
        //UniTensor& combineBond(const std::vector<int>& combined_labels);
        //std::string printRawElem(bool print=true)const;
        //static std::string profile(bool print = true);
        //std::vector<_Swap> exSwap(const UniTensor& Tb)const;
        //void addGate(const std::vector<_Swap>& swaps);
        //Complex trace()const;
        //UniTensor& partialTrace(int la, int lb);
        //Matrix getRawElem()const;
        //UniTensor& assign(const std::vector<Bond>& _bond);
        //UniTensor& operator*= (Real a);
        //UniTensor& operator*= (Complex a);
        //UniTensor& operator*= (const UniTensor& Tb);
        //friend UniTensor operator*(const UniTensor& Ta, const Complex& a);
        //friend UniTensor operator*(const Complex& a, const UniTensor& Ta);
        //friend UniTensor operator*(const UniTensor& Ta, Real a);
        //friend UniTensor operator*(Real a, const UniTensor& Ta);
        //friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
        //UniTensor& operator=(const UniTensor& UniT);
        //UniTensor& operator+= (const UniTensor& Tb);
        //friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
        //void setRawElem(const Block& blk);

        //bool similar(const UniTensor& Tb)const;
        //bool elemCmp(const UniTensor& UniT)const;

        //void clear();
        //Real operator[](size_t idx) const;
        //Complex operator()(size_t idx) const;

        //friend class Node;
        //friend class Network;

    };

  template <typename uni10_type>
    uni10_int32 UniTensor<uni10_type>::COUNTER = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::ELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMTEN = 0;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const UniTensor<uni10_type>& UniT){
      std::cout << "===============================" << std::endl;
      if(!(*(UniT.status) & UniT.HAVEBOND)){
        if(UniT.U_elem->__ongpu)
          std::cout<<"\nScalar: " << UniT.U_elem->__elem[0]<<", onGPU";
        else
          std::cout<<"\nScalar: " << UniT.U_elem->__elem[0];
        std::cout<<"\n\n";
        return os;
      }

      uni10_uint64 row = 0;
      uni10_uint64 col = 0;

      std::vector<Bond>bonds = (*UniT.bonds);
      for(uni10_uint64 i = 0; i < bonds.size(); i++)
        if(bonds[i].type() == BD_IN)
          row++;
        else
          col++;
      uni10_uint64 layer = std::max(row, col);
      uni10_uint64 nmlen = (*UniT.name).length() + 2;
      uni10_uint64 star = 12 + (14 - nmlen) / 2;
      for(uni10_uint64 s = 0; s < star; s++)
        std::cout << "*";
      if((*UniT.name).length() > 0)
        std::cout << " " << (*UniT.name) << " ";
      for(uni10_uint64 s = 0; s < star; s++)
        std::cout <<"*";
      std::cout<<std::endl;

      if(UniT.U_elem->__uni10_typeid == 1)
        std::cout << "REAL" << std::endl;
      else if(UniT.U_elem->__uni10_typeid == 2)
        std::cout << "COMPLEX" << std::endl;

      if(UniT.U_elem->__ongpu)
        std::cout<<"\n                 onGPU";
      std::cout << "\n             ____________\n";
      std::cout << "            |            |\n";
      uni10_uint64 llab = 0;
      uni10_uint64 rlab = 0;
      char buf[128];
      for(uni10_uint64 l = 0; l < layer; l++){
        if(l < row && l < col){
          llab = (*UniT.labels)[l];
          rlab = (*UniT.labels)[row + l];
          sprintf(buf, "    %5ld___|%-4d    %4d|___%-5ld\n", llab, bonds[l].dim(), bonds[row + l].dim(), rlab);
          std::cout<<buf;
        }
        else if(l < row){
          llab = (*UniT.labels)[l];
          sprintf(buf, "    %5d___|%-4d    %4s|\n", llab, bonds[l].dim(), "");
          std::cout<<buf;
        }
        else if(l < col){
          rlab = (*UniT.labels)[row + l];
          sprintf(buf, "    %5s   |%4s    %4d|___%-5ld\n", "", "", bonds[row + l].dim(), rlab);
          std::cout << buf;
        }
        std::cout << "            |            |   \n";
      }
      std::cout << "            |____________|\n";

      std::cout << "\n================BONDS===============\n";
      for(uni10_uint64 b = 0; b < bonds.size(); b++)
        std::cout << bonds[b];

      std::cout<<"\n===============BLOCKS===============\n";
      uni10_bool printElem = true;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = UniT.blocks->begin();
      for (; it != UniT.blocks->end(); it++ ){
        std::cout << "--- " << it->first << ": ";// << Rnum << " x " << Cnum << " = " << Rnum * Cnum << " ---\n\n";
        if(!((*(UniT.status) & UniT.HAVEELEM) && printElem))
          std::cout<<it->second;
        else
          std::cout<<it->second.row() << " x "<<it->second.col()<<": "<<it->second.elemNum()<<std::endl<<std::endl;
      }
      std::cout << "Total elemNum: "<<(*UniT.U_elemNum)<<std::endl;
      std::cout << "===============================" << std::endl;
      os << "\n";
      return os;
    }

};

//UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast = false);
//UniTensor contract(rflag tp, UniTensor& Ta, UniTensor& Tb, bool fast = false);
//UniTensor contract(cflag tp, UniTensor& Ta, UniTensor& Tb, bool fast = false);
//UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
//UniTensor otimes(rflag tp, const UniTensor& Ta, const UniTensor& Tb);
//UniTensor otimes(cflag tp, const UniTensor& Ta, const UniTensor& Tb);

#endif
