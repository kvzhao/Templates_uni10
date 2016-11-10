#ifndef __UNI10_ELEM_CPU_H__
#define __UNI10_ELEM_CPU_H__

#include <iostream>

#include "uni10/uni10_tools.h"

namespace uni10{

  
  template<typename uni10_type>
    class Block;
  template<typename uni10_type>
    class Matrix;

  template<typename uni10_type>
    class uni10_elem_lapack_cpu{
      
      public:

        explicit uni10_elem_lapack_cpu();

        explicit uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum);

        explicit uni10_elem_lapack_cpu(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum);

        explicit uni10_elem_lapack_cpu(const uni10_elem_lapack_cpu& _elem);

        ~uni10_elem_lapack_cpu();

        inline bool empty() const{ return elem == NULL; };

        inline uni10_type_id typeID() const { return __uni10_id; };

        inline uni10_type* getElem() const { return elem; };

        inline uni10_uint64 elemNum() const { return __elemNum; };

        void setElem(const uni10_type* src, bool src_ongpu = false);

        void assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum);

        void print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag) const;

        friend class Block<uni10_type>;

        friend class Matrix<uni10_type>;

      private:

        uni10_type_id __uni10_id;  

        uni10_uint64 __elemNum;

        uni10_type* elem;

        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, const uni10_type* src=NULL);

    };

}

#endif
