// Have to write something.
#ifndef __UNI10_ELEM_MAGMA_GPU_H__
#define __UNI10_ELEM_MAGMA_GPU_H__

#include <iostream>

#include "uni10/uni10_tools.h"

namespace uni10{
  
  template<typename uni10_type>

    class uni10_elem_gpu{
      
      public:

        uni10_elem_gpu();

        uni10_elem_gpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag = false);

        uni10_elem_gpu(uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag = false);

        uni10_elem_gpu(const uni10_elem_gpu& _elem);

        ~uni10_elem_gpu();

        inline uni10_type* getElem() const { /*Developping*/ };

        inline uni10_uint64 elemNum() const { return __elemNum; };

        inline bool isdiag() const { return __isdiag; };

        void setElem(uni10_type* src);

        void print_elem_gpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum) const;

      private:

        bool __isdiag;

        uni10_uint64 __elemNum;

        uni10_type* elem;

        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag, uni10_type* src=NULL);

    };

}

#endif
