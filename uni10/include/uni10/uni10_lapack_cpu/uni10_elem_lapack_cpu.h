#ifndef __UNI10_ELEM_LAPACK_CPU_H__
#define __UNI10_ELEM_LAPACK_CPU_H__

#include <iostream>

#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_cpu.h"

namespace uni10{

  template<typename uni10_type>
    class uni10_elem_lapack_cpu{
      
      public:

        explicit uni10_elem_lapack_cpu();

        explicit uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false, uni10_bool _ongpu = false);

        explicit uni10_elem_lapack_cpu(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false, uni10_bool _ongpu = false);

        explicit uni10_elem_lapack_cpu(const uni10_elem_lapack_cpu& _elem);

        
        uni10_elem_lapack_cpu& operator=(const uni10_elem_lapack_cpu& _m){
          //std::cout << "QQ!!! QQ!!QQ!!!!!\n\n\n";
          __uni10_typeid = _m.__uni10_typeid;
          __ongpu        = _m.__ongpu;
          this->init(1, _m.__elemNum, false, _m.__elem);
          return *this;
        }

        ~uni10_elem_lapack_cpu();

        inline bool empty() const{ return __elem == NULL; }

        void set_zeros();

        void setElem(const uni10_type* src, bool src_ongpu = false);

        void assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum);

        void clear();

        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src=NULL);

        void copy(uni10_uint64 begin_idx, const uni10_elem_lapack_cpu<uni10_type>& src, uni10_uint64 len);

        void resize(uni10_uint64 _row, uni10_uint64 _col, uni10_uint64& _Rnum, uni10_uint64& _Cnum, uni10_bool& _isdiag, uni10_const_bool& _fixHead = true);

        void print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag) const;

        uni10_type_id __uni10_typeid;  

        uni10_bool __ongpu;

        uni10_uint64 __elemNum;

        uni10_type* __elem;

    };

}

#endif
