#include "uni10/uni10_type.h"
#include "uni10/uni10_elem/uni10_elem_cpu.h"

namespace uni10{

  template<typename uni10_type>
    uni10_elem_cpu<uni10_type>::uni10_elem_cpu(): __isdiag(false), __elemNum(0), elem(NULL){};

  template<typename uni10_type>
    uni10_elem_cpu<uni10_type>::uni10_elem_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag): elem(NULL){

      init(_Rnum, _Cnum, _isdiag, NULL);

    }

  template<typename uni10_type>
    uni10_elem_cpu<uni10_type>::uni10_elem_cpu(uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag): elem(NULL){

      init(_Rnum, _Cnum, _isdiag, src);

    };

  template<typename uni10_type>
    uni10_elem_cpu<uni10_type>::uni10_elem_cpu(const uni10_elem_cpu& _elem): __isdiag(_elem.isdiag()), __elemNum(_elem.elemNum()){

      init(1, __elemNum, false, _elem.elem);

    };

  template<typename uni10_type>
    uni10_elem_cpu<uni10_type>::~uni10_elem_cpu(){

      uni10_elem_free_cpu(elem, __elemNum * sizeof(uni10_type));

    };

  template<typename uni10_type>
    void uni10_elem_cpu<uni10_type>::setElem(uni10_type* src){

      uni10_error_msg( elem == NULL, "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      uni10_error_msg( src  == NULL, "The source ptr is NULL.");

      uni10_elem_copy_cpu( elem, src, __elemNum * sizeof(uni10_type) );

    };

  template<typename uni10_type>
    void uni10_elem_cpu<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag, uni10_type* src){

      __isdiag = _isdiag;

      __elemNum = __isdiag ? fmax(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        elem = (uni10_type*)uni10_elem_alloc_cpu( memsize );

        if(src != NULL)
          uni10_elem_copy_cpu( elem, src, memsize );
        else
          uni10_elemBzero_cpu( elem, memsize );

      }

    };

  template<typename uni10_type>
    void uni10_elem_cpu<uni10_type>::print_elem_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum) const{

      for(int i = 0; i < (int)_Rnum; i++){

        for(int j = 0; j < (int)_Cnum; j++){

          if(__isdiag){

            if(i == j)
              std::cout << elem[i*_Cnum + j] << " ";
            else
              std::cout << 0.000 << " ";

          }else{

            std::cout << elem[i*_Cnum + j] << " ";

          }

        }

        std::cout << std::endl;

      }

    }

  template class uni10_elem_cpu<uni10_double64>;
  template class uni10_elem_cpu<uni10_complex128>;

} /* namespace uni10 */
