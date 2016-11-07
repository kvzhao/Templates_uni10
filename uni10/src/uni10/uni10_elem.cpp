#include "uni10/uni10_type.h"
#include "uni10/uni10_elem.h"

namespace uni10{

  template<typename uni10_type>
    uni10_elem<uni10_type>::uni10_elem(): __isdiag(false), __elemNum(0), elem(NULL){};

  template<typename uni10_type>
    uni10_elem<uni10_type>::uni10_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag): elem(NULL){

#ifdef CPU 
      init_cpu(_Rnum, _Cnum, _isdiag, NULL);
#endif

    }

  template<typename uni10_type>
    uni10_elem<uni10_type>::uni10_elem(uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag): elem(NULL){

#ifdef CPU 
      init_cpu(_Rnum, _Cnum, _isdiag, src);
#endif

    };

  template<typename uni10_type>
    uni10_elem<uni10_type>::uni10_elem(const uni10_elem& _elem): __isdiag(_elem.isdiag()), __elemNum(_elem.elemNum()){

#ifdef CPU 
      init_cpu(1, __elemNum, false, _elem.elem);
#endif

    };

  template<typename uni10_type>
    uni10_elem<uni10_type>::~uni10_elem(){

#ifdef CPU 
      uni10_elem_free(elem, __elemNum * sizeof(uni10_type));
#endif

    };

  template<typename uni10_type>
    void uni10_elem<uni10_type>::setElem(uni10_type* src){

      uni10_error_msg( elem == NULL, "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      uni10_error_msg( src  == NULL, "The source ptr is NULL.");

#ifdef CPU 
      uni10_elem_copy( elem, src, __elemNum * sizeof(uni10_type) );
#endif

    };

  template<typename uni10_type>
    void uni10_elem<uni10_type>::init_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag, uni10_type* src){

      __isdiag = _isdiag;

      __elemNum = __isdiag ? fmax(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        elem = (uni10_type*)uni10_elem_alloc( memsize );

        if(src != NULL)
          uni10_elem_copy( elem, src, memsize );
        else
          uni10_elemBzero( elem, memsize );

      }

    };

  template<typename uni10_type>
    void uni10_elem<uni10_type>::print_elem_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum) const{

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

  template class uni10_elem<uni10_double64>;
  template class uni10_elem<uni10_complex128>;

} /* namespace uni10 */
