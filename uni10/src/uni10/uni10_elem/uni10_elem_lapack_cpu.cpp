#include "uni10/uni10_type.h"
#include "uni10/uni10_elem/uni10_elem_lapack_cpu.h"

namespace uni10{

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(): __uni10_id(UNI10_TYPE_ID(uni10_type)), __elemNum(0), elem(NULL){};

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum): __uni10_id(UNI10_TYPE_ID(uni10_type)), elem(NULL){

      init(_Rnum, _Cnum, NULL);

    }

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum): __uni10_id(UNI10_TYPE_ID(uni10_type)), elem(NULL){

      init(_Rnum, _Cnum, src);

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(const uni10_elem_lapack_cpu& _elem): __uni10_id(_elem.__uni10_id), __elemNum(_elem.__elemNum){

      init(1, __elemNum, _elem.elem);

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::~uni10_elem_lapack_cpu(){

      uni10_elem_free_cpu(elem, __elemNum * sizeof(uni10_type));

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::setElem(const uni10_type* src, bool src_ongpu){

      uni10_error_msg( src_ongpu, " The source pointer is on the device. Please install MAGMA or CUDAONLY gpu version instead.");
      uni10_error_msg( elem == NULL, "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      uni10_error_msg( src  == NULL, "The source ptr is NULL.");

      uni10_elem_copy_cpu( elem, src, __elemNum * sizeof(uni10_type) );

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, const uni10_type* src){

      __elemNum =  _Rnum * _Cnum ;

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
    void uni10_elem_lapack_cpu<uni10_type>::assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum){

      __elemNum = _Rnum * _Cnum;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        elem = (uni10_type*)uni10_elem_alloc_cpu( memsize );
        uni10_elemBzero_cpu( elem, memsize );

      }

    }


  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag) const{

      if ( _Rnum == 1 ) {
        fprintf(stdout, "[ " );
      }
      else {
        fprintf(stdout, "[\n" );
      }

      if(_isdiag){

        for( int i = 0; i < (int)_Rnum; ++i ) {
          for( int j = 0; j < (int)_Cnum; ++j ) {
            if ( i != j) {
              if(__uni10_id == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(elem[ i ]);
            }
          }
          if ( _Rnum > 1 ) 
            fprintf(stdout, "\n" );
          else 
            fprintf(stdout, " " );
        }
        fprintf(stdout, "];\n" );

      }
      else{

        for( int i = 0; i < (int)_Rnum; ++i ) {
          for( int j = 0; j < (int)_Cnum; ++j ) {
            if ( elem[ i * _Cnum + j] == 0.) {
              if(typeID() == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(elem[ i * _Cnum + j ]);
            }
          }
          if ( _Rnum > 1 ) 
            fprintf(stdout, "\n" );
          else 
            fprintf(stdout, " " );
        }
        fprintf(stdout, "];\n" );

      }

    }

  template class uni10_elem_lapack_cpu<uni10_double64>;
  template class uni10_elem_lapack_cpu<uni10_complex128>;

} /* namespace uni10 */
