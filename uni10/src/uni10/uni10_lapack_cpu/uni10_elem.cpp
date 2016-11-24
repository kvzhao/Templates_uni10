#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"

namespace uni10{

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __elemNum(0), __ongpu(false), __elem(NULL){};

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_bool _ongpu): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(_ongpu),__elem(NULL){

      init(_Rnum, _Cnum, _isdiag, NULL);

    }

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_bool _ongpu): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(_ongpu), __elem(NULL){

      init(_Rnum, _Cnum, _isdiag, src);

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(const uni10_elem_lapack_cpu& _elem): __uni10_typeid(_elem.__uni10_typeid), __elemNum(_elem.__elemNum), __ongpu(_elem.__ongpu){

      init(1, __elemNum, false, _elem.__elem);

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::~uni10_elem_lapack_cpu(){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free_cpu(__elem, __elemNum * sizeof(uni10_type));

      __elem    = NULL;
      __elemNum = 0;

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::setElem(const uni10_type* src, bool src_ongpu){

      uni10_error_msg( src_ongpu, "%s", " The source pointer is on the device. Please install MAGMA or CUDAONLY gpu version instead.");
      uni10_error_msg( __elem == NULL, "%s", "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      uni10_error_msg( src  == NULL, "%s", "The source ptr is NULL.");

      uni10_elem_copy_cpu( __elem, src, __elemNum * sizeof(uni10_type) );

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src){

      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      if(__elem != NULL)
        uni10_elem_free_cpu(__elem, __elemNum*sizeof(uni10_type));

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc_cpu( memsize );
        if(src != NULL){
          uni10_elem_copy_cpu( __elem, src, memsize );
        }
        else{
          uni10_elemBzero_cpu( __elem, memsize );
        }

      }

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum){

      __elemNum = _Rnum * _Cnum;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc_cpu( memsize );
        uni10_elemBzero_cpu( __elem, memsize );

      }

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::clear(){
      
      __elemNum = 0;

      if(__elem != NULL)
        uni10_elem_free_cpu(__elem, __elemNum * sizeof(uni10_type));

      __elem = NULL;

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
              if(__uni10_typeid == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(__elem[ i ]);
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
            if ( __elem[ i * _Cnum + j] == 0.) {
              if(__uni10_typeid == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(__elem[ i * _Cnum + j ]);
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

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::resize(uni10_uint64 _row, uni10_uint64 _col, uni10_uint64& Rnum, uni10_uint64& Cnum, uni10_bool& isdiag, uni10_const_bool& _fixHead){

      if(_fixHead){

        if(isdiag){

          uni10_uint64 _elemNum = _row < _col ? _row : _col;

          if(_elemNum > __elemNum){

            uni10_type* _elem = (uni10_type*)uni10_elem_alloc_cpu(_elemNum * sizeof(uni10_type));

            uni10_elemBzero_cpu(_elem, _elemNum * sizeof(uni10_type));

            uni10_elem_copy_cpu(_elem, __elem, __elemNum * sizeof(uni10_type));

            if(__elem != NULL)
              uni10_elem_free_cpu( __elem, __elemNum * sizeof(uni10_type) );

            __elem = _elem;

          }
          else
            shrinkWithoutFree( (__elemNum - _elemNum) * sizeof(uni10_type) );

          __elemNum = _elemNum;

          Rnum = _row;
          Cnum = _col;

        }
        else{

          uni10_uint64 _elemNum = _row * _col;

          if(_col == Cnum){

            if(_row > Rnum){

              uni10_type* _elem = (uni10_type*)uni10_elem_alloc_cpu( _elemNum * sizeof(uni10_type) );

              uni10_elemBzero_cpu( _elem, _elemNum * sizeof(uni10_type) );

              uni10_elem_copy_cpu(_elem, __elem, __elemNum * sizeof(uni10_type) );

              if(__elem != NULL)
                uni10_elem_free_cpu( __elem, __elemNum * sizeof(uni10_type) );

              __elem = _elem;

            }
            else
              shrinkWithoutFree( (__elemNum - _elemNum) * sizeof(uni10_type) );

            __elemNum = _elemNum;

            Rnum = _row;
          }
          else{

            uni10_uint64 data_row = _row < Rnum ? _row : Rnum;
            uni10_uint64 data_col = _col < Cnum ? _col : Cnum;

            uni10_type* _elem = (uni10_type*)uni10_elem_alloc_cpu( _elemNum * sizeof(uni10_type) );

            uni10_elemBzero_cpu( _elem, _elemNum * sizeof(uni10_type) );

            for(int r = 0; r < (int)data_row; r++)
              uni10_elem_copy_cpu( &(_elem[r * _col]), &(__elem[r * Cnum]), data_col * sizeof(uni10_type) );

            if(__elem != NULL)
              uni10_elem_free_cpu(__elem, __elemNum * sizeof(uni10_type) );

            __elem    = _elem;
            __elemNum = _elemNum;

            Rnum = _row;
            Cnum = _col;

          }

        }

      }else{

        uni10_error_msg(true, "%s", "Resize fixTail is developping !!!");

      }

    }

  template class uni10_elem_lapack_cpu<uni10_double64>;
  template class uni10_elem_lapack_cpu<uni10_complex128>;

} /* namespace uni10 */
