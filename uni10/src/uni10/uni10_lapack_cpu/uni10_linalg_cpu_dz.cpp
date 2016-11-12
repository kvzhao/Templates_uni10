#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/uni10_tools_cpu.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_z.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_dz.h"

#ifdef MKL
#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(std::complex<double>* Y, double* X, size_t N){
      for(size_t i = 0; i < N; i++)
        Y[i] += X[i];
    }

    void vectorSub(std::complex<double>* Y, double* X, size_t N){
      for(size_t i = 0; i < N; i++)
        Y[i] -= X[i];
    }

    void matrixMul(double* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C){

      int size_A = M * K;
      std::complex<double>* CA = (std::complex<double>*)malloc(size_A*sizeof(std::complex<double>));
      uni10_elem_cast_cpu(CA, A, size_A);
      std::complex<double> alpha = 1.0, beta = 0.0;
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, CA, &K, &beta, C, &N);
    }

    void matrixMul(std::complex<double>* A, double* B, int M, int N, int K, std::complex<double>* C){

      int size_B = K * N;
      std::complex<double>* CB = (std::complex<double>*)malloc(size_B*sizeof(std::complex<double>));
      std::complex<double> alpha = 1.0, beta = 0.0;
      uni10_elem_cast_cpu(CB, B, size_B);
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, CB, &N, A, &K, &beta, C, &N);
    }

    void vectorScal(double a, std::complex<double>* X, size_t N){
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zdscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void diagRowMul(std::complex<double>* mat, double* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++)
        vectorScal(diag[i], &(mat[i * N]), N);
    }

    void diagColMul(std::complex<double>* mat, double* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++){
        size_t ridx = i * N;
        for(size_t j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }
    }

    void vectorExp(double a, std::complex<double>* X, size_t N){
      for(size_t i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);
    }

    void eigDecompose(double* Kij_ori, int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      std::complex<double> *Kij = (std::complex<double>*) malloc(N * N * sizeof(std::complex<double>));
      uni10_elem_cast_cpu(Kij, Kij_ori, N * N);
      eigDecompose(Kij, N, Eig, EigVec);
      free(Kij);
    }

    void eigSyDecompose(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec){
      //eigDecompose(Kij, N, Eig, EigVec, ongpu);
      memcpy(EigVec, Kij, N * N * sizeof(std::complex<double>));
      int ldA = N;
      int lwork = -1;
      std::complex<double> worktest;
      double* rwork = (double*) malloc((3*N+1) * sizeof(double));
      int info;
      zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, rwork, &info);

      uni10_lapack_error_msg(info != 0, "Error in Lapack function 'zheev': Lapack INFO = ", info);

      lwork = (int)worktest.real();
      std::complex<double>* work= (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
      zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, rwork, &info);

      uni10_lapack_error_msg(info != 0, "Error in Lapack function 'zheev': Lapack INFO = ", info);

      free(work);
      free(rwork);
    }

    void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, M * N * sizeof(std::complex<double>));
      int min = std::min(M, N);
      int ldA = N, ldu = N, ldvT = min;
      int lwork = -1;
      std::complex<double> worktest;
      int info;
      double *rwork = (double*) malloc(std::max(1, 5 * min) * sizeof(double));
      zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, rwork, &info);

      uni10_lapack_error_msg(info != 0, "Error in Lapack function 'zgesvd': Lapack INFO = ", info);

      lwork = (int)(worktest.real());
      std::complex<double> *work = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, rwork, &info);

      uni10_lapack_error_msg(info != 0, "Error in Lapack function 'zgesvd': Lapack INFO = ", info);

      free(rwork);
      free(work);
      free(Mij);
    }

  };/* namespace uni10_linalg */

};/* namespace uni10 */

