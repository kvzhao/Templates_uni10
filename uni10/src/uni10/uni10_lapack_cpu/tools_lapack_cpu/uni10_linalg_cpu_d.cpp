#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_cpu.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_cpu_d.h"

#ifdef MKL
#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(double a, double* X, int incx, double* Y, int incy, size_t N){   // Y = aX + Y
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &incx, Y + offset, &incy);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorAdd(double* Y, double* X, size_t N){   // Y = Y + X
      double a = 1.0;
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorSub(double* Y, double* X, size_t N){   // Y = Y + X
      double a = -1.0;
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorMul(double* Y, double* X, size_t N){ // Y = Y * X, element-wise multiplication;
      for(size_t i = 0; i < N; i++)
        Y[i] *= X[i];
    }

    void vectorScal(double a, double* X, size_t N){
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        dscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorExp(double a, double* X, size_t N){
      for(size_t i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);
    }

    double vectorSum(double* X, size_t N, int inc){
      double sum = 0;
      size_t idx = 0;
      for(size_t i = 0; i < N; i++){
        sum += X[idx];
        idx += inc;
      }
      return sum;
    }

    double vectorNorm(double* X, size_t N, int inc){
      double norm2 = 0;
      double tmp = 0;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        tmp = dnrm2(&chunk, X + offset, &inc);
        norm2 += tmp * tmp;
        offset += chunk;
        left -= INT_MAX;
      }
      return sqrt(norm2);
    }

    void matrixDot(double* A, double* B, int M, int N, int K, double* C){
      double alpha = 1, beta = 0;
      dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
    }

    void diagRowMul(double* mat, double* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++)
        vectorScal(diag[i], &(mat[i * N]), N);
    }

    void diagColMul(double *mat, double* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++){
        size_t ridx = i * N;
        for(size_t j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }
    }

    void setTranspose(double* A, size_t M, size_t N, double* AT){
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          AT[j * M + i] = A[i * N + j];
    }

    void setTranspose(double* A, size_t M, size_t N){
      size_t memsize = M * N * sizeof(double);
      double *AT = (double*)malloc(memsize);
      setTranspose(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setDagger(double* A, size_t M, size_t N, double* AT){
      setTranspose(A, M, N, AT);
    }

    void setDagger(double* A, size_t M, size_t N){
      setTranspose(A, M, N);
    }

    void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT){

      double* Mij = (double*)malloc(M * N * sizeof(double));
      memcpy(Mij, Mij_ori, M * N * sizeof(double));
      int min = std::min(M, N);
      int ldA = N, ldu = N, ldvT = min;
      int lwork = -1;
      double worktest;
      int info;

      dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesvd': Lapack INFO = ", info);

      lwork = (int)worktest;
      double *work = (double*)malloc(lwork*sizeof(double));
      dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesvd': Lapack INFO = ", info);

      free(work);
      free(Mij);
    }

    // lapack is builded by fortran which is load by column, so we use 
    // dorgqr -> lq
    // dorglq -> qr
    // dorgrq -> ql 
    // dorgql -> rq
    void matrixQR(double* Mij_ori, int M, int N, double* Q, double* R){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQR()");

      double* Mij = (double*)malloc(N*M*sizeof(double));
      memcpy(Mij, Mij_ori, N*M*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      int lda = N;
      int lwork = -1;
      double worktestdge;
      double worktestdor;
      int info;
      int K = N;
      dgelqf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
      dorglq(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgelqf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
      //getQ
      uni10_getUpTri_cpu(Mij, R, M, N);
      lwork = (int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorglq(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(double));
      //getR
      //double alpha = 1, beta = 0;
      //dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(Mij);
      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixRQ(double* Mij_ori, int M, int N, double* Q, double* R){

      uni10_error_msg(N < M, "%s", "N must be larger than M in matrixRQ()");

      double* Mij = (double*)malloc(M*N*sizeof(double));
      memcpy(Mij, Mij_ori, M*N*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      int lda = N;
      int lwork = -1;
      double worktestdge;
      double worktestdor;
      int info;
      int K = M;
      dgeqlf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
      dorgql(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgeqlf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
      //getR
      //printf("\n\n");
      //for(int i = 0; i < M; i++){
      //  for(int j = 0; j < N; j++){
      //    printf("%.4f ", Mij[i*N+j]);
      //  }
      //  printf("\n\n");
      //}
      uni10_getUpTri_cpu(Mij, R, M, N);
      ///getQ
      lwork = (int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgql(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(double));
      //double alpha = 1, beta = 0;
      //dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

      free(Mij);
      free(tau);
      free(workdge);
      free(workdor);

    }

    void matrixLQ(double* Mij_ori, int M, int N, double* Q, double* L){

      uni10_error_msg(N < M, "%s","N must be larger than M in matrixLQ()");

      double* Mij = (double*)malloc(M*N*sizeof(double));
      memcpy(Mij, Mij_ori, M*N*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      int lda = N;
      int lwork = -1;
      double worktestdge;
      double worktestdor;
      int info;
      int K = M;
      dgeqrf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
      dorgqr(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgeqrf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
      //getL
      uni10_getDnTri_cpu(Mij, L, M, N);
      //getQ
      lwork = (int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgqr(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(double));
      //getR
      //double alpha = 1, beta = 0;
      //dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(Mij);
      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixQL(double* Mij_ori, int M, int N, double* Q, double* L){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQL()");

      double* Mij = (double*)malloc(N*M*sizeof(double));
      memcpy(Mij, Mij_ori, N*M*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      int lda = N;
      int lwork = -1;
      double worktestdge;
      double worktestdor;
      int info;
      int K = N;
      dgerqf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
      dorgrq(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgerqf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
      //getR
      uni10_getDnTri_cpu(Mij, L, M, N);
      //getQ
      lwork = (int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgrq(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(double));
      //double alpha = 1, beta = 0;
      //dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(Mij);
      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixInv(double* A, int N){

      int *ipiv = (int*)malloc((N+1)*sizeof(int));
      int info;
      dgetrf(&N, &N, A, &N, ipiv, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetrf': Lapack INFO = ", info);

      int lwork = -1;
      double worktest = 0.;
      dgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetri': Lapack INFO = ", info);

      lwork = (int)worktest;
      double *work = (double*)malloc(lwork * sizeof(double));
      dgetri(&N, A, &N, ipiv, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetri': Lapack INFO = ", info);

      free(ipiv);
      free(work);
    }

    double matrixDet(double* A, int N){

      int *ipiv = (int*)malloc((N+1)*sizeof(int));
      int lwork = 64 * N;
      double *work = (double*)malloc(lwork * sizeof(double));
      int info;
      dgetrf(&N,&N,A,&N,ipiv,&info);
      uni10_error_msg( info != 0, "%s %d", "Error in Lapack function 'dgetrf': Lapack INFO = ", info );
      double det = 1;
      int neg = 0;
      for (int i = 0; i < N; i++) {
        det *= A[i * N + i];
        if (ipiv[i] != (i+1)) neg = !neg;
      }
      free(ipiv);
      free(work);
      return neg?-det:det;

    }

    void setIdentity(double* elem, size_t M, size_t N){
      size_t min;
      if(M < N)
        min = M;
      else  
        min = N;
      memset(elem, 0, M * N * sizeof(double));
      for(size_t i = 0; i < min; i++)
        elem[i * N + i] = 1;
    }

    void eigSyDecompose(double* Kij, int N, double* Eig, double* EigVec){

      memcpy(EigVec, Kij, N * N * sizeof(double));
      int ldA = N;
      int lwork = -1;
      double worktest;
      int info;
      dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dsyev': Lapack INFO = ", info);

      lwork = (int)worktest;
      double* work= (double*)malloc(sizeof(double)*lwork);
      dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dsyev': Lapack INFO = ", info);

      free(work);
    }


  } /* namespace uni10_linalg */

} /* namespace uni10 */
