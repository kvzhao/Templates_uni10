#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_cpu.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_cpu_z.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_cpu_dz.h"

#ifdef MKL
#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(std::complex<double> a, std::complex<double>* X, int incx, std::complex<double>* Y, int incy, size_t N){   // Y = aX + Y
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &incx, Y + offset, &incy);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    // Blas
    void vectorAdd(std::complex<double>* Y, std::complex<double>* X, size_t N){

      std::complex<double> a = 1.0;
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorSub(std::complex<double>* Y, std::complex<double>* X, size_t N){

      std::complex<double> a = -1.0;
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorMul(std::complex<double>* Y, std::complex<double>* X, size_t N){ 
      for(size_t i = 0; i < N; i++)
        Y[i] *= X[i];
    }

    void vectorScal(std::complex<double> a, std::complex<double>* X, size_t N){
      int inc = 1;
      int64_t left = N;
      size_t offset = 0;
      int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorExp(std::complex<double> a, std::complex<double>* X, size_t N){
      for(size_t i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);
    }

    std::complex<double> vectorSum(std::complex<double>* X, size_t N, int inc){
      std::complex<double> sum = 0.0;
      size_t idx = 0;
      for(size_t i = 0; i < N; i++){
        sum += X[idx];
        idx += inc;
      }
      return sum;
    }

    double vectorNorm(std::complex<double>* X, size_t N, int inc){
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
        tmp = dznrm2(&chunk, X + offset, &inc);
        norm2 += tmp * tmp;
        offset += chunk;
        left -= INT_MAX;
      }
      return sqrt(norm2);
    }

    void matrixDot(std::complex<double>* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C){
      std::complex<double> alpha = 1.0, beta = 0.0;
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
    }

    void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++)
        vectorScal(diag[i], &(mat[i * N]), N);
    }

    void diagColMul(std::complex<double> *mat, std::complex<double>* diag, size_t M, size_t N){
      for(size_t i = 0; i < M; i++){
        size_t ridx = i * N;
        for(size_t j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }
    }

    void setTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT){
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          AT[j * M + i] = A[i * N + j];
    }

    void setTranspose(std::complex<double>* A, size_t M, size_t N){
      size_t memsize = M * N * sizeof(std::complex<double>);
      std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
      setTranspose(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setDagger(std::complex<double>* A, size_t M, size_t N, std::complex<double> *AT){
      for(size_t i = 0; i < M; i++)
        for(size_t j = 0; j < N; j++)
          AT[j * M + i] = std::conj(A[i * N + j]);
    }

    void setDagger(std::complex<double>* A, size_t M, size_t N){
      size_t memsize = M * N * sizeof(std::complex<double>);
      std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
      setDagger(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setConjugate(std::complex<double> *A, size_t N, std::complex<double> *A_conj){
      for(size_t i = 0; i < N; i++)
        A_conj[i] = std::conj(A[i]);
    }

    void setConjugate(std::complex<double> *A, size_t N){
      for(size_t i = 0; i < N; i++)
        A[i] = std::conj(A[i]);
    }

    //LAPACK
    //
    void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      int min = std::min(M, N);
      double* S = (double*)malloc(min * sizeof(double));
      matrixSVD(Mij_ori, M, N, U, S, vT);
      uni10_elem_cast_cpu(S_ori, S, min);
      free(S);

    }

    void matrixQR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R){

      std::complex<double>* Mij = (std::complex<double>*)malloc(N*M*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = N;
      zgelqf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zunglq(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgelqf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri_cpu(Mij, R, M, N);
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zunglq(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);
    }

    void matrixRQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* R, std::complex<double>* Q){

      std::complex<double>* Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = M;
      zgeqlf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zungql(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqlf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri_cpu(Mij, R, M, N);
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungql(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixLQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* L, std::complex<double>* Q){

      std::complex<double>* Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = M;
      zgeqrf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zungqr(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqrf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri_cpu(Mij, L, M, N);
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungqr(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixQL(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L){

      std::complex<double>* Mij = (std::complex<double>*)malloc(N*M*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = N;
      zgerqf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zungrq(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgerqf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri_cpu(Mij, L, M, N);
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungrq(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //getR
      //std::complex<double> alpha (1.0, 0.0), beta (1.0, 1.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, L, &N);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);
    }

    void matrixQDR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      std::complex<double>* Mij = (std::complex<double>*)malloc(N*M*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = N;
      zgelqf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zunglq(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgelqf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri_cpu(Mij, R, M, N);
      uni10_getDiag_cpu(R, D, M, N, N);
      for(int i = 0; i < N; i++)
        for(int j = 0; j < N-i; j++)
          R[i*N+i+j] /= D[i];
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zunglq(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);
    }


    void matrixLDQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* L, std::complex<double>* D, std::complex<double>* Q){

      std::complex<double>* Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      int lda = N;
      int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      int info;
      int K = M;
      zgeqrf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
      zungqr(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqrf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri_cpu(Mij, L, M, N);
      uni10_getDiag_cpu(L, D, M, N, M);
      for(int i = 0; i < M; i++)
        for(int j = 0; j < M-i; j++)
          L[(i+j)*M+i] /= D[i];
      //getQ
      lwork = (int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungqr(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
      memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(Mij);
      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixQDRCPIVOT(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(M != N, "%s", "M must be equalt to N");

      std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>)); // Q(M x M): orthogonal basis
      setTranspose(Mij_ori, M, N, Mij);                      //column major + square matrix  // TP(M x M)
      int max = M > N ? M : N;
      int min = M < N ? M : N;
      int lda = max;
      int lwork = N;
      int lrwork = 2 * N;
      std::complex<double>* work = (std::complex<double>*) malloc(lwork * sizeof(std::complex<double>));
      double* rwork = (double*) malloc(lrwork * sizeof(double));
      int* jpvt = (int*)malloc(N * sizeof(int));  //column vectors
      memset(jpvt, 0, N * sizeof(int));
      std::complex<double>* tau = (std::complex<double>*)malloc(min * sizeof(std::complex<double>));
      int info;
      zgeqpf(&M, &N, Mij, &lda, jpvt, tau, work, rwork, &info);
      uni10_error_msg(info != 0, "Lapack Info = %d", info);
      for(int i = 0; i < M; i++)
        D[i] = Mij[i * N + i];                               // D
      std::complex<double>* T = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>));
      for(int i = 0; i < M; i++)
        for(int j = i; j < N; j++)
          if(i == j)
            T[i * N + j] = 1;
          else
            T[i * N + j] = Mij[j * N + i] / D[i];
      for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
          R[i * N + (jpvt[j]-1)] = T[i * N + j];              // R 
      zungqr(&M, &N, &N, Mij, &lda, tau, work, &lwork, &info);
      uni10_error_msg(info != 0, "Lapack Info = %d", info);
      setTranspose(Mij, M, N, Q);                             // Q
      free(Mij);
      free(work);
      free(jpvt);
      free(tau);
      free(T);

    }


    void matrixInv(std::complex<double>* A, int N){
      //if(diag){
      //  for(int i = 0; i < N; i++)
      //    A[i] = std::abs(A[i]) == 0 ? 0.0 : 1.0/A[i];
      //  return;
      //}
      int *ipiv = (int*)malloc((N+1) * sizeof(int));
      int info;
      zgetrf(&N, &N, A, &N, ipiv, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetrf': Lapack INFO = ", info);

      int lwork = -1;
      std::complex<double> worktest;
      zgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetri': Lapack INFO = ", info);

      lwork = (int)(worktest.real());
      std::complex<double> *work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
      zgetri(&N, A, &N, ipiv, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetri': Lapack INFO = ", info);

      free(ipiv);
      free(work);
    }

    std::complex<double> matrixDet(std::complex<double>* A, int N){

      int *ipiv = (int*)malloc((N+1)*sizeof(int));
      int lwork = 64 * N;
      std::complex<double> *work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
      int info;
      zgetrf(&N,&N,A,&N,ipiv,&info);
      uni10_error_msg( info != 0, "%s %d", "Error in Lapack function 'zgetrf': Lapack INFO = ", info );
      std::complex<double> det = 1;
      int neg = 0;
      for (int i = 0; i < N; i++) {
        det *= A[i * N + i];
        if (ipiv[i] != (i+1)) neg = !neg;
      }
      free(ipiv);
      free(work);
      return neg?-det:det;

    }

    void eigDecompose(std::complex<double>* Kij, int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      size_t memsize = N * N * sizeof(std::complex<double>);
      std::complex<double> *A = (std::complex<double>*) malloc(memsize);
      memcpy(A, Kij, memsize);
      int ldA = N;
      int ldvl = 1;
      int ldvr = N;
      int lwork = -1;
      double *rwork = (double*) malloc(2 * N * sizeof(double));
      std::complex<double> worktest;
      int info;
      zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, &worktest, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s, %d", "Error in Lapack function 'zgeev': Lapack INFO = ", info);

      lwork = (int)worktest.real();
      std::complex<double>* work = (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
      zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, work, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s, %d", "Error in Lapack function 'zgeev': Lapack INFO = ", info);

      free(work);
      free(rwork);
      free(A);
    }

    void setIdentity(std::complex<double>* elem, size_t M, size_t N){
      size_t min;
      if(M < N)
        min = M;
      else
        min = N;
      memset(elem, 0, M * N * sizeof(std::complex<double>));
      for(size_t i = 0; i < min; i++)
        elem[i * N + i] = 1.0;
    }




  } /* namespace uni10_linalg */

}
