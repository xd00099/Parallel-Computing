#include <immintrin.h>
const char* dgemm_desc = "Simple blocked dgemm.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))

// 4-way multiply-add
#define microkernel(As, Cs, Br) { \
    Ar = _mm256_load_pd(As); \
    Cr = _mm256_load_pd(Cs); \
    _mm256_store_pd(Cs, _mm256_fmadd_pd(Ar, Br, Cr)); \
}

// For 4-aligned matrix multiplication:
// block multiplication for 4-padding
static void do_block_4(int lda, int M, int N, int K, double* A, double* B, double* C) {
    if(M == 32){
        for (int j = 0; j < N; ++j){
            for (int k = 0; k < K; ++k){
                __m256d Ar, Cr, Br = _mm256_set1_pd(B[k + j * lda]);

                double *As = A + k * lda, *Cs = C + j * lda;
                microkernel(As, Cs, Br);
                microkernel(As + 4, Cs + 4, Br);
                microkernel(As + 8, Cs + 8, Br);
                microkernel(As + 12, Cs + 12, Br);
                microkernel(As + 16, Cs + 16, Br);
                microkernel(As + 20, Cs + 20, Br);
                microkernel(As + 24, Cs + 24, Br);
                microkernel(As + 28, Cs + 28, Br);
            }
        }
    }
    else{
        for (int j = 0; j < N; ++j){
            for (int k = 0; k < K; ++k){
                __m256d Ar, Cr, Br = _mm256_set1_pd(B[k + j * lda]);

                double *As = A + k * lda, *Cs = C + j * lda;
                for(int i = 0; i < M; i += 4, As += 4, Cs += 4)
                    microkernel(As, Cs, Br);
            }
        }
    }
    
}

// copy original matrix to padded and aligned matrix
void copy_aligned_4(double *out, double *in, int out_size, int in_size){
    for(int j = 0; j < out_size; j++)
        for(int i = 0; i < out_size; i++){
            if(i < in_size && j < in_size) out[j * out_size + i] = in[j * in_size + i];
            else out[j * out_size + i] = 0;
        }
}

// main function
void square_dgemm_4(int lda, double * A, double *B, double *C){
    int LDA = lda;
    lda = ((lda - 1) / 4 + 2) * 4;
    double *A_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
    double *B_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
    double *C_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);

    // magic malloc
    _mm_malloc(lda * lda * sizeof(double), 32);

    copy_aligned_4(A_aligned, A, lda, LDA);
    copy_aligned_4(B_aligned, B, lda, LDA);
    copy_aligned_4(C_aligned, C, lda, LDA);
    double *As, *Bs, *Cs;
    for (int j = 0; j < lda; j += BLOCK_SIZE) 
        for (int k = 0; k < lda; k += BLOCK_SIZE) 
            for (int i = 0; i < lda; i += BLOCK_SIZE) {
                int M = min(BLOCK_SIZE, lda - i);
                int N = min(BLOCK_SIZE, lda - j);
                int K = min(BLOCK_SIZE, lda - k);
                As = A_aligned + (i + k * lda);
                Bs = B_aligned + (k + j * lda);
                Cs = C_aligned + (i + j * lda);
                do_block_4(lda, M, N, K, As, Bs, Cs);
            }
    
    for(int j = 0; j < LDA; j++)
        for(int i = 0; i < LDA; i++)
            C[i + j * LDA] = C_aligned[i + j * lda];
    
    _mm_free(A_aligned);
    _mm_free(B_aligned);
    _mm_free(C_aligned);
}

// For 32-aligned matrix multiplication: 
// block multiplication for 32-padding: 4 * 8 micro kernel
static void do_block_32(int lda, double* A, double* B, double* C) {
    double *As, *Cs = C, *Bs = B;
    for (int j = 0; j < BLOCK_SIZE; ++j, Cs += BLOCK_SIZE){
        As = A;
        for (int k = 0; k < BLOCK_SIZE; ++k, As += BLOCK_SIZE, Bs++){
            __m256d Ar, Cr, Br = _mm256_set1_pd(*Bs);
            microkernel(As, Cs, Br);
            microkernel(As + 4, Cs + 4, Br);
            microkernel(As + 8, Cs + 8, Br);
            microkernel(As + 12, Cs + 12, Br);
            microkernel(As + 16, Cs + 16, Br);
            microkernel(As + 20, Cs + 20, Br);
            microkernel(As + 24, Cs + 24, Br);
            microkernel(As + 28, Cs + 28, Br);
        }
    }
}

// copy and repack
void copy_aligned_and_repack(double *out, double *in, int out_size, int in_size){
    int cnt = 0;
    for (int i = 0; i < out_size; i += BLOCK_SIZE)
        for(int j = 0; j < out_size; j += BLOCK_SIZE)
            for(int ii = 0; ii < BLOCK_SIZE; ii++)
                for(int jj = 0; jj < BLOCK_SIZE; jj++){
                    if(i + ii < in_size && j + jj < in_size) out[cnt++] = in[(i + ii) * in_size + (j + jj)];
                    else out[cnt++] = 0;
                }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int lda, double* A, double* B, double* C) {

    // check if 4-aligned needed
    if((lda &31) == 1){
        square_dgemm_4(lda, A, B, C);
        return;
    }


    int LDA = lda;
    lda = ((lda - 1) / 32 + 1) * 32;
    double *A_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
    double *B_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
    double *C_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
    
    // magic malloc
    if(lda > 32){
        _mm_malloc(lda * lda * sizeof(double), 32);
    }

    copy_aligned_and_repack(A_aligned, A, lda, LDA);
    copy_aligned_and_repack(B_aligned, B, lda, LDA);
    copy_aligned_and_repack(C_aligned, C, lda, LDA);

    double *As, *Bs, *Cs;
    for (int j = 0; j < lda; j += BLOCK_SIZE) {
        for (int k = 0; k < lda; k += BLOCK_SIZE) {
            for (int i = 0; i < lda; i += BLOCK_SIZE) {
                // C[i,j] += A[i,k] * B[k,j]
                As = A_aligned + (i * BLOCK_SIZE + k * lda);
                Bs = B_aligned + (k * BLOCK_SIZE + j * lda);
                Cs = C_aligned + (i * BLOCK_SIZE + j * lda);
                do_block_32(lda, As, Bs, Cs);
            }
        }
    }

    //repack back to original C
    int cnt=0;
    for (int i = 0; i < lda; i += BLOCK_SIZE)
        for(int j = 0; j < lda; j += BLOCK_SIZE)
            for(int ii = 0; ii < BLOCK_SIZE; ii++)
                for(int jj = 0; jj < BLOCK_SIZE; jj++){
                    if(i+ii<LDA && j+jj<LDA) C[(i + ii) * LDA + (j + jj)] = C_aligned[cnt];
                    cnt++;
                }
    
    _mm_free(A_aligned);
    _mm_free(B_aligned);
    _mm_free(C_aligned);
}
