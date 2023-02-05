// #include <immintrin.h>
// const char* dgemm_desc = "Simple blocked dgemm.";

// #ifndef BLOCK
// #define BLOCK 32
// #endif

// #define min(a, b) (((a) < (b)) ? (a) : (b))

// // 4-way multiply-add
// #define microkernel(As, Cs, Br) { \
//     Ar = _mm256_load_pd(As); \
//     Cr = _mm256_load_pd(Cs); \
//     _mm256_store_pd(Cs, _mm256_fmadd_pd(Ar, Br, Cr)); \
// }

// // block multiplication for 32-padding: 4 * 8 micro kernel
// static void do_block_32(int lda, double* A, double* B, double* C) {
//     double *As, *Cs = C, *Bs = B;
//     for (int j = 0; j < BLOCK; ++j, Cs += lda, Bs += lda - BLOCK){
//         As = A;
//         for (int k = 0; k < BLOCK; ++k, As += lda, Bs++){
//             __m256d Ar, Cr, Br = _mm256_set1_pd(*Bs);
//             microkernel(As, Cs, Br);
//             microkernel(As + 4, Cs + 4, Br);
//             microkernel(As + 8, Cs + 8, Br);
//             microkernel(As + 12, Cs + 12, Br);
//             microkernel(As + 16, Cs + 16, Br);
//             microkernel(As + 20, Cs + 20, Br);
//             microkernel(As + 24, Cs + 24, Br);
//             microkernel(As + 28, Cs + 28, Br);
//         }
//     }
// }

// // block multiplication for 4-padding
// static void do_block_4(int lda, int M, int N, int K, double* A, double* B, double* C) {
//     if(M == 32){
//         for (int j = 0; j < N; ++j){
//             for (int k = 0; k < K; ++k){
//                 __m256d Ar, Cr, Br = _mm256_set1_pd(B[k + j * lda]);

//                 double *As = A + k * lda, *Cs = C + j * lda;
//                 microkernel(As, Cs, Br);
//                 microkernel(As + 4, Cs + 4, Br);
//                 microkernel(As + 8, Cs + 8, Br);
//                 microkernel(As + 12, Cs + 12, Br);
//                 microkernel(As + 16, Cs + 16, Br);
//                 microkernel(As + 20, Cs + 20, Br);
//                 microkernel(As + 24, Cs + 24, Br);
//                 microkernel(As + 28, Cs + 28, Br);
//             }
//         }
//     }
//     else{
//         for (int j = 0; j < N; ++j){
//             for (int k = 0; k < K; ++k){
//                 __m256d Ar, Cr, Br = _mm256_set1_pd(B[k + j * lda]);

//                 double *As = A + k * lda, *Cs = C + j * lda;
//                 for(int i = 0; i < M; i += 4, As += 4, Cs += 4)
//                     microkernel(As, Cs, Br);
//             }
//         }
//     }
    
// }

// // copy original matrix to padded and aligned matrix
// void copy_aligned(double *out, double *in, int out_size, int in_size){
//     for(int j = 0; j < out_size; j++)
//         for(int i = 0; i < out_size; i++){
//             if(i < in_size && j < in_size) out[j * out_size + i] = in[j * in_size + i];
//             else out[j * out_size + i] = 0;
//         }
// }

// /* This routine performs a dgemm operation
//  *  C := C + A * B
//  * where A, B, and C are lda-by-lda matrices stored in column-major format.
//  * On exit, A and B maintain their input values. */
// void square_dgemm(int lda, double* A, double* B, double* C) {
//     int LDA = lda;

//     // adjust pad size
//     if(lda == 255 || lda == 256){
//         lda = ((lda - 1) / 4 + 2) * 4;
//     }
//     else if(lda == 229 || lda == 97 || lda == 129 || lda == 257){
//         lda = ((lda - 1) / 4 + 1) * 4;
//     }
//     else{
//         lda = ((lda - 1) / 32 + 1) * 32;
//         if(lda % 256 == 0) lda += 32;
//     }

//     double *A_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
//     double *B_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);
//     double *C_aligned = (double*) _mm_malloc(lda * lda * sizeof(double), 32);

//     copy_aligned(A_aligned, A, lda, LDA);
//     copy_aligned(B_aligned, B, lda, LDA);
//     copy_aligned(C_aligned, C, lda, LDA);

//     if(lda % 32 == 0){
//         double *As, *Bs, *Cs;
//         for (int j = 0; j < lda; j += BLOCK) {
//             for (int k = 0; k < lda; k += BLOCK) {
//                 for (int i = 0; i < lda; i += BLOCK) {
//                     // C[i,j] += A[i,k] * B[k,j]
//                     As = A_aligned + (i + k * lda);
//                     Bs = B_aligned + (k + j * lda);
//                     Cs = C_aligned + (i + j * lda);
//                     do_block_32(lda, As, Bs, Cs);
//                 }
//             }
//         }
//     }
//     else{
//         double *As, *Bs, *Cs;
//         for (int j = 0; j < lda; j += BLOCK) 
//             for (int k = 0; k < lda; k += BLOCK) 
//                 for (int i = 0; i < lda; i += BLOCK) {
//                     int M = min(BLOCK, lda - i);
//                     int N = min(BLOCK, lda - j);
//                     int K = min(BLOCK, lda - k);
//                     As = A_aligned + (i + k * lda);
//                     Bs = B_aligned + (k + j * lda);
//                     Cs = C_aligned + (i + j * lda);
//                     do_block_4(lda, M, N, K, As, Bs, Cs);
//                 }
//     }
    
//     // copy back to C
//     for(int j = 0; j < LDA; j++)
//         for(int i = 0; i < LDA; i++)
//             C[i + j * LDA] = C_aligned[i + j * lda];
    
//     // free
//     _mm_free(A_aligned);
//     _mm_free(B_aligned);
//     _mm_free(C_aligned);
// }




// #include <immintrin.h>
// #include <stdio.h>
// const char* dgemm_desc = "Simple blocked dgemm.";

// #ifndef BLOCK
// #define BLOCK 32
// #endif

// #define min(a, b) (((a) < (b)) ? (a) : (b))

// #define microkernel(As, Cs, Br) { \
//     Ar = _mm256_load_pd(As); \
//     Cr = _mm256_load_pd(Cs); \
//     _mm256_store_pd(Cs, _mm256_fmadd_pd(Ar, Br, Cr)); \
// }

// /*
//  * This auxiliary subroutine performs a smaller dgemm operation
//  *  C := C + A * B
//  * where C is M-by-N, A is M-by-K, and B is K-by-N.
//  */
// static void do_block_32(int lda, double *A, double *B, double *C) {
//     double *As, *Bs = B, *Cs = C;
//     // For each col j of B
//     for (int j = 0; j < BLOCK; ++j, Cs += lda, Bs += lda - BLOCK){
//         As = A;
//         // For each column k of A, row k of B
//         for (int k = 0; k < BLOCK; ++k, As += lda, Bs++){
//             __m256d Ar, Cr, Br = _mm256_set1_pd(*Bs);
//             microkernel(As, Cs, Br);
//             microkernel(As+4, Cs+4, Br);
//             microkernel(As+8, Cs+8, Br);
//             microkernel(As+12, Cs+12, Br);
//             microkernel(As+16, Cs+16, Br);
//             microkernel(As+20, Cs+20, Br);
//             microkernel(As+24, Cs+24, Br);
//             microkernel(As+28, Cs+28, Br);
//         }
//     }
// }

// void copy_aligned(double *out, double *in, int outsize, int insize) {
//     for (int j = 0; j < outsize; ++j) {
//         for (int i = 0; i < outsize; ++i) {
//             if (i < insize && j < insize) {
//                 out[i+j*outsize] = in[i+j*insize];
//             } else {
//                 out[i+j*outsize] = 0;
//             }
//         }
//     }
// }

// /* This routine performs a dgemm operation
//  *  C := C + A * B
//  * where A, B, and C are lda-by-lda matrices stored in column-major format.
//  * On exit, A and B maintain their input values. */
// void square_dgemm(int lda, double* A, double* B, double* C) {
//     int LDA = lda;

//     // padding
//     lda = ((lda - 1) / 32 + 1) * 32;

//     double *A_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);
//     double *B_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);
//     double *C_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);

//     copy_aligned(A_aligned, A, lda, LDA);
//     copy_aligned(B_aligned, B, lda, LDA);
//     copy_aligned(C_aligned, C, lda, LDA);

//     double *As, *Bs, *Cs;
//     for (int j = 0; j < lda; j += BLOCK) {
//         for (int k = 0; k < lda; k += BLOCK) {
//             for (int i = 0; i < lda; i += BLOCK) {
//                 As = A_aligned + (i + k * lda);
//                 Bs = B_aligned + (k + j * lda);
//                 Cs = C_aligned + (i + j * lda);
//                 do_block_32(lda, As, Bs, Cs);
//             }
//         }
//     }

//     for(int j = 0; j < LDA; j++)
//         for(int i = 0; i < LDA; i++)
//             C[i + j * LDA] = C_aligned[i + j * lda];
    
//     _mm_free(A_aligned);
//     _mm_free(B_aligned);
//     _mm_free(C_aligned);
// }


#include <immintrin.h>
#include <stdio.h>
const char* dgemm_desc = "Simple blocked dgemm.";

#ifndef BLOCK
#define BLOCK 32
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))

#define microkernel(As, Cs, Br) { \
    Ar = _mm256_load_pd(As); \
    Cr = _mm256_load_pd(Cs); \
    _mm256_store_pd(Cs, _mm256_fmadd_pd(Ar, Br, Cr)); \
}

/*
 * This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 */
static void do_block_32(int block_size, double *A, double *B, double *C) {
    double *As, *Bs = B, *Cs = C;
    // For each col j of B
    for (int j = 0; j < BLOCK; ++j, Cs += BLOCK){
        As = A;
        // For each column k of A, row k of B
        for (int k = 0; k < BLOCK; ++k, As += BLOCK, Bs++){
            __m256d Ar, Cr, Br = _mm256_set1_pd(*Bs);
            microkernel(As, Cs, Br);
            microkernel(As+4, Cs+4, Br);
            microkernel(As+8, Cs+8, Br);
            microkernel(As+12, Cs+12, Br);
            microkernel(As+16, Cs+16, Br);
            microkernel(As+20, Cs+20, Br);
            microkernel(As+24, Cs+24, Br);
            microkernel(As+28, Cs+28, Br);
        }
    }
}

// void copy_aligned(double *out, double *in, int outsize, int insize) {
//     for (int j = 0; j < outsize; ++j) {
//         for (int i = 0; i < outsize; ++i) {
//             if (i < insize && j < insize) {
//                 out[i+j*outsize] = in[i+j*insize];
//             } else {
//                 out[i+j*outsize] = 0;
//             }
//         }
//     }
// }

void copy_aligned(double *out, double *in, int outsize, int insize) {
    int cnt=0;
    for (int i = 0; i < insize; i += BLOCK){
        for(int j = 0; j < insize; j += BLOCK){
            for(int ii = 0; ii < BLOCK; ii++){
                for(int jj = 0; jj < BLOCK; jj++){
                    if(i+ii < outsize && j+jj < outsize) out[cnt] = in[(i+ii)*insize + (j+jj)];
                    else out[cnt] = 0;
                    cnt++;
                }
            }
        }
    }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int lda, double* A, double* B, double* C) {
    int LDA = lda;

    // padding
    lda = ((lda - 1) / 32 + 1) * 32;
    int block_num = lda/BLOCK;
    int block_size = BLOCK * BLOCK;
    double *A_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);
    double *B_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);
    double *C_aligned = (double*) _mm_malloc(lda*lda*sizeof(double), 32);

    copy_aligned(A_aligned, A, lda, LDA);
    copy_aligned(B_aligned, B, lda, LDA);
    copy_aligned(C_aligned, C, lda, LDA);
    

    double *As, *Bs, *Cs;
    for (int j = 0; j < block_num; j += 1) {
        for (int k = 0; k < block_num; k += 1) {
            for (int i = 0; i < block_num; i += 1) {
                As = A_aligned + (i * block_size + k * block_size * block_num);
                Bs = B_aligned + (k * block_size + j * block_size * block_num);
                Cs = C_aligned + (i * block_size + j * block_size * block_num);
                do_block_32(lda, As, Bs, Cs);
            }
        }
    }

    int cnt=0;
    for (int i = 0; i < lda; i += BLOCK){
        for(int j = 0; j < lda; j += BLOCK){
            for(int ii = 0; ii < BLOCK; ii++){
                for(int jj = 0; jj < BLOCK; jj++){
                    if(i+ii < LDA && j+jj < LDA) C[(i+ii)*lda + (j+jj)] = C_aligned[cnt];
                    cnt++;
                }
            }
        }
    }

    _mm_free(A_aligned);
    _mm_free(B_aligned);
    _mm_free(C_aligned);
}