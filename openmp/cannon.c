/*
 * =====================================================================================
 *
 *       Filename:  cannon.c
 *
 *    Description:  Cannon algorithm
 *
 *        Version:  1.0
 *       Compiler:  gcc
 *
 *         Author:  Jorge Ramirez <jorgeramirez1990@gmail.com> 
 *
 * =====================================================================================
 */

#include <stdio.h>
#include "../util/matrix.h"
#include "omp.h"

#define P_SQRT 2 
#define P (P_SQRT * P_SQRT) //number of processes
#define N 4 //matrix size
#define BLOCK_SZ (N / P_SQRT) //block size

/**
 * shift the given matrix left.
 *
 * @param m: the matrix to shift.
 * @param incstep: a value > 0 indicates that it is a first shift, otherwise is
 *                 normal shift.
 */
void shift_matrix_left(Matrix *m, int incstep) {
    int i, j, k, s, step = BLOCK_SZ;
    Matrix aux;
    
    create_matrix(&aux, 1, m->ncol);
    for(k = 0, s = 0; k < N; k += BLOCK_SZ, s++){
        for(i = k; i < (k + BLOCK_SZ); i++){
            if(incstep > 0){
                step = s * BLOCK_SZ;
            }
            for(j = 0; j < m->ncol; j++){
                aux.data[0][j] = m->data[i][(j + step) % m->ncol];
            }
            for(j = 0; j < m->ncol; j++){
                m->data[i][j] = aux.data[0][j];
            }
        }
    }
}


/**
 * shift the given matrix up.
 *
 * @param m: the matrix to shift.
 * @param incstep: a value > 0 indicates that it is a first shift, otherwise is
 *                 normal shift.
 */
void shift_matrix_up(Matrix *m, int incstep) {
    int i, j, k, s, step = BLOCK_SZ;
    Matrix aux;
    
    create_matrix(&aux, 1, m->nrow);
    for(k = 0, s = 0; k < N; k += BLOCK_SZ, s++){
        for(i = k; i < (k + BLOCK_SZ); i++){
            if(incstep > 0){
                step = s * BLOCK_SZ;
            }
            for(j = 0; j < m->nrow; j++){
                aux.data[0][j] = m->data[(j + step) % m->nrow][i];
            }
            for(j = 0; j < m->nrow; j++){
                m->data[j][i] = aux.data[0][j];
            }
        }
    }
}

/**
 * multiply the corresponding submatrices of A and B.
 */
void process_mult(Matrix *A, Matrix *B, Matrix *C) {
    int r, c, id, k, 
        rbegin, rend, cbegin, cend, // block delimiters
        l, m;
    Matrix sa, sb, sc;

    #pragma omp parallel default(none) private(l, m, r, c, k, rbegin, rend, cbegin, cend, id, sa, sb, sc) shared(A, B, C) num_threads(P)
    {
        id = omp_get_thread_num();
        rbegin = (id / P_SQRT) * BLOCK_SZ;
        rend = rbegin + BLOCK_SZ;

        cbegin = (id % P_SQRT) * BLOCK_SZ;
        cend = cbegin + BLOCK_SZ;

        create_matrix(&sa, BLOCK_SZ, BLOCK_SZ);
        create_matrix(&sb, BLOCK_SZ, BLOCK_SZ);
        create_matrix(&sc, BLOCK_SZ, BLOCK_SZ);

        //copy the blocks for this process
        for(r = rbegin, l = 0; r < rend; r++, l++){
            for(c = cbegin, m = 0; c < cend; c++, m++){
                sa.data[l][m] = A->data[r][c];
                sb.data[l][m] = B->data[r][c];
                sc.data[l][m] = C->data[r][c];
            }
        }
        
        //do multiplication of the block
        for(r = 0; r < BLOCK_SZ; r++){
            for(c = 0; c < BLOCK_SZ; c++){
                for(k = 0; k < BLOCK_SZ; k++){
                    sc.data[r][c] += sa.data[r][k] * sb.data[k][c];
                }
            }
        }

        //put results back to C
        for(r = rbegin, l = 0; r < rend; r++, l++){
            for(c = cbegin, m = 0; c < cend; c++, m++){
                C->data[r][c] = sc.data[l][m];
            }
        }
    }
}

int main() {
    Matrix A, B, C;

    create_matrix(&A, N, N);
    create_matrix(&B, N, N);
    create_matrix(&C, N, N);

    populate_matrix(&A);
    populate_matrix(&B);

    printf("Matrices generadas\n");
    print_matrix(&A, 'A');
    printf("\n\n");
    print_matrix(&B, 'B');
    
    shift_matrix_left(&A, 1);
    shift_matrix_up(&B, 1);

    int i;
    for(i = 0; i < P_SQRT; i++){
        process_mult(&A, &B, &C);
        shift_matrix_left(&A, 0);
        shift_matrix_up(&B, 0);
    }

    printf("\nResultado\n\n");
    print_matrix(&C, 'C');
    return 0;
}

