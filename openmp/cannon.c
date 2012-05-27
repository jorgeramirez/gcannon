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
        
        matrix_product(&sc, &sa, &sb);

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
    
    shift_matrix_left(&A, BLOCK_SZ, 1);
    shift_matrix_up(&B, BLOCK_SZ, 1);

    int i;
    for(i = 0; i < P_SQRT; i++){
        process_mult(&A, &B, &C);
        shift_matrix_left(&A, BLOCK_SZ, 0);
        shift_matrix_up(&B, BLOCK_SZ, 0);
    }

    printf("\nResultado\n\n");
    print_matrix(&C, 'C');
    return 0;
}

