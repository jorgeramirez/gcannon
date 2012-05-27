/*
 * =====================================================================================
 *
 *       Filename:  gcannon.c
 *
 *    Description:  Generalized Cannon's Algorithm
 *
 *        Version:  1.0
 *       Compiler:  gcc
 *
 *         Author:  Jorge Ramirez <jorgeramirez1990@gmail.com>
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <math.h>
#include "../util/matrix.h"
#include "omp.h"

#define P 2  //physical processes array P x Q
#define Q 2
#define M 4  //virtual processes array M x N
#define N 4
#define K 4

#define A_RN 8 //matrix A of A_RN x A_CN
#define A_CN 8 
#define B_RN 8 //matrix B of B_RN x B_CN
#define B_CN 8


#define BLOCK_SZ (A_RN / M) //block size for A and B, since we asume their are square matrices for simplicity.

#define FALSE 0
#define TRUE 1

int lcm(a, b) {
    int n;
    for(n = 1;; n++) {
        if(n % a == 0 && n % b == 0)
            return n;
    }
}


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
    for(k = 0, s = 0; k < A_CN; k += BLOCK_SZ, s++){
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
    for(k = 0, s = 0; k < A_RN; k += BLOCK_SZ, s++){
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


void rsync_submatrix(Matrix *mat, Matrix *sub, int id, int push) {
    int r_pq, c_pq, // row and column of PxQ processes matrices
        r_mn, c_mn, // row and column of MxN block matrices
        ri, rj, // row and column of input matrices, which are square for simplicity.
        roffset = A_RN / M, //row offset
        coffset = A_CN / N, //column offset
        rbegin, rend, cbegin, cend,
        r, s;

    r_pq = (id / P);
    c_pq = (id % Q);
    
    int step = 0, subr, subc;

    for(r_mn = r_pq; r_mn < M; r_mn += P){
        
        rbegin = r_mn * roffset;
        rend = rbegin + roffset;

        for(c_mn = c_pq; c_mn < N; c_mn += Q, step++){
            
            cbegin = c_mn * coffset;
            cend = cbegin + coffset;

            for(ri = rbegin, r = 0; ri < rend; ri++, r++){
                for(rj = cbegin, s = 0; rj < cend; rj++, s++){
                    
                    subr = (step / roffset);
                    subc = (step % coffset);

                    if(push){
                        mat->data[ri][rj] = sub->data[r + (subr * roffset)][s + (subc * coffset)];
                    }else{
                        sub->data[r + (subr * roffset)][s + (subc * coffset)] = mat->data[ri][rj];
                    }
                }
            }
        }
    }
}



void rsync_local(Matrix *mat, Matrix *sub, int row, int col, int push) {
    int i, j, r, c,
        offset = mat->nrow / sub->nrow,
        ibegin = row * offset,
        jbegin = col * offset;

    for(i = ibegin, r = 0; i < (ibegin + offset); i++, r++){
        for(j = jbegin, c = 0; j < (jbegin + offset); j++, c++){
            if(push){
                mat->data[i][j] = sub->data[r][c];
            }else{
                sub->data[r][c] = mat->data[i][j];
            }
        }
    }
}


void matrix_product(Matrix *c, Matrix *a, Matrix *b){
    int r, s, k;

    for(r = 0; r < a->nrow; r++){
        for(s = 0; s < b->ncol; s++){
            for(k = 0; k < a->ncol; k++){
                c->data[r][s] += a->data[r][k] * b->data[k][s];
            }
        }
    }
}


int main() {
    int id, sarn, sacn, sbcn, sbrn, LCM, t, i, j, l, jp, ip;
    Matrix A, B, C, // global matrices
           sa, sb, sc, //local matrices
           subsa, subsb, subsc; //local submatrices

    create_matrix(&A, A_RN, A_CN);
    create_matrix(&B, B_RN, B_CN);
    create_matrix(&C, A_RN, B_CN);

    populate_matrix(&A);
    populate_matrix(&B);

    printf("\nMatrices generadas\n");
    print_matrix(&A, 'A');
    printf("\n");
    print_matrix(&B, 'B');

    //submatrices sizes
    sarn = (A_RN / M) * (M / P);
    sacn = (A_CN / K) * (K / Q);
    sbrn = (B_RN / K) * (K / P);
    sbcn = (B_CN / N) * (N / Q);
    
    int local_block = sarn / (A_RN / M);
    
    LCM = lcm(P, Q);
    shift_matrix_left(&A, 1);
    shift_matrix_up(&B, 1);

    #pragma omp parallel default(none) shared(A, B, C, sarn, sacn, sbrn, sbcn, LCM, local_block) \
                                       private(sa, sb, sc, id, t, i, j, l, jp, ip, subsa, subsb, subsc) num_threads(P * Q)
    {
        id = omp_get_thread_num();

        create_matrix(&sa, sarn, sacn);
        create_matrix(&sb, sbrn, sbcn);
        create_matrix(&sc, sarn, sbcn);
        
        for(t = 0; t < LCM; t++){
            
            create_matrix(&subsa, local_block, local_block);
            create_matrix(&subsb, local_block, local_block);
            create_matrix(&subsc, local_block, local_block);
            
            rsync_submatrix(&A, &sa, id, FALSE);
            rsync_submatrix(&B, &sb, id, FALSE);
            
            for(i = 0; i < (M / P); i++){
                for(j = 0; j < (N / Q); j++){
                    for(l = 0; l < (K / LCM); l++){
                        jp = (j % K + l * LCM / Q) % (K / Q);
                        ip = (i % K + l * LCM / P) % (K / P);


                        rsync_local(&sc, &subsc, i, j, FALSE);
                        rsync_local(&sa, &subsa, i, jp, FALSE);
                        rsync_local(&sb, &subsb, ip, j, FALSE);
                        
                        matrix_product(&subsc, &subsa, &subsb);

                        rsync_local(&sc, &subsc, i, j, TRUE);
                    }
                }
            }
            
            rsync_submatrix(&C, &sc, id, TRUE);
            
            #pragma omp barrier

            #pragma omp single
            {
                shift_matrix_left(&A, 0);
                shift_matrix_up(&B, 0);
            }
        }
    }
    
    printf("\nResultado\n");
    print_matrix(&C, 'C');
    return 0;
}
