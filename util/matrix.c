/*
 * =====================================================================================
 *
 *       Filename:  matrix.c
 *
 *    Description:  implements functions in matrix.h
 *
 *        Version:  1.0
 *       Compiler:  gcc
 *
 *         Author:  Jorge Ramirez <jorgeramirez1990@gmail.com>, 
 *
 * =====================================================================================
 */
#include "matrix.h"
#include <stdlib.h>

void create_matrix(Matrix *m, int nrow, int ncol) { 
    int i;

    m->nrow = nrow; 
    m->ncol = ncol;
    m->data = malloc(nrow *  sizeof(int*));
    for(i = 0; i < ncol; i++){
        m->data[i] = calloc(ncol, sizeof(int));
    }
}

void populate_matrix(Matrix *m) {
    int i, j;
    for(i = 0; i < m->nrow; i++){
        for(j = 0; j < m->ncol; j++){
            m->data[i][j] = rand() % 10;
        }
    }
}

void print_matrix(Matrix *m, char iden) {
    int i, j;

    for(i = 0; i < m->nrow; i++){
        for(j = 0; j < m->ncol; j++){
            printf("%c[%d][%d] = %d  ", iden, i, j , m->data[i][j]);
        }
        printf("\n");
    }
}


//void rsync_submatrix(Matrix *mat, Matrix *sub, int id, int push, int P, int Q, int M, int SIZE) {
//    int r_pq, c_pq, // row and column of PxQ processes matrices
//        r_mn, c_mn, // row and column of MxN block matrices
//        ri, rj, // row and column of input matrices, which are square for simplicity.
//        roffset = SIZE / M, //row offset
//        coffset = SIZE / M, //column offset
//        rbegin, rend, cbegin, cend,
//        r, s;
//
//    r_pq = (id / P);
//    c_pq = (id % Q);
//    
//    int step = 0, subr, subc;
//
//    for(r_mn = r_pq; r_mn < M; r_mn += P){
//        
//        rbegin = r_mn * roffset;
//        rend = rbegin + roffset;
//
//        for(c_mn = c_pq; c_mn < M; c_mn += Q, step++){
//            
//            cbegin = c_mn * coffset;
//            cend = cbegin + coffset;
//
//            for(ri = rbegin, r = 0; ri < rend; ri++, r++){
//                for(rj = cbegin, s = 0; rj < cend; rj++, s++){
//                    
//                    subr = (step / roffset);
//                    subc = (step % coffset);
//
//                    if(push){
//                        mat->data[ri][rj] = sub->data[r + (subr * roffset)][s + (subc * coffset)];
//                    }else{
//                        sub->data[r + (subr * roffset)][s + (subc * coffset)] = mat->data[ri][rj];
//                    }
//                }
//            }
//        }
//    }
//}
//
//
//
//void rsync_local(Matrix *mat, Matrix *sub, int row, int col, int push) {
//    int i, j, r, c,
//        offset = mat->nrow / sub->nrow,
//        ibegin = row * offset,
//        jbegin = col * offset;
//
//    for(i = ibegin, r = 0; i < (ibegin + offset); i++, r++){
//        for(j = jbegin, c = 0; j < (jbegin + offset); j++, c++){
//            if(push){
//                mat->data[i][j] = sub->data[r][c];
//            }else{
//                sub->data[r][c] = mat->data[i][j];
//            }
//        }
//    }
//}
//
//
//void matrix_product(Matrix *c, Matrix *a, Matrix *b){
//    int r, s, k;
//
//    for(r = 0; r < a->nrow; r++){
//        for(s = 0; s < b->ncol; s++){
//            for(k = 0; k < a->ncol; k++){
//                c->data[r][s] += a->data[r][k] * b->data[k][s];
//            }
//        }
//    }
//}
