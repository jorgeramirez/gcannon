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
