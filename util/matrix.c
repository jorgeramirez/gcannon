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
#include <stdlib.h>
#include "matrix.h"

void create_matrix(Matrix *m, int nrow, int ncol) { 
    int i;

    m->nrow = nrow; 
    m->ncol = ncol;
    m->data = malloc(nrow *  sizeof(int*));
    for(i = 0; i < ncol; i++){
        m->data[i] = calloc(ncol, sizeof(int));
    }
}

/*int main() {
    Matrix m;
    create_matrix(&m, 4, 4);
    for(int i = 0; i < m.nrow; i++){
        for(int j = 0; j < m.ncol; j++){
            printf("c[%d][%d] = %d", i, j , m.data[i][j]);
        }
        printf("\n");
    }
    return 0;
}*/
