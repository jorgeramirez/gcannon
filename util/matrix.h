/*
 * =====================================================================================
 *
 *       Filename:  matrix.h
 *
 *    Description:  Utility functions for matrices
 *
 *        Version:  1.0
 *       Compiler:  gcc
 *
 *         Author:  Jorge Ramirez <jorgeramirez1990@gmail.com>
 *
 * =====================================================================================
 */
#include <stdio.h>

#define TRUE 1
#define FALSE 0

struct struct_matrix {
    int nrow;
    int ncol;
    int **data;
};

typedef struct struct_matrix Matrix; 

void create_matrix(Matrix *m, int nrow, int ncol);
void populate_matrix(Matrix *m);
void print_matrix(Matrix *m, char iden);
void shift_matrix_left(Matrix *m, int block_sz, int initial);
void shift_matrix_up(Matrix *m, int block_sz, int initial);
void matrix_product(Matrix *c, Matrix *a, Matrix *b);
int* create_array_as_matrix(int r, int c);
void populate_array_as_matrix(int *arr, int r, int c);
int array_as_matrix_equals(int *a, int *b, int r, int c);
