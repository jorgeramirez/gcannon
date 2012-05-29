#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "../util/matrix.h"


#define N 320
#define BLOCK_SZ 10

void imprimematriz(int *matriz, int n) {
	int i;
	int j;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			printf("%d\t", matriz[i*n+j]);
		}
		printf("\n");
	}
}


void copy_block(int *src, int *dest, int n) {
    int i;
    for(i = 0; i < n*n; i++){
        dest[i] = src[i];
    }
}


void copy_blocks(int **src, int **dest, int m, int n){
    int i, j;
    for(i = 0; i < m; i++){
        for(j = 0; j < n*n; j++){
            dest[i][j] = src[i][j];
        }
    }
}


int** create_block(int blocks_pp, int block_sz) {
    int i, j, k, **a;

    a = malloc(sizeof(int *)* blocks_pp);
    for(i = 0; i < blocks_pp; i++){
        a[i] = malloc(sizeof(int) * block_sz * block_sz);
    }
    return a;
}


void shift_left(int **buf, int bpp, int bz){
    int i, j, k, step, r, s, pos;

    step = (int)sqrt(bpp);
    int aux[step];
    for(i = 0, s = 0; i < bpp; i += step, s++){    
        for(k = 0; k < bz*bz; k++){
            for(j = i, r = 0; j < (i + step); j++, r++){
                
                pos = j + 1;
                if(pos >= (step + i)){
                    pos = s * step;
                }
                aux[r] = buf[pos][k];           
            }
            for(j = i, r = 0; j < i + step; j++, r++){
                buf[j][k] = aux[r];           
            }
        }  
    }
}

void shift_up(int **buf, int bpp, int bz){

    int i, j, k, step, r, s, pos;
    step = (int)sqrt(bpp);
    int aux[bpp];
    
    for(k = 0; k < bz*bz; k++){
        for(i = 0; i < bpp; i++){
            aux[i] = buf[(i+step) % bpp][k];
        }

        for(i = 0; i < bpp; i++){
            buf[i][k] = aux[i];   
        }
    }
}


void print_blocks(int **buf, int bpp, int bz){
    int i, j, k;
    for(i = 0; i < bpp; i++){
        printf("bloque %d\n", i);
        for(j = 0; j < bz*bz; j++){
            printf("%d ", buf[i][j]);
        }
        printf("\n");
    }
}

/* This function performs a serial matrix-matrix multiplication c = a*b */
MatrixMultiply(int n, int *a, int *b, int *c)
{
	int i, j, k;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				c[i * n + j] += a[i * n + k] * b[k * n + j];
}


MatrixMatrixMultiply(int n, int **a, int **b, int **c, int *c_final, 
                     int blocks_per_process, int block_sz, MPI_Comm comm) {
	int i;
	int npes, dims[2], periods[2];
	int rank, my2drank, mycoords[2];
	int uprank, downrank, leftrank, rightrank, coords[2];
	int shiftsource, shiftdest;
	MPI_Status status;
	MPI_Comm comm_2d;

	/* Get the communicator related information */
	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &rank);

	/* Set up the Cartesian topology */
	dims[0] = dims[1] = sqrt(npes);

	/* Set the periods for wraparound connections */
	periods[0] = periods[1] = 1;

	/* Create the Cartesian topology, with rank reordering */
	MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

	/* Get the rank and coordinates with respect to the new topology */
	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	/* Compute ranks of the up and left shifts */
	MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    //initial shift
    int **asend_buffer = create_block(blocks_per_process, block_sz);
    int **bsend_buffer = create_block(blocks_per_process, block_sz);
    int **bcopy = create_block(blocks_per_process, block_sz);
    
    copy_blocks(&b[0], &bcopy[0], blocks_per_process, block_sz);
    
    int j;
    int limit = (int)sqrt(blocks_per_process);

    int asend[block_sz * block_sz], arecv[block_sz * block_sz];
    int bsend[block_sz * block_sz], brecv[block_sz * block_sz];
    
    copy_blocks(&a[0], &asend_buffer[0], blocks_per_process, block_sz);
    
    for(i = 0; i < limit; i++){
        
        copy_blocks(&bcopy[0], &bsend_buffer[0], blocks_per_process, block_sz);
        int rdisp = (mycoords[0] + i * sqrt(npes));

        if((mycoords[1] - rdisp) < 0){
            shift_left(&asend_buffer[0], blocks_per_process, block_sz);
        }
        
        for(j = 0; j < limit; j++){

            int cdisp = (mycoords[1] + j * sqrt(npes));
            
            if((mycoords[0] - cdisp) < 0){
                shift_up(&bsend_buffer[0], blocks_per_process, block_sz);
            }

            copy_block(&asend_buffer[i * limit + j][0], &asend[0], block_sz);
            copy_block(&bsend_buffer[i * limit + j][0], &bsend[0], block_sz);

            MPI_Cart_shift(comm_2d, 1, -rdisp, &shiftsource, &shiftdest); //A
            MPI_Sendrecv(asend, block_sz*block_sz, MPI_INT, shiftdest, 1, arecv, block_sz*block_sz, MPI_INT, shiftsource, 1, comm_2d, &status); //a

            MPI_Cart_shift(comm_2d, 0, -cdisp,  &shiftsource, &shiftdest); // B
            MPI_Sendrecv(bsend, block_sz*block_sz, MPI_INT, shiftdest, 1, brecv, block_sz*block_sz, MPI_INT, shiftsource, 1, comm_2d, &status); //b
            
            copy_block(&arecv[0], &a[i * limit + j][0], block_sz);
            copy_block(&brecv[0], &b[i * limit + j][0], block_sz);
        }
    }

    //main loop    
    int sqrtp = (int)sqrt(npes);
    int m = n / block_sz;
    int max = m / sqrtp;
    int t, l, j_prima, i_prima;

    for (t = 0; t < sqrtp; t++){
        for (i = 0; i < max; i++){
            for (j = 0; j < max; j++){
                for (l = 0; l < max; l++){
                    j_prima = (j % m + l) % (m / sqrtp);
                    i_prima = (i % m + l) %(m / sqrtp);
                    MatrixMultiply(block_sz, &a[i * max + j_prima][0], &b[i_prima * max + j][0], &c[i * max + j][0]);    
                }
            }
        }

        // Do shift. First A, then B
        if(mycoords[1] == 0){ 
            shift_left(&a[0], blocks_per_process, block_sz);
        }
        /* Shift matrix a left by one */
        for(i = 0; i < blocks_per_process; i++){
            MPI_Sendrecv_replace(a[i], block_sz*block_sz, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, &status);
        }

        if(mycoords[0] == 0){
            shift_up(&b[0], blocks_per_process, block_sz);
        }

        /* Shift matrix b up by one */
        for(i = 0; i < blocks_per_process; i++){
            MPI_Sendrecv_replace(b[i], block_sz*block_sz, MPI_INT, uprank, 1, downrank, 1, comm_2d, &status);
        }
    }
        
	MPI_Comm_free(&comm_2d); /* Free up communicator */

    //reduce
	int proc_num, block_num, pos;
	proc_num = rank;

    int col_proc = proc_num % (int)sqrt(npes);
	int row_proc = (proc_num - col_proc) / sqrt(npes);

	for(block_num = 0; block_num < blocks_per_process; block_num++)
	{
       	int col_block = block_num % (int)sqrt(blocks_per_process);
		int row_block = (block_num - col_block) / sqrt(blocks_per_process);

		for(pos = 0; pos < block_sz * block_sz; pos++)
		{
			int col_pos = pos % block_sz;
			int row_pos = (pos - col_pos) / block_sz;

			int big_row = row_block * sqrt(npes) + row_proc;
			int big_col = col_block * sqrt(npes) + col_proc;
			c_final[big_row * N * block_sz + row_pos * N + big_col * block_sz + col_pos] = c[block_num][pos];
		}
	}

	if(rank != 0)
	{
		MPI_Reduce(c_final, c_final, npes * blocks_per_process * block_sz * block_sz, MPI_INT, MPI_SUM, 0, comm);
	}
	else
	{
		MPI_Reduce(MPI_IN_PLACE, c_final, npes * blocks_per_process * block_sz * block_sz, MPI_INT, MPI_SUM, 0, comm);
	}

}


int main(int argc, char *argv[]) {

    int *a = create_array_as_matrix(N, N);
    int *b = create_array_as_matrix(N, N);
    int *c = create_array_as_matrix(N, N);

    populate_array_as_matrix(&a[0], N, N);
    populate_array_as_matrix(&b[0], N, N);

    int block_sz = BLOCK_SZ;
	MPI_Init(&argc,&argv);
	int nro_procesos;
	MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);
    
    int P = sqrt(nro_procesos);
	int M = (N / block_sz);
    int blocks_per_process = (M / P) * (M / P);
    
    int*** bloque_a = malloc(sizeof(int**) * nro_procesos);
    int*** bloque_b = malloc(sizeof(int**) * nro_procesos);
    int*** bloque_c = malloc(sizeof(int**) * nro_procesos);
    int i; 
    for(i = 0; i < nro_procesos; i++){
        bloque_a[i] = create_block(blocks_per_process, block_sz);
        bloque_b[i] = create_block(blocks_per_process, block_sz);
        bloque_c[i] = create_block(blocks_per_process, block_sz);
    }

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int r_pq, c_pq, // row and column of PxQ processes matrices
        r_mn, c_mn, // row and column of MxN block matrices
        ri, rj, // row and column of input matrices, which are square for simplicity.
        offset = N / M, //row offset
        rbegin, rend, cbegin, cend,
        r, s;

    r_pq = (rank / P);
    c_pq = (rank % P);
    
    int step = 0;

    for(r_mn = r_pq; r_mn < M; r_mn += P){
        
        rbegin = r_mn * offset;
        rend = rbegin + offset;

        for(c_mn = c_pq; c_mn < M; c_mn += P, step++){
            
            cbegin = c_mn * offset;
            cend = cbegin + offset;

            for(ri = rbegin, r = 0; ri < rend; ri++, r++){
                for(rj = cbegin, s = 0; rj < cend; rj++, s++){
                    bloque_a[rank][step][r * block_sz + s] = a[ri * N + rj];
                    bloque_b[rank][step][r * block_sz + s] = b[ri * N + rj];
                    bloque_c[rank][step][r * block_sz + s] = 0;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    
    MatrixMatrixMultiply(N, &bloque_a[rank][0], &bloque_b[rank][0], &bloque_c[rank][0], &c[0], blocks_per_process, block_sz, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    
    if(rank == 0){
        // Serial multiplication. In order to check the results.
        int *d = create_array_as_matrix(N, N);
        MatrixMultiply(N, a, b, d);
        
        int equal = array_as_matrix_equals(&d[0], &c[0], N, N);
        if(equal){
            printf("\nSon iguales\n");
        }

        printf("\nTiempo: %.4f segundos\n", (end - start));
    }

    MPI_Finalize();
	return 0;
}

