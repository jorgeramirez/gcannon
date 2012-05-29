#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "../util/matrix.h"

#define N 320

void imprimematriz(int* matriz, int n);

MatrixMatrixMultiply(int n, int *a, int *b, int *c, int *c_grande, MPI_Comm comm)
{
	int i;
	int nlocal;
	int npes, dims[2], periods[2];
	int myrank, my2drank, mycoords[2];
	int uprank, downrank, leftrank, rightrank, coords[2];
	int shiftsource, shiftdest;
	MPI_Status status;
	MPI_Comm comm_2d;

	/* Get the communicator related information */
	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &myrank);

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

	/* Determine the dimension of the local matrix block */
	nlocal = n/dims[0];

	/* Perform the initial matrix alignment. First for A and then for B */
	MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	/* Get into the main computation loop */
	for (i=0; i<dims[0]; i++)
	{
		MatrixMultiply(nlocal, a, b, c); /*c=c+a*b*/
		/* Shift matrix a left by one */
		MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, &status);
		/* Shift matrix b up by one */
		MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_INT, uprank, 1, downrank, 1, comm_2d, &status);
	}

	/* Restore the original distribution of a and b */
	MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_INT,shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Comm_free(&comm_2d); /* Free up communicator */
    
    int ind;
	int n_grande = n;
	int ind_col;
	for(ind = 0; ind < nlocal; ind++)
	{
		for(ind_col = 0; ind_col < nlocal; ind_col++)
		{
			int fila_grande = mycoords[0]*nlocal+ind;
			int columna_grande = mycoords[1]*nlocal+ind_col;
			c_grande[fila_grande*n_grande+columna_grande] = c[ind*nlocal+ind_col];
		}
	}

	if(myrank != 0)
	{
		MPI_Reduce(c_grande, c_grande, n_grande*n_grande, MPI_INT, MPI_SUM, 0, comm);
	}
	else
	{
		MPI_Reduce(MPI_IN_PLACE, c_grande, n_grande*n_grande, MPI_INT, MPI_SUM, 0, comm);
	}
}

/* This function performs a serial matrix-matrix multiplication c = a*b */
MatrixMultiply(int n, int *a, int *b, int *c)
{
	int i, j, k;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			for (k=0; k<n; k++)
				c[i*n+j] += a[i*n+k]*b[k*n+j];
}

void imprimematriz(int *matriz, int n)
{
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

int main(int argc, char *argv[])
{

    int *a = create_array_as_matrix(N, N);
    int *b = create_array_as_matrix(N, N);
    int *c = create_array_as_matrix(N, N);
    
    populate_array_as_matrix(&a[0], N, N);
    populate_array_as_matrix(&b[0], N, N);

    MPI_Init(&argc,&argv);
	int nro_procesos;
	MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);

	int n_chica = N/sqrt(nro_procesos);

	int max_fila_bloque = n_chica;
	int max_columna_bloque = n_chica;

	int fila_bloque;
	int columna_bloque;
	int fila;
	int columna;

	int n_bloques = N/n_chica;

	int bloque_a[nro_procesos][n_chica*n_chica];
	int bloque_b[nro_procesos][n_chica*n_chica];
	int bloque_c[nro_procesos][n_chica*n_chica];

	int mi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mi_rank);

	columna_bloque = mi_rank % n_bloques;
	fila_bloque = (mi_rank - columna_bloque) / n_bloques;

	int indice_bloque = 0;
	
	for(fila = fila_bloque * n_chica; fila < fila_bloque * n_chica + n_chica; fila++)
	{
		for(columna = columna_bloque * n_chica; columna < columna_bloque * n_chica + n_chica; columna++)
		{
			bloque_a[mi_rank][indice_bloque] = a[fila * N + columna];
			bloque_b[mi_rank][indice_bloque] = b[fila * N + columna];
			bloque_c[mi_rank][indice_bloque] = 0;
			indice_bloque++;
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
	
    MatrixMatrixMultiply(N, bloque_a[mi_rank], bloque_b[mi_rank], bloque_c[mi_rank], &c[0], MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    
    if(mi_rank == 0){
        // Serial multiplication. In order to check te results.
        int *d = create_array_as_matrix(N, N);
        MatrixMultiply(N, a, b, d);
//        imprimematriz(&d[0], N);
        printf("\n\n");
//        imprimematriz(&c[0], N);
        
        int equal = array_as_matrix_equals(&d[0], &c[0], N, N);
        if(equal){
            printf("\nSon iguales\n");
        }

        printf("\nTiempo: %.4f segundos\n", (end - start));
    }

	MPI_Finalize();
	return 0;
}

