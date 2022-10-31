#include "mpi.h"
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>

void initmatrixbyrows(double* m, int n)
{
    for (int i = 0; i < n * n; i++)
        m[i] = (double)rand() / double(1000);
}

void initmatrixbycols(double* m, int n)
{
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            m[j + i * n] = (double)rand() / double(1000);
}

void printmatrix(double* m, int r, int c)
{
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            fprintf(stdout, "%f \t", m[j + i * c]);
            fflush(stdout);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    fprintf(stdout, "\n\n");
    fflush(stdout);
}

void P_multiplication(double* bufA, double* bufB, double* bufC, int* fragmentation, int* elems2send, int* col_indexes, int n, int numprocs, int myid)
{
    MPI_Status status;
    int i, j, k, dest = 0, source = 0, shift = 0;
    double temp = 0.0;


    for (int p = 0; p < numprocs; p++)
    {
        shift = (myid + p) % numprocs;
        for (i = 0; i < fragmentation[myid]; i++)
            for (j = 0; j < fragmentation[shift]; j++)
            {
                for (k = 0; k < n; k++) temp += bufA[k + i * n] * bufB[k + j * n];
                bufC[j + i * n + col_indexes[shift]] = temp;
                temp = 0.0;
            }

        dest = (myid == 0) ? (numprocs - 1) : (myid - 1);
        source = (myid + 1) % numprocs;
        MPI_Sendrecv_replace(bufB, elems2send[0], MPI_DOUBLE, dest, 0, source, 0, MPI_COMM_WORLD, &status);


    }

}

void S_multiplication(double* A, double* B, double* Cs, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                Cs[j + i * n] += A[k + i * n] * B[k + j * n];
}

void init_and_distribute(double*& A, double*& B, double*& Cs, double*& Cp, double*& bufA, double*& bufB, double*& bufC, int*& fragmentation, int*& elems2send, int*& start_indexes, int*& col_indexes, int n, int numprocs, int myid)
{
    A = new double[n * n];
    B = new double[n * n];
    Cs = new double[n * n];
    Cp = new double[n * n];
    fragmentation = new int[numprocs];
    elems2send = new int[numprocs];
    start_indexes = new int[numprocs];
    col_indexes = new int[numprocs];

    if (myid == 0)
    {
        initmatrixbyrows(A, n);
        initmatrixbycols(B, n);
    }

    for (int i = 0; i < numprocs; i++)
        fragmentation[i] = 0;

    for (int i = 0; i < n;)
        for (int j = 0; (j < numprocs) && (i < n); j++, i++)
            fragmentation[j]++;

    for (int i = 0; i < numprocs; i++)
        elems2send[i] = fragmentation[i] * n;

    start_indexes[0] = 0;
    for (int i = 1; i < numprocs; i++)
        start_indexes[i] = start_indexes[i - 1] + elems2send[i - 1];

    for (int i = 0; i < numprocs; i++)
        col_indexes[i] = start_indexes[i] / n;

    bufA = new double[elems2send[0]];
    bufB = new double[elems2send[0]];
    bufC = new double[elems2send[0]];

    for (int i = 0; i < elems2send[0]; i++)
        bufC[i] = 0;

    for (int i = 0; i < n * n; i++)
        Cs[i] = 0;

    MPI_Scatterv(A, elems2send, start_indexes, MPI_DOUBLE, bufA, elems2send[myid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, elems2send, start_indexes, MPI_DOUBLE, bufB, elems2send[myid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void deinit(double*& A, double*& B, double*& Cs, double*& Cp, double*& bufA, double*& bufB, double*& bufC, int*& fragmentation, int*& elems2send, int*& start_indexes, int*& col_indexes)
{
    delete[] A; A = NULL;
    delete[] B; B = NULL;
    delete[] Cs; Cs = NULL;
    delete[] Cp; Cp = NULL;
    delete[] bufA; bufA = NULL;
    delete[] bufB; bufB = NULL;
    delete[] bufC; bufC = NULL;
    delete[] fragmentation; fragmentation = NULL;
    delete[] elems2send; elems2send = NULL;
    delete[] start_indexes; start_indexes = NULL;
    delete[] col_indexes; col_indexes = NULL;
}


int main(int argc, char* argv[])
{
    double* A = NULL;
    double* B = NULL;
    double* Cs = NULL;
    double* Cp = NULL;
    double* bufA = NULL;
    double* bufB = NULL;
    double* bufC = NULL;
    int n = 1000;
    int myid, numprocs=9;
    int* fragmentation = NULL;
    int* elems2send = NULL;
    int* start_indexes = NULL;
    int* col_indexes = NULL;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    double t1 = .0, t2 = .0, p_dt, s_dt, ac;

    if (myid == 0)
        scanf_s("%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


    init_and_distribute(A, B, Cs, Cp, bufA, bufB, bufC, fragmentation, elems2send, start_indexes, col_indexes, n, numprocs, myid);

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    P_multiplication(bufA, bufB, bufC, fragmentation, elems2send, col_indexes, n, numprocs, myid);
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();
    p_dt = t2 - t1;

    MPI_Gatherv(bufC, elems2send[myid], MPI_DOUBLE, Cp, elems2send, start_indexes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    t1 = MPI_Wtime();
    S_multiplication(A, B, Cs, n);
    t2 = MPI_Wtime();
    s_dt = t2 - t1;
    ac = s_dt / p_dt;


    if (myid == 0)
    {
        int flag = 1;
        for (int i = 0; i < n * n; i++)
        {
            if (Cp[i] != Cs[i]) { flag = 0; break; }
        }

        if (flag)
        {
            fprintf(stdout, "OK\n\n parallel %f\n serial %f\n speedup %f\n", p_dt, s_dt, ac);
            fflush(stdout);
        }
        else
        {
            fprintf(stdout, "NOT OK");
            fflush(stdout);
        }
    }

    deinit(A, B, Cs, Cp, bufA, bufB, bufC, fragmentation, elems2send, start_indexes, col_indexes);
    MPI_Finalize();
    return 0;
}
