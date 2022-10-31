#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
void DummyDataInitialization(double* pAMatrix, double* pBMatrix, int Size) {
	int i, j; 
	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++) {
			pAMatrix[i * Size + j] = 1;
			pBMatrix[i * Size + j] = 1;
		}
}
void RandomDataInitialization(double* pAMatrix, double* pBMatrix,
	int Size) {
	int i, j; 
	srand(unsigned(clock()));
	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++) {
			pAMatrix[i * Size + j] = rand() / double(1000);
			pBMatrix[i * Size + j] = rand() / double(1000);
		}
}
void ProcessInitialization(double*& pAMatrix, double*& pBMatrix,
	double*& pCMatrix, int& Size) {
	do {
		printf("\nEnter the size of matrices: ");
		scanf_s("%d", &Size);
		printf("\nChosen matrices' size = %d\n", Size);
		if (Size <= 0)
			printf("\nSize of objects must be greater than 0!\n");
	} while (Size <= 0);
	pAMatrix = new double[Size * Size];
	pBMatrix = new double[Size * Size];
	pCMatrix = new double[Size * Size];
	DummyDataInitialization(pAMatrix, pBMatrix, Size);
	for (int i = 0; i < Size * Size; i++) {
		pCMatrix[i] = 0;
	}
}
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	int i, j; 
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}
void SerialResultCalculation(double* pAMatrix, double* pBMatrix,
	double* pCMatrix, int Size) {
		int i, j, k;
	for (i = 0; i < Size; i++) {
		for (j = 0; j < Size; j++)
			for (k = 0; k < Size; k++)
				pCMatrix[i * Size + j] += pAMatrix[i * Size + k] * pBMatrix[k * Size + j];
	}
}
void ProcessTermination(double* pAMatrix, double* pBMatrix,
	double* pCMatrix) {
	delete[] pAMatrix;
	delete[] pBMatrix;
	delete[] pCMatrix;
}
void main() {
	double* pAMatrix; 
	double* pBMatrix; 
	double* pCMatrix; 
	int Size; 
	time_t start, finish;
	double duration;
	printf("Serial matrix multiplication program\n");
	ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, Size);
	//printf("Initial A Matrix \n");
	//PrintMatrix(pAMatrix, Size, Size);
	//printf("Initial B Matrix \n");
	//PrintMatrix(pBMatrix, Size, Size);
	start = clock();
	SerialResultCalculation(pAMatrix, pBMatrix, pCMatrix, Size);
	finish = clock();
	duration = (finish - start) / double(CLOCKS_PER_SEC);

	//printf("\n Result Matrix: \n");
	//PrintMatrix(pCMatrix, Size, Size);
	printf("\n Time of execution: %f\n", duration);
	ProcessTermination(pAMatrix, pBMatrix, pCMatrix);
}