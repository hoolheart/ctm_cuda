/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "ctm_cuda.h"

#define BLOCK_SIZE 256

CtmCell *dCells;
CtmLink *dLinks;
float *dPosIn;
float *dPosOut;
float *dIn;
float *dOut;
int numCells;

//__device__ float min(float d1, float d2) {
//	if (d1<d2)
//		return d1;
//	else
//		return d2;
//}

__device__ float mid(float d1, float d2, float d3) {
	if (d1<=d2) {
		if (d2<=d3)
			return d2;
		else {
			if (d1<=d3)
				return d3;
			else
				return d1;
		}
	}
	else {
		if (d1<=d3)
			return d1;
		else {
			if (d2<=d3)
				return d3;
			else
				return d2;
		}
	}
}

__global__ void updateAccess(CtmCell *c, int *a, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n) {}

    __syncthreads();
}

__global__ void calPosFlow(CtmCell *ListCell, float *CellPosIn, float *CellPosOut, float dt, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i<n) {
		switch (ListCell[i].type) {
		case CELL_TYPE_INPUT:
			CellPosIn[i] = 0;
			ListCell[i].length += ListCell[i].rate*dt;
			CellPosOut[i] = ListCell[i].length;
			break;
		case CELL_TYPE_OUTPUT:
			CellPosIn[i] = ListCell[i].rate*dt;
			CellPosOut[i] = 0;
			break;
		case CELL_TYPE_NORMAL:
			break;
		default:
			break;
		}
    }

    __syncthreads();
}

__global__ void calFlow(
		CtmCell *ListCell,
		CtmLink *ListLink,
		float *CellPosIn,
		float *CellPosOut,
		float *CellIn,
		float *CellOut,
		int n) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i<n)
		CellOut[i] = 0;
	if (i<n) {
	}

    __syncthreads();
}

__global__ void updateCells(
		CtmCell *ListCell,
		float *CellIn,
		float *CellOut,
		int n) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i<n)
		ListCell[i].length += CellIn[i]-CellOut[i];

    __syncthreads();
}

void deleteCudaEnv() {
	cudaFree(dCells);
	cudaFree(dLinks);
	cudaFree(dPosIn);
	cudaFree(dPosOut);
	cudaFree(dIn);
	cudaFree(dOut);
}

void createCudaEnv(
		CtmCell *hCells,
		CtmLink *hLinks,
		float *hPosIn,
		float *hPosOut,
		float *hIn,
		float *hOut,
		int n) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // Create CUDA variables
    numCells = n;
    err = cudaMalloc((void **)&dCells, n*sizeof(CtmCell));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dLinks, n*sizeof(CtmLink));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dPosIn, n*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dPosOut, n*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dIn, n*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dOut, n*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy host data to device
    err = cudaMemcpy(dCells, hCells, n*sizeof(CtmCell), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dLinks, hLinks, n*sizeof(CtmLink), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dPosIn, hPosIn, n*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dPosOut, hPosOut, n*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dIn, hIn, n*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dOut, hOut, n*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void updateCudaAcc(CtmCell *hCells) {
	// Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // size of vector of accesses
    size_t size = numCells * sizeof(int);

    // Allocate the host input vector A
    int *h_A = (int *)malloc(size);

    // Verify that allocations succeeded
    if (h_A == NULL)
    {
        fprintf(stderr, "Failed to allocate host vectors!\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the host access data
    for (int i = 0; i < numCells; ++i) {
    }

    // Allocate the device access data
    int *d_A = NULL;
    err = cudaMalloc((void **)&d_A, size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the host access data to the device
    err = cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector A from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    //update access data
    int blocksPerGrid=(numCells + BLOCK_SIZE - 1) / BLOCK_SIZE;
    updateAccess<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,d_A,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // free all temporary data
    cudaFree(d_A);
    free(h_A);
}

void loadCudaLen(CtmCell *hCells) {
	// Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // Copy the device cell data to the host
    err = cudaMemcpy(hCells, dCells, numCells*sizeof(CtmCell), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy cell data from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

bool simCuda(float dt) {
	// Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // calculate possible maximum flows
    int blocksPerGrid=(numCells + BLOCK_SIZE - 1) / BLOCK_SIZE;
    calPosFlow<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dPosIn,dPosOut,dt,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to calculate possible maximum flows (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // calculate actual flows
    calFlow<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dLinks,dPosIn,dPosOut,dIn,dOut,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to calculate actual flows (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // update cell lengths
    updateCells<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dIn,dOut,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to update the cell lengths (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

	return true;
}
