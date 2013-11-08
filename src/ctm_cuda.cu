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
int numLinks;

//__device__ float min(float d1, float d2) {
//    if (d1<d2)
//        return d1;
//    else
//        return d2;
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

__global__ void updateAccess(CtmLink *ListLink, int begin, int end, bool acc) {
    int i = begin + blockDim.x * blockIdx.x + threadIdx.x;

    if (i <= end) {
    	ListLink[i].access = acc;
    }

    __syncthreads();
}

__global__ void calPosFlow(CtmCell *ListCell,
		float *CellPosIn,
		float *CellPosOut,
        float *CellIn,
        float *CellOut,
		float dt, float w_vf, int n) {
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
        default:
            float maxFlow = ListCell[i].rate*dt;
            CellPosIn[i] = min(maxFlow,w_vf*(ListCell[i].cap-ListCell[i].length));
            CellPosOut[i] = min(maxFlow,ListCell[i].length);
            break;
        }
        CellIn[i] = 0;
        CellOut[i] = 0;
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

    if (i<n) {
        if(ListLink[i].access) {
        	int c1 = ListLink[i].cells[0], c2 = ListLink[i].cells[1], c3 = ListLink[i].cells[2];
        	float b1 = ListLink[i].ratio; float b2 = 1-b1;
            float f1,f2;
            switch (ListLink[i].type) {
            case LINK_TYPE_MERGE:
                if (CellPosIn[c3]>=CellPosOut[c1]+CellPosOut[c2]) {
                    f1 = CellPosOut[c1];
                    f2 = CellPosOut[c2];
                }
                else {
                    f1 = mid(CellPosOut[c1],CellPosIn[c3]-CellPosOut[c2],b1*CellPosIn[c3]);
                    f2 = mid(CellPosOut[c2],CellPosIn[c3]-CellPosOut[c1],b2*CellPosIn[c3]);
                }
				CellOut[c1] = f1;
				CellOut[c2] = f2;
				CellIn[c3] = f1+f2;
                break;
            case LINK_TYPE_DIVERGE:
            	f1 = min(CellPosOut[c1],min(CellPosIn[c2]/b1,CellPosIn[c3]/b2));
            	CellOut[c1] = f1;
            	CellIn[c2] = f1*b1;
            	CellIn[c3] = f1*b2;
                break;
            case LINK_TYPE_DIRECT:
            default:
            	f1 = min(CellPosOut[c1],CellPosIn[c2]);
            	CellOut[c1] = f1;
            	CellIn[c2] = f1;
                break;
            }
        }
    }

    __syncthreads();
}

__global__ void updateCells(
        CtmCell *ListCell,
        float *CellIn,
        float *CellOut,
        float dt,
        float veh_len,
        float vf,
        int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<n) {
        ListCell[i].length += CellIn[i]-CellOut[i];
        ListCell[i].delay += dt*(ListCell[i].length-CellOut[i]*ListCell[i].cap*veh_len/vf);
    }

    __syncthreads();
}

__global__ void clearCellDelay(
        CtmCell *ListCell,
        int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<n)
        ListCell[i].delay = 0;

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
		CtmCell *hCells, int n_cell,
		CtmLink *hLinks, int n_link) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // Create CUDA variables
    numCells = n_cell;
    numLinks = n_link;
    err = cudaMalloc((void **)&dCells, n_cell*sizeof(CtmCell));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dLinks, n_link*sizeof(CtmLink));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dPosIn, n_cell*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dPosOut, n_cell*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dIn, n_cell*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&dOut, n_cell*sizeof(float));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy host data to device
    err = cudaMemcpy(dCells, hCells, n_cell*sizeof(CtmCell), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dLinks, hLinks, n_link*sizeof(CtmLink), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)     {
        fprintf(stderr, "Failed to copy host data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void updateCudaAcc(int begin,int end,bool acc) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // check input parameters
    begin = (begin>0)?begin:0; begin = (begin<numLinks)?begin:(numLinks-1);
    end = (end>begin)?end:begin; end = (end<numLinks)?end:(numLinks-1);

    //update access data
    int blocksPerGrid=(end-begin+1 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    updateAccess<<<blocksPerGrid, BLOCK_SIZE>>>(dLinks,begin,end,acc);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to update accesses (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

bool simCuda(float w_vf,float veh_len,float vf,float dt) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // calculate possible maximum flows
    int blocksPerGrid=(numCells + BLOCK_SIZE - 1) / BLOCK_SIZE;
    calPosFlow<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dPosIn,dPosOut,dIn,dOut,dt,w_vf,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to calculate possible maximum flows (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // calculate actual flows
    calFlow<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dLinks,dPosIn,dPosOut,dIn,dOut,numLinks);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to calculate actual flows (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // update cell lengths
    updateCells<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,dIn,dOut,dt,veh_len,vf,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to update the cell lengths (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return true;
}

void clearCudaDelay() {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    // clear delay of all cells
    int blocksPerGrid=(numCells + BLOCK_SIZE - 1) / BLOCK_SIZE;
    clearCellDelay<<<blocksPerGrid, BLOCK_SIZE>>>(dCells,numCells);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to clear cell delays (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void loadCudaCells(CtmCell *hCells) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
    err = cudaMemcpy(dCells, hCells, numCells*sizeof(CtmCell), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to load cell data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void loadCudaCells(int begin,int len,CtmCell *hCells) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
    err = cudaMemcpy(dCells+begin, hCells+begin, len*sizeof(CtmCell), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to load cell data to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void loadCudaLinks(CtmLink *hLinks) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(dLinks, hLinks, numLinks*sizeof(CtmLink), cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to load link data to device (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}

void loadCudaLinks(int begin,int len,CtmLink *hLinks) {
	cudaError_t err = cudaSuccess;
	err = cudaMemcpy(dLinks+begin, hLinks+begin, len*sizeof(CtmLink), cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to load link data to device (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}

void saveCudaCells(CtmCell *hCells) {
	cudaError_t err = cudaSuccess;
    err = cudaMemcpy(hCells, dCells, numCells*sizeof(CtmCell), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to save cell data from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void saveCudaLinks(CtmLink *hLinks) {
	cudaError_t err = cudaSuccess;
    err = cudaMemcpy(hLinks, dLinks, numLinks*sizeof(CtmCell), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to save link data from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}
