/*
 * ctm_cuda.h
 *
 *  Created on: 17 avr. 2013
 *      Author: Edward Chou
 */

#ifndef CTM_CUDA_H_
#define CTM_CUDA_H_

#define CELL_TYPE_NORMAL 0
#define CELL_TYPE_INPUT 1
#define CELL_TYPE_OUTPUT 2
#define CELL_TYPE_SWITCH 3

typedef struct Cell {
	int type;
	float capacity;
	float length;
	float rate;
	int access;
}CELL;

//#define LINK_TYPE_STRAIGHT 0
//#define LINK_TYPE_MERGE 1
//#define LINK_TYPE_DIVERGE 2

typedef struct Link {
	float p;
	int c1;
	float p1;
	int c2;
	float p2;
}LINK;

void deleteCudaEnv();
void createCudaEnv(
		CELL *hCells,
		LINK *hLinks,
		float *hPosIn,
		float *hPosOut,
		float *hIn,
		float *hOut,
		int n);
void updateCudaAcc(CELL *hCells);
void loadCudaLen(CELL *hCells);
bool simCuda(float dt);

#endif /* CTM_CUDA_H_ */
