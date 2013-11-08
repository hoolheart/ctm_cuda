/*
 * ctm_cuda.h
 *
 *  Created on: 17 avr. 2013
 *      Author: Edward Chou
 */

#ifndef CTM_CUDA_H_
#define CTM_CUDA_H_
#include <vector>
#include <string>
using namespace std;

#define CELL_TYPE_NORMAL 0
#define CELL_TYPE_INPUT 1
#define CELL_TYPE_OUTPUT 2

typedef struct ctm_cell {
	int type;
	float rate;
	float cap;
	float length;
	float delay;
}CtmCell;

#define LINK_TYPE_DIRECT 0
#define LINK_TYPE_MERGE 1
#define LINK_TYPE_DIVERGE 2

typedef struct ctm_link {
	int type;
	int cells[3];
	float ratio;
	bool access;
}CtmLink;

class CtmInfo {
public:
	bool is_valid;
	bool is_sim_on;
	bool is_synchronized;
	float vf;
	float w_vf;
	float veh_len;
	float cell_cap;
};

#define LANE_TYPE_NORMAL 0
#define LANE_TYPE_ENTRY 1
#define LANE_TYPE_EXIT 2

class CtmLane {
public:
	string id;
	int type;
	float cap;
	float sat_rate;
	float in_rate;
	float out_ratio;
	int in_cell;
	int out_cell;
	int o_cell;
	int d_cell;
	int in_link;
	int out_link;
	int getHeadCell() const {return (type==LANE_TYPE_EXIT)?out_cell:o_cell;}
	int getTailCell() const {return (type==LANE_TYPE_EXIT)?out_cell:d_cell;}
};

class CtmInnerCell {
public:
	float cap;
	float sat_rate;
	int index;
	CtmInnerCell() {
		cap = 0; sat_rate = 0; index = -1;
	}
	CtmInnerCell(float _cap,float _rate) {
		cap = _cap;
		sat_rate = _rate;
		index = -1;
	}
};

typedef struct ctm_link_info {
	int type;
	float ratio;
	int c1s;
	int c1i;
	int c2s;
	int c2i;
	int c3s;
	int c3i;
}CtmLinkInfo;

class CtmPhase {
public:
	vector<CtmLinkInfo> info;
	int head_link;
	int tail_link;
};

class CtmIntersection {
public:
	string id;
	vector<CtmLane *> in_lanes;
	vector<CtmLane *> out_lanes;
	vector<CtmInnerCell *> inner_cells;
	vector<CtmPhase *> phases;
	int cur_phase;
	vector<int> in_cells;
	vector<int> out_cells;
	~CtmIntersection() {
		for (int i=0;i<(int)inner_cells.size();i++)
			delete inner_cells[i];
		inner_cells.clear();
		for (int i=0;i<(int)phases.size();i++)
			delete phases[i];
		phases.clear();
	}
};

// create and destroy the CUDA environment for CTM simulation
void createCudaEnv(
		CtmCell *hCells, int n_cell,
		CtmLink *hLinks, int n_link);
void deleteCudaEnv();

// data import and export
void loadCudaCells(CtmCell *hCells);
void loadCudaCells(int begin,int len,CtmCell *hCells);
void loadCudaLinks(CtmLink *hLinks);
void loadCudaLinks(int begin,int len,CtmLink *hLinks);
void saveCudaCells(CtmCell *hCells);
void saveCudaLinks(CtmLink *hLinks);

// simulation
bool simCuda(float w_vf,float veh_len,float vf,float dt);
void clearCudaDelay();
void updateCudaAcc(int begin,int end,bool acc);

#endif /* CTM_CUDA_H_ */
