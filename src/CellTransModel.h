/*
 * CellTransModel.h
 *
 *  Created on: 26 mars 2013
 *      Author: Edward
 */

#ifndef CELLTRANSMODEL_H_
#define CELLTRANSMODEL_H_
#include <string>
#include <sstream>
using namespace std;

/*#define CELL_TYPE_NORMAL 0
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
}LINK;*/

#include "ctm_cuda.h"

class CellTransModel {
private:
	CELL *ListCell;
	LINK *ListLink;
	float *CellPosIn;
	float *CellPosOut;
	float *CellIn;
	float *CellOut;
	int n;
	string err;
	bool isSimOn;
public:
	CellTransModel();
	CellTransModel(int _n);
	virtual ~CellTransModel();
	string getErr() const {return err;}
	void cleanErr() {err="";}
	string getErrClean() {
		string t = err;
		err = "";
		return t;
	}
	void resetModel();
	bool initialModel(int _n, float _cap=10, float _s=0.5);
	bool setCell(int i, int _type, float _cap, float _s);
	bool setLink(int i, float _p, int _c1, float _p1=1, int _c2=0, float _p2=0);
	bool startSim(float const _len[], int const _acc[]);
	bool sim(float dt, int steps=1);
	bool changeAccess(int i, int _acc);
	bool changeAccesses(int const _acc[]);
	bool getCurrentLengths(float _len[]);
	void stopSim();
	bool resumeSim();
private:
	void calPosFlows(float dt);
	void calRealFlows();
	bool updateCells();
	float min(float d1, float d2);
	float mid(float d1, float d2, float d3);
};

#endif /* CELLTRANSMODEL_H_ */
