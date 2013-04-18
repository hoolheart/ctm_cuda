/*
 * CellTransModel.cpp
 *
 *  Created on: 26 mars 2013
 *      Author: Edward
 */

#include "CellTransModel.h"

CellTransModel::CellTransModel() {
	// initialize the parameters
	ListCell = NULL;
	ListLink = NULL;
	CellPosIn = NULL;
	CellPosOut = NULL;
	CellIn = NULL;
	CellOut = NULL;
	err = "";
	n = 0;
	isSimOn = false;
}

CellTransModel::CellTransModel(int _n) {
	// initialize the parameters
	ListCell = NULL;
	ListLink = NULL;
	CellPosIn = NULL;
	CellPosOut = NULL;
	CellIn = NULL;
	CellOut = NULL;
	err = "";
	n = 0;
	isSimOn = false;
	// initialize model with n cells and m links
	initialModel(_n);
}

CellTransModel::~CellTransModel() {
	if(ListCell!=NULL)
		delete ListCell;
	if(ListLink!=NULL)
		delete ListLink;
	if(CellPosIn!=NULL)
		delete CellPosIn;
	if(CellPosOut!=NULL)
		delete CellPosOut;
	if(CellIn!=NULL)
		delete CellIn;
	if(CellOut!=NULL)
		delete CellOut;
	if(isSimOn) {
		deleteCudaEnv();
	}
}

void CellTransModel::resetModel() {
	if(ListCell!=NULL) {
		delete ListCell;
		ListCell = NULL;
	}
	if(ListLink!=NULL) {
		delete ListLink;
		ListLink = NULL;
	}
	if(CellPosIn!=NULL) {
		delete CellPosIn;
		CellPosIn = NULL;
	}
	if(CellPosOut!=NULL) {
		delete CellPosOut;
		CellPosOut = NULL;
	}
	if(CellIn!=NULL) {
		delete CellIn;
		CellIn = NULL;
	}
	if(CellOut!=NULL) {
		delete CellOut;
		CellOut = NULL;
	}
	err = "";
	n = 0;
	if(isSimOn) {
		deleteCudaEnv();
	}
	isSimOn = false;
}

bool CellTransModel::initialModel(int _n, float _cap, float _s) {
	if (_n<=0) {
		err = "Non-positive number of cells!";
		return false;
	}
	resetModel();
	n = _n;
	ListCell = new CELL[n];
	ListLink = new LINK[n];
	CellPosIn = new float[n];
	CellPosOut = new float[n];
	CellIn = new float[n];
	CellOut = new float[n];
	for (int i=0;i<n;i++) {
		ListCell[i].type = CELL_TYPE_NORMAL;
		ListCell[i].capacity = _cap;
		ListCell[i].length = 0;
		ListCell[i].rate = _s;
		ListCell[i].access = 1;
		CellPosIn[i] = 0;
		CellPosOut[i] = 0;
		CellIn[i] = 0;
		CellOut[i] = 0;
		ListLink[i].p = 1;
		ListLink[i].c1 = (i-1)>=0?(i-1):0;
		ListLink[i].c2 = 0;
		ListLink[i].p1 = 1;
		ListLink[i].p2 = 0;
	}
	return true;
}

bool CellTransModel::setCell(int i, int _type, float _cap, float _s) {
	if(n==0) {
		err = "No cell in the model yet!";
		return false;
	}
	if (i>=n || i<0) {
		err = "Wrong index of cell!";
		return false;
	}
	if (_s<=0) {
		err = "Negative maximal flow rate!";
		return false;
	}
	switch (_type) {
	case CELL_TYPE_INPUT:
		ListCell[i].type = _type;
		ListCell[i].rate = _s;
		break;
	case CELL_TYPE_OUTPUT:
		ListCell[i].type = _type;
		ListCell[i].rate = _s;
		break;
	case CELL_TYPE_SWITCH:
		if (_cap<=0) {
			err = "Negative capacity of cell!";
			return false;
		}
		ListCell[i].type = _type;
		ListCell[i].capacity = _cap;
		ListCell[i].rate = _s;
		break;
	case CELL_TYPE_NORMAL:
	default:
		if (_cap<=0) {
			err = "Negative capacity of cell!";
			return false;
		}
		ListCell[i].type = CELL_TYPE_NORMAL;
		ListCell[i].capacity = _cap;
		ListCell[i].rate = _s;
		break;
	}
	return true;
}

bool CellTransModel::setLink(int i, float _p, int _c1, float _p1, int _c2, float _p2) {
	if(n==0) {
		err = "No link in the model yet!";
		return false;
	}
	if (i>=n || i<0) {
		err = "Wrong index of link!";
		return false;
	}
	if (_p<=0 || _p>1 || _p1<=0 || _p1>1 || _p2<0 || _p2>1) {
		err = "Wrong proportion!";
		return false;
	}
	if (_c1<0 || _c1>=n || _c2<0 || _c2>=n) {
		err = "Wrong index of input cell(s)!";
		return false;
	}
	ListLink[i].p = _p;
	ListLink[i].c1 = _c1;
	ListLink[i].p1 = _p1;
	ListLink[i].c2 = _c2;
	ListLink[i].p2 = _p2;
	return true;
}

bool CellTransModel::startSim(float const _len[], int const _acc[]) {
	if (n==0) {
		err = "No cell and/or link in the model yet!";
		return false;
	}
	if (isSimOn) {
		err = "The Simulation was started!";
		return false;
	}
	for (int i=0;i<n;i++) {
		switch (ListCell[i].type) {
		case CELL_TYPE_INPUT:
		case CELL_TYPE_OUTPUT:
			ListCell[i].length = 0;
			break;
		case CELL_TYPE_SWITCH:
		case CELL_TYPE_NORMAL:
		default:
			if (_len[i]<0) {
				ostringstream oss;
				oss << "Negative length of Cell " << i << "!";
				err = oss.str();
				return false;
			}
			else if (_len[i]>ListCell[i].capacity) {
				ostringstream oss;
				oss << "Over the capacity of Cell " << i << "!";
				err = oss.str();
				return false;
			}
			ListCell[i].length = _len[i];
			break;
		}
		if (ListCell[i].type==CELL_TYPE_SWITCH) {
			if(_acc[i]!=0)
				ListCell[i].access = 1;
			else
				ListCell[i].access = 0;
		}
	}
	createCudaEnv(ListCell,ListLink,CellPosIn,CellPosOut,CellIn,CellOut,n);
	isSimOn = true;
	return true;
}

bool CellTransModel::changeAccess(int i, int _acc) {
	if (i>=n || i<0) {
		err = "Wrong index of cell!";
		return false;
	}
	else if (ListCell[i].type != CELL_TYPE_SWITCH) {
		err = "This cell is not controlled by signals!";
		return false;
	}
	if (_acc!=0)
		ListCell[i].access = 1;
	else
		ListCell[i].access = 0;
	if(isSimOn)
		updateCudaAcc(ListCell);
	return true;
}

bool CellTransModel::changeAccesses(int const _acc[]) {
	if(n==0) {
		err = "No cell in the model yet!";
		return false;
	}
	for (int i=0;i<n;i++) {
		if (ListCell[i].type == CELL_TYPE_SWITCH) {
			if (_acc[i]!=0)
				ListCell[i].access = 1;
			else
				ListCell[i].access = 0;
		}
	}
	if(isSimOn)
		updateCudaAcc(ListCell);
	return true;
}

bool CellTransModel::getCurrentLengths(float _len[]) {
	if(n==0) {
		err = "No cell in the model yet!";
		return false;
	}
	if (isSimOn)
		loadCudaLen(ListCell);
	for (int i=0;i<n;i++) {
		_len[i] = ListCell[i].length;
	}
	return true;
}

bool CellTransModel::sim(float dt, int steps) {
	if (!isSimOn) {
		err = "The Simulation has not been started!";
		return false;
	}
	if (dt<=0) {
		err = "Non-positive interval!";
		return false;
	}
	if (steps<=0) {
		err = "Non-positive steps!";
		return false;
	}
	for (int i=0;i<steps;i++) {
//		calPosFlows(dt);
//		calRealFlows();
//		if(!updateCells())
//			return false;
		if(!simCuda(dt))
			return false;
	}
	return true;
}

void CellTransModel::stopSim() {
	loadCudaLen(ListCell);
	deleteCudaEnv();
	isSimOn = false;
}

bool CellTransModel::resumeSim() {
	if (n==0) {
		err = "No cell and/or link in the model yet!";
		return false;
	}
	if (isSimOn) {
		err = "The Simulation was started!";
		return false;
	}
	createCudaEnv(ListCell,ListLink,CellPosIn,CellPosOut,CellIn,CellOut,n);
	isSimOn = true;
	return true;
}

void CellTransModel::calPosFlows(float dt) {
	for (int i=0;i<n;i++) {
		switch (ListCell[i].type) {
		case CELL_TYPE_INPUT:
			CellPosIn[i] = 0;
			CellPosOut[i] = ListCell[i].rate*dt;
			break;
		case CELL_TYPE_OUTPUT:
			CellPosIn[i] = ListCell[i].rate*dt;
			CellPosOut[i] = 0;
			break;
		case CELL_TYPE_NORMAL:
		case CELL_TYPE_SWITCH:
		default:
			float maxFlow = ListCell[i].rate*dt;
			CellPosIn[i] = min(maxFlow,ListCell[i].capacity-ListCell[i].length);
			CellPosOut[i] = min(maxFlow,ListCell[i].length)*ListCell[i].access;
			break;
		}
	}
}

float CellTransModel::min(float d1, float d2) {
	if (d1<d2)
		return d1;
	else
		return d2;
}

void CellTransModel::calRealFlows() {
	for (int i=0;i<n;i++) {
		CellOut[i] = 0;
	}
	for (int i=0;i<n;i++) {
		if (ListCell[i].type==CELL_TYPE_INPUT)
			continue;
		if (ListLink[i].p==1) {
			float f = min(CellPosIn[i],CellPosOut[ListLink[i].c1]*ListLink[i].p1);
			CellIn[i] = f;
			CellOut[ListLink[i].c1] += f;
		}
		else {
			float r = CellPosIn[i];
			float s1 = CellPosOut[ListLink[i].c1]*ListLink[i].p1;
			float s2 = CellPosOut[ListLink[i].c2]*ListLink[i].p2;
			float f1, f2;
			if (r>=(s1+s2)) {
				f1 = s1; f2 = s2;
			}
			else {
				f1 = mid(s1, r-s2, r*ListLink[i].p);
				f2 = mid(s2, r-s1, r*(1-ListLink[i].p));
			}
			CellIn[i] = f1+f2;
			CellOut[ListLink[i].c1] += f1;
			CellOut[ListLink[i].c2] += f2;
		}
	}
}

float CellTransModel::mid(float d1, float d2, float d3) {
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

bool CellTransModel::updateCells() {
	for (int i=0;i<n;i++) {
		if (CellIn[i]>CellPosIn[i]+1e-6 || CellOut[i]>CellPosOut[i]+1e-6) {
			ostringstream oss;
			oss << "Simulation failed with unknown reason in Cell " << i << "!";
			err = oss.str();
			return false;
		}
		ListCell[i].length += CellIn[i]-CellOut[i];
	}
	return true;
}
