//============================================================================
// Name        : ctm_c.cpp
// Author      : Edward
// Version     :
// Copyright   : GPL
// Description : Implementation of CTM by C++
//============================================================================

#include <iostream>
#include <math.h>
using namespace std;
#include "CellTransModel.h"

void initialModel(CellTransModel &mdl);
void startSim(CellTransModel &mdl, float lens[]);
void simOneCycle(CellTransModel &mdl, float c, float g1, float g2, float lens[]);
void setAccesses(CellTransModel &mdl, int acc[]);
void printLens(float lens[]);
void updateLens(CellTransModel &mdl, float lens[]);

int main() {
	cout << "Hello World!" << endl; // prints Hello World!

	// testing the CTM in CPU
	CellTransModel model;
	initialModel(model);
	float lens[]= {31,10,40,45,10,42,48,51};
	startSim(model,lens);
	printLens(lens);
	float c = 90;
	simOneCycle(model,c,45,45,lens);
	printLens(lens);
	simOneCycle(model,c,45,45,lens);
	printLens(lens);
	simOneCycle(model,c,55,45,lens);
	printLens(lens);

//	string str;
//	cin >> str;
	return 0;
}

void initialModel(CellTransModel &mdl) {
	float cap = 20, s = 0.6;
	if (!mdl.initialModel(40,cap,s))
		cout << mdl.getErr() << endl;
	for (int i=0;i<8;i++) {
		if (!mdl.setCell(i*4,CELL_TYPE_INPUT,cap,s/2)) cout << mdl.getErr() << endl;
		if (!mdl.setCell(i*4+3,CELL_TYPE_SWITCH,cap,s)) cout << mdl.getErr() << endl;
		if (!mdl.setCell(32+i,CELL_TYPE_OUTPUT,cap,s*2)) cout << mdl.getErr() << endl;
		if (!mdl.setLink(32+i,1,i*4+3,1,0,0)) cout << mdl.getErr() << endl;
	}
	if (!mdl.setLink(5,0.9,23,0.9,4,1)) cout << mdl.getErr() << endl;
	if (!mdl.setLink(17,0.9,3,0.9,16,1)) cout << mdl.getErr() << endl;
	if (!mdl.setCell(4,CELL_TYPE_INPUT,cap,s/20)) cout << mdl.getErr() << endl;
	if (!mdl.setCell(16,CELL_TYPE_INPUT,cap,s/20)) cout << mdl.getErr() << endl;
	if (!mdl.setLink(32,1,3,0.1,0,0)) cout << mdl.getErr() << endl;
	if (!mdl.setLink(37,1,23,0.1,0,0)) cout << mdl.getErr() << endl;
}

void startSim(CellTransModel &mdl, float lens[]) {
	float len_cell[40] = {0};
	float cap = 20;
	for (int j=0;j<8;j++) {
		float l = lens[j];
		for (int i=3;i>0;i--) {
			float t = l<cap?l:cap;
			len_cell[4*j+i] = t;
			l -= t;
		}
		len_cell[4*j] = l;
	}
	int acc_cell[40];
	for (int i=0;i<40;i++)
		acc_cell[i] = 1;
	acc_cell[11] = 0; acc_cell[15] = 0; acc_cell[27] = 0; acc_cell[31] = 0;
	if (!mdl.startSim(len_cell,acc_cell))
		cout << mdl.getErr() << endl;
}

void simOneCycle(CellTransModel &mdl, float c, float g1, float g2, float lens[]) {
	int s1[] = {1,1,0,0,1,1,0,0};
	int *s2;
	int s2p[] = {0,0,1,1,1,1,0,0};
	int s2n[] = {1,1,0,0,0,0,1,1};
	int s3[] = {0,0,1,1,0,0,1,1};
	float t1,t2,t3;
	if (g1<g2) {
		s2 = s2p;
		t1 = g1;
		t2 = g2-g1;
		t3 = c-g2;
	}
	else {
		s2 = s2n;
		t1 = g2;
		t2 = g1-g2;
		t3 = c-g1;
	}

	setAccesses(mdl,s1);
	if (t1>=1) {
		int sp = (int)floor(t1);
		if (!mdl.sim(1.0,sp)) cout << mdl.getErr() << endl;
		t1 -= sp;
	}
	if (t1>0)
		mdl.sim(t1);

	setAccesses(mdl,s2);
	if (t2>=1) {
		int sp = (int)floor(t2);
		if (!mdl.sim(1.0,sp)) cout << mdl.getErr() << endl;
		t2 -= sp;
	}
	if (t2>0)
		mdl.sim(t2);

	setAccesses(mdl,s3);
	if (t3>=1) {
		int sp = (int)floor(t3);
		mdl.sim(1.0,sp);
		t3 -= sp;
	}
	if (t3>0)
		mdl.sim(t3);

	updateLens(mdl,lens);
}

void setAccesses(CellTransModel &mdl, int acc[]) {
	for (int i=0;i<8;i++)
		if (!mdl.changeAccess(i*4+3,acc[i]))
			cout << mdl.getErr();
}

void printLens(float lens[]) {
	cout << "queue lengths = ";
	for (int i=0;i<7;i++)
		cout << lens[i] << "\t";
	cout << lens[7] << endl;
}

void updateLens(CellTransModel &mdl, float lens[]) {
	float len_cell[40];
	mdl.getCurrentLengths(len_cell);
	for (int i=0;i<8;i++)
		lens[i] = len_cell[i*4+1]+len_cell[i*4+2]+len_cell[i*4+3];
}
