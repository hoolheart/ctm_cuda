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
void startSim(CellTransModel &mdl, double lens[]);
void printLens(vector<double> lens);
void simOneCycle(CellTransModel &mdl, double c, double g1, double g2, vector<double> &lens);

int main() {
	cout << "Hello World!" << endl; // prints Hello World!

	// testing the CTM in CPU
	CellTransModel model;
	initialModel(model);
	double lens[]= {31,10,40,45,10,42,48,51};
	startSim(model,lens);
	vector<double> lengths;
	model.readLanes(lengths);
	printLens(lengths);
	double c = 90;
	simOneCycle(model,c,45,45,lengths);
	printLens(lengths);
	simOneCycle(model,c,45,45,lengths);
	printLens(lengths);
	simOneCycle(model,c,55,45,lengths);
	printLens(lengths);

	string str;
	cin >> str;
	return 0;
}

void initialModel(CellTransModel &mdl) {
	// lanes
	// entry lanes for intersection 1
	cout<<mdl.addLane("L1",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<mdl.addLane("L2",LANE_TYPE_NORMAL,60,0.6,0.03,0.1);
	cout<<mdl.addLane("L3",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<mdl.addLane("L4",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<endl;
	// entry lanes for intersection 2
	cout<<mdl.addLane("L5",LANE_TYPE_NORMAL,60,0.6,0.03,0.1);
	cout<<mdl.addLane("L6",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<mdl.addLane("L7",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<mdl.addLane("L8",LANE_TYPE_ENTRY,60,0.6,0.3,0);
	cout<<endl;
	// exit lanes
	mdl.addLane("LE1",LANE_TYPE_EXIT,0,0,0,1);
	mdl.addLane("LE2",LANE_TYPE_EXIT,0,0,0,1);
	mdl.addLane("LE3",LANE_TYPE_EXIT,0,0,0,1);
	mdl.addLane("LE4",LANE_TYPE_EXIT,0,0,0,1);
	mdl.addLane("LE5",LANE_TYPE_EXIT,0,0,0,1);
	mdl.addLane("LE6",LANE_TYPE_EXIT,0,0,0,1);
	// intersections
	vector<string> in_lanes, out_lanes;
	double inner_cells[2][2] = {{2,0.6},{2,0.6}};
	double phase1[4][8] = {{LINK_TYPE_DIRECT, 1, 1,0, 0,0,0,0},
			{LINK_TYPE_DIRECT, 1, 0,0, 2,0,0,0},
			{LINK_TYPE_DIRECT, 1, 1,1, 0,1,0,0},
			{LINK_TYPE_DIRECT, 1, 0,1, 2,1,0,0}};
	double phase2[4][8] = {{LINK_TYPE_DIRECT, 1, 1,2, 0,0,0,0},
			{LINK_TYPE_DIRECT, 1, 0,0, 2,2,0,0},
			{LINK_TYPE_DIRECT, 1, 1,3, 0,1,0,0},
			{LINK_TYPE_DIRECT, 1, 0,1, 2,3,0,0}};
	// intersection 1
	in_lanes.push_back("L1");in_lanes.push_back("L2");
	in_lanes.push_back("L3");in_lanes.push_back("L4");
	out_lanes.push_back("L5");out_lanes.push_back("LE1");
	out_lanes.push_back("LE2");out_lanes.push_back("LE3");
	cout<<mdl.addIntersection("I1",in_lanes,out_lanes,2,inner_cells);
	cout<<mdl.addPhase("I1",4,phase1);cout<<mdl.addPhase("I1",4,phase2);
	cout<<endl;
	in_lanes.clear(); out_lanes.clear();
	// intersection 2
	in_lanes.push_back("L5");in_lanes.push_back("L6");
	in_lanes.push_back("L7");in_lanes.push_back("L8");
	out_lanes.push_back("LE4");out_lanes.push_back("L2");
	out_lanes.push_back("LE5");out_lanes.push_back("LE6");
	cout<<mdl.addIntersection("I2",in_lanes,out_lanes,2,inner_cells);
	cout<<mdl.addPhase("I2",4,phase1);cout<<mdl.addPhase("I2",4,phase2);
	cout<<endl;
	in_lanes.clear(); out_lanes.clear();
	// build CTM
	if(mdl.buildCTM())
		cout << "Building successfully." <<endl;
	else
		cout << "Failed to build CTM." << endl;
}

void startSim(CellTransModel &mdl, double lens[]) {
	vector<double> lanes;
	for(int i=0;i<8;i++)
		lanes.push_back(lens[i]);
	for(int i=0;i<6;i++)
		lanes.push_back(0);
	vector<int> ints;
	ints.push_back(0);ints.push_back(0);
	mdl.startSim(lanes,ints);
}

void printLens(vector<double> lens) {
	cout << "queue lengths = ";
	for (int i=0;i<7;i++)
			cout << lens[i] << "\t";
	cout << lens[7] << endl;
}

void simOneCycle(CellTransModel &mdl,
		double c, double g1, double g2,
		vector<double> &lens) {
	string s1,s2;
	double t1,t2,t3;
	if (g1<g2) {
		s1 = "I1";
		s2 = "I2";
		t1 = g1;
		t2 = g2-g1;
		t3 = c-g2;
	}
	else {
		s1 = "I2";
		s2 = "I1";
		t1 = g2;
		t2 = g1-g2;
		t3 = c-g1;
	}

	// step 1
	if (t1>=1) {
		int sp = (int)floor(t1);
		if (!mdl.sim(1.0,sp))
			cout << "Simulation Error!" << endl;
		t1 -= sp;
	}
	if (t1>1e-6)
		mdl.sim(t1);
	mdl.switchIntersection(s1);

	// step 2
	if (t2>=1) {
		int sp = (int)floor(t2);
		if (!mdl.sim(1.0,sp))
			cout << "Simulation Error!" << endl;
		t2 -= sp;
	}
	if (t2>0)
		mdl.sim(t2);
	mdl.switchIntersection(s2);

	// step 3
	if (t3>=1) {
		int sp = (int)floor(t3);
		if (!mdl.sim(1.0,sp))
			cout << "Simulation Error!" << endl;
		t3 -= sp;
	}
	if (t3>0)
		mdl.sim(t3);
	mdl.switchIntersection("I1");mdl.switchIntersection("I2");

	mdl.readLanes(lens);
}
