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
#include "ctm_cuda.h"
using namespace std;

class CellTransModel {
private:
    vector<CtmCell *> list_cells;
    vector<CtmLink *> list_links;
    vector<double> list_pos_in,list_pos_out,list_in,list_out;
    CtmInfo info;
    vector<CtmLane *> list_lanes;
    vector<CtmIntersection *> list_ints;
public:
    CellTransModel();
    virtual ~CellTransModel();
    bool sim(double dt, int steps=1);
    void resetSystem(double vf, double w, double veh_len, double pos_dt);
    int addLane(const string &id, int type, double cap,
            double sat_rate, double in_rate, double out_ratio);
    CtmLane * getLaneById(string id) {
    	for (int i=0;i<(int)list_lanes.size();i++)
    		if (list_lanes[i]->id==id) return list_lanes[i];
    	return (CtmLane *)0;
    }
    int getLaneIndexById(string id) {
    	for (int i=0;i<(int)list_lanes.size();i++)
    		if (list_lanes[i]->id==id) return i;
    	return -1;
    }
    int addIntersection(string id,const vector<string> &in_lanes,
    		const vector<string> &out_lanes,
            int n_inner,double inner_cells[][2]);
    CtmIntersection * getIntersectionById(string id) {
    	for (int i=0;i<(int)list_ints.size();i++)
    		if (list_ints[i]->id==id) return list_ints[i];
    	return (CtmIntersection *)0;
    }
    int getIntersectionIndexById(string id) {
    	for (int i=0;i<(int)list_ints.size();i++)
    		if (list_ints[i]->id==id) return i;
    	return -1;
    }
    bool addPhase(string id,int n_links,double info[][8]);
    bool buildCTM();
    bool checkCells();
    bool checkPhases();
    bool startSim();
    bool setLaneQueue(int i,double x);
    bool setLaneQueue(string id,double x);
    bool setIntersectionPhase(int i,int p);
    bool setIntersectionPhase(string id,int p);
    bool startSim(const vector<double> &x,const vector<int> &p);
    bool stopSim();
    bool cleanAllCells();
    bool addInputs(const vector<double> &in);
    bool modifyLaneInRate(string id,double r);
    bool modifyLaneSatRate(string id,double r);
    bool modifyLaneOutRatio(string id,double r);
    bool switchIntersection(string id);
    bool readCells(vector<double> &tar);
    bool readLanes(vector<double> &tar);
    bool readPhases(vector<int> &tar);
    bool readLaneDelays(vector<double> &tar);
    double readTotalDelay();
    void resetDelay();
private:
    void calPosFlows(double dt);
    void calRealFlows();
    bool updateCells(double dt);
    double min(double d1, double d2);
    double mid(double d1, double d2, double d3);
    void modelingLane(CtmLane *l);
    void modelingIntersection(CtmIntersection *l);
};

#endif /* CELLTRANSMODEL_H_ */
