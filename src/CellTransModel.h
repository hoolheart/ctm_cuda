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
    vector<CtmCell> list_cells;
    vector<CtmLink> list_links;
//    vector<float> list_pos_in,list_pos_out,list_in,list_out;
    CtmInfo info;
    vector<CtmLane *> list_lanes;
    vector<CtmIntersection *> list_ints;
public:
    CellTransModel();
    virtual ~CellTransModel();
    bool sim(float dt, int steps=1);
    void resetSystem(float vf, float w, float veh_len, float pos_dt);
    int addLane(const string &id, int type, float cap,
            float sat_rate, float in_rate, float out_ratio);
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
            int n_inner,float inner_cells[][2]);
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
    bool addPhase(string id,int n_links,float info[][8]);
    bool buildCTM();
    bool checkCells();
    bool checkPhases();
    bool startSim();
    bool setLaneQueue(int i,float x);
    bool setLaneQueue(string id,float x);
    bool setIntersectionPhase(int i,int p);
    bool setIntersectionPhase(string id,int p);
    bool startSim(const vector<float> &x,const vector<int> &p);
    bool stopSim();
    bool cleanAllCells();
    bool addInputs(const vector<float> &in);
    bool modifyLaneInRate(string id,float r);
    bool modifyLaneSatRate(string id,float r);
    bool modifyLaneOutRatio(string id,float r);
    bool switchIntersection(string id);
    bool readCells(vector<float> &tar);
    bool readLanes(vector<float> &tar);
    bool readPhases(vector<int> &tar);
    bool readLaneDelays(vector<float> &tar);
    float readTotalDelay();
    void resetDelay();
private:
//    void calPosFlows(float dt);
//    void calRealFlows();
//    bool updateCells(float dt);
//    float min(float d1, float d2);
//    float mid(float d1, float d2, float d3);
    void modelingLane(CtmLane *l);
    void modelingIntersection(CtmIntersection *l);
};

#endif /* CELLTRANSMODEL_H_ */
