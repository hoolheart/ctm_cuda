/*
 * CellTransModel.cpp
 *
 *  Created on: 26 mars 2013
 *      Author: Edward
 */

#include "CellTransModel.h"
#include <iostream>
#include <limits>
#include <cmath>
using namespace std;

CellTransModel::CellTransModel() {
    // initialize the parameters
	resetSystem(10,10,5,1);
}

CellTransModel::~CellTransModel() {
	for (int i=0;i<(int)list_cells.size();i++)
		delete list_cells[i];
	list_cells.clear();
	for (int i=0;i<(int)list_links.size();i++)
		delete list_links[i];
	list_links.clear();
	list_pos_in.clear();
	list_pos_out.clear();
	list_in.clear();
	list_out.clear();

	for (int i=0;i<(int)list_lanes.size();i++)
		delete list_lanes[i];
	list_lanes.clear();

	for (int i=0;i<(int)list_ints.size();i++)
		delete list_ints[i];
	list_ints.clear();
}

bool CellTransModel::sim(double dt, int steps) {
    if (!info.is_sim_on) {
//        err = "The Simulation has not been started!";
        return false;
    }
    if (dt<=0) {
//        err = "Non-positive interval!";
        return false;
    }
    if (steps<=0) {
//        err = "Non-positive steps!";
        return false;
    }
    for (int i=0;i<steps;i++) {
        calPosFlows(dt);
        calRealFlows();
        if(!updateCells(dt))
            return false;
    }
    return true;
}

void CellTransModel::calPosFlows(double dt) {
	int n_cells = list_cells.size();
    for (int i=0;i<n_cells;i++) {
        list_in[i] = 0;
        list_out[i] = 0;
        switch (list_cells[i]->type) {
        case CELL_TYPE_INPUT:
            list_pos_in[i] = 0;
            list_cells[i]->length += list_cells[i]->rate*dt;
            list_pos_out[i] = list_cells[i]->length;
            break;
        case CELL_TYPE_OUTPUT:
            list_pos_in[i] = std::numeric_limits<double>::max();
            list_pos_out[i] = 0;
            break;
        case CELL_TYPE_NORMAL:
        default:
            double maxFlow = list_cells[i]->rate*dt;
            list_pos_in[i] = min(maxFlow,info.w_vf*(list_cells[i]->cap-list_cells[i]->length));
            list_pos_out[i] = min(maxFlow,list_cells[i]->length);
            break;
        }
    }
}

double CellTransModel::min(double d1, double d2) {
    if (d1<d2)
        return d1;
    else
        return d2;
}

void CellTransModel::calRealFlows() {
	int n_links = list_links.size();
    for (int i=0;i<n_links;i++) {
        if(list_links[i]->access) {
        	int c1 = list_links[i]->cells[0], c2 = list_links[i]->cells[1], c3 = list_links[i]->cells[2];
        	double b1 = list_links[i]->ratio; double b2 = 1-b1;
            double f1,f2;
            switch (list_links[i]->type) {
            case LINK_TYPE_MERGE:
                if (list_pos_in[c3]>=list_pos_out[c1]+list_pos_out[c2]) {
                    f1 = list_pos_out[c1];
                    f2 = list_pos_out[c2];
                }
                else {
                    f1 = mid(list_pos_out[c1],list_pos_in[c3]-list_pos_out[c2],b1*list_pos_in[c3]);
                    f2 = mid(list_pos_out[c2],list_pos_in[c3]-list_pos_out[c1],b2*list_pos_in[c3]);
                }
				list_out[c1] = f1;
				list_out[c2] = f2;
				list_in[c3] = f1+f2;
                break;
            case LINK_TYPE_DIVERGE:
            	f1 = min(list_pos_out[c1],min(list_pos_in[c2]/b1,list_pos_in[c3]/b2));
            	list_out[c1] = f1;
            	list_in[c2] = f1*b1;
            	list_in[c3] = f1*b2;
                break;
            case LINK_TYPE_DIRECT:
            default:
            	f1 = min(list_pos_out[c1],list_pos_in[c2]);
            	list_out[c1] = f1;
            	list_in[c2] = f1;
                break;
            }
        }
    }
}

double CellTransModel::mid(double d1, double d2, double d3) {
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

bool CellTransModel::updateCells(double dt) {
	int n_cells = list_cells.size();
    for (int i=0;i<n_cells;i++) {
        if (list_in[i]>list_pos_in[i]+1e-6 || list_out[i]>list_pos_out[i]+1e-6) {
//            ostringstream oss;
//            oss << "Simulation failed with unknown reason in Cell " << i << "!";
//            err = oss.str();
            return false;
        }
        list_cells[i]->length += list_in[i]-list_out[i];
        list_cells[i]->delay += dt*(list_cells[i]->length-list_out[i]*list_cells[i]->cap*info.veh_len/info.vf);
    }
    return true;
}

void CellTransModel::resetSystem(double vf,
		double w,
		double veh_len,
		double pos_dt) {
	for (int i=0;i<(int)list_cells.size();i++)
		delete list_cells[i];
	list_cells.clear();
	for (int i=0;i<(int)list_links.size();i++)
		delete list_links[i];
	list_links.clear();
	list_pos_in.clear();
	list_pos_out.clear();
	list_in.clear();
	list_out.clear();

	info.is_valid = false;
	info.is_sim_on = false;
	info.vf = vf;
	info.w_vf = w/vf;
	info.veh_len = veh_len;
	info.cell_cap = vf*pos_dt/veh_len;

	for (int i=0;i<(int)list_lanes.size();i++)
		delete list_lanes[i];
	list_lanes.clear();

	for (int i=0;i<(int)list_ints.size();i++)
		delete list_ints[i];
	list_ints.clear();
}

int CellTransModel::addLane(const string &id,
		int type,
		double cap,
		double sat_rate,
		double in_rate,
		double out_ratio) {
	if (info.is_sim_on) {
		return -1;
	}
	if (getLaneIndexById(id)+1) {
		return -1;
	}

	CtmLane * l = new CtmLane();
	l->id = id;
	l->type = type;
	bool isValid = false;
	switch (type) {
	case LANE_TYPE_NORMAL:
		if (cap<=0)
			break;
		else if (sat_rate <= 0)
			break;
		else if (in_rate < 0)
			break;
		else if (out_ratio<0 || out_ratio>1)
			break;
		isValid = true;
		l->cap = cap;
		l->sat_rate = sat_rate;
		l->in_rate = in_rate;
		l->out_ratio = out_ratio;
		break;
	case LANE_TYPE_ENTRY:
		if (cap<=0)
			break;
		else if (sat_rate <= 0)
			break;
		else if (in_rate < 0)
			break;
		isValid = true;
		l->cap = cap;
		l->sat_rate = sat_rate;
		l->in_rate = in_rate;
		l->out_ratio = 0;
		break;
	case LANE_TYPE_EXIT:
		isValid = true;
		l->cap = 0;
		l->sat_rate = std::numeric_limits<double>::max();
		l->in_rate = 0;
		l->out_ratio = 1;
		break;
	default:
		break;
	}

	if (isValid) {
		list_lanes.push_back(l);
		info.is_valid = false;
		return list_lanes.size()-1;
	}
	else {
		delete l;
		return -1;
	}
}

int CellTransModel::addIntersection(string id,const vector<string> &in_lanes,
		const vector<string> &out_lanes,
        int n_inner,double inner_cells[][2]) {
	if (info.is_sim_on) {
			return -1;
	}
	if (getIntersectionIndexById(id)+1) {
		return -1;
	}

	CtmIntersection * l = new CtmIntersection();
	l->id = id;
	for (int i=0;i<(int)in_lanes.size();i++) {
		int tmp = getLaneIndexById(in_lanes[i]);
		if (tmp>=0)
			l->in_lanes.push_back(list_lanes[tmp]);
	}
	for (int i=0;i<(int)out_lanes.size();i++) {
		int tmp = getLaneIndexById(out_lanes[i]);
		if (tmp>=0)
			l->out_lanes.push_back(list_lanes[tmp]);
	}
	for (int i=0;i<n_inner;i++) {
		CtmInnerCell * c = new CtmInnerCell(inner_cells[i][0],inner_cells[i][1]);
		l->inner_cells.push_back(c);
	}

	list_ints.push_back(l);
	info.is_valid = false;
	return list_ints.size()-1;
}

bool CellTransModel::addPhase(string id,int n_links,double info[][8]) {
	int index = getIntersectionIndexById(id);
	if (index<0 || index>=(int)list_ints.size()) {
		return false;
	}
	CtmPhase *p = new CtmPhase();
	p->info.resize(n_links);
	for (int i=0;i<n_links;i++) {
		p->info[i].type = (int)info[i][0];
		p->info[i].ratio = info[i][1];
		p->info[i].c1s = (int)info[i][2];
		p->info[i].c1i = (int)info[i][3];
		p->info[i].c2s = (int)info[i][4];
		p->info[i].c2i = (int)info[i][5];
		p->info[i].c3s = (int)info[i][6];
		p->info[i].c3i = (int)info[i][7];
	}
	list_ints[index]->phases.push_back(p);
	return true;
}

bool CellTransModel::buildCTM() {
	if (info.is_sim_on) {
		return false;
	}
	else if (info.is_valid) {
		return true;
	}
	else if (list_lanes.size()==0 || list_ints.size()==0) {
		return false;
	}

	// clear the current CTM
	for (int i=0;i<(int)list_cells.size();i++)
		delete list_cells[i];
	list_cells.clear();
	for (int i=0;i<(int)list_links.size();i++)
		delete list_links[i];
	list_links.clear();
	list_pos_in.clear();
	list_pos_out.clear();
	list_in.clear();
	list_out.clear();

	// build new CTM
	for (int i=0;i<(int)list_lanes.size();i++) {
		modelingLane(list_lanes[i]);
	}
	for (int i=0;i<(int)list_ints.size();i++) {
		modelingIntersection(list_ints[i]);
	}
	list_pos_in = vector<double>(list_cells.size(),0);
	list_pos_out = vector<double>(list_cells.size(),0);
	list_in = vector<double>(list_cells.size(),0);
	list_out = vector<double>(list_cells.size(),0);
	info.is_valid = true;
	return true;
}

void CellTransModel::modelingLane(CtmLane *l) {
	int n_cells = (int)list_cells.size();
	int n_links = (int)list_links.size();
	int n;
	CtmCell *c;
	CtmLink *k;
	switch (l->type) {
	case LANE_TYPE_NORMAL:
		n = (l->cap<3*info.cell_cap)?3:(int)round(l->cap/info.cell_cap);
		for (int i=0;i<n;i++) {
			c = new CtmCell;
			c->type = CELL_TYPE_NORMAL;
			c->cap = l->cap/n;
			c->rate = l->sat_rate;
			c->delay = 0;
			c->length = 0;
			list_cells.push_back(c);
		}
		c = new CtmCell;
		c->type = CELL_TYPE_INPUT;
		c->cap = numeric_limits<double>::max(); c->rate = l->in_rate;
		c->delay = 0; c->length = 0;
		list_cells.push_back(c);
		c = new CtmCell;
		c->type = CELL_TYPE_OUTPUT;
		c->cap = numeric_limits<double>::max(); c->rate = numeric_limits<double>::max();
		c->delay = 0; c->length = 0;
		list_cells.push_back(c);
		k = new CtmLink;
		k->type = LINK_TYPE_DIVERGE; k->cells[0] = n_cells; k->cells[1] = n_cells+1; k->cells[2] = n_cells+n+1;
		k->ratio = 1-l->out_ratio; k->access = true;
		list_links.push_back(k);
		k = new CtmLink;
		k->type = LINK_TYPE_MERGE; k->cells[0] = n_cells+1; k->cells[1] = n_cells+n; k->cells[2] = n_cells+2;
		k->ratio = l->sat_rate/(l->sat_rate+l->in_rate); k->access = true;
		list_links.push_back(k);
		for (int i=2;i<n-1;i++) {
			k = new CtmLink;
			k->type = LINK_TYPE_DIRECT; k->cells[0] = n_cells+i; k->cells[1] = n_cells+i+1; k->cells[2] = -1;
			k->ratio = 1; k->access = true;
			list_links.push_back(k);
		}
		l->o_cell = n_cells; l->d_cell = n_cells+n-1;
		l->in_cell = n_cells+n; l->out_cell = n_cells+n+1;
		l->in_link = n_links+1; l->out_link = n_links;
		break;
	case LANE_TYPE_ENTRY:
		n = (int)round(l->cap/info.cell_cap);
		c = new CtmCell;
		c->type = CELL_TYPE_INPUT;
		c->cap = numeric_limits<double>::max(); c->rate = l->in_rate;
		c->delay = 0; c->length = 0;
		list_cells.push_back(c);
		for (int i=0;i<n;i++) {
			c = new CtmCell;
			c->type = CELL_TYPE_NORMAL;
			c->cap = l->cap/n;
			c->rate = l->sat_rate;
			c->delay = 0;
			c->length = 0;
			list_cells.push_back(c);
		}
		for (int i=0;i<n;i++) {
			k = new CtmLink;
			k->type = LINK_TYPE_DIRECT; k->cells[0] = n_cells+i; k->cells[1] = n_cells+i+1; k->cells[2] = -1;
			k->ratio = 1; k->access = true;
			list_links.push_back(k);
		}
		l->o_cell = n_cells+1; l->d_cell = n_cells+n;
		l->in_cell = n_cells; l->out_cell = -1;
		l->in_link = n_links; l->out_link = -1;
		break;
	case LANE_TYPE_EXIT:
		c = new CtmCell;
		c->type = CELL_TYPE_OUTPUT;
		c->cap = numeric_limits<double>::max(); c->rate = numeric_limits<double>::max();
		c->delay = 0; c->length = 0;
		list_cells.push_back(c);
		l->o_cell = -1; l->d_cell = -1;
		l->in_cell = -1; l->out_cell = n_cells;
		l->in_link = -1; l->out_link = -1;
		break;
	}
}

void CellTransModel::modelingIntersection(CtmIntersection *l) {
	// clear
	l->in_cells.clear();
	l->out_cells.clear();
	l->cur_phase = -1;

	// modeling
	for (int i=0;i<(int)l->in_lanes.size();i++) {
		l->in_cells.push_back(l->in_lanes[i]->getTailCell());
	}
	for (int i=0;i<(int)l->out_lanes.size();i++) {
		l->out_cells.push_back(l->out_lanes[i]->getHeadCell());
	}
	CtmCell *c;
	for (int i=0;i<(int)l->inner_cells.size();i++) {
		c = new CtmCell;
		c->type = CELL_TYPE_NORMAL;
		c->cap = l->inner_cells[i]->cap;
		c->rate = l->inner_cells[i]->sat_rate;
		c->delay = 0;
		c->length = 0;
		list_cells.push_back(c);
		l->inner_cells[i]->index = (int)list_cells.size()-1;
	}
	CtmLink *k;
	for (int i=0;i<(int)l->phases.size();i++) {
		CtmPhase *p = l->phases[i];
		if(p->info.empty()) {
			p->head_link = -1;
			p->tail_link = -1;
			continue;
		}
		for (int j=0;j<(int)p->info.size();j++) {
			k = new CtmLink;
			k->access = false;
			k->type = p->info[j].type;
			switch (p->info[j].c1s) {
			case 0:
				k->cells[0] = l->inner_cells[p->info[j].c1i]->index;
				break;
			case 1:
				k->cells[0] = l->in_cells[p->info[j].c1i];
				break;
			case 2:
			default:
				k->cells[0] = l->out_cells[p->info[j].c1i];
				break;
			}
			switch (p->info[j].c2s) {
			case 0:
				k->cells[1] = l->inner_cells[p->info[j].c2i]->index;
				break;
			case 1:
				k->cells[1] = l->in_cells[p->info[j].c2i];
				break;
			case 2:
			default:
				k->cells[1] = l->out_cells[p->info[j].c2i];
				break;
			}
			switch (p->info[j].type) {
			case LINK_TYPE_MERGE:
			case LINK_TYPE_DIVERGE:
				switch (p->info[j].c3s) {
				case 0:
					k->cells[2] = l->inner_cells[p->info[j].c3i]->index;
					break;
				case 1:
					k->cells[2] = l->in_cells[p->info[j].c3i];
					break;
				case 2:
				default:
					k->cells[2] = l->out_cells[p->info[j].c3i];
					break;
				}
				k->ratio = p->info[j].ratio;
				break;
			case LINK_TYPE_DIRECT:
			default:
				k->type = LINK_TYPE_DIRECT;
				k->cells[2] = -1;
				k->ratio = 1;
				break;
			}
			list_links.push_back(k);
		}
		p->tail_link = (int)list_links.size()-1;
		p->head_link = (int)(list_links.size()-p->info.size());
	}
	return;
}

bool CellTransModel::checkCells() {
	if(info.is_valid) {
		for(int i=0;i<(int)list_cells.size();i++) {
			if(list_cells[i]->length<0 || list_cells[i]->length>list_cells[i]->cap)
				return false;
		}
		return true;
	}
	else
		return false;
}

bool CellTransModel::checkPhases() {
	if(info.is_valid) {
		for(int i=0;i<(int)list_ints.size();i++) {
			if(list_ints[i]->cur_phase<0 ||
					list_ints[i]->cur_phase>=(int)list_ints[i]->phases.size())
				return false;
		}
		return true;
	}
	else
		return false;
}

bool CellTransModel::startSim() {
	if(!info.is_valid)
		return false;
	if(info.is_sim_on)
		return true;
	if(checkCells()&&checkPhases()) {
		info.is_sim_on = true;
		return true;
	}
	else
		return false;
}

bool CellTransModel::setLaneQueue(int i,double x) {
	if(!info.is_valid)
		return false;
	if(info.is_sim_on)
		return false;
	if(i<0 || i>=(int)list_lanes.size())
		return false;
	if(list_lanes[i]->type==LANE_TYPE_EXIT)
		return true;

	CtmLane * l = list_lanes[i];
	x = (x>0)?x:0; x = (x<l->cap)?x:l->cap;
	//int n = l->d_cell-l->o_cell+1;
	CtmCell *c;
	for (int j=l->d_cell;j>=l->o_cell;j--) {
		c = list_cells[j];
		double clen = (x>c->cap)?c->cap:x;
		c->length = clen;
		x -= clen;
	}
	return true;
}

bool CellTransModel::setLaneQueue(string id,double x) {
	int i = getLaneIndexById(id);
	if (i<0)
		return false;
	else
		return setLaneQueue(i,x);
}

bool CellTransModel::setIntersectionPhase(int i,int p) {
	if(!info.is_valid)
		return false;
//	if(info.is_sim_on)
//		return false;
	if(i<0 || i>=(int)list_ints.size())
		return false;

	CtmIntersection *l = list_ints[i];
	p = (p>0)?p:0;
	p = (p<(int)l->phases.size())?p:((int)l->phases.size()-1);
	bool acc; CtmPhase *f;
	for(int j=0;j<(int)l->phases.size();j++) {
		acc = (j==p);
		f = l->phases[j];
		for(int k=f->head_link;k<=f->tail_link;k++) {
			list_links[k]->access = acc;
		}
	}
	l->cur_phase = p;
	return true;
}

bool CellTransModel::setIntersectionPhase(string id, int p) {
	int i = getIntersectionIndexById(id);
	if (i<0)
		return false;
	else
		return setIntersectionPhase(i,p);
}

bool CellTransModel::startSim(const vector<double> &x,const vector<int> &p) {
	if(!info.is_valid)
		return false;
	if(info.is_sim_on)
		return false;
	if(x.size()!=list_lanes.size() || p.size()!=list_ints.size())
		return false;

	for(int i=0;i<(int)list_lanes.size();i++) {
		setLaneQueue(i,x[i]);
	}
	for(int i=0;i<(int)list_ints.size();i++) {
		setIntersectionPhase(i,p[i]);
	}
	info.is_sim_on = true;
	return true;
}

bool CellTransModel::stopSim() {
	if(!info.is_valid)
		return false;
	info.is_sim_on = false;
	return true;
}

bool CellTransModel::cleanAllCells() {
	if(!info.is_valid)
		return false;

	for(int i=0;i<(int)list_cells.size();i++) {
		list_cells[i]->length = 0;
		list_cells[i]->delay = 0;
	}
	return true;
}

bool CellTransModel::modifyLaneInRate(string id,double r) {
	if(!info.is_valid)
		return false;
	int i = getLaneIndexById(id);
	if(i<0 || i>=(int)list_lanes.size())
		return false;
	if(list_lanes[i]->type==LANE_TYPE_EXIT)
		return false;

	CtmLane * l = list_lanes[i];
	r = (r>0)?r:0;
	l->in_rate = r;
	list_cells[l->in_cell]->rate = r;
	return true;
}

bool CellTransModel::modifyLaneSatRate(string id,double r) {
	if(!info.is_valid)
		return false;
	int i = getLaneIndexById(id);
	if(i<0 || i>=(int)list_lanes.size())
		return false;
	if(list_lanes[i]->type==LANE_TYPE_EXIT)
		return false;

	CtmLane * l = list_lanes[i];
	r = (r>0)?r:0;
	l->sat_rate = r;
	for(int j=l->o_cell;j<=l->d_cell;j++)
		list_cells[j]->rate = r;
	return true;
}

bool CellTransModel::modifyLaneOutRatio(string id,double r) {
	if(!info.is_valid)
		return false;
	int i = getLaneIndexById(id);
	if(i<0 || i>=(int)list_lanes.size())
		return false;
	if(list_lanes[i]->type!=LANE_TYPE_NORMAL)
		return false;

	CtmLane * l = list_lanes[i];
	r = (r>0)?r:0; r = (r<1)?r:1;
	l->out_ratio = r;
	list_links[l->out_link]->ratio = 1-r;
	return true;
}

bool CellTransModel::switchIntersection(string id) {
	if(!info.is_valid)
		return false;
	int i = getIntersectionIndexById(id);
	if(i<0 || i>=(int)list_ints.size())
		return false;

	CtmIntersection *l = list_ints[i];
	if(l->phases.size()>0) {
		int p = l->cur_phase;
		CtmPhase *f;
		f = l->phases[p];
		for(int k=f->head_link;k<=f->tail_link;k++) {
			list_links[k]->access = false;
		}
		p++;
		p = (p<(int)l->phases.size())?p:0;
		f = l->phases[p];
		for(int k=f->head_link;k<=f->tail_link;k++) {
			list_links[k]->access = true;
		}
		l->cur_phase = p;
	}
	return true;
}

bool CellTransModel::readCells(vector<double> &tar) {
	if(!info.is_valid)
		return false;
	tar.resize(list_cells.size());
	for(int i=0;i<(int)list_cells.size();i++)
		tar[i] = list_cells[i]->length;
	return true;
}

bool CellTransModel::readLanes(vector<double> &tar) {
	if(!info.is_valid)
		return false;
	tar.resize(list_lanes.size());
	for(int i=0;i<(int)list_lanes.size();i++) {
		tar[i] = 0;
		if(list_lanes[i]->type==LANE_TYPE_EXIT)
			continue;
		for(int j=list_lanes[i]->o_cell;j<=list_lanes[i]->d_cell;j++)
			tar[i] += list_cells[j]->length;
	}
	return true;
}

bool CellTransModel::readPhases(vector<int> &tar) {
	if(!info.is_valid)
		return false;
	tar.resize(list_ints.size());
	for(int i=0;i<(int)list_ints.size();i++)
		tar[i] = list_ints[i]->cur_phase;
	return true;
}

bool CellTransModel::readLaneDelays(vector<double> &tar) {
	if(!info.is_valid)
		return false;
	tar.resize(list_lanes.size());
	for(int i=0;i<(int)list_lanes.size();i++) {
		tar[i] = 0;
		if(list_lanes[i]->type==LANE_TYPE_EXIT)
			continue;
		for(int j=list_lanes[i]->o_cell;j<=list_lanes[i]->d_cell;j++)
			tar[i] += list_cells[j]->delay;
	}
	return true;
}

double CellTransModel::readTotalDelay() {
	vector<double> delays;
	if(readLaneDelays(delays)) {
		double total=0;
		for(int i=0;i<(int)delays.size();i++)
			total += delays[i];
		return total;
	}
	else
		return 0;
}

void CellTransModel::resetDelay() {
	for (int i=0;i<(int)list_cells.size();i++)
		list_cells[i]->delay = 0;
}
