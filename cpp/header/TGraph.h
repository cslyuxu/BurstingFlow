#pragma once
/*
 * TGraph.h for temporal network
 *
 *  Created on: October 8, 2023
 *      Author: Lyu
 */

#ifndef TGRAPH_H_
#define TGRAPH_H_




#include "GlobalDefinition.h"
#include <vector>

//#include <tr1/unordered_map>
/*For VS*/
#include <unordered_map>
#include <iostream>
#include <queue>
#include <fstream>
#include <string> 
#include<stdlib.h>
#include<time.h>
#include<math.h>
using namespace std;


template<class WeightType, class TimestampType>
class TGRAPH {
public:
	/**
	 * Adj definition.
	 */
	
	//template<class WeightType>
	class AdjElement {
	public:
		//    VertexID v;
		EdgeID eid;
		WeightType eweight;

		/**
		 * used and initialized in EL
		 */
		bool isVisited;

		AdjElement() {
			//      v =
			//      eid = elabel = -1;
		}

		// VertexID v : can be removed?
		AdjElement(EdgeID eid, const WeightType& _eweight)
			://v(v),
			eid(eid),
			eweight(_eweight) {
		}

		~AdjElement() {

		}
	};

	/*
	 * data structures for basic operations
	 */
	typedef unordered_map<VertexID, int> VLabels;
	//typedef AdjElement ADJELE;
	typedef unordered_map<TimestampType, AdjElement> TOUTEDGE;
	typedef unordered_map<VertexID, TOUTEDGE> AdjList;
	typedef unordered_map<VertexID, AdjList> OutEdge;			//for directed graph, use adjacent link
	typedef unordered_map<TimestampType, bool> AdjListBool;
	typedef unordered_map<VertexID, AdjListBool> TINEDGE;
	typedef unordered_map<VertexID, TINEDGE> InVertex;		//record the inVertex	
	typedef unordered_map<VertexID, int> MapMatrix;				//record the map for the matrix: make the vertexID within [0. integer]
	typedef unordered_map<VertexID, set<TimestampType>> TempVertexTimestamp;   //to index the VertexTimestamp
	typedef unordered_map<VertexID, unordered_map<int, TimestampType>> VertexNumTimestamp; //: <-1, # of timestamp> 
	typedef unordered_map<VertexID, unordered_map<TimestampType, int>> VertexTimestampNum; 
 
	GraphID graphId;	//for multiple graph
	int Vcnt, Ecnt;		//vertex count && edge count
	VLabels _vlabels;  // vertex label
	OutEdge _outEdges;  // out edge list
	InVertex _inVertex;  // in vertex
	int diameter;
	VertexID e;
	MapMatrix _map;		//for the array matrix
	VertexID* matrix;
	set<TimestampType> TimeStampSet;						//record the timestamps of vertices
	unordered_map<TimestampType, TimestampType> func;
	VertexNumTimestamp _OutVNT, _InVNT, _VNT;
	VertexTimestampNum _OutVTN, _InVTN, _VTN;
	TempVertexTimestamp _OutVT, _InVT, _VT;
	TimestampType Tmax, Tmin;
	unordered_map<int, VertexID> fastVisitVertex;
	int tempVcount;

	/*
	 * basic operations
	 */
	VLabels& getVLabel();
	OutEdge& getOutEdge();
	InVertex& getInVertex();
	int getVcnt();
	int getEcnt();
	TGRAPH();
	~TGRAPH();
	void reset();	
	bool isVertex(VertexID v);
	bool isEdge(VertexID s, VertexID d, TimestampType t);
	bool isTimeStamp(const TimestampType& timestamp);
	int& getVLabel(VertexID v);
	void insertVertex(VertexID v, const int& label);	
	void insertEdge(VertexID s, VertexID d, const TimestampType& el, const WeightType& wl);
	void initMap();
	void IndexTimestamp(TimestampType, TimestampType);
	int maxdegree();
	int avgdegree();
	int stddegree();
	VertexID getRV();
	MapMatrix& getMap();



	/*
	 * Do *NOT* use the below data structures
	 */
	MapLabelCnt _outLabelCnt;
	MapLabelCnt _inLabelCnt;
	AdjListBool _vVisited;

};







/**
 * implementations
 */

#ifndef DEFAULT_VERTEX_NUMBER
#define DEFAULT_VERTEX_NUMBER   (256)
#endif


template<class WeightType, class TimestampType>
typename TGRAPH<WeightType, TimestampType>::VLabels& TGRAPH<WeightType,
	TimestampType>::getVLabel() {
	return _vlabels;
}

template<class WeightType, class TimestampType>
typename TGRAPH<WeightType, TimestampType>::OutEdge& TGRAPH<WeightType,
	TimestampType>::getOutEdge() {
	return _outEdges;
}

template<class WeightType, class TimestampType>
typename TGRAPH<WeightType, TimestampType>::InVertex& TGRAPH<WeightType,
	TimestampType>::getInVertex() {
	return _inVertex;
}

template<class WeightType, class TimestampType>
int TGRAPH<WeightType, TimestampType>::getVcnt() {
	return Vcnt;
}

template<class WeightType, class TimestampType>
int TGRAPH<WeightType, TimestampType>::getEcnt() {
	return Ecnt;
}



template<class WeightType, class TimestampType>
TGRAPH<WeightType, TimestampType>::TGRAPH()
	: Vcnt(0),
	Ecnt(0),
	graphId(INVALID_GRAPH_ID) {
		this->matrix = NULL;
}

template<class WeightType, class TimestampType>
TGRAPH<WeightType, TimestampType>::~TGRAPH() {
	getVLabel().clear();
	getOutEdge().clear();
	getInVertex().clear();
	if(this->matrix != NULL)
		delete[] this->matrix;
	TimeStampSet.clear();
	func.clear();
}



template<class WeightType, class TimestampType>
void TGRAPH<WeightType, TimestampType>::reset() {
	getVLabel().clear();
	getOutEdge().clear();
	getInVertex().clear();
	Vcnt = Ecnt = 0;
}






/******************************************/
template<class WeightType, class TimestampType>
bool TGRAPH<WeightType, TimestampType>::isVertex(VertexID v) {
	return (getVLabel().find(v) != getVLabel().end());
}


template<class WeightType, class TimestampType>
bool TGRAPH<WeightType, TimestampType>::isEdge(VertexID s, VertexID d, TimestampType t) {
	ASSERT(isVertex(s));
	ASSERT(isVertex(d));
	//ASSERT(isTimeStamp(t));
	if(getOutEdge()[s].find(d) != getOutEdge()[s].end())
		return (getOutEdge()[s][d].find(t) != getOutEdge()[s][d].end());
	return false;
}


template<class WeightType, class TimestampType>
bool TGRAPH<WeightType, TimestampType>::isTimeStamp(const TimestampType& timestamp) {
	return (TimeStampSet.find(timestamp) != TimeStampSet.end());
}



template<class WeightType, class TimestampType>
int& TGRAPH<WeightType, TimestampType>::getVLabel(VertexID v) {
	ASSERT(isVertex(v));
	return getVLabel()[v];
}



template<class WeightType, class TimestampType>
void TGRAPH<WeightType, TimestampType>::insertVertex(VertexID v,
	const int& label) {
	if (!isVertex(v)) {
		Vcnt++;
		fastVisitVertex[tempVcount] = v;
		tempVcount++;
	}
	getVLabel()[v] = label;
}



template<class WeightType, class TimestampType>
void TGRAPH<WeightType, TimestampType>::insertEdge(VertexID s, VertexID d,
	const TimestampType& el, const WeightType& wl) {
	//  getOutEdge()[s][d] = AdjElement<TimestampType>(s, Ecnt, el);
	if(!isEdge(s, d, el)){
		Ecnt++;
		getOutEdge()[s][d][el] = AdjElement(Ecnt, wl);			//The number of edge 'Ecnt' is used as an ID for the edge.
		getInVertex()[d][s][el] = true;
		if (TimeStampSet.find(el) == TimeStampSet.end())
			TimeStampSet.insert(el);
		if(_VT[s].find(el) == _VT[s].end())
			_VT[s].insert(el);
		if(_VT[d].find(el) == _VT[d].end())
			_VT[d].insert(el);
		if(_OutVT[s].find(el) == _OutVT[s].end())
			_OutVT[s].insert(el);
		if(_InVT[d].find(el) == _InVT[d].end())
			_InVT[d].insert(el);		
	}
	//else{    //for checking the repreated edges
	//	cout<<"Edge <"<<s<<", "<<d << ", "<<el<<", "<< wl <<"> is repeated!"<<endl;
	//}
}


//initial the mapping from VertexID to the matrix number*/
template<class WeightType, class TimestampType>
void TGRAPH<WeightType, TimestampType>::initMap() {
	int size = getVcnt();
	this->matrix = new VertexID[size];
	int i = 0;
	for (typename VLabels::iterator it = getVLabel().begin();
		it != getVLabel().end(); it++) {
		VertexID s = it->first;
		getMap()[s] = i;
		this->matrix[i] = s;
		//cout << "The number of vetex " << s << " is " << i << endl;	
		i++;
		if (i == size)
			break;
	}
}


//Index the timestamp for each vertex
template<class WeightType, class TimestampType>
void TGRAPH<WeightType, TimestampType>::IndexTimestamp(TimestampType max, TimestampType min) {
	VertexID tempV;
	int count;
	TimestampType tempcount;

	//make the timestamp standard
	tempcount = 1;	
	this->Tmin = tempcount;	
	for(auto it = TimeStampSet.begin();it != TimeStampSet.end();it++){
		func[*it] = tempcount;
		tempcount++;
	}
	this->Tmax = tempcount-1;

	for(auto it = getVLabel().begin();it != getVLabel().end(); it++){
		tempV = it->first;		
		//For out-edge timestamp
		count = 0;
		for(auto it1= _OutVT[tempV].begin(); it1 != _OutVT[tempV].end(); it1++){
			_OutVNT[tempV][count] = func[*it1];
			_OutVTN[tempV][func[*it1]] = count;
			//if(_OutVNT[tempV][count]<_OutVNT[tempV][count-1])
			//	cout<<"wrong!"<<endl;
			count++;
		}
		_OutVNT[tempV][-1] = count;  


		//For in-edge timestamp
		count = 0;
		for(auto it2= _InVT[tempV].begin(); it2 != _InVT[tempV].end(); it2++){
			_InVNT[tempV][count] = func[*it2];
			_InVTN[tempV][func[*it2]] = count;
			count++;
		}
		_InVNT[tempV][-1] = count;  


		//For edge timestamp
		count = 0;
		for(auto it3= _VT[tempV].begin(); it3 != _VT[tempV].end(); it3++){
			_VNT[tempV][count] = func[*it3];
			_VTN[tempV][func[*it3]] = count;
			count++;
			//if(func[*it3]<5)
				//cout <<"V: "<<tempV<<"		T:"<< func[*it3]<<endl;
		}
		_VNT[tempV][-1] = count;  
		
		//cout <<"V: "<<tempV<<"		Outdegree:"<< _OutVNT[tempV][-1]<<"		Indegree:"<< _InVNT[tempV][-1]<<"		Degree:"<< _VNT[tempV][-1]<<endl;
	}



	_OutVT.clear();
	_InVT.clear();
	_VT.clear();
}


template<class WeightType, class TimestampType>
int TGRAPH<WeightType, TimestampType>::maxdegree() {
	return 0;
}

template<class WeightType, class TimestampType>
int TGRAPH<WeightType, TimestampType>::avgdegree() {
	double avg = getEcnt()*2/getVcnt();
	return avg;
}

template<class WeightType, class TimestampType>
int TGRAPH<WeightType, TimestampType>::stddegree() {
	double std, avg, temp, total;
	avg = avgdegree();
	std = 0;
	total = 0;
	VertexID u, v;
	for(auto it = getVLabel().begin();it != getVLabel().end();it++){
		u = it->first;
		temp = 0;
		for(auto it1 = getOutEdge()[u].begin();it1 != getOutEdge()[u].end();it1++){
			v = it1->first;
			for(auto it2 = getOutEdge()[u][v].begin();it2 != getOutEdge()[u][v].end();it2++)
				temp++;
		}

		for(auto it1 = getInVertex()[u].begin();it1 != getInVertex()[u].end();it1++){
			v = it1->first;
			for(auto it2 = getInVertex()[u][v].begin();it2 != getInVertex()[u][v].end();it2++)
				temp++;
		}
		temp = temp - avg;
		total = total + (temp * temp);
	}
	total = total / getVcnt();
	return sqrt(total);
}

template<class WeightType, class TimestampType>
typename TGRAPH<WeightType, TimestampType>::MapMatrix& TGRAPH<WeightType, TimestampType>::getMap() {
	return _map;
}

template<class WeightType, class TimestampType>
VertexID TGRAPH<WeightType, TimestampType>::getRV() {
	VertexID num = rand()%getVcnt();
	//VertexID count = 0;
	//for(auto it = getVLabel().begin();it!=getVLabel().end();it++){
	//	if(count==num)
	//		return it->first;
	//	count++;
	//}
	//return -1;
	return fastVisitVertex[num];
}

#endif /* TGRAPH_H_ */
