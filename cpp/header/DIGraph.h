#pragma once
/*
 * DIGraph.h for transformed network
 *
 *  Created on: October 8, 2023
 *      Author: Lyu
 */

#ifndef DIGRAPH_H_
#define DIGRAPH_H_


#include <thread>
#include "TGraph.h"
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
using namespace std;


template<class WeightType, class TimestampType>
class DIGRAPH {
public:
	/**
	 * Adj definition.
	 */
	
	//template<class WeightType>
	class AdjElement {
	public:
		//    VertexID v;
		EdgeID eid;
		WeightType eweight; // WEIGHT_MAX = +infy
		WeightType flow;
		WeightType residual; // have not been used

		/**
		 * used and initialized in EL
		 */
		bool isVisited;

		AdjElement() {
			//      v =
			//      eid = elabel = -1;
		}

		// VertexID v : can be removed?
		AdjElement(EdgeID eid, const WeightType& _eweight, const WeightType& _flow, const WeightType& _residual)
			://v(v),
			eid(eid),
			eweight(_eweight), flow(_flow), residual(_residual) {
		}

		~AdjElement() {

		}
	};

	class VERTEX {
	public:
		VertexID vid;
		TimestampType timestamp;

		VERTEX(){
			//
		}

		VERTEX(VertexID _vid, const TimestampType& _timestamp)
			:
			vid(_vid),
			timestamp(_timestamp) {

		}

		void setvalue(VertexID _vid, const TimestampType& _timestamp){
			vid = _vid;
			timestamp = _timestamp;
		}

		~VERTEX() {

		}
	};

	/*
	 * data structures for basic operations
	 */
	typedef unordered_map<VertexID, int> VLabels;
	typedef unordered_map<VertexID, long> TLabels;
	//typedef AdjElement ADJELE;
	typedef unordered_map<VertexID, AdjElement> AdjList;
	typedef unordered_map<VertexID, AdjList> OutEdge;			//for directed graph, use adjacent link
	typedef unordered_map<VertexID, bool> AdjListBool;
	typedef unordered_map<VertexID, AdjListBool> InVertex;		//record the inVertex	
	typedef unordered_map<VertexID, int> MapMatrix;				//record the map for the matrix: make the vertexID within [0. integer]
	typedef unordered_map<VertexID, unordered_map<VertexID, TimestampType>> VVT;
	typedef unordered_map<VertexID, unordered_map<TimestampType, VertexID>> VTV;
	typedef unordered_map<VertexID, unordered_map<int, VertexID>> VNV; //for Dinic: record the number of out-going edges
	typedef unordered_map<VertexID, VertexID> NextV;
	typedef unordered_map<VertexID, TimestampType> NextT;

	VertexID s, t;
	VERTEX source, sink;
	GraphID graphId;	//for multiple graph
	int Vcnt, Ecnt;		//vertex count && edge count
	VLabels _vlabels;  // vertex label
	TLabels _tlabels;  // <temporal vertex, earliest timestamp>
	OutEdge _outEdges;  // out edge list
	InVertex _inVertex;  // in vertex
	int diameter;
	VertexID e;
	MapMatrix _map;		//for the array matrix
	VertexID* matrix;
	//unordered_set<TimestampType> TimeStampSet;						//record the timestamps of vertices
	//VVT _VtoTV;
	VTV _TVtoV;
	VNV OutEdgeNum;		// used for Dinic
	NextV nextVertex; //used for TransMAPSE to record the next vertex for each temporal vertex
	NextT nextTime; //used for TransMAPSE to record the timestamp of nextVertex
	queue<VertexID> nodes;
	WeightType maxflow;
	/*
	 * basic operations
	 */
	VLabels& getVLabel();
	OutEdge& getOutEdge();
	InVertex& getInVertex();
	int getVcnt();
	int getEcnt();
	DIGRAPH();
	~DIGRAPH();
	void reset();	
	void copy(DIGRAPH<WeightType, TimestampType> *);
	bool isVertex(VertexID v);
	bool isEdge(VertexID s, VertexID d);
	//bool isTimeStamp(const TimestampType& timestamp);
	int& getVLabel(VertexID v);
	bool insertVertex(VertexID tv, const TimestampType& timestamp);	
	void insertEdge(VertexID s, VertexID d, const WeightType& wl);
	void initMap();
	void eraseVertex(VertexID s, VertexID tv, const TimestampType& timestamp);
	void removeAllInEdges(VertexID s);
	void removeEdge(VertexID s, VertexID d, bool verify = true);
	//bool isEdgeHasFlow(VertexID s, VertexID d);
	MapMatrix& getMap();

	void clear(queue<int>& q) {
		queue<int> empty;
		swap(empty, q);
	};

	/*
	void setVLabel(VertexID v, WeightType& label);
	void setELabel(VertexID s, VertexID d, TimestampType& el);
	TimestampType& getELabel(VertexID s, VertexID d);
	void removeEdge(VertexID s, VertexID d, bool verify = true);
	void removeAllOutEdges(VertexID s);
	void removeAllInEdges(VertexID s);
	void removeVertex(VertexID s);
	void eraseVertex(VertexID s);
	int getDegree(VertexID v);
	int getOutDegree(VertexID v);
	int getInDegree(VertexID v);
	void prinDIGRAPH(ostream& out);
	void getDiameter();
	void getDNeighbor(VertexID, int, unordered_set<VertexID>&, DIGRAPH*);
	void createBall(int);
	void countBall(int);
	void getDNeighbor_without_label_check_of_query(VertexID, int, unordered_set<VertexID>&
	);
	void ConstructInducedGraph(unordered_set<VertexID>&, DIGRAPH*);
	void initVL();
	void initVL_Random();
	void initVertexLabel(unordered_map<VertexID, WeightType>&);	
	void initEL(VertexID x);
	void initEdgeVisited();
	void constLabelCnt();
	void setVertexVisited();
	*/

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
typename DIGRAPH<WeightType, TimestampType>::VLabels& DIGRAPH<WeightType,
	TimestampType>::getVLabel() {
	return _vlabels;
}

template<class WeightType, class TimestampType>
typename DIGRAPH<WeightType, TimestampType>::OutEdge& DIGRAPH<WeightType,
	TimestampType>::getOutEdge() {
	return _outEdges;
}

template<class WeightType, class TimestampType>
typename DIGRAPH<WeightType, TimestampType>::InVertex& DIGRAPH<WeightType,
	TimestampType>::getInVertex() {
	return _inVertex;
}

template<class WeightType, class TimestampType>
int DIGRAPH<WeightType, TimestampType>::getVcnt() {
	return Vcnt;
}

template<class WeightType, class TimestampType>
int DIGRAPH<WeightType, TimestampType>::getEcnt() {
	return Ecnt;
}



template<class WeightType, class TimestampType>
DIGRAPH<WeightType, TimestampType>::DIGRAPH()
	: Vcnt(0),
	Ecnt(0),
	maxflow(0),
	graphId(INVALID_GRAPH_ID) {
		this->matrix = NULL;
}

template<class WeightType, class TimestampType>
DIGRAPH<WeightType, TimestampType>::~DIGRAPH() {
	getVLabel().clear();
	getOutEdge().clear();
	getInVertex().clear();
	if(this->matrix != NULL)
		delete[] this->matrix;
	//TimeStampSet.clear();
	//_VtoTV.clear();
	_TVtoV.clear();
	OutEdgeNum.clear();
	_tlabels.clear();
	nextVertex.clear();
	nextTime.clear();
	//nodes.clear();
	clear(nodes);
	
}



template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::reset() {
	getVLabel().clear();
	getOutEdge().clear();
	getInVertex().clear();
	Vcnt = Ecnt = 0;
	maxflow = 0;
	//TimeStampSet.clear();						
	//_VtoTV.clear();
	_TVtoV.clear();
	OutEdgeNum.clear();
	_tlabels.clear();
	//nodes.clear();
	clear(nodes);
	nextVertex.clear();
	nextTime.clear();
}





template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::copy(DIGRAPH<WeightType, TimestampType> *B) {
	this->Vcnt = B->Vcnt;
	this->Ecnt = B->Ecnt;	
	/*
	DIGRAPH<WeightType, TimestampType> *A = this;
			thread t1([A, B]{A->getVLabel() = B->getVLabel();});



			thread t2([A, B]{A->getOutEdge() = B->getOutEdge();});



			thread t3([A, B]{A->getInVertex() = B->getInVertex();});



			thread t4([A, B]{A->_TVtoV = B->_TVtoV;});


			thread t5([A, B]{A->OutEdgeNum = B->OutEdgeNum;});


			thread t6([A, B]{A->nodes = B->nodes;});



			t1.join();
			t2.join();
			t3.join();
			t4.join();
			t5.join();
			t6.join();
	*/
	
	this->getVLabel() = B->getVLabel();
	this->getOutEdge() = B->getOutEdge();
	this->getInVertex() = B->getInVertex();
	//B->TimeStampSet = this->TimeStampSet;		
	//B->_VtoTV = this->_VtoTV;
	this->_TVtoV = B->_TVtoV;
	this->OutEdgeNum = B->OutEdgeNum;
	this->nodes = B->nodes;
	this->maxflow = B->maxflow;
	this->nextVertex = B->nextVertex;
	this->nextTime = B->nextTime;
	this->_tlabels = B->_tlabels;
}







/******************************************/
template<class WeightType, class TimestampType>
bool DIGRAPH<WeightType, TimestampType>::isVertex(VertexID v) {
	return (getVLabel().find(v) != getVLabel().end());
}


template<class WeightType, class TimestampType>
bool DIGRAPH<WeightType, TimestampType>::isEdge(VertexID s, VertexID d) {
	ASSERT(isVertex(s));
	ASSERT(isVertex(d));
	//ASSERT(isTimeStamp(t));
	return (getOutEdge()[s].find(d) != getOutEdge()[s].end());
}

/*
template<class WeightType, class TimestampType>
bool DIGRAPH<WeightType, TimestampType>::isTimeStamp(const TimestampType& timestamp) {
	return (TimeStampSet.find(timestamp) != TimeStampSet.end());
}*/



template<class WeightType, class TimestampType>
int& DIGRAPH<WeightType, TimestampType>::getVLabel(VertexID v) {
	ASSERT(isVertex(v));
	return getVLabel()[v];
}



template<class WeightType, class TimestampType>
bool DIGRAPH<WeightType, TimestampType>::insertVertex(VertexID tv, const TimestampType& timestamp) {

	if(_TVtoV[tv].find(timestamp)==_TVtoV[tv].end())
		Vcnt++;
	else{
		//cout<<"Vertex <"<<tv<<","<<timestamp<<"> exists!"<<endl;
		//cout<<"hello"<<endl;
		return false;
	}

	if(_tlabels.find(tv) == _tlabels.end()){
		_tlabels[tv] = timestamp;
	}


	VertexID ps;
	if(!nodes.empty()){
		ps = nodes.front();
		nodes.pop();
		getVLabel()[ps] = 0;
		//_tlabels[tv] = 0;
		//_VtoTV[v][tv] = timestamp;
		_TVtoV[tv][timestamp] = ps;
		OutEdgeNum[ps][-1] = 0;
	}else{
		getVLabel()[Vcnt-1] = 0;
		//_tlabels[tv] = 0;
		//_VtoTV[v][tv] = timestamp;
		_TVtoV[tv][timestamp] = Vcnt-1;
		OutEdgeNum[Vcnt-1][-1] = 0;
	}

	/*
	if (!isVertex(v)) {
		Vcnt++;
	}
	getVLabel()[v] = 0;
	_tlabels[tv] = 0;
	//_VtoTV[v][tv] = timestamp;
	_TVtoV[tv][timestamp] = v;
	OutEdgeNum[v][-1] = 0;
	*/
	return true;
}



template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::insertEdge(VertexID s, VertexID d, const WeightType& wl) {

	int neighbor_position;
	
	//  getOutEdge()[s][d] = AdjElement<TimestampType>(s, Ecnt, el);
	if(!isEdge(s, d)){
		Ecnt++;
		getOutEdge()[s][d] = AdjElement(Ecnt, wl, 0, 0);			//The number of edge 'Ecnt' is used as an ID for the edge.
		getInVertex()[d][s] = true;

		//assume no bi-directed edges
		getOutEdge()[d][s] = AdjElement(Ecnt, 0, 0, 0);

		neighbor_position = OutEdgeNum[s][-1];
  		OutEdgeNum[s][neighbor_position] = d;  
  		OutEdgeNum[s][-1]++;
  		neighbor_position = OutEdgeNum[d][-1];
  		OutEdgeNum[d][neighbor_position] = s;  
  		OutEdgeNum[d][-1]++;
	}	

	

}


/*initial the mapping from VertexID to the matrix number*/
template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::initMap() {
	int size = getVcnt();
	this->matrix = new VertexID[size];
	int i = 0;
	for (typename VLabels::iterator it = getVLabel().begin();
		it != getVLabel().end(); it++) {
		//VertexID s = it->first;
		//getMap()[s] = i;
		//this->matrix[i] = s;
		//cout << "The number of vetex " << s << " is " << i << endl;	
		i++;
		if (i == size)
			break;
	}
}

template<class WeightType, class TimestampType>
typename DIGRAPH<WeightType, TimestampType>::MapMatrix& DIGRAPH<WeightType, TimestampType>::getMap() {
	return _map;
}

/*
template<class WeightType, class TimestampType>
bool DIGRAPH<WeightType, TimestampType>::isEdgeHasFlow(VertexID s, VertexID d) {
	ASSERT(isVertex(s));
	ASSERT(isVertex(d));
	ASSERT(isEdge(s, d));
	if((getOutEdge()[s][d].eweight > 0)||(getOutEdge()[s][d].flow < getOutEdge()[s][d].eweight)) 
		return true;
	return false;
}*/

template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::eraseVertex(VertexID s, VertexID tv, const TimestampType& timestamp) {
	ASSERT(isVertex(s));

	/*
	getVLabel().erase(s);
	getOutEdge().erase(s);
	getInVertex().erase(s);
	Vcnt--;
	*/

	//ReduceSizeOfVertex
	Vcnt--;
	VertexID ps;
	int pos1;
	for(auto it = getOutEdge()[s].begin(); it!= getOutEdge()[s].end();it++){
		ps = it->first;
  		pos1 = OutEdgeNum[ps][-1];
		for(int i = 0; i<pos1;i++){
			if(OutEdgeNum[ps][i] == s){
				OutEdgeNum[ps][i] = OutEdgeNum[ps][pos1-1];
				OutEdgeNum[ps][-1]--;
				OutEdgeNum[ps].erase(pos1-1);
				getOutEdge()[ps].erase(s);
				getInVertex()[ps].erase(s);
				Ecnt--;
				break;
			}
		}

	}

	/*
	for(auto it = getInVertex()[s].begin(); it!= getInVertex()[s].end();it++){
		ps = it->first;
  		pos1 = OutEdgeNum[ps][-1];
		for(int i = 0; i<pos1;i++){
			if(OutEdgeNum[ps][i] == s){
				OutEdgeNum[ps][i] = OutEdgeNum[ps][pos1-1];
				OutEdgeNum[ps][-1]--;
				OutEdgeNum[ps].erase(pos1-1);
				getInVertex()[ps].erase(s);
				break;
			}
		}

	}*/

	getVLabel().erase(s);
	getOutEdge().erase(s);
	getInVertex().erase(s);
	OutEdgeNum.erase(s);
	nodes.push(s);	
	nextVertex.erase(s);
	nextTime.erase(s);
	_TVtoV[tv].erase(timestamp);
	//TimeStampSet.clear();		

}


template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::removeAllInEdges(VertexID s) {
	ASSERT(isVertex(s));

	int pos1, pos2;

	for (typename AdjListBool::iterator it = getInVertex()[s].begin();
		it != getInVertex()[s].end(); it++) {
		VertexID ps = it->first;
		getOutEdge()[ps].erase(s);		
		getInVertex()[ps].erase(s);
		Ecnt--;

		pos1 = OutEdgeNum[ps][-1];
		for(int i = 0; i< pos1;i++){
			if(OutEdgeNum[ps][i] == s){
				OutEdgeNum[ps][i] = OutEdgeNum[ps][pos1-1];
				OutEdgeNum[ps][-1]--;
				OutEdgeNum[ps].erase(pos1-1);
				break;
			}


		}
		//pos2 = OutEdgeNum[ps][s]; //###here
  		//OutEdgeNum[ps][pos2] = OutEdgeNum[ps][pos1];  
		//OutEdgeNum[ps].erase(pos1);
  		//OutEdgeNum[ps][-1]--;
	}
	//getInVertex().erase(s);
	//getOutEdge().erase(s);
	//OutEdgeNum.erase(s);
	OutEdgeNum[s][-1] = 0;
	//nextVertex.erase(s);
}


template<class WeightType, class TimestampType>
void DIGRAPH<WeightType, TimestampType>::removeEdge(VertexID s, VertexID d,
	bool verify) {		
	if (verify){
		ASSERT(isEdge(s, d));		
	}
	Ecnt--;


	int pos1, pos2;
  	pos1 = OutEdgeNum[s][-1];
	for(int i = 0; i<pos1;i++){
		if(OutEdgeNum[s][i] == d){
			OutEdgeNum[s][i] = OutEdgeNum[s][pos1-1];
			OutEdgeNum[s][-1]--;
			OutEdgeNum[s].erase(pos1-1);
			break;
		}
	}

	pos1 = OutEdgeNum[d][-1];
	for(int i = 0; i<pos1;i++){
		if(OutEdgeNum[d][i] == s){
			OutEdgeNum[d][i] = OutEdgeNum[d][pos1-1];
			OutEdgeNum[d][-1]--;
			OutEdgeNum[d].erase(pos1-1);
			break;
		}
	}

	//   remove out edge
	getOutEdge()[s].erase(d);

	// remove in vertex
	getInVertex()[d].erase(s);




	// check isolated
	//if (getDegree(s) == 0) {
	//	eraseVertex(s);
	//}
	//if (getDegree(d) == 0) {
	//	eraseVertex(d);
	//}
}


#endif /* DIGRAPH_H_ */