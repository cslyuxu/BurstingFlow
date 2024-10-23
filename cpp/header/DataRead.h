#pragma once
/*
 * data_read_with_label.h
 *
 *  Created on: October 8, 2023
 */

#ifndef DATAREAD_H_
#define DATAREAD_H_

//#include "TGraph.h"
#include "DIGraph.h"

//#include <fstream.h>		For VS2017
#include <fstream>
#include <iostream>
#include <string> 
using namespace std;


//#define     DEFAULTMSGSIZE        2048
//#define     DEFAULTRANDOM         32

template<class WeightType, class TimestampType>
class DataRead {
public:
	// variables
	TGRAPH<WeightType, TimestampType> graph;
	TGRAPH<WeightType, TimestampType>* graph_ptr;
	unordered_map<VertexID, VertexID> query_set;
	int querynum;
	//void Mapping();
	void ReadGraph(string);
	void ReadQuery(string);
	DataRead();
	~DataRead();
	//private:

};



template<class WeightType, class TimestampType>
DataRead<WeightType, TimestampType>::DataRead() {

}


template<class WeightType, class TimestampType>
void DataRead<WeightType, TimestampType>::ReadGraph(string FileName){
	cout<<"Reading graph data!"<<endl;

	ifstream OpenFile(FileName);
	char str[100];
	string temp;
	//int counter = 0;
	VertexID src = 0, dst = 0;
	WeightType weight=0;
	TimestampType timestamp = 0;
	TimestampType Tmax = 0, Tmin = WEIGHT_MAX;
	//unordered_map<VertexID, WeightType> VLabels;
	graph.tempVcount = 0; //for createing the fastVisitVertex for visiting random vertex to generate random queries
	
	while (OpenFile >> src >> dst >> timestamp >> weight) {		
		//cout << src <<" " << dst <<" " << weight <<" " << timestamp <<endl;
		
		//Initial the graph
		graph.insertVertex(src, 0);
		graph.insertVertex(dst, 0);
		if(src == dst)
			continue;		

		graph.insertEdge(src, dst, timestamp, weight);

		if(Tmax < timestamp)
			Tmax = timestamp;
		if(Tmin > timestamp)
			Tmin = timestamp;		
		/*
		while(OpenFile >> temp){
			if(temp == "#"){
				OpenFile >> temp;
				src = stol(temp);
				OpenFile >> temp;
				label = stol(temp);
				VLabels[src] = label;
				graph.insertVertex(src, 0);
				continue;
			}else{
				dst = stol(temp);
				if(src == dst){
					if (!graph.isEdge(src, dst))				//to make sure the graph is not a multi-directed graph
					graph.insertEdge(src, dst, 0);
				}else{
					graph.insertVertex(dst, 0);
					if (!graph.isEdge(src, dst))				//to make sure the graph is not a multi-directed graph
						graph.insertEdge(src, dst, 0);
				}
			}

		}*/	

	}
	//cout<<"     Done!"<<endl;
	cout <<"\nVCnt: " << graph.getVcnt()<<endl;
	cout <<"ECnt: " << graph.getEcnt()<<endl;
	cout <<"Timestamp: " << graph.TimeStampSet.size()<<endl;
	graph.IndexTimestamp(Tmax, Tmin);
	//cout <<"The maximum degree: "<< graph.maxdegree() <<endl;
	cout <<"The average degree: "<< graph.avgdegree() <<endl;
	cout <<"The stddev degree: "<< graph.stddegree() <<endl;

	
	//Mapping();

	/*Set the labels for the vertices*/
	//graph.initVertexLabel(VLabels);
	//for(auto it = graph.TimeStampSet.begin();it != graph.TimeStampSet.end();it++){
	//		cout << *it <<endl;
	//}

	this->graph_ptr = &this->graph;
	//VLabels.clear();
	OpenFile.close();
}



template<class WeightType, class TimestampType>
void DataRead<WeightType, TimestampType>::ReadQuery(string FileName){
	cout<<"Reading query data!"<<endl;

	ifstream OpenFile(FileName);
	int count = 0;
	VertexID src = 0, dst = 0;	
	while (OpenFile >> src >> dst) {		
		//cout << src <<" " << dst << endl;
		query_set[src] = dst;
		count++;
	}
	querynum = count;
	//cout <<"        Done: " << count<<" queries!"<<endl;
	OpenFile.close();
}


/*
template<class WeightType, class TimestampType>
void DataRead<WeightType, TimestampType>::Mapping() {
	//initial the mapping from VertexID to the matrix number
	graph.initMap();
}*/




template<class WeightType, class TimestampType>
DataRead<WeightType, TimestampType>::~DataRead() {}



#endif 
