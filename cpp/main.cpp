#define DETAIL

#include <cstdio>
#include <fstream>
#include <iostream>
#include <ctime>    //or #include <time.h>




#include "header/DataRead.h"
#include "header/BFlow.h"
#include<string>

using namespace std;



int main(int argc, char *argv[]) {

  //cout<<"program start..."<<endl;
   if (argc != 4) {  
     cout<<"Please refer to ReadMe for the usage!"<<endl;
     return 0;
   }

  clock_t startTime, endTime; 
  startTime = clock();
  
  //load the query and flow network data   
	cout << "\n*************************************" << endl;
	cout << "********   Loading Data    **********" << endl;
	cout << "*************************************" << endl;
  DataRead<double, Timestamp> Graph;

  string QueryName = argv[1];
  cout << "Query: "<<QueryName<<endl;
  Graph.ReadQuery(QueryName+ ".qry");

  
	//Loading network data
  string GraphName = argv[2];	
	cout << "\nGraph: "<<GraphName<<endl;
  Graph.ReadGraph(GraphName+ ".gra");
	cout << endl;
   
  //Timestamp d = stoi(argv[3]);    // delta of the bursting flow
  int k = stoi(argv[3]);    // type of experiments
  
  endTime = clock();
  
  cout<<"The query and flow network have been loaded in "<< (double)(endTime - startTime) / CLOCKS_PER_SEC <<" seconds..."<<endl;


  //Querying   
	cout << "\n*************************************" << endl;
	cout << "***   Bursting Flow Computation   ***" << endl;
	cout << "*************************************" << endl;
    
  //srand(0);
  //d = d*(Graph.graph_ptr->Tmax-Graph.graph_ptr->Tmin)/10;
  //VertexID s, t;
  //for(auto it = Graph.query_set.begin(); it != Graph.query_set.end();it++){
  //s = it->first; 
  //t = it->second;


	int pos = GraphName.find_last_of('/');
	string temp(GraphName.substr(pos+1));
	//replace(temp.begin(), temp.end(), '/', '-');

 
  BFlow<double, Timestamp> test(Graph.graph_ptr, Graph.query_set, Graph.querynum, temp);

  
  
  //test.QueryGenerator(16, 3); //   (number of query, reachable with path at least length)
  if(k == 1)
    test.Overall(k);    // k is used here
  if(k == 2)
    test.VaryDelta(k);
  if(k == 3)
    test.VaryNum(k);
  if(k == 4)
    test.TransTime(k);
  if(k == 5)
    test.VarySize(k);

  return 0;  
}