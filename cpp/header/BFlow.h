#pragma once

#include <map>
#include <vector>

#include "DIGraph.h"
#include <time.h>

template <class WeightType, class TimestampType>
class BFlow
{
public:
  TGRAPH<WeightType, TimestampType> *Tgraph_ptr;
  int delta;
  VertexID source, sink;
  unordered_map<VertexID, VertexID> queryset;
  int querynum;
  string OutputName;
  bool FlagOfPruning;
  int QueryGeneratorLength;
  unordered_map<int, double> ALG1, ALG2, ALG3; // used for recording runtime for VaryNum
  unordered_map<int, int> InCase, DeCase; //used for recording times of incremental Maxflow Computation for VaryNum

  BFlow(TGRAPH<WeightType, TimestampType> *, unordered_map<VertexID, VertexID>, int, string);
  ~BFlow();
  void BFlowBase();
  void OBFlowBase(ofstream &, int, int);
  void MAPE(ofstream &, int, int, bool);
  void MAPSE(ofstream &, int, int, int, bool);
  void Overall(int);
  void VaryDelta(int);
  void VaryNum(int);
  void VarySize(int);
  void TransTime(int);
  void Transformation(TGRAPH<WeightType, TimestampType> *, DIGRAPH<WeightType, TimestampType> *, TimestampType, TimestampType);
  bool Observation2(TGRAPH<WeightType, TimestampType> *, TimestampType, TimestampType, TimestampType, WeightType, WeightType);
  void TransformationMAPE(TGRAPH<WeightType, TimestampType> *, DIGRAPH<WeightType, TimestampType> *, TimestampType, TimestampType, bool);
  double TransformationMAPSE(TGRAPH<WeightType, TimestampType> *, DIGRAPH<WeightType, TimestampType> *, TimestampType, TimestampType, TimestampType);
  bool ValidTimestamp(const TimestampType &, const TimestampType &, const TimestampType &);
  void QueryGenerator(int, int);
  bool Reach(VertexID, TimestampType, TimestampType, int, clock_t, unordered_map<VertexID, unordered_map<TimestampType, bool>> &);
  WeightType Dinic(DIGRAPH<WeightType, TimestampType> *, VertexID, TimestampType, VertexID, TimestampType);
  WeightType ReverseDinic(DIGRAPH<WeightType, TimestampType> *, VertexID, TimestampType, VertexID, TimestampType);
  bool BFS(DIGRAPH<WeightType, TimestampType> *, VertexID, TimestampType, VertexID, TimestampType, unordered_map<VertexID, int> &);
  bool ReverseBFS(DIGRAPH<WeightType, TimestampType> *, VertexID, TimestampType, VertexID, TimestampType, unordered_map<VertexID, int> &);
  WeightType SendFlow(DIGRAPH<WeightType, TimestampType> *, VertexID, WeightType, VertexID, unordered_map<VertexID, int> &, unordered_map<VertexID, int> &);
  // WeightType ReverseSendFlow(DIGRAPH<WeightType, TimestampType> *, VertexID, WeightType, VertexID, unordered_map<VertexID, int> &, unordered_map<VertexID, int> &);
  void TestFlow(ofstream &, TimestampType, TimestampType, int, int);
};

template <class WeightType, class TimestampType>
BFlow<WeightType, TimestampType>::BFlow(TGRAPH<WeightType, TimestampType> *graph_ptr, unordered_map<VertexID, VertexID> query_set, int query_num, string GraphName)
{
  this->Tgraph_ptr = graph_ptr;
  // this-> source = s;
  // this-> sink = t;
  // this-> delta = _d;
  this->delta = 3 * (graph_ptr->Tmax - graph_ptr->Tmin) / 100; // default
  this->queryset = query_set;
  this->querynum = query_num;
  this->OutputName = GraphName;
  FlagOfPruning = false;
  // for(auto it = queryset.begin(); it!= queryset.end();it++){
  //   cout<<it->first<<" "<<it->second<<endl;
  // }
  // cout<<"There are "<<querynum<<" queries!"<<endl;
}

template <class WeightType, class TimestampType>
BFlow<WeightType, TimestampType>::~BFlow(){}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::Overall(int k)
{
  // source = 0, sink = 6; //for verifying correctness.gra
  // source = 150340, sink = 295606; //CTU-13

  int num = 0;
  string OutFileName = "result/overall-" + OutputName;
  ofstream OutFile(OutFileName);
  cout << OutFileName << endl;

  // cout<<"Start the overall evaluation:    delta is : 0."<< this-> delta <<"|T| = "<< this-> delta*(this->Tgraph_ptr->Tmax-this->Tgraph_ptr->Tmin)/100 <<"!"<< endl;
  cout << "Start the overall evaluation:    delta is : " << this->delta << "!" << endl;

  OutFile << "# TAP, TAPE, TAPSE" << endl;

  for (auto it = queryset.begin(); it != queryset.end(); it++)
  {

    this->source = it->first;
    this->sink = it->second;
    cout << "\n***************  Q" << ++num << "  *****************" << endl;
    cout << "Source: " << this->source << "         Sink: " << this->sink << endl;
    cout << "\n";

    // if (num != 7)
    //   continue;

    // type: overall = 1, delta = 2
    // BFlowBase();
    // TestFlow(OutFile, 232, 1127, num, k);
    // TestFlow(OutFile, 297, 1127, num, k);
    OBFlowBase(OutFile, num, 1);
    MAPE(OutFile, num, 1, true);
    // MAPSE: 0 = Dinice, 1 = ReverseDinic
    MAPSE(OutFile, num, 1, 0, true);
    // MAPSE(OutFile, num, k, 1);

    cout << "\n\n"
         << endl;
  }

  OutFile.close();
}


template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::TransTime(int k)
{
  // source = 0, sink = 6; //for verifying correctness.gra
  // source = 150340, sink = 295606; //CTU-13

  int num = 0;
  string OutFileName = "result/transTime-" + OutputName;
  ofstream OutFile(OutFileName);
  cout << OutFileName << endl;

  // cout<<"Start the overall evaluation:    delta is : 0."<< this-> delta <<"|T| = "<< this-> delta*(this->Tgraph_ptr->Tmax-this->Tgraph_ptr->Tmin)/100 <<"!"<< endl;
  cout << "Start the overall evaluation:    delta is : " << this->delta << "!" << endl;

  OutFile << "# TAP, TAPE, TAPSE, transformation" << endl;

  for (auto it = queryset.begin(); it != queryset.end(); it++)
  {

    this->source = it->first;
    this->sink = it->second;
    cout << "\n***************  Q" << ++num << "  *****************" << endl;
    cout << "Source: " << this->source << "         Sink: " << this->sink << endl;
    cout << "\n";

    // if (num != 7)
    //   continue;

    // type: overall = 1, delta = 2
    // BFlowBase();
    // TestFlow(OutFile, 232, 1127, num, k);
    // TestFlow(OutFile, 297, 1127, num, k);
    OBFlowBase(OutFile, num, 3);
    MAPE(OutFile, num, 3, true);
    // MAPSE: 0 = Dinice, 1 = ReverseDinic
    MAPSE(OutFile, num, 3, 0, true);
    // MAPSE(OutFile, num, k, 1);

    cout << "\n\n"
         << endl;
  }

  OutFile.close();
}


template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::VaryDelta(int k)
{
  for (int i = 1; i <= 3; i++)
  {
    if (i == 1)
      this->delta = 3 * (this->Tgraph_ptr->Tmax - this->Tgraph_ptr->Tmin) / 100;
    if (i == 2)
      this->delta = 6 * (this->Tgraph_ptr->Tmax - this->Tgraph_ptr->Tmin) / 100; // default
    if (i == 3)
      this->delta = 9 * (this->Tgraph_ptr->Tmax - this->Tgraph_ptr->Tmin) / 100;

    int num = 0;
    string OutFileName = "result/delta-" + OutputName + "-" + to_string(i);
    ofstream OutFile(OutFileName);
    cout << OutFileName << endl;

    // cout<<"Start the overall evaluation:    delta is : 0."<< this-> delta <<"|T| = "<< this-> delta*(this->Tgraph_ptr->Tmax-this->Tgraph_ptr->Tmin)/100 <<"!"<< endl;
    cout << "Start the <<varying delta>> evaluation:    delta is : " << this->delta << "!" << endl;

    OutFile << "# TAP, TAPE, TAPSE" << endl;

    for (auto it = queryset.begin(); it != queryset.end(); it++)
    {

      this->source = it->first;
      this->sink = it->second;
      cout << "\n***************  Q" << ++num << "  *****************" << endl;
      cout << "Source: " << this->source << "         Sink: " << this->sink << endl;
      cout << "\n";

      // if (num != 7)
      //   continue;

      // type: overall = 1
      // BFlowBase();
      // TestFlow(OutFile, 232, 1127, num, k);
      // TestFlow(OutFile, 297, 1127, num, k);
      OBFlowBase(OutFile, num, 1);
      MAPE(OutFile, num, 1, true);
      // MAPSE: 0 = Dinice, 1 = ReverseDinic
      MAPSE(OutFile, num, 1, 0, true);
      // MAPSE(OutFile, num, k, 1);

      cout << "\n\n"
           << endl;
    }
    OutFile.close();
  }
 
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::VaryNum(int k)
{
  int num = 0;
  int number;
  double time;

  string OutFileName = "result/num-" + OutputName + "-1";
  ofstream OutFile(OutFileName);
  cout << OutFileName << endl;


  //double average, min, max;


  // cout<<"Start the overall evaluation:    delta is : 0."<< this-> delta <<"|T| = "<< this-> delta*(this->Tgraph_ptr->Tmax-this->Tgraph_ptr->Tmin)/100 <<"!"<< endl;
  cout << "Start the <<varying # of incremental MaxFLow computation>> evaluation!" << endl;

  for (auto it = queryset.begin(); it != queryset.end(); it++)
  {

    this->source = it->first;
    this->sink = it->second;
    cout << "\n***************  Q" << ++num << "  *****************" << endl;
    cout << "Source: " << this->source << "         Sink: " << this->sink << endl;
    cout << "\n";

    // if (num != 7)
    //   continue;

    // type: overall = 1
    // BFlowBase();
    // TestFlow(OutFile, 232, 1127, num, k);
    // TestFlow(OutFile, 297, 1127, num, k);
    OBFlowBase(OutFile, num, 2);
    MAPE(OutFile, num, 2, false);
    // MAPSE: 0 = Dinice, 1 = ReverseDinic
    MAPSE(OutFile, num, 2, 0, false);
    // MAPSE(OutFile, num, k, 1);

    cout << "\n\n"
         << endl;
  }



  OutFile << "#  \% of TAPE";
  map<int, int> times;
  map<int, bool> flag;
  for(auto it = ALG2.begin();it != ALG2.end();it++){
    flag[InCase[it->first]] = true;
  }

  map<int, double> average, max, min; 
  for(auto it = ALG2.begin();it != ALG2.end();it++){
    number = it->first;  
    time = it->second;
    cout << InCase[number] <<"  "<<ALG1[number]/time<<endl;
    //OutFile <<endl;
    //OutFile << InCase[number] <<"  "<<ALG1[number]/time;
    if(flag[InCase[number]]){
      flag[InCase[number]] = false;
      average[InCase[number]] = ALG1[number]/time;
      times[InCase[number]] = 1;
      max[InCase[number]] = ALG1[number]/time;
      min[InCase[number]] = ALG1[number]/time;
    }else{
      average[InCase[number]] += ALG1[number]/time;
      if(ALG1[number]/time >  max[InCase[number]])  
        max[InCase[number]] = ALG1[number]/time;
      if(ALG1[number]/time < min[InCase[number]])
        min[InCase[number]] = ALG1[number]/time;
      times[InCase[number]]++;
    }
    //if(tempmap.find(InCase[number]) == tempmap.end())
    //  tempmap[InCase[number]] = time/ALG1[number];
    //else if(time/ALG1[number]>tempmap[InCase[number]])
    //  tempmap[InCase[number]] = time/ALG1[number];
  }

  for(auto it = average.begin(); it != average.end(); it++){
    OutFile <<endl;
    OutFile << it->first <<" "<<it->second/times[it->first]<<"  "<<min[it->first]<<"  "<<max[it->first];
  }
  OutFile.close();

  string OutFileName2 = "result/num-" + OutputName + "-2";
  ofstream OutFile2(OutFileName2);
  cout << OutFileName2 << endl;
  min.clear();
  max.clear();
  average.clear();
  times.clear();
  flag.clear();

  for(auto it = ALG2.begin();it != ALG2.end();it++){
    flag[DeCase[it->first]] = true;
  }



  OutFile2 << "#  \% of TAPSE";
  for(auto it = ALG2.begin();it != ALG2.end();it++){
    number = it->first;  
    time = it->second;
    cout << DeCase[number] <<"  "<<time/ALG3[number]<<endl;
    //OutFile2 << endl;
    //OutFile2 << DeCase[number] <<"  "<<time/ALG3[number];
    if(flag[DeCase[number]]){
      flag[DeCase[number]] = false;
      average[DeCase[number]] = time/ALG3[number];
      times[DeCase[number]] = 1;
      max[DeCase[number]] = time/ALG3[number];
      min[DeCase[number]] = time/ALG3[number];
    }else{
      average[DeCase[number]] += time/ALG3[number];
      if(time/ALG3[number] >  max[DeCase[number]])  
        max[DeCase[number]] = time/ALG3[number];
      if(time/ALG3[number] < min[DeCase[number]])
        min[DeCase[number]] = time/ALG3[number];
      times[DeCase[number]]++;
    }
    //if(tempmap.find(DeCase[number]) == tempmap.end())
    //  tempmap[DeCase[number]] = ALG3[number]/time;
    //else if(ALG3[number]/time>tempmap[DeCase[number]])
    //  tempmap[DeCase[number]] = ALG3[number]/time;
  }

  for(auto it = average.begin(); it != average.end(); it++){
    OutFile2 <<endl;
    OutFile2 << it->first <<" "<<it->second/times[it->first]<<"  "<<min[it->first]<<"  "<<max[it->first];
  }
  //for(auto it = tempmap.begin(); it != tempmap.end(); it++){
   // OutFile2 <<endl;
  //  OutFile2 << it->first <<" "<<it->second;
  //}
  OutFile2.close();
}



template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::VarySize(int k)
{
  FlagOfPruning = true;
  for (int i = 1; i <= 3; i++)
  {
    int num = 0;
    string OutFileName = "result/size-" + OutputName + "-" + to_string(i)+"-";
    ofstream OutFile(OutFileName);
    cout << OutFileName << endl;

    // cout<<"Start the overall evaluation:    delta is : 0."<< this-> delta <<"|T| = "<< this-> delta*(this->Tgraph_ptr->Tmax-this->Tgraph_ptr->Tmin)/100 <<"!"<< endl;
    cout << "Start the <<varying size>> evaluation!" << endl;

    OutFile << "# size, runtime" << endl;

    for (auto it = queryset.begin(); it != queryset.end(); it++)
    {

      this->source = it->first;
      this->sink = it->second;
      cout << "\n***************  Q" << ++num << "  *****************" << endl;
      cout << "Source: " << this->source << "         Sink: " << this->sink << endl;
      cout << "\n";

      // if (num != 7)
      //   continue;

      // type: overall = 1
      // BFlowBase();
      // TestFlow(OutFile, 232, 1127, num, k);
      // TestFlow(OutFile, 297, 1127, num, k);
      if(i == 1)
        OBFlowBase(OutFile, num, 4);
      if(i == 2)
        MAPE(OutFile, num, 4, false);
      // MAPSE: 0 = Dinice, 1 = ReverseDinic
      if(i == 3)
        MAPSE(OutFile, num, 4, 0, false);
      // MAPSE(OutFile, num, k, 1);

      cout << "\n\n"
           << endl;
    }
    OutFile.close();
  }
 
}



template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::BFlowBase()
{
  clock_t startTime, endTime;
  startTime = clock();
  DIGRAPH<WeightType, TimestampType> *Graph_ptr = new DIGRAPH<WeightType, TimestampType>;
  TimestampType ts, te, ts_temp, te_temp;
  WeightType density = -1, tempdensity = 0;

  int source_num = Tgraph_ptr->_OutVNT[this->source][-1];
  if (source_num == 0)
  {
    cout << "The source has no out-going edges!" << endl;
    return;
  }

  int sink_num = Tgraph_ptr->_InVNT[this->sink][-1];
  if (sink_num == 0)
  {
    cout << "The sink has no in-coming edges!" << endl;
    return;
  }

  // return;

  // Line 4 of Base
  TimestampType bound1 = Tgraph_ptr->_OutVNT[this->source][0] - this->delta;
  TimestampType bound2 = Tgraph_ptr->_OutVNT[this->source][Tgraph_ptr->_OutVNT[this->source][-1] - 1];
  // Line 5 of Base
  TimestampType bound3 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1] + this->delta;

  for (ts = Tgraph_ptr->Tmin; ts != Tgraph_ptr->Tmax; ts++)
  {
    // for(auto it1 =Tgraph_ptr->TimeStampSet.begin();it1 != Tgraph_ptr->TimeStampSet.end();it1++){
    // auto temp = it1;
    // ts = *it1;
    // pruning for baseline of source's timestamp
    if (ts < bound1)
      continue;
    if (ts > bound2)
      break;

    for (te = ts + this->delta; te != Tgraph_ptr->Tmax; te++)
    {
      // for(auto it2 =(++temp);it2 != Tgraph_ptr->TimeStampSet.end();it2++){
      // te = *it2;
      // pruning for baseline of sink's timestamp
      // if((te-ts)<this->delta)
      //   continue;
      if (te > bound3)
        break;

      Graph_ptr->reset();
      // cout<<"Transforming for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) !"<<endl;
      Transformation(this->Tgraph_ptr, Graph_ptr, ts, te);

      tempdensity = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
      tempdensity = tempdensity / (te - ts);
      // if(tempdensity != 0)
      cout << "MaxFlow for (<" << this->source << ", " << ts << ">, <" << this->sink << ", " << te << ">) is " << tempdensity << " !" << endl;
      if (tempdensity > density)
      {
        density = tempdensity;
        ts_temp = ts;
        te_temp = te;
      }
      // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;
    }
  }
  endTime = clock();
  cout << "The runtime of Base is " << (double)(endTime - startTime) / CLOCKS_PER_SEC << " seconds!" << endl;
  cout << "The density of BFlow is " << density << " !" << endl;

  delete Graph_ptr;
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::Transformation(TGRAPH<WeightType, TimestampType> *T_ptr, DIGRAPH<WeightType, TimestampType> *ptr, TimestampType start, TimestampType end)
{
  // duplicate source and sink
  VertexID tempV;

  // vertices for s at start
  ptr->insertVertex(this->source, start);
  tempV = ptr->_TVtoV[this->source][start];

  for (int i = 0; i < T_ptr->_OutVNT[this->source][-1]; i++)
  {
    // cout<<"hello!"<<endl;

    // time interval
    if ((T_ptr->_OutVNT[this->source][i] <= start))
      continue;
    if ((T_ptr->_OutVNT[this->source][i] >= end))
      break;

    ptr->insertVertex(this->source, T_ptr->_OutVNT[this->source][i]);
    // insert edges with infy capacities for the same tv
    ptr->insertEdge(tempV, ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]], WEIGHT_MAX);
    ptr->nextVertex[tempV] = ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]];
    ptr->nextTime[tempV] = T_ptr->_OutVNT[this->source][i];

    tempV = ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]];
  }

  // vertices for s at end
  ptr->insertVertex(this->source, end);
  ptr->insertEdge(tempV, ptr->_TVtoV[this->source][end], WEIGHT_MAX);
  ptr->nextVertex[tempV] = ptr->_TVtoV[this->source][end];
  ptr->nextTime[tempV] = end;

  // vertices for t at start
  ptr->insertVertex(this->sink, start);
  tempV = ptr->_TVtoV[this->sink][start];

  for (int i = 0; i < T_ptr->_InVNT[this->sink][-1]; i++)
  {

    // time interval
    if ((T_ptr->_InVNT[this->sink][i] <= start))
      continue;
    if ((T_ptr->_InVNT[this->sink][i] >= end))
      break;

    ptr->insertVertex(this->sink, T_ptr->_InVNT[this->sink][i]);

    // insert edges with infy capacities for the same tv
    ptr->insertEdge(tempV, ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]], WEIGHT_MAX);
    ptr->nextVertex[tempV] = ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]];
    ptr->nextTime[tempV] = T_ptr->_InVNT[this->sink][i];
    // cout<<"The "<<VCnt<<" node is <"<<this->sink <<", "<< T_ptr->_InVNT[this->sink][i] <<"> !"<<endl;
    tempV = ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]];
  }

  // vertices for t at end
  ptr->insertVertex(this->sink, end);
  ptr->insertEdge(tempV, ptr->_TVtoV[this->sink][end], WEIGHT_MAX);
  ptr->nextVertex[tempV] = ptr->_TVtoV[this->sink][end];
  ptr->nextTime[tempV] = end;

  // duplicate the rest nodes
  VertexID x;
  // starting from the source, do the BFS
  x = this->source;
  unordered_map<VertexID, int> map_hop;
  queue<VertexID> nodes;

  nodes.push(x);
  map_hop[x] = 0;
  map_hop[this->sink] = -1;

  // begin BFS
  while (!nodes.empty())
  {
    VertexID v = nodes.front();
    nodes.pop();
    // for each out node u of v
    for (auto it1 = T_ptr->getOutEdge()[v].begin(); it1 != T_ptr->getOutEdge()[v].end(); it1++)
    {
      VertexID u = it1->first;

      /*
      // u is visited
      if(map_hop.find(u) != map_hop.end()){
        for(auto it2 = T_ptr->getOutEdge()[v][u].begin(); it2 != T_ptr->getOutEdge()[v][u].end(); it2++){
            TimestampType _t = it2->first;
            //prune edges not within the time interval [start, end]
            if(!ValidTimestamp(_t, start, end))
              continue;

            //cout<<"yet!"<<endl;
            ptr->insertEdge(ptr->_TVtoV[v][_t], ptr->_TVtoV[u][_t], it2->second.eweight);
            //cout<<"The edge from "<< ptr->_TVtoV[v][_t] << " to " << ptr->_TVtoV[u][_t] <<" with timestamp "<< _t << " has capacity " << it2->second.eweight <<endl;
            //cout<<"yet2!"<<endl;
        }
        continue;
      }*/

      for (auto it3 = T_ptr->getOutEdge()[v][u].begin(); it3 != T_ptr->getOutEdge()[v][u].end(); it3++)
      {
        TimestampType t1 = T_ptr->func[it3->first];
        // prune edges not within the time interval [start, end]
        if (!ValidTimestamp(t1, start, end))
          continue;

        // add u into queue for next iteration
        if (map_hop.find(u) == map_hop.end())
        {

          if (u != this->sink)
            nodes.push(u);

          // vertices for u at start
          ptr->insertVertex(u, start);
          tempV = ptr->_TVtoV[u][start];

          for (int i = 0; i < T_ptr->_VNT[u][-1]; i++)
          {
            // time interval
            if ((T_ptr->_VNT[u][i] <= start))
              continue;
            if ((T_ptr->_VNT[u][i] >= end))
              break;

            TimestampType _t = T_ptr->_VNT[u][i];
            // create nodes for timestamps within the time interval [start, end]
            // if(!ValidTimestamp(_t, start, end))
            //  continue;
            ptr->insertVertex(u, _t);
            // insert edges with infy capacities for the same tv
            ptr->insertEdge(tempV, ptr->_TVtoV[u][_t], WEIGHT_MAX);
            // cout<<"The "<<VCnt<<" node is <"<<u <<", "<< _t <<"> !"<<endl;
            ptr->nextVertex[tempV] = ptr->_TVtoV[u][_t];
            ptr->nextTime[tempV] = _t;
            tempV = ptr->_TVtoV[u][_t];
          }

          // vertices for u at end
          ptr->insertVertex(u, end);
          ptr->insertEdge(tempV, ptr->_TVtoV[u][end], WEIGHT_MAX);
          ptr->nextVertex[tempV] = ptr->_TVtoV[u][end];
          ptr->nextTime[tempV] = end;

          // index design may locate at here
          map_hop[u] = map_hop[v] + 1;
        }

        ptr->insertEdge(ptr->_TVtoV[v][t1], ptr->_TVtoV[u][t1], it3->second.eweight);
        // cout<<"The edge from "<< ptr->_TVtoV[v][_t] << " to " << ptr->_TVtoV[u][_t] <<" with timestamp "<< _t << " has capacity " << it3->second.eweight <<endl;
      }
    }
  }
}

template <class WeightType, class TimestampType>
bool BFlow<WeightType, TimestampType>::ValidTimestamp(const TimestampType &t, const TimestampType &start, const TimestampType &end)
{
  if (t < start)
    return false;
  if (t > end)
    return false;
  return true;
}

template <class WeightType, class TimestampType>
WeightType BFlow<WeightType, TimestampType>::Dinic(DIGRAPH<WeightType, TimestampType> *ptr, VertexID s, TimestampType ts, VertexID t, TimestampType te)
{
  if (s == t)
    return -1;

  WeightType total = 0;
  WeightType flow;

  unordered_map<VertexID, int> level;
  unordered_map<VertexID, int> start;

  while (BFS(ptr, s, ts, t, te, level) == true)
  {
    start.clear();
    // initialization of start
    for (auto it = ptr->getVLabel().begin(); it != ptr->getVLabel().end(); it++)
      start[it->first] = 0;
    while (flow = SendFlow(ptr, s, WEIGHT_MAX, t, level, start))
      total += flow;
  }

  level.clear();
  return total;
}

// Finds if more flow can be sent from s to t.
// Also assigns levels to nodes.
template <class WeightType, class TimestampType>
bool BFlow<WeightType, TimestampType>::BFS(DIGRAPH<WeightType, TimestampType> *ptr, VertexID s, TimestampType ts, VertexID t, TimestampType te, unordered_map<VertexID, int> &level)
{
  queue<VertexID> nodes;

  for (auto it = ptr->getVLabel().begin(); it != ptr->getVLabel().end(); it++)
    level[it->first] = -1;

  level[s] = 0; // level of source
  nodes.push(s);

  // begin BFS
  while (!nodes.empty())
  {
    VertexID v = nodes.front();
    nodes.pop();
    // for each out node u of v
    for (auto it1 = ptr->getOutEdge()[v].begin(); it1 != ptr->getOutEdge()[v].end(); it1++)
    {
      VertexID u = it1->first;

      // cout<<"("<<v<<","<<u<<") "<<ptr->getOutEdge()[v][u].flow <<" "<<  ptr->getOutEdge()[v][u].eweight<<endl;
      if ((level[u] < 0) && (ptr->getOutEdge()[v][u].flow < ptr->getOutEdge()[v][u].eweight))
      {
        level[u] = level[v] + 1;
        // add u into queue for next iteration
        nodes.push(u);
      }
    }
  }
  return level[t] < 0 ? false : true;
}

// A DFS based function to send flow after BFS has
// figured out that there is a possible flow and
// constructed levels. This function called multiple
// times for a single call of BFS.
// flow : Current flow send by parent function call
// u : Current vertex
// t : Sink
// https://www.zhihu.com/question/263856769/answer/889283306
template <class WeightType, class TimestampType>
WeightType BFlow<WeightType, TimestampType>::SendFlow(DIGRAPH<WeightType, TimestampType> *ptr, VertexID u, WeightType flow, VertexID t, unordered_map<VertexID, int> &level, unordered_map<VertexID, int> &start)
{
  // Sink reached
  if (u == t)
    return flow;

  // Traverse all adjacent edges one -by - one.
  for (; start[u] < ptr->OutEdgeNum[u][-1]; start[u]++)
  {
    VertexID v = ptr->OutEdgeNum[u][start[u]];

    if ((level[v] == level[u] + 1) && (ptr->getOutEdge()[u][v].flow < ptr->getOutEdge()[u][v].eweight))
    {
      // find minimum flow from u to t
      WeightType curr_flow = min(flow, ptr->getOutEdge()[u][v].eweight - ptr->getOutEdge()[u][v].flow);

      WeightType temp_flow = SendFlow(ptr, v, curr_flow, t, level, start);

      // flow is greater than zero
      if (temp_flow > 0)
      {
        // add flow to current edge
        ptr->getOutEdge()[u][v].flow += temp_flow;

        // subtract flow from reverse edge of current edge
        ptr->getOutEdge()[v][u].flow -= temp_flow;
        return temp_flow;
      }
    }
  }

  return 0;
}

template <class WeightType, class TimestampType>
WeightType BFlow<WeightType, TimestampType>::ReverseDinic(DIGRAPH<WeightType, TimestampType> *ptr, VertexID s, TimestampType ts, VertexID t, TimestampType te)
{
  if (s == t)
    return -1;

  WeightType total = 0;
  WeightType flow;

  unordered_map<VertexID, int> level;
  unordered_map<VertexID, int> start;

  while (ReverseBFS(ptr, s, ts, t, te, level) == true)
  {
    start.clear();
    // initialization of start
    for (auto it = ptr->getVLabel().begin(); it != ptr->getVLabel().end(); it++)
      start[it->first] = 0;
    while (flow = SendFlow(ptr, t, WEIGHT_MAX, s, level, start))
      total += flow;
  }

  level.clear();
  // ptr->maxflow = total;
  return total;
}

template <class WeightType, class TimestampType>
bool BFlow<WeightType, TimestampType>::ReverseBFS(DIGRAPH<WeightType, TimestampType> *ptr, VertexID s, TimestampType ts, VertexID t, TimestampType te, unordered_map<VertexID, int> &level)
{
  queue<VertexID> nodes;

  for (auto it = ptr->getVLabel().begin(); it != ptr->getVLabel().end(); it++)
    level[it->first] = -1;

  level[t] = 0; // level of source
  nodes.push(t);

  // begin BFS
  while (!nodes.empty())
  {
    VertexID v = nodes.front();
    nodes.pop();
    // for each out node u of v
    for (auto it1 = ptr->getInVertex()[v].begin(); it1 != ptr->getInVertex()[v].end(); it1++)
    {
      VertexID u = it1->first;

      if ((level[u] < 0) && (ptr->getOutEdge()[u][v].flow < ptr->getOutEdge()[u][v].eweight))
      {
        level[u] = level[v] + 1;
        // add u into queue for next iteration
        nodes.push(u);
      }
    }
  }
  return level[s] < 0 ? false : true;
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::MAPE(ofstream &OutFile, int num, int type, bool prune)
{
  double avg_vertex, avg_edge;
  double _runtime = 0, _density = 0, _avg_vertex = 0, _avg_edge = 0, transTime = 0;
  int countnum, _countnum = 0, incase = 0;
  clock_t startTime, endTime, tempstart, tempend, TransformationStart, TransformationEnd;

  DIGRAPH<WeightType, TimestampType> *Graph_ptr = new DIGRAPH<WeightType, TimestampType>;
  TimestampType ts, te, ts_temp, te_temp, te_previous, te_last;
  WeightType density, tempdensity, tempflow, lastflow;

  int source_num = Tgraph_ptr->_OutVNT[this->source][-1];
  if (source_num == 0)
  {
    cout << "The source has no out-going edges!" << endl;
    return;
  }

  int sink_num = Tgraph_ptr->_InVNT[this->sink][-1];
  if (sink_num == 0)
  {
    cout << "The sink has no in-coming edges!" << endl;
    return;
  }

  // return;

  // Line 4 of Base
  TimestampType bound1 = Tgraph_ptr->_OutVNT[this->source][0];
  if (bound1 < 0)
    bound1 = 0;
  TimestampType bound2 = Tgraph_ptr->_OutVNT[this->source][Tgraph_ptr->_OutVNT[this->source][-1] - 1];
  // Line 5 of Base
  TimestampType bound3 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1] + this->delta;
  TimestampType bound4 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1];

  for (int rp = 0; rp < REPEAT; rp++)
  {
    Graph_ptr->reset();
    startTime = clock();
    countnum = 0;
    avg_vertex = 0, avg_edge = 0;
    density = -1;
    for (int j = 0; j < source_num; j++)
    {
      ts = Tgraph_ptr->_OutVNT[this->source][j];
      // for(ts = bound1; ts <= bound4; ts++){
      // for(auto it1 =Tgraph_ptr->TimeStampSet.begin();it1 != Tgraph_ptr->TimeStampSet.end();it1++){
      //   auto temp = it1;
      //   ts = *it1;
      // pruning for baseline of source's timestamp
      if (ts > bound2)
        break;

      te = ts + this->delta;
      if (te > Tgraph_ptr->Tmax)
        break;
      te_previous = te;
      int count = -1; // find the minimum timestamp of t
      for (int i = 0; i < sink_num; i++)
      {
        if (Tgraph_ptr->_InVNT[this->sink][i] > te)
        {
          count = i;
          break;
        }
      }

      Graph_ptr->reset();
      TransformationStart = clock();
      Transformation(Tgraph_ptr, Graph_ptr, ts, te);
      TransformationEnd = clock();
      transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
      // cout<<Graph_ptr->getEcnt()<<" "<<Graph_ptr->getVcnt()<<endl;
      tempstart = clock();
      tempflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
      Graph_ptr->maxflow = tempflow;
      tempend = clock();

      countnum++;
      avg_edge = avg_edge + Graph_ptr->getEcnt();
      avg_vertex = avg_vertex + Graph_ptr->getVcnt();

      lastflow = tempflow; // used for Observation2
      te_last = te;        // used for Observation2
      tempdensity = tempflow / (te - ts);
      // if(tempdensity != 0)
      // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !  runtime:"<<(double)(tempend - tempstart)  / CLOCKS_PER_SEC << endl;
      // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
      if (tempdensity > density)
      {
        density = tempdensity;
        ts_temp = ts;
        te_temp = te;
      }
      // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;

      if (count == -1)
        continue;

      // Observation 1

      for (int i = count; i < sink_num; i++)
      {
        te = Tgraph_ptr->_InVNT[this->sink][i];
        
        TransformationStart = clock();
        tempstart = clock();
        TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_previous, te, true);
        TransformationEnd = clock();
        transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
        // cout<<Graph_ptr->getEcnt()<<" "<<Graph_ptr->getVcnt()<<endl;
        // Observation 2
        if (Observation2(this->Tgraph_ptr, ts, te_last, te, lastflow, density)&&prune)
        {
          te_previous = te;
          continue;
        }  

        if(FlagOfPruning)
          tempstart = clock();

        tempflow += Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
        tempend = clock();
        if(type == 4){
          OutFile << Graph_ptr->Vcnt<<" "<< (double)(tempend - tempstart) / CLOCKS_PER_SEC<<endl;
        }
        incase++;
        Graph_ptr->maxflow = tempflow;
        // lastflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
        // tempflow += lastflow;
        // cout<<"FLow: "<<lastflow<<" ts: "<<ts<<" te:"<<te<<endl;
        // cout<< Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[this->source][3]][Graph_ptr->_TVtoV[2][3]].eweight<< " " << Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[2][3]][Graph_ptr->_TVtoV[2][4]].eweight<< " " << Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[2][4]][Graph_ptr->_TVtoV[3][4]].eweight<< " " << Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[3][4]][Graph_ptr->_TVtoV[4][4]].eweight<< " " << Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[4][4]][Graph_ptr->_TVtoV[6][4]].eweight<< " " << endl;
        // cout<< Graph_ptr->getOutEdge()[Graph_ptr->_TVtoV[6][3]][Graph_ptr->_TVtoV[6][4]].flow << endl;

        countnum++;
        avg_edge = avg_edge + Graph_ptr->getEcnt();
        avg_vertex = avg_vertex + Graph_ptr->getVcnt();

        lastflow = tempflow;
        te_last = te;
        tempdensity = tempflow / (te - ts);
        // if(tempdensity != 0)
        // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
        if (tempdensity > density)
        {
          density = tempdensity;
          ts_temp = ts;
          te_temp = te;
        }

        te_previous = te;
        // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;
      }
    }
    endTime = clock();
    _avg_vertex += avg_vertex / countnum;
    _avg_edge += avg_edge / countnum;
    _countnum += countnum;
    _runtime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
    _density += density;
  }

  // cout<<"The runtime of MAPE is " << (double)(endTime - startTime)  / CLOCKS_PER_SEC <<" seconds!"<<endl;
  // cout<<"The density of BFlow is " << density << " !"<<endl;
  // cout<<"Avg_Edge: "<<avg_edge /countnum << "   Avg_Vertex: "<<avg_vertex/countnum << "   # of Dinic: " << countnum<<endl;
  cout << "The runtime of MAPE is " << _runtime / REPEAT << " seconds!" << endl;
  cout << "The density of BFlow is " << _density / REPEAT << "during [" << ts_temp << "," << te_temp << "] !" << endl;
  cout << "Avg_Edge: " << _avg_edge / REPEAT << "   Avg_Vertex: " << _avg_vertex / REPEAT << "   # of Dinic: " << _countnum / REPEAT << endl;

  if (type == 1)
    OutFile << "  " << _runtime / REPEAT << "  ";
  if (type == 2)
  {
    ALG2[num] = _runtime / REPEAT;
    InCase[num] = incase / REPEAT;
  }
  if (type == 3)
    OutFile << "  " << _runtime / REPEAT << "  " << transTime / REPEAT << "  ";

  delete Graph_ptr;
}

template <class WeightType, class TimestampType>
bool BFlow<WeightType, TimestampType>::Observation2(TGRAPH<WeightType, TimestampType> *T_ptr, TimestampType ts, TimestampType te_last, TimestampType te, WeightType lastflow, WeightType density)
{
  WeightType tempflow = lastflow;
  VertexID InNeighbor;
  TimestampType temptime;
  for (auto it = T_ptr->getInVertex()[this->sink].begin(); it != T_ptr->getInVertex()[this->sink].end(); it++)
  {
    InNeighbor = it->first;
    for (auto it2 = T_ptr->getOutEdge()[InNeighbor][this->sink].begin(); it2 != T_ptr->getOutEdge()[InNeighbor][this->sink].end(); it2++)
    {
      temptime = T_ptr->func[it2->first];
      if ((temptime < te_last) || (temptime > te))
        continue;
      // tempflow += T_ptr->getOutEdge()[InNeighbor][this->sink][temptime].eweight;
      tempflow += it2->second.eweight; // better description
    }
  }

  // if((te - ts)<= 0)
  //   cout<<" the time interval is wrong! "<<endl;
  if (tempflow <= (density * (te - ts)))
  {
    return true;
  }
  else
    return false;
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::TransformationMAPE(TGRAPH<WeightType, TimestampType> *T_ptr, DIGRAPH<WeightType, TimestampType> *ptr, TimestampType start, TimestampType end, bool flagOfZero)
{
  // duplicate source and sink
  // int VCnt = ptr->getVcnt(); //VertexID for vertices in transfomred networks
  bool flag, check;
  VertexID tempV;

  if (start == end)
    return;

  // vertices for s
  flag = false;
  if (ptr->_TVtoV[this->source].find(start) == ptr->_TVtoV[this->source].end())
  {
    ptr->insertVertex(this->source, start);
    // VCnt++;
    cout << "case 1 that should not happen!" << endl;
  }

  tempV = ptr->_TVtoV[this->source][start];

  for (int i = 0; i < T_ptr->_OutVNT[this->source][-1]; i++)
  {
    // cout<<"hello!"<<endl;

    // time interval
    if ((T_ptr->_OutVNT[this->source][i] <= start))
      continue;
    if ((T_ptr->_OutVNT[this->source][i] >= end))
      break;

    ptr->insertVertex(this->source, T_ptr->_OutVNT[this->source][i]);
    // insert edges with infy capacities for the same tv
    ptr->insertEdge(tempV, ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]], WEIGHT_MAX);
    ptr->nextVertex[tempV] = ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]];
    ptr->nextTime[tempV] = T_ptr->_OutVNT[this->source][i];
    // cout<<"The "<<VCnt<<" node is <"<<this->source <<", "<< T_ptr->_OutVNT[this->source][i] <<"> !"<<endl;
    // VCnt++;
    tempV = ptr->_TVtoV[this->source][T_ptr->_OutVNT[this->source][i]];
  }

  // vertices for s at end
  ptr->insertVertex(this->source, end);
  ptr->insertEdge(tempV, ptr->_TVtoV[this->source][end], WEIGHT_MAX);
  ptr->nextVertex[tempV] = ptr->_TVtoV[this->source][end];
  ptr->nextTime[tempV] = end;
  // VCnt++;

  // vertices for t at start, for the sink, the flow should be also added
  if (ptr->_TVtoV[this->sink].find(start) == ptr->_TVtoV[this->sink].end())
  {
    ptr->insertVertex(this->sink, start);
    // VCnt++;
    cout << "case 2 that should not happen!" << endl;
  }

  tempV = ptr->_TVtoV[this->sink][start];

  for (int i = 0; i < T_ptr->_InVNT[this->sink][-1]; i++)
  {

    // time interval
    if ((T_ptr->_InVNT[this->sink][i] <= start))
      continue;
    if ((T_ptr->_InVNT[this->sink][i] >= end))
      break;

    ptr->insertVertex(this->sink, T_ptr->_InVNT[this->sink][i]);
    // insert edges with infy capacities for the same tv
    ptr->insertEdge(tempV, ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]], WEIGHT_MAX);
    ptr->nextVertex[tempV] = ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]];
    ptr->nextTime[tempV] = T_ptr->_InVNT[this->sink][i];
    if (flagOfZero)
    {
      ptr->getOutEdge()[tempV][ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]]].flow = ptr->maxflow;
      ptr->getOutEdge()[ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]]][tempV].flow = -ptr->maxflow;
    }
    tempV = ptr->_TVtoV[this->sink][T_ptr->_InVNT[this->sink][i]];
  }

  // vertices for t at end
  // if (ptr->_TVtoV[this->sink].find(end) != ptr->_TVtoV[this->sink].end())
  //  cout << this->source << ":" << start << ", " << this->sink << ":" << end << endl;
  ptr->insertVertex(this->sink, end); // ###here
  ptr->insertEdge(tempV, ptr->_TVtoV[this->sink][end], WEIGHT_MAX);
  ptr->nextVertex[tempV] = ptr->_TVtoV[this->sink][end];
  ptr->nextTime[tempV] = end;
  ptr->getOutEdge()[tempV][ptr->_TVtoV[this->sink][end]].flow = ptr->maxflow;
  ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][tempV].flow = -ptr->maxflow;

  // duplicate the rest nodes
  VertexID x;
  // starting from the source, do the BFS
  x = this->source;
  unordered_map<VertexID, int> map_hop;
  queue<VertexID> nodes;

  nodes.push(x);
  map_hop[x] = 0; // whether was visited
  map_hop[this->sink] = -1;

  for (auto it = ptr->_tlabels.begin(); it != ptr->_tlabels.end(); it++)
  {
    VertexID u = it->first;

    if ((u == this->source) || (u == this->sink))
      continue;

    // map_hop[u] = -1;
    tempV = ptr->_TVtoV[u][start];
    for (int i = 0; i < T_ptr->_VNT[u][-1]; i++)
    {
      TimestampType _t = T_ptr->_VNT[u][i];
      // time interval
      if ((_t <= start))
        continue;
      if ((_t >= end))
        break;

      // create nodes for timestamps within the time interval [start, end]
      // if(!ValidTimestamp(_t, start, end))
      //  continue;
      ptr->insertVertex(u, _t); // ###here
      // insert edges with infy capacities for the same tv
      ptr->insertEdge(tempV, ptr->_TVtoV[u][_t], WEIGHT_MAX);
      ptr->nextVertex[tempV] = ptr->_TVtoV[u][_t];
      ptr->nextTime[tempV] = _t;
      // cout<<"The "<<VCnt<<" node is <"<<u <<", "<< _t <<"> !"<<endl;
      // VCnt++;
      tempV = ptr->_TVtoV[u][_t];
    }

    // vertices for u at end
    ptr->insertVertex(u, end); // ###here
    ptr->insertEdge(tempV, ptr->_TVtoV[u][end], WEIGHT_MAX);
    ptr->nextVertex[tempV] = ptr->_TVtoV[u][end];
    ptr->nextTime[tempV] = end;
    nodes.push(u);
    // VCnt++;
  }

  // begin BFS
  while (!nodes.empty())
  {
    VertexID v = nodes.front();
    nodes.pop();
    map_hop[v] = -1;
    //  for each out node u of v
    for (auto it1 = T_ptr->getOutEdge()[v].begin(); it1 != T_ptr->getOutEdge()[v].end(); it1++)
    {
      VertexID u = it1->first;

      for (auto it3 = T_ptr->getOutEdge()[v][u].begin(); it3 != T_ptr->getOutEdge()[v][u].end(); it3++)
      {
        TimestampType __t = it3->first;
        // prune edges not within the time interval [start, end]

        if (T_ptr->func[__t] <= start)
          continue;

        if (T_ptr->func[__t] > end)
          continue;

        if (map_hop.find(u) == map_hop.end())
        {
          // add u into queue for next iteration
          // if(u != this->sink)
          nodes.push(u);

          // vertices for u at start
          // if (ptr->_TVtoV[u].find(start) == ptr->_TVtoV[u].end())
          //{
          ptr->insertVertex(u, start);
          // VCnt++;
          //}

          tempV = ptr->_TVtoV[u][start];

          for (int i = 0; i < T_ptr->_VNT[u][-1]; i++)
          {
            // time interval
            if ((T_ptr->_VNT[u][i] <= start))
              continue;
            if ((T_ptr->_VNT[u][i] >= end))
              break;

            TimestampType _t = T_ptr->_VNT[u][i];
            // create nodes for timestamps within the time interval [start, end]
            // if(!ValidTimestamp(_t, start, end))
            //  continue;
            ptr->insertVertex(u, _t); // ###here
            // insert edges with infy capacities for the same tv
            ptr->insertEdge(tempV, ptr->_TVtoV[u][T_ptr->_VNT[u][i]], WEIGHT_MAX);
            ptr->nextVertex[tempV] = ptr->_TVtoV[u][T_ptr->_VNT[u][i]];
            ptr->nextTime[tempV] = T_ptr->_VNT[u][i];
            tempV = ptr->_TVtoV[u][T_ptr->_VNT[u][i]];
          }
          // vertices for u at end
          ptr->insertVertex(u, end); // ###here
          ptr->insertEdge(tempV, ptr->_TVtoV[u][end], WEIGHT_MAX);
          ptr->nextVertex[tempV] = ptr->_TVtoV[u][end];
          ptr->nextTime[tempV] = end;

          // index design may locate at here
          map_hop[u] = map_hop[v] + 1;
        }
        ptr->insertEdge(ptr->_TVtoV[v][T_ptr->func[__t]], ptr->_TVtoV[u][T_ptr->func[__t]], it3->second.eweight);
        // cout<<"The edge from "<< v << " to " << u <<" with timestamp "<< _t << " has capacity " << it3->second.eweight <<endl;
      }
    }
  }
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::OBFlowBase(ofstream &OutFile, int num, int type)
{
  double avg_vertex, avg_edge;
  double _runtime = 0, _density = 0, _avg_vertex = 0, _avg_edge = 0, transTime = 0;
  int countnum, _countnum = 0;
  clock_t startTime, endTime, tempstart, tempend, TransformationStart, TransformationEnd;

  DIGRAPH<WeightType, TimestampType> *Graph_ptr = new DIGRAPH<WeightType, TimestampType>;
  TimestampType ts, te, ts_temp, te_temp, te_previous, te_last;
  WeightType density, tempdensity, tempflow, lastflow;

  int source_num = Tgraph_ptr->_OutVNT[this->source][-1];
  if (source_num == 0)
  {
    cout << "The source has no out-going edges!" << endl;
    return;
  }

  int sink_num = Tgraph_ptr->_InVNT[this->sink][-1];
  if (sink_num == 0)
  {
    cout << "The sink has no in-coming edges!" << endl;
    return;
  }

  // return;

  // Line 4 of Base

  TimestampType bound1 = Tgraph_ptr->_OutVNT[this->source][0];
  if (bound1 < 0)
    bound1 = 0;
  TimestampType bound2 = Tgraph_ptr->_OutVNT[this->source][Tgraph_ptr->_OutVNT[this->source][-1] - 1];
  // Line 5 of Base
  TimestampType bound3 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1] + this->delta;
  TimestampType bound4 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1];

  for (int rp = 0; rp < REPEAT; rp++)
  {
    Graph_ptr->reset();
    startTime = clock();
    countnum = 0;
    avg_vertex = 0, avg_edge = 0;
    density = -1;
    for (int j = 0; j < source_num; j++)
    {
      ts = Tgraph_ptr->_OutVNT[this->source][j];
      // for(ts = bound1; ts <= bound4; ts++){
      // for(auto it1 =Tgraph_ptr->TimeStampSet.begin();it1 != Tgraph_ptr->TimeStampSet.end();it1++){
      //   auto temp = it1;
      //   ts = *it1;
      // pruning for baseline of source's timestamp
      if (ts > bound2)
        break;

      te = ts + this->delta;
      if (te > Tgraph_ptr->Tmax)
        break;
      te_previous = te;
      int count = -1; // find the minimum timestamp of t
      for (int i = 0; i < sink_num; i++)
      {
        if (Tgraph_ptr->_InVNT[this->sink][i] > te)
        {
          count = i;
          break;
        }
      }

      Graph_ptr->reset();
      TransformationStart = clock();
      tempstart = clock();
      Transformation(Tgraph_ptr, Graph_ptr, ts, te);
      TransformationEnd = clock();

      if(FlagOfPruning)
        tempstart = clock();

      transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
      // cout<<Graph_ptr->getEcnt()<<" "<<Graph_ptr->getVcnt()<<endl;
      tempflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
      tempend = clock();
      if(type == 4){
        OutFile << Graph_ptr->Vcnt<<" "<< (double)(tempend - tempstart) / CLOCKS_PER_SEC<<endl;
      }
      countnum++;
      avg_edge = avg_edge + Graph_ptr->getEcnt();
      avg_vertex = avg_vertex + Graph_ptr->getVcnt();

      lastflow = tempflow; // used for Observation2
      te_last = te;        // used for Observation2
      tempdensity = tempflow / (te - ts);
      // if(tempdensity != 0)
      //   cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !  runtime:"<<(double)(tempend - tempstart)  / CLOCKS_PER_SEC << endl;
      if (tempdensity > density)
      {
        density = tempdensity;
        ts_temp = ts;
        te_temp = te;
      }
      // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;

      // cout<<"hello!"<<endl;

      if (count == -1)
        continue;

      // Observation 1

      for (int i = count; i < sink_num; i++)
      {
        te = Tgraph_ptr->_InVNT[this->sink][i];
        Graph_ptr->reset();        
        TransformationStart = clock();
        tempstart = clock();
        Transformation(Tgraph_ptr, Graph_ptr, ts, te);
        TransformationEnd = clock();

        if(FlagOfPruning)
          tempstart = clock();

        transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
        // cout<<Graph_ptr->getEcnt()<<" "<<Graph_ptr->getVcnt()<<endl;
        tempflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
        tempend = clock();
        if(type == 4){
          OutFile << Graph_ptr->Vcnt<<" "<< (double)(tempend - tempstart) / CLOCKS_PER_SEC<<endl;
        }
        countnum++;
        avg_edge = avg_edge + Graph_ptr->getEcnt();
        avg_vertex = avg_vertex + Graph_ptr->getVcnt();

        lastflow = tempflow;
        te_last = te;
        tempdensity = tempflow / (te - ts);
        // if(tempdensity != 0)
        //  cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !  runtime:"<<(double)(tempend - tempstart)  / CLOCKS_PER_SEC << endl;
        if (tempdensity > density)
        {
          density = tempdensity;
          ts_temp = ts;
          te_temp = te;
        }

        te_previous = te;
        // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;
      }
    }
    endTime = clock();
    _avg_vertex += avg_vertex / countnum;
    _avg_edge += avg_edge / countnum;
    _countnum += countnum;
    _runtime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
    _density += density;
  }

  cout << "The runtime of OBase is " << _runtime / REPEAT << " seconds!" << endl;
  cout << "The density of BFlow is " << _density / REPEAT << "during [" << ts_temp << "," << te_temp << "] !" << endl;
  cout << "Avg_Edge: " << _avg_edge / REPEAT << "   Avg_Vertex: " << _avg_vertex / REPEAT << "   # of Dinic: " << _countnum / REPEAT << endl;

  if (type == 1)
    OutFile << num << "  " << source_num << "  " << sink_num << "  " << _runtime / REPEAT;
  if (type == 2)
  {
    ALG1[num] = _runtime / REPEAT;
  }
  if (type == 3)
    OutFile << num << "  " << source_num << "  " << sink_num << "  " << _runtime / REPEAT << "  " << transTime / REPEAT;


  delete Graph_ptr;
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::MAPSE(ofstream &OutFile, int num, int type, int DinicType, bool prune)
{
  double avg_vertex, avg_edge, flowvaluetest;
  double _runtime = 0, _density = 0, _avg_vertex = 0, _avg_edge = 0, transTime = 0;
  int countnum, _countnum = 0, decase = 0;
  clock_t startTime, endTime, copystart, copyend, trans_end, trans_start, t1, t2, tempstart, tempend, TransformationStart, TransformationEnd;
  double copytime = 0, transtime = 0, dinictime = 0;
  int Timeofcopy = 0, count;
  startTime = clock();
  TimestampType ts, te, ts_temp, te_temp, te_previous, te_last, ts_previous, te_mark, te_mark_previous; // te_previous for incremental transformation, te_last for Observation 2
  WeightType density, tempdensity, tempflow, lastflow, tempvalue, flow_mark;
  bool flag;

  int source_num = Tgraph_ptr->_OutVNT[this->source][-1];
  if (source_num == 0)
  {
    cout << "The source has no out-going edges!" << endl;
    return;
  }

  int sink_num = Tgraph_ptr->_InVNT[this->sink][-1];
  if (sink_num == 0)
  {
    cout << "The sink has no in-coming edges!" << endl;
    return;
  }

  // cout<<source_num<<"     "<<sink_num<<endl;

  // Line 4 of Base
  TimestampType bound1 = Tgraph_ptr->_OutVNT[this->source][0];
  if (bound1 < 0)
    bound1 = 0;
  TimestampType bound2 = Tgraph_ptr->_OutVNT[this->source][Tgraph_ptr->_OutVNT[this->source][-1] - 1];
  // Line 5 of Base
  TimestampType bound3 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1] + this->delta;
  TimestampType bound4 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1];

  bool FlagOfFirst;

  for (int rp = 0; rp < REPEAT; rp++)
  {
    DIGRAPH<WeightType, TimestampType> *Graph_ptr = new DIGRAPH<WeightType, TimestampType>;
    DIGRAPH<WeightType, TimestampType> *Temp_ptr = new DIGRAPH<WeightType, TimestampType>;
    Graph_ptr->reset();
    Temp_ptr->reset();
    startTime = clock();
    countnum = 0;
    avg_vertex = 0, avg_edge = 0;
    density = -1;
    FlagOfFirst = true;

    for (int j = 0; j < source_num; j++)
    {
      ts = Tgraph_ptr->_OutVNT[this->source][j];
      // cout<<"j: "<<j<<"        ts: "<<ts<<endl;
      //  for(ts = bound1; ts <= bound4; ts++){
      //  for(auto it1 =Tgraph_ptr->TimeStampSet.begin();it1 != Tgraph_ptr->TimeStampSet.end();it1++){
      //    auto temp = it1;
      //    ts = *it1;
      //  pruning for baseline of source's timestamp
      if (ts > bound2)
        break;

      // cout<<"begin: "<<ts<<" "<<te<<endl;

      // First Dinic
      if (FlagOfFirst)
      {
        te = ts + this->delta;
        if (te > Tgraph_ptr->Tmax)
          break;
        TransformationStart = clock();
        Transformation(Tgraph_ptr, Graph_ptr, ts, te);
        TransformationEnd = clock();
        transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
        tempstart = clock();
        tempflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
        tempend = clock();
        Graph_ptr->maxflow = tempflow;

        countnum++;
        avg_edge = avg_edge + Graph_ptr->getEcnt();
        avg_vertex = avg_vertex + Graph_ptr->getVcnt();
        // cout<< avg_edge <<" "<<avg_vertex<<endl;

        lastflow = tempflow;  // used for Observation2
        flow_mark = tempflow; // used for the second case of flow value in the incremental ts
        te_last = te;         // used for Observation2
        tempdensity = tempflow / (te - ts);
        // if(tempdensity != 0)
        //      cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
        if (tempdensity > density)
        {
          density = tempdensity;
          ts_temp = ts;
          te_temp = te;
        }
        // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;

        // Below are incremental Dinic
        flag = false; // record whether copy is needed
        count = -1;   // find the minimum timestamp of t

        te_previous = te;
        te_mark_previous = te;
        ts_previous = ts;
        te_mark = te;

        for (int i = 0; i < sink_num; i++)
        {
          // cout<<"te: "<< Tgraph_ptr->_InVNT[this->sink][i]<<endl;
          if (Tgraph_ptr->_InVNT[this->sink][i] > te)
          {
            // cout<<"hello1"<<endl;
            te_mark = te;
            count = i;
            break;
          }
        }

        if (count == -1)
        {
          FlagOfFirst = false;
          continue;
        }

        // cout<<"count: "<<count<<endl;
        // cout<<"<<"<<ts<<", "<<te<<">>"<<endl;
        //       if(j+1==source_num)
        //         continue;
        if ((j + 1) != source_num)
        {
          for (; count < sink_num; count++)
          {
            te = Tgraph_ptr->_InVNT[this->sink][count];
            if (te >= (Tgraph_ptr->_OutVNT[this->source][j + 1] + this->delta))
            { // find the timestamp of te that is just larger than the next ts + delta
              // cout << "hello2" << endl;
              // te_mark = te;
              flag = true;
              // cout<<"1: "<<ts<<"   "<<te_mark_previous<<"    "<<te<<endl;
              // cout<<"hi"<<endl;
              break;
            }
            // cout << "case 0:" << ts << "~" << te << endl;
            //  cout<<"hello0"<<endl;
            TransformationStart = clock();
            TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_previous, te, true);
            TransformationEnd = clock();
            transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
            te_mark_previous = te;
            te_previous = te;
            // cout<< Temp_ptr->Vcnt <<" "<<Graph_ptr->Vcnt<<endl;
            // cout<<"Max for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
            // Observation 2
            if (Observation2(this->Tgraph_ptr, ts, te_last, te, lastflow, density)&&prune)
            {
              continue;
            }

            tempstart = clock();
            tempvalue = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
            tempend = clock();
            tempflow += tempvalue;
            Graph_ptr->maxflow = tempflow;


            countnum++;
            avg_edge = avg_edge + Graph_ptr->getEcnt();
            avg_vertex = avg_vertex + Graph_ptr->getVcnt();

            lastflow = tempflow;
            flow_mark = tempflow;
            te_last = te;
            tempdensity = tempflow / (te - ts);
            // if (tempdensity != 0)
            //  cout << "MaxFlow for (<" << this->source << ", " << ts << ">, <" << this->sink << ", " << te << ">) is " << tempdensity << " !" << endl;
            if (tempdensity > density)
            {
              density = tempdensity;
              ts_temp = ts;
              te_temp = te;
            }
            // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;
          }
          te_mark = te;
        }

        // cout<<te_mark<<endl;
        FlagOfFirst = false;
      }
      else
      { // incremental computation for ts
        // ############## corner case  ######################  ts_next > te_current

        // corner case: te_max < ts
        if (Tgraph_ptr->_InVNT[this->sink][sink_num - 1] < ts)
        {
          break;
        }

        te = te_mark;
        if (te < ts + this->delta)
          te = ts + this->delta;
        // call once of TransformationMAPE because there is a <<break>> during FlagOfFirst = true;
        // cout<<"hello0"<<endl;
        // cout<<"MAPSE: "<<ts_previous<<"   "<<ts<<"   "<<te_mark_previous<<"    "<<te_mark<<endl;

        trans_start = clock();

        // For different Dinic
        // cout<<"ts1: "<< ts_previous<<" "<<ts<<endl;
        // cout<<"te1: "<< te_mark_previous<<" "<<te<<endl;

        if (te_mark_previous >= ts)
        {
          TransformationStart = clock();
          tempstart = clock();
          flowvaluetest = TransformationMAPSE(this->Tgraph_ptr, Graph_ptr, ts_previous, ts, te);
          TransformationEnd = clock();
          transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
          if (flowvaluetest > 0)
          {
            if (DinicType == 0)
            {
              tempstart = clock();
              tempvalue = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->sink][te_mark_previous], te_mark_previous, Graph_ptr->_TVtoV[this->source][-1], -1);
              tempend = clock();
              
              // cout << "tempflow-: " << tempvalue << "at [" << te << ", -1]" << endl;
              flow_mark -= tempvalue;
            }
            if (DinicType == 1)
            {
              tempvalue = ReverseDinic(Graph_ptr, Graph_ptr->_TVtoV[this->sink][te_mark_previous], te_mark_previous, Graph_ptr->_TVtoV[this->source][-1], -1);
              // cout<<"tempflow-: "<<tempvalue<<endl;
              flow_mark -= tempvalue;
            }
          }
          TransformationStart = clock();
          TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_mark_previous, te, true);
          TransformationEnd = clock();
          transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
        }
        else
        { // corner case
          TransformationStart = clock();
          TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_mark_previous, te, false);
          TransformationMAPSE(this->Tgraph_ptr, Graph_ptr, ts_previous, ts, te);
          TransformationEnd = clock();
          transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
          flow_mark = 0;
        }


        if(FlagOfPruning)
          tempstart = clock();

        trans_end = clock();        
        tempvalue = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
        tempend = clock();
        if(type == 4){
          OutFile << Graph_ptr->Vcnt<<" "<< (double)(tempend - tempstart) / CLOCKS_PER_SEC<<endl;
        }
        decase++;
        // cout << "tempflow+: " << tempvalue << "at [" << ts << ", " << te << "]" << endl;
        flow_mark += tempvalue;
        tempflow = flow_mark;
        Graph_ptr->maxflow = tempflow;
        tempdensity = tempflow / (te - ts);
        if ((ts + this->delta) > te)
        {
          tempdensity = tempflow / this->delta;
        }
        // if (tempdensity != 0)
        //   cout << "MaxFlow for (<" << this->source << ", " << ts << ">, <" << this->sink << ", " << te << ">) is " << tempdensity << " !" << endl;
        //  cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
        if (tempdensity > density)
        {
          density = tempdensity;
          ts_temp = ts;
          te_temp = te;
          if ((ts + this->delta) > te)
          {
            ts_temp = te - this->delta;
          }
        }
        transtime += (double)(trans_end - trans_start) / CLOCKS_PER_SEC;
        countnum++;
        avg_edge = avg_edge + Graph_ptr->getEcnt();
        avg_vertex = avg_vertex + Graph_ptr->getVcnt();
        // cout<< avg_edge <<" "<<avg_vertex<<endl;

        lastflow = tempflow;
        te_last = te;
        te_previous = te;

        /*

        if ((ts + this->delta) > te)
        {
          te = ts + this->delta;

          if (te > this->Tgraph_ptr->Tmax) // corner case: [ts, te_mark] use MAPE to [ts, ts+delta]
            te = this->Tgraph_ptr->Tmax;

          TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_previous, te);
          tempvalue = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
          tempflow += tempvalue;
          Graph_ptr->maxflow = tempflow;

          countnum++;
          avg_edge = avg_edge + Graph_ptr->getEcnt();
          avg_vertex = avg_vertex + Graph_ptr->getVcnt();

          lastflow = tempflow;
          flow_mark = tempflow;
          te_last = te;
          tempdensity = tempflow / this->delta;
          // if(tempdensity != 0)
          // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
          // if(tempdensity != 0)
          // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
          if (tempdensity > density)
          {
            density = tempdensity;
            te_temp = te;
            ts_temp = te - this->delta;
          }
        }*/

        // Below are incremental Dinic
        flag = false; // record whether copy is needed
        count = -1;   // find the minimum timestamp of t
        for (int i = 0; i < sink_num; i++)
        {
          if (Tgraph_ptr->_InVNT[this->sink][i] > te)
          {
            count = i;
            break;
          }
        }

        te_previous = te;
        te_mark_previous = te;
        ts_previous = ts;
        te_mark = te;

        if (count == -1)
        {
          continue;
        }

        // cout<<"hihi"<<endl;
        if ((j + 1) != source_num)
        {
          for (; count < sink_num; count++)
          {

            te = Tgraph_ptr->_InVNT[this->sink][count];
            if (te >= (Tgraph_ptr->_OutVNT[this->source][j + 1] + this->delta))
            { // find the timestamp of te that is just larger than the next ts + delta
              te_mark = te;
              flag = true;
              break;
            }
            // cout<<"case 2:"<<ts<<"~"<<te<<endl;
            // cout<<"hello2"<<endl;
            TransformationStart = clock();
            TransformationMAPE(this->Tgraph_ptr, Graph_ptr, te_previous, te, true);
            TransformationEnd = clock();
            transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;
            te_previous = te;
            te_mark_previous = te;
            // cout<< Temp_ptr->Vcnt <<" "<<Graph_ptr->Vcnt<<endl;
            // cout<<"Max for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;

            // Observation 2
            if (Observation2(this->Tgraph_ptr, ts, te_last, te, lastflow, density)&&prune)
            {
              continue;
            }

            tempstart = clock();
            tempvalue = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
            tempend = clock();
 
            tempflow += tempvalue;
            Graph_ptr->maxflow = tempflow;

            countnum++;
            avg_edge = avg_edge + Graph_ptr->getEcnt();
            avg_vertex = avg_vertex + Graph_ptr->getVcnt();

            lastflow = tempflow;
            flow_mark = tempflow;
            te_last = te;
            tempdensity = tempflow / (te - ts);
            // if(tempdensity != 0)
            // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
            if (tempdensity > density)
            {
              density = tempdensity;
              ts_temp = ts;
              te_temp = te;
            }
            // cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te)<<" !"<<endl;
            te_mark = te;
          }
        }
      }

      // Below continue incremental Dinic
      if (flag)
      { // the case of network copy
        copystart = clock();
        Temp_ptr->reset();
        Temp_ptr->copy(Graph_ptr);
        copyend = clock();
        copytime += (double)(copyend - copystart) / CLOCKS_PER_SEC;
        Timeofcopy++;
      }
      else
      {
        Temp_ptr = Graph_ptr;
      }

      // Observation 1
      for (; count < sink_num; count++)
      {
        te = Tgraph_ptr->_InVNT[this->sink][count];
        // cout<<"case 3:"<<ts<<"~"<<te_previous<<"~"<<te<<endl;
        //  cout<<"hello3"<<endl;
        TransformationStart = clock();
        TransformationMAPE(this->Tgraph_ptr, Temp_ptr, te_previous, te, true);
        TransformationEnd = clock();
        transTime += (double)(TransformationEnd - TransformationStart) / CLOCKS_PER_SEC;

        // cout<< Temp_ptr->Vcnt <<" "<<Graph_ptr->Vcnt<<endl;
        // cout<<"Max for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;

        // Observation 2
        if (Observation2(this->Tgraph_ptr, ts, te_last, te, lastflow, density))
        {
          te_previous = te;
          continue;
        }

        tempstart = clock();
        tempvalue = Dinic(Temp_ptr, Temp_ptr->_TVtoV[this->source][ts], ts, Temp_ptr->_TVtoV[this->sink][te], te);
        tempend = clock();

        tempflow += tempvalue;
        Temp_ptr->maxflow = tempflow;

        countnum++;
        avg_edge = avg_edge + Temp_ptr->getEcnt();
        avg_vertex = avg_vertex + Temp_ptr->getVcnt();

        lastflow = tempflow;
        te_last = te;
        tempdensity = tempflow / (te - ts);
        // if(tempdensity != 0)
        //  cout<<"MaxFlow for (<" << this->source <<", " << ts << ">, <" << this->sink<<", " <<te <<">) is "<< tempdensity <<" !"<<endl;
        if (tempdensity > density)
        {
          density = tempdensity;
          ts_temp = ts;
          te_temp = te;
        }

        te_previous = te;
      }
    }
    endTime = clock();
    if (Temp_ptr != Graph_ptr)
      delete Temp_ptr;
    delete Graph_ptr;
    _avg_vertex += avg_vertex / countnum;
    _avg_edge += avg_edge / countnum;
    _countnum += countnum;
    _runtime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
    _density += density;

    // cout << "The runtime of TransformationMAPSE is " << transtime << " seconds!" << endl;
    // cout << "The runtime of latter dinic is " << dinictime << " seconds!" << endl;
    // cout << Timeofcopy << " copy take " << copytime << " seconds!" << endl;
  }
  if (DinicType == 0)
  {
    cout << "The runtime of MAPSE is " << _runtime / REPEAT << " seconds!" << endl;
    cout << "The density of BFlow is " << _density / REPEAT << "during [" << ts_temp << "," << te_temp << "] !" << endl;
    cout << "Avg_Edge: " << _avg_edge / REPEAT << "   Avg_Vertex: " << _avg_vertex / REPEAT << "   # of Dinic: " << _countnum / REPEAT << endl;
  }

  if (DinicType == 1)
  {
    cout << "The runtime of MAPSE* is " << _runtime / REPEAT << " seconds!" << endl;
    cout << "The density of BFlow is " << _density / REPEAT << "during [" << ts_temp << "," << te_temp << "] !" << endl;
    cout << "Avg_Edge: " << _avg_edge / REPEAT << "   Avg_Vertex: " << _avg_vertex / REPEAT << "   # of Dinic: " << _countnum / REPEAT << endl;
  }

  // Output overall results
  if (type == 1)
  {
    if (DinicType == 0)
      OutFile << _runtime / REPEAT << endl;
    if (DinicType == 1)
      OutFile << _runtime / REPEAT << endl;
  }
  if (type == 2)
  {
    ALG3[num] = _runtime / REPEAT;
    DeCase[num] = decase / REPEAT;
  }
  if (type == 3)
  {
    if (DinicType == 0)
      OutFile << _runtime / REPEAT << "  " << transTime / REPEAT << endl;
    if (DinicType == 1)
      OutFile << _runtime / REPEAT << "  " << transTime / REPEAT << endl;
  }

}

template <class WeightType, class TimestampType>
double BFlow<WeightType, TimestampType>::TransformationMAPSE(TGRAPH<WeightType, TimestampType> *T_ptr, DIGRAPH<WeightType, TimestampType> *ptr, TimestampType start, TimestampType end, TimestampType t_end)
{
  double flowvalue = 0;
  // duplicate source and sink
  // int VCnt = ptr->getVcnt(); //VertexID for vertices in transfomred networks
  int count, i;
  VertexID ps, next, NextV, CurrentV;
  WeightType tempW;
  TimestampType temptime, mark, CurrentT, NextT;

  // virtual node <source, -1>
  if (ptr->_TVtoV[this->source].find(-1) != ptr->_TVtoV[this->source].end())
  {
    ptr->eraseVertex(ptr->_TVtoV[this->source][-1], this->source, -1);
  }
  ptr->insertVertex(this->source, -1);

  /*
  if (ptr->_TVtoV[this->source].find(-1) == ptr->_TVtoV[this->source].end())
  {
    ptr->insertVertex(this->source, -1);
    // VCnt++;
  }
  else
  {
    ptr->removeAllInEdges(ptr->_TVtoV[this->source][-1]);
  }*/

  // for source
  CurrentV = ptr->_TVtoV[this->source][start];
  CurrentT = start;

  /*
  if(ptr->_tlabels[this->source] < start){
    cout<<"hi1"<<endl;
    CurrentV = ptr->_TVtoV[this->source][ptr->_tlabels[this->source]];
    CurrentT = ptr->_tlabels[this->source];
  }*/

  // cout<<"start: "<<start<<"       end: "<<end<<"      t_end: "<<t_end<<endl;
  while (ptr->nextTime[CurrentV] < end)
  {
    // cout<<CurrentV << "   "<<CurrentT<<endl;
    NextV = ptr->nextVertex[CurrentV];
    NextT = ptr->nextTime[CurrentV];
    // cout<<"erase the source at "<< CurrentT<<endl;
    ptr->eraseVertex(CurrentV, this->source, CurrentT);
    CurrentV = NextV;
    CurrentT = NextT;
  }

  // cout <<"hello1"<<endl;
  // NextT >= end
  NextT = ptr->nextTime[CurrentV];
  if (NextT > end)
  {
    NextV = ptr->nextVertex[CurrentV];
    ptr->insertVertex(this->source, end);
    ptr->insertEdge(ptr->_TVtoV[this->source][end], NextV, WEIGHT_MAX);
    ptr->getOutEdge()[ptr->_TVtoV[this->source][end]][NextV].flow = ptr->getOutEdge()[CurrentV][NextV].flow;
    ptr->getOutEdge()[NextV][ptr->_TVtoV[this->source][end]].flow = -ptr->getOutEdge()[CurrentV][NextV].flow;
    ptr->nextVertex[ptr->_TVtoV[this->source][end]] = NextV;
    ptr->nextTime[ptr->_TVtoV[this->source][end]] = NextT;
    // cout<<"erase the source at "<< CurrentT<<endl;
    ptr->eraseVertex(CurrentV, this->source, CurrentT);
  }
  else
  {
    ptr->eraseVertex(CurrentV, this->source, CurrentT);
  }
  ptr->_tlabels[this->source] = end;

  // cout <<"hello2"<<endl;

  // for sink
  CurrentV = ptr->_TVtoV[this->sink][start];
  CurrentT = start;

  /*
  if(ptr->_tlabels[this->sink] < start){
    cout<<"hi2"<<endl;
    CurrentV = ptr->_TVtoV[this->sink][ptr->_tlabels[this->sink]];
    CurrentT = ptr->_tlabels[this->sink];
  }*/

  while (ptr->nextTime[CurrentV] < end)
  {
    NextV = ptr->nextVertex[CurrentV];
    NextT = ptr->nextTime[CurrentV];
    // cout<<"erase the sink at "<< CurrentT<<endl;
    ptr->eraseVertex(CurrentV, this->sink, CurrentT);
    CurrentV = NextV;
    CurrentT = NextT;
  }

  NextT = ptr->nextTime[CurrentV];
  if (NextT > end)
  {
    NextV = ptr->nextVertex[CurrentV];
    ptr->insertVertex(this->sink, end);
    ptr->insertEdge(ptr->_TVtoV[this->sink][end], NextV, WEIGHT_MAX);
    ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][NextV].flow = ptr->getOutEdge()[CurrentV][NextV].flow;
    ptr->getOutEdge()[NextV][ptr->_TVtoV[this->sink][end]].flow = ptr->getOutEdge()[CurrentV][NextV].flow;
    ptr->nextVertex[ptr->_TVtoV[this->sink][end]] = NextV;
    ptr->nextTime[ptr->_TVtoV[this->sink][end]] = NextT;
    ptr->insertEdge(ptr->_TVtoV[this->sink][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][NextV].flow);
    flowvalue += ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][NextV].flow;
    // cout<<"erase the sink at "<< CurrentT<<endl;
    // cout<<"flow of sink from "<< CurrentT<<" to " <<NextT<<" is "<<ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][NextV].flow<<endl;
    ptr->eraseVertex(CurrentV, this->sink, CurrentT);
  }
  else
  {
    NextV = ptr->nextVertex[CurrentV];
    ptr->insertEdge(NextV, ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[CurrentV][NextV].flow);
    flowvalue += ptr->getOutEdge()[CurrentV][NextV].flow;
    // cout<<"erase (else) the sink at "<< CurrentT<<endl;
    // cout<<"flow of sink from "<< CurrentT<<" to " <<NextT<<" is "<<ptr->getOutEdge()[CurrentV][NextV].flow<<endl;
    ptr->eraseVertex(CurrentV, this->sink, CurrentT);
  }
  ptr->_tlabels[this->sink] = end;

  for (auto it = ptr->_tlabels.begin(); it != ptr->_tlabels.end(); it++)
  {
    ps = it->first;
    if ((ps == this->source) || (ps == this->sink))
      continue;

    if (ptr->_tlabels[ps] > end)
      continue;

    // special case: ps is inserted into _tlabels during TransformationMAPE, add one vertex for ps at start
    if (ptr->_tlabels[ps] > start)
    {
      ptr->insertVertex(ps, start);
      // cout <<"hello1"<<endl;
      ptr->insertEdge(ptr->_TVtoV[ps][start], ptr->_TVtoV[ps][ptr->_tlabels[ps]], WEIGHT_MAX);
      ptr->nextVertex[ptr->_TVtoV[ps][start]] = ptr->_TVtoV[ps][ptr->_tlabels[ps]];
      ptr->nextTime[ptr->_TVtoV[ps][start]] = ptr->_tlabels[ps];
      ptr->_tlabels[ps] = start;
      continue;
    }

    // ptr->_tlabels[ps] <= start, start to delete vertices
    CurrentT = ptr->_tlabels[ps];
    CurrentV = ptr->_TVtoV[ps][CurrentT];

    while (ptr->nextTime[CurrentV] < end)
    {
      NextV = ptr->nextVertex[CurrentV];
      NextT = ptr->nextTime[CurrentV];
      // cout<<"erase "<<ps<<" at "<< CurrentT<<":   "<<CurrentV<<endl;
      ptr->eraseVertex(CurrentV, ps, CurrentT);
      ptr->_tlabels[ps] = NextT;
      CurrentV = NextV;
      CurrentT = NextT;
    }

    NextT = ptr->nextTime[CurrentV];
    if (NextT > end)
    {
      NextV = ptr->nextVertex[CurrentV];
      ptr->insertVertex(ps, end);
      // cout <<"hello2"<<endl;
      ptr->insertEdge(ptr->_TVtoV[ps][end], NextV, WEIGHT_MAX);
      ptr->getOutEdge()[ptr->_TVtoV[ps][end]][NextV].flow = ptr->getOutEdge()[CurrentV][NextV].flow;
      ptr->getOutEdge()[NextV][ptr->_TVtoV[ps][end]].flow = -ptr->getOutEdge()[CurrentV][NextV].flow;
      ptr->nextVertex[ptr->_TVtoV[ps][end]] = NextV;
      ptr->nextTime[ptr->_TVtoV[ps][end]] = NextT;
      // cout <<"hello3"<<endl;
      ptr->insertEdge(ptr->_TVtoV[ps][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[ps][end]][NextV].flow);
      flowvalue += ptr->getOutEdge()[ptr->_TVtoV[ps][end]][NextV].flow;
      // cout<<"erase "<<ps<<" at "<< CurrentT<<endl;
      ptr->eraseVertex(CurrentV, ps, CurrentT);
      // ptr->_tlabels[ps] =NextT;
      ptr->_tlabels[ps] = end;
    }
    else
    {
      NextV = ptr->nextVertex[CurrentV];
      // cout <<"hello4"<<endl;
      ptr->insertEdge(NextV, ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[CurrentV][NextV].flow);
      flowvalue += ptr->getOutEdge()[CurrentV][NextV].flow;
      // cout<<"erase (else) "<<ps<<" at "<< CurrentT<<endl;
      ptr->eraseVertex(CurrentV, ps, CurrentT);
      ptr->_tlabels[ps] = end;
    }

    //  ### here
    /*
    if (ptr->_TVtoV[ps].size() == 2){
      if((ptr->getOutEdge()[ptr->_TVtoV[ps][start]].size()==1)&&(ptr->getOutEdge()[ptr->_TVtoV[ps][end]].size()==1)){
        ptr->eraseVertex(ptr->_TVtoV[ps][start], ps, start);
        ptr->eraseVertex(ptr->_TVtoV[ps][end], ps, end);
      ptr->_tlabels.erase(ps);
      }
    }*/
  }

  return flowvalue;

  /*
  //for source
  if (ptr->_TVtoV[this->source].find(end) == ptr->_TVtoV[this->source].end())
  {

    next = ptr->nextVertex[ptr->_TVtoV[this->source][start]];

    ptr->insertVertex(this->source, end);
    ptr->insertEdge(ptr->_TVtoV[this->source][end], next, WEIGHT_MAX);
    ptr->getOutEdge()[ptr->_TVtoV[this->source][end]][next].flow = ptr->getOutEdge()[ptr->_TVtoV[this->source][start]][next].flow;
    ptr->removeEdge(ptr->_TVtoV[this->source][start], next, true);
    ptr->nextVertex[ptr->_TVtoV[this->source][end]] = next;
    ptr->nextVertex.erase(ptr->_TVtoV[this->source][start]);
    ptr->eraseVertex(ptr->_TVtoV[this->source][start], this->source, start);
  }
  else
  {
    ptr->removeEdge(ptr->_TVtoV[this->source][start], ptr->_TVtoV[this->source][end], true);
    ptr->eraseVertex(ptr->_TVtoV[this->source][start], this->source, start);
    ptr->nextVertex.erase(ptr->_TVtoV[this->source][start]);
    //}
  }


  // for sink
  if (ptr->_TVtoV[this->sink].find(end) == ptr->_TVtoV[this->sink].end())
  {

    next = ptr->nextVertex[ptr->_TVtoV[this->sink][start]];

    ptr->insertVertex(this->sink, end);
    ptr->insertEdge(ptr->_TVtoV[this->sink][end], next, WEIGHT_MAX);
    ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][next].flow = ptr->getOutEdge()[ptr->_TVtoV[this->sink][start]][next].flow;
    ptr->insertEdge(ptr->_TVtoV[this->sink][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[this->sink][end]][next].flow);
    ptr->removeEdge(ptr->_TVtoV[this->sink][start], next, true);
    ptr->nextVertex[ptr->_TVtoV[this->sink][end]] = next;
    ptr->nextVertex.erase(ptr->_TVtoV[this->sink][start]);
    ptr->eraseVertex(ptr->_TVtoV[this->sink][start], this->sink, start);
  }
  else
  {
    ptr->insertEdge(ptr->_TVtoV[this->sink][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[this->sink][start]][ptr->_TVtoV[this->sink][end]].flow);
    ptr->removeEdge(ptr->_TVtoV[this->sink][start], ptr->_TVtoV[this->sink][end], true);
    ptr->nextVertex.erase(ptr->_TVtoV[this->sink][start]);
    ptr->eraseVertex(ptr->_TVtoV[this->sink][start], this->sink, start);
  }
  */

  /*
   for (auto it = ptr->_tlabels.begin(); it != ptr->_tlabels.end(); it++)
   {
     ps = it->first;
     if ((ps == this->source) || (ps == this->sink))
       continue;

     // special case: ps is inserted into _tlabels during TransformationMAPE
     if (ptr->_tlabels[ps] > start)
       continue;

     if (ptr->_TVtoV[ps].size() == 0)
     {
       ptr->_tlabels.erase(ps);
       continue;
     }

     if (ptr->_TVtoV[ps].find(end) == ptr->_TVtoV[ps].end())
     {

       next = ptr->nextVertex[ptr->_TVtoV[ps][start]];

       ptr->insertVertex(ps, end);
       ptr->insertEdge(ptr->_TVtoV[ps][end], next, WEIGHT_MAX);
       ptr->getOutEdge()[ptr->_TVtoV[ps][end]][next].flow = ptr->getOutEdge()[ptr->_TVtoV[ps][start]][next].flow;
       ptr->insertEdge(ptr->_TVtoV[ps][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[ps][end]][next].flow);
       ptr->removeEdge(ptr->_TVtoV[ps][start], next, true);
       ptr->nextVertex[ptr->_TVtoV[ps][end]] = next;
       ptr->nextVertex.erase(ptr->_TVtoV[ps][start]);
       ptr->eraseVertex(ptr->_TVtoV[ps][start], ps, start);
     }
     else
     {
       ptr->insertEdge(ptr->_TVtoV[ps][end], ptr->_TVtoV[this->source][-1], ptr->getOutEdge()[ptr->_TVtoV[ps][start]][ptr->_TVtoV[ps][end]].flow);
       ptr->removeEdge(ptr->_TVtoV[ps][start], ptr->_TVtoV[ps][end], true);
       ptr->nextVertex.erase(ptr->_TVtoV[ps][start]);
       ptr->eraseVertex(ptr->_TVtoV[ps][start], ps, start);
     }
   } */
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::QueryGenerator(int num, int length)
{
  int count = 0, threshold = 50;
  bool flag = false;
  clock_t currentTime;
  string OutFileName = "query/" + OutputName + ".qry";
  ofstream OutFile(OutFileName);
  cout << OutFileName << endl;
  unordered_map<VertexID, unordered_map<TimestampType, bool>> flagarray;
  // srand((int)time(0));
  srand(1);
  while (1)
  {
    this->source = this->Tgraph_ptr->getRV();
    // for(auto it2 = this->Tgraph_ptr->getVLabel().begin(); it2 != this->Tgraph_ptr->getVLabel().end(); it2++){
    // this->sink = it2->first;
    // continue;
    // for(int i = 0; i < threshold; i++){
    this->sink = this->Tgraph_ptr->getRV();
    if (this->source == this->sink)
      continue;
    cout << this->source << "  " << this->sink << ":" << endl;
    // while(this->source == this->sink){
    //   this->sink = this->Tgraph_ptr->getRV();
    // }

    // this->sink = 79045;
    // this->source = 12398;
    currentTime = clock();
    flagarray.clear();
    if (Reach(0, -1, 0, 0, currentTime, flagarray))
    {
      if (this->QueryGeneratorLength >= length)
      {
        cout << "Reachable from " << this->source << " to " << this->sink << " !" << endl;
        if (flag)
        {
          OutFile << endl;
          OutFile << this->source << "  " << this->sink;
        }
        else
        {
          OutFile << this->source << "  " << this->sink;
          flag = true;
        }
        count++;

        if (count == num)
        {
          OutFile.close();
          return;
        }

        // break;
      }
    }
    // return;
    //}
    //}
    // return;
  }
}

template <class WeightType, class TimestampType>
bool BFlow<WeightType, TimestampType>::Reach(VertexID current_node, TimestampType current_time, TimestampType maxT, int length, clock_t startTime, unordered_map<VertexID, unordered_map<TimestampType, bool>> &flag)
{
  VertexID tempV;
  TimestampType tempT;

  if (current_time == -1)
  {
    // if((this->source !=12398)||(this->sink!=79045))
    // return false;
    // cout<<"hello4"<<endl;
    // this->sink = 79045;
    // this->source = 12398;

    clock_t tempTime;
    TimestampType max = 0;
    for (auto it = this->Tgraph_ptr->getInVertex()[this->sink].begin(); it != this->Tgraph_ptr->getInVertex()[this->sink].end(); it++)
    {
      tempV = it->first;
      if (tempV == this->source)
      {
        this->QueryGeneratorLength = 1;
        return true;
      }
      for (auto it1 = this->Tgraph_ptr->getOutEdge()[tempV][this->sink].begin(); it1 != this->Tgraph_ptr->getOutEdge()[tempV][this->sink].end(); it1++)
      {
        tempT = it1->first;
        if (tempT > max)
          max = tempT;
      }
    }

    for (auto it = this->Tgraph_ptr->getOutEdge()[this->source].begin(); it != this->Tgraph_ptr->getOutEdge()[this->source].end(); it++)
    {
      tempV = it->first;
      if (tempV == this->sink)
      {
        this->QueryGeneratorLength = 1;
        return true;
      }
      for (auto it1 = this->Tgraph_ptr->getOutEdge()[this->source][tempV].begin(); it1 != this->Tgraph_ptr->getOutEdge()[this->source][tempV].end(); it1++)
      {
        tempT = it1->first;
        if (tempT > max)
          continue;

        flag[tempV][tempT] = true;
        if (Reach(tempV, tempT, max, 1, startTime, flag))
        {
          // cout<<tempV<<", "<<tempT<<endl;
          return true;
        }
        flag[tempV].erase(tempT);
        // control the time
        tempTime = clock();
        if (((double)(tempTime - startTime) / CLOCKS_PER_SEC) > 5)
          return false;
      }
    }
  }
  else
  {
    clock_t tempTime;
    // control the time
    tempTime = clock();
    if (((double)(tempTime - startTime) / CLOCKS_PER_SEC) > 5)
      return false;
    for (auto it = this->Tgraph_ptr->getOutEdge()[current_node].begin(); it != this->Tgraph_ptr->getOutEdge()[current_node].end(); it++)
    {
      tempV = it->first;

      if (tempV == this->source)
        continue;

      // cout<<"V: "<<tempV<<endl;
      for (auto it1 = this->Tgraph_ptr->getOutEdge()[current_node][tempV].begin(); it1 != this->Tgraph_ptr->getOutEdge()[current_node][tempV].end(); it1++)
      {
        tempT = it1->first;
        // cout<<"T: "<<tempT<<endl;
        if (tempT < current_time)
          continue;
        if (tempT > maxT)
          continue;
        if (tempV == this->sink)
        {
          // cout<<current_node<<", "<<current_time<<endl;
          this->QueryGeneratorLength = length;
          return true;
        }

        if (flag[tempV].find(tempT) != flag[tempV].end())
          return false;

        if (Reach(tempV, tempT, maxT, length + 1, startTime, flag))
        {
          // cout<<current_node<<", "<<current_time<<endl;
          return true;
        }
      }
    }
  }
  return false;
}

template <class WeightType, class TimestampType>
void BFlow<WeightType, TimestampType>::TestFlow(ofstream &OutFile, TimestampType t1, TimestampType t2, int num, int type)
{
  double avg_vertex, avg_edge;
  double _runtime = 0, _density = 0, _avg_vertex = 0, _avg_edge = 0;
  int countnum, _countnum = 0;
  clock_t startTime, endTime, tempstart, tempend;

  DIGRAPH<WeightType, TimestampType> *Graph_ptr = new DIGRAPH<WeightType, TimestampType>;
  TimestampType ts, te, ts_temp, te_temp, te_previous, te_last;
  WeightType density, tempdensity, tempflow, lastflow;

  int source_num = Tgraph_ptr->_OutVNT[this->source][-1];
  if (source_num == 0)
  {
    cout << "The source has no out-going edges!" << endl;
    return;
  }

  int sink_num = Tgraph_ptr->_InVNT[this->sink][-1];
  if (sink_num == 0)
  {
    cout << "The sink has no in-coming edges!" << endl;
    return;
  }

  // return;

  // Line 4 of Base

  TimestampType bound1 = Tgraph_ptr->_OutVNT[this->source][0];
  if (bound1 < 0)
    bound1 = 0;
  TimestampType bound2 = Tgraph_ptr->_OutVNT[this->source][Tgraph_ptr->_OutVNT[this->source][-1] - 1];
  // Line 5 of Base
  TimestampType bound3 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1] + this->delta;
  TimestampType bound4 = Tgraph_ptr->_InVNT[this->sink][sink_num - 1];

  Graph_ptr->reset();
  startTime = clock();
  countnum = 0;
  avg_vertex = 0, avg_edge = 0;
  density = -1;

  ts = t1;
  te = t2;

  Transformation(Tgraph_ptr, Graph_ptr, ts, te);
  // cout<<Graph_ptr->getEcnt()<<" "<<Graph_ptr->getVcnt()<<endl;
  tempstart = clock();
  tempflow = Dinic(Graph_ptr, Graph_ptr->_TVtoV[this->source][ts], ts, Graph_ptr->_TVtoV[this->sink][te], te);
  tempend = clock();
  countnum++;
  avg_edge = avg_edge + Graph_ptr->getEcnt();
  avg_vertex = avg_vertex + Graph_ptr->getVcnt();

  lastflow = tempflow; // used for Observation2
  te_last = te;        // used for Observation2
  tempdensity = tempflow / (te - ts);
  if (tempdensity > density)
  {
    density = tempdensity;
    ts_temp = ts;
    te_temp = te;
  }

  endTime = clock();
  _avg_vertex += avg_vertex / countnum;
  _avg_edge += avg_edge / countnum;
  _countnum += countnum;
  _runtime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
  _density += density;

  cout << "The runtime of Test is " << _runtime / 1 << " seconds!" << endl;
  cout << "The density of BFlow is " << _density / 1 << "during [" << ts_temp << "," << te_temp << "] !" << endl;
  cout << "Avg_Edge: " << _avg_edge / 1 << "   Avg_Vertex: " << _avg_vertex / 1 << "   # of Dinic: " << _countnum / 1 << endl;

  delete Graph_ptr;
}