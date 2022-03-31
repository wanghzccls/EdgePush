#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
//#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include "priority.h"
#include <unordered_map>
#include <unordered_set>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <queue>
#include <cmath>
#include <random>
#include <ctime>
//#include <set>
#include <queue>


typedef unsigned int uint;


class SimStruct{
	public:
		Graph g;//Class Graph
		Random R;//Class Random
		uint vert;//# of vertice
		uint nedge;//# of edges
		string filelabel;
		double alpha;
		double eps;
		double avg_time;
		uint seed;
		uint *H[2];
		uint *U[2];
		uint* C[2];//finalreserve
		double *finalReserve;
		//uint residue_count;
		uint finalreserve_count;
		
		SimStruct(string name, string file_label, double epsilon, double para_alpha) {
			filelabel = file_label;
			g.inputGraph(name,file_label,epsilon);
			R = Random();
			vert = g.n;
			nedge = g.m;
			alpha=para_alpha;
			eps = epsilon;
			avg_time = 0;
			seed=(uint)time(0);
			finalreserve_count=0;
		}

		~SimStruct() {
		
		}      
  
		virtual void query(uint u){ }  
};



class EdgePushstruct:public SimStruct{
	public:
	Priority PQ;
	double avg_ppr_time;
	uint candidate_start;
	uint candidate_end;
	uint *H[2];//mark the node expense, vert
	uint *U[2];//mark the node income(residue), vert
	uint* HC[2];//NodeResidue
	uint* UC[2];//resetPriority
	double *NodeResidue;
	double *residue;//node income
	double *expense;
	bool* isInArray;
	uint *candidate_set;
	uint NodeResidue_count;
	uint residue_count;
	uint expense_count;
	uint resetPQ_count;
	uint tmpind;
	bool globalflag;
		
		
	EdgePushstruct(string name, string file_label, double epsilon, double para_alpha):
	SimStruct(name, file_label, epsilon, para_alpha) {
		PQ.initPriority(vert,g.initsort,g.outPL);
		avg_ppr_time=0;
		U[0] = new uint[vert];
		U[1] = new uint[vert];
		H[0] = new uint[nedge]();
		H[1] = new uint[nedge]();
		C[0] = new uint[vert];
		C[1] = new uint[vert];
		UC[0] = new uint[vert];
		UC[1] = new uint[vert];
		HC[0] = new uint[vert];
		HC[1] = new uint[vert];
		NodeResidue = new double[vert];
		residue = new double[vert];
		expense = new double[nedge]();//reset as all zero
		finalReserve = new double[vert];
		isInArray = new bool[vert];
		candidate_set=new uint[vert+1];
		finalreserve_count=0;
		NodeResidue_count=0;
		residue_count=0;
		expense_count=0;
		resetPQ_count=0;
		candidate_start=0;
		candidate_end=0;
		globalflag=false;

		for(uint i = 0; i < vert; i++){
			isInArray[i] = false;
			NodeResidue[i] = -1;
			residue[i]=0;
			finalReserve[i]=0;
			U[0][i] = 0;
			U[1][i] = 0;
			C[0][i] = 0;
			C[1][i] = 0;
			UC[0][i] = 0;
			UC[1][i] = 0;
			HC[0][i] = 0;
			HC[1][i] = 0;
		}
	}

	~EdgePushstruct() {
		delete[] H[0];
		delete[] H[1];
		delete[] U[0];
		delete[] U[1];
		delete[] C[0];
		delete[] C[1];
		delete[] UC[0];
		delete[] UC[1];
		delete[] HC[0];
		delete[] HC[1];
		delete[] NodeResidue;
		delete[] residue;
		delete[] expense;
		delete[] finalReserve;
		delete[] isInArray;
		delete[] candidate_set;
	}      
	
	uint count_candidate_size(uint start,uint end){
		if(end>=start){
			return (end-start);
		}
		else{
			return (end+vert+1-start);
		}
	}
		
	void query(uint u){
		//cout<<"alpha="<<alpha<<endl;
		for(uint j=0;j<NodeResidue_count;j++){
			NodeResidue[HC[0][j]]=-1;
			HC[1][HC[0][j]]=0;
		}
		NodeResidue_count=0;
		
		for(uint j=0;j<finalreserve_count;j++){	
			finalReserve[C[0][j]]=0;
			C[1][C[0][j]]=0;
		}
		finalreserve_count=0;

		for(uint j=0;j<residue_count;j++){
			residue[U[0][j]]=0;
			U[1][U[0][j]]=0;
		}
			residue_count=0;
	
		for(uint j=0;j<expense_count;j++){
			expense[H[0][j]]=0;
			H[1][H[0][j]]=0;
		}
		expense_count=0;	
	
		for(uint j=0;j<resetPQ_count;j++){
			PQ.resetPriority(g.initsort[UC[0][j]],UC[0][j],g.getOutSize(UC[0][j]));
			UC[1][UC[0][j]]=0;
		}
		resetPQ_count=0;

		for(uint j=candidate_start;j<candidate_end;j++){
			isInArray[candidate_set[j]]=false;
		}
		candidate_end=candidate_start;

		globalflag=false;
		uint tempLevel = 0;

		residue[u]=1;
		U[0][residue_count++]=u;
		U[1][u]=1;

		finalReserve[u]=alpha;
		C[0][finalreserve_count++]=u;
		C[1][u]=1;	

		double Gku=PQ.findMin(u)-(1-alpha)/g.getOutVertWeight(u);
		if(Gku<0){
			candidate_set[candidate_end]=u;
			candidate_end=(candidate_end+1)%(vert+1);
			isInArray[u]=true;
		}
		uint global_switchSize=vert+1;

		while(candidate_end!=candidate_start){
			uint tempi=candidate_set[candidate_start];
			candidate_start=(candidate_start+1)%(vert+1);
			isInArray[tempi]=false;
			double outVertWeight=g.getOutVertWeight(tempi);

			if(NodeResidue[tempi]<0){
				uint tempi_OutSize=g.getOutSize(tempi);
				uint tempi_switchSize=1;
				for(uint itej=0;itej<tempi_OutSize;itej++){
					if(itej>=tempi_switchSize){
						NodeResidue[tempi]=0;	
						if(HC[1][tempi]==0){
							HC[1][tempi]=1;
							HC[0][NodeResidue_count++]=tempi;
						}
						uint start_pos=g.outPL[tempi];
						uint end_pos=g.outPL[tempi+1];
						for(uint itek=start_pos;itek<end_pos;itek++){
							uint newNode=g.outEL[itek];
							double newEdgeWeight=g.outWEL[itek];
							double incre_newNode=(1-alpha)*residue[tempi]*newEdgeWeight/outVertWeight-expense[itek];
							if(incre_newNode==0){
								continue;
							}
							finalReserve[newNode]+=alpha*incre_newNode;
							if(C[1][newNode]==0){
								C[1][newNode]=1;
								C[0][finalreserve_count++]=newNode;
							}
							double outVertWeight_newNode=g.getOutVertWeight(newNode);
							if(NodeResidue[newNode]>=0){
								NodeResidue[newNode]+=incre_newNode;
								if(HC[1][newNode]==0){
									HC[1][newNode]=1;
									HC[0][NodeResidue_count++]=newNode;
								}
								if((isInArray[newNode]==false)&&(NodeResidue[newNode]>eps*outVertWeight_newNode)){
									candidate_set[candidate_end]=newNode;
									candidate_end=(candidate_end+1)%(vert+1);
									isInArray[newNode]=true;
								}
							}
							else{
								residue[newNode]+=incre_newNode;
								if(U[1][newNode]==0){
									U[1][newNode]=1;
									U[0][residue_count++]=newNode;
								}

								if(isInArray[newNode]==false){
									double newGKey_newNode=PQ.findMin(newNode)-(1-alpha)*residue[newNode]/outVertWeight_newNode;
									if(newGKey_newNode<0){
										candidate_set[candidate_end]=newNode;
										candidate_end=(candidate_end+1)%(vert+1);
										isInArray[newNode]=true;
									}
								}
							}
						}
						break;
					}

					uint tempj_index=PQ.findMin_index(tempi);
					uint tempj_PL=g.outPL[tempi]+tempj_index;
					uint tempj=g.getOutVert(tempi,tempj_index);	    
	
					double outEdgeWeight_tempj=g.getOutEdgeWeight(tempi,tempj_index);
					double incre;
					incre=(1-alpha)*residue[tempi]*outEdgeWeight_tempj/outVertWeight-expense[tempj_PL];
					finalReserve[tempj]+=alpha*incre;
					if(C[1][tempj]==0){
						C[1][tempj]=1;
						C[0][finalreserve_count++]=tempj;
					}
			
					expense[tempj_PL]+=incre;
					if(H[1][tempj_PL]==0){
						H[1][tempj_PL]=1;
						H[0][expense_count++]=tempj_PL;
					}
					if(UC[1][tempi]==0){
						UC[1][tempi]=1;
						UC[0][resetPQ_count++]=tempi;
					}

					double newKey=PQ.findMin(tempi)+incre/outEdgeWeight_tempj;	
					PQ.increaseKey(tempi,newKey);  
				
					double outVertWeight_tempj=g.getOutVertWeight(tempj);
					if(NodeResidue[tempj]<0){
						residue[tempj]+=incre;
						if(U[1][tempj]==0){
							U[1][tempj]=1;
							U[0][residue_count++]=tempj;
						}
						if(isInArray[tempj]==false){
							double newGKey_j=PQ.findMin(tempj)-(1-alpha)*residue[tempj]/outVertWeight_tempj;
							if(newGKey_j<0){
								candidate_set[candidate_end]=tempj;
								candidate_end=(candidate_end+1)%(vert+1);
								isInArray[tempj]=true;
							}
						}
					}
					else{
						NodeResidue[tempj]+=incre;
						if(HC[1][tempj]==0){
							HC[1][tempj]=1;
							HC[0][NodeResidue_count++]=tempj;
						}
						if(isInArray[tempj]==false){
							if(NodeResidue[tempj]>eps*outVertWeight_tempj){
								candidate_set[candidate_end]=tempj;
								candidate_end=(candidate_end+1)%(vert+1);
								isInArray[tempj]=true;
							}
						}
					}
	  
					double newGKey_i=PQ.findMin(tempi)-(1-alpha)*residue[tempi]/outVertWeight;
					if(newGKey_i>=0){
						break;
					}
				}
			}
			else{
				double tempR=NodeResidue[tempi];
				NodeResidue[tempi]=0;
				if(HC[1][tempi]==0){
					HC[1][tempi]=1;
					HC[0][NodeResidue_count++]=tempi;
				}
				double partincre=(1-alpha)*tempR/outVertWeight;
				uint start_pos=g.outPL[tempi];
				uint end_pos=g.outPL[tempi+1];
				for(uint itek=start_pos;itek<end_pos;itek++){
					uint newNode=g.outEL[itek];
					double newEdgeWeight=g.outWEL[itek];
					double incre=partincre*newEdgeWeight;		
					finalReserve[newNode]+=alpha*incre;
					if(C[1][newNode]==0){
						C[1][newNode]=1;
						C[0][finalreserve_count++]=newNode;
					}
					double outVertWeight_newNode=g.getOutVertWeight(newNode);
					if(NodeResidue[newNode]>=0){
						NodeResidue[newNode]+=incre;
						if(HC[1][newNode]==0){
							HC[1][newNode]=1;
							HC[0][NodeResidue_count++]=newNode;
						}
						if(isInArray[newNode]==false){
							if(NodeResidue[newNode]>eps*outVertWeight_newNode){
								candidate_set[candidate_end]=newNode;
								candidate_end=(candidate_end+1)%(vert+1);
								isInArray[newNode]=true;
							}
						}
					}
					else{
						residue[newNode]+=incre;
						if(U[1][newNode]==0){
							U[1][newNode]=1;
							U[0][residue_count++]=newNode;
						}
						if(isInArray[newNode]==false){
							double newGKey_newNode=PQ.findMin(newNode)-(1-alpha)*residue[newNode]/outVertWeight_newNode;
							if(newGKey_newNode<0){
								candidate_set[candidate_end]=newNode;
								candidate_end=(candidate_end+1)%(vert+1);
								isInArray[newNode]=true;
							}
						}
					}
				}
			}
		}
		tempLevel++;
	}
  
};




class powermethod:public SimStruct{
	public:
	uint *candidate_set[2];
	uint candidate_count[2];
	double *residue[2];

	powermethod(string name, string file_label, double epsilon, double para_alpha): 
	SimStruct(name, file_label,epsilon,para_alpha){
		candidate_count[0]=0;
		candidate_count[1]=0;
		H[0] = new uint[vert];
		H[1] = new uint[vert];
		U[0] = new uint[vert];
		U[1] = new uint[vert];
		C[0] = new uint[vert];
		C[1] = new uint[vert];
		candidate_set[0] = new uint[vert];
		candidate_set[1] = new uint[vert];
		residue[0]=new double[vert];
		residue[1]=new double[vert];
		finalReserve = new double[vert];
		finalreserve_count=0;
		for(uint i = 0; i < vert; i++){
			residue[0][i]=0;
			residue[1][i]=0;
			finalReserve[i]=0;
			H[0][i] = 0;
			H[1][i] = 0;
			U[0][i] = 0;
			U[1][i] = 0;
			C[0][i] = 0;
			C[1][i] = 0;
			candidate_set[0][i]=0;
			candidate_set[1][i]=0;
		}
	}

	~powermethod() {
		delete[] H[0];
		delete[] H[1];
		delete[] U[0];
		delete[] U[1];
		delete[] C[0];
		delete[] C[1];
		delete[] residue[0];
		delete[] residue[1];
		delete[] finalReserve;
		delete[] candidate_set[0];
		delete[] candidate_set[1];
	}      
	
	void query(uint u){
		for(uint j = 0; j < finalreserve_count; j++){
			finalReserve[H[0][j]] = 0;
			H[1][H[0][j]] = 0;
		}
		finalreserve_count=0;
		uint tempLevel = 0;
	
		double w_i=alpha;
		double Y_i=1;

		residue[0][u] = 1;
		candidate_set[0][0]=u;
		candidate_count[0]=1;
		candidate_count[1]=0;
		uint L=100;
		//cout<<"L="<<L<<endl;
		
		while(tempLevel<=L){
			uint tempLevelID=tempLevel%2;
			uint newLevelID=(tempLevel+1)%2;
			uint candidateCnt=candidate_count[tempLevelID];
			//cout<<"Iteration "<<tempLevel<<": candidateCnt="<<candidateCnt<<endl;
			if(candidateCnt==0){
				//cout<<"candidateCnt=0 tempLevel="<<tempLevel<<endl;
				break;
			}
			candidate_count[tempLevelID]=0;
		
			for(uint j = 0; j < candidateCnt; j++){
				uint tempNode = candidate_set[tempLevelID][j];
				double tempR = residue[tempLevelID][tempNode];
				U[tempLevelID][tempNode]=0;
				residue[tempLevelID][tempNode] = 0;
				if(H[1][tempNode] == 0){
					H[0][finalreserve_count++] = tempNode;
					H[1][tempNode] = 1;
				}
				finalReserve[tempNode]+=alpha*tempR;
				if(tempLevel==L){
					continue;
				}

				uint outSize = g.getOutSize(tempNode);
				double outVertWt = g.getOutVertWeight(tempNode);
				double incre = tempR* (1-alpha)/outVertWt;
		
				for(uint k = 0; k < outSize; k++){
					uint newNode = g.getOutVert(tempNode, k);
					residue[newLevelID][newNode] += incre*g.getOutEdgeWeight(tempNode,k);
					if(U[newLevelID][newNode] == 0){
						U[newLevelID][newNode] = 1;
						candidate_set[newLevelID][candidate_count[newLevelID]++]=newNode;
					}
				}
			}
			tempLevel++;
		}
	}
};



#endif
