#ifndef PRIORITY_H
#define PRIORITY_H

#include <random>
#include <algorithm>
//#include <stack>
#include "Graph.h"

using namespace std;

typedef unsigned int uint;

class Priority{
public:
	pair<double,uint> **prior;
	uint vert;
	uint *length;

	Priority(){

	}

	~Priority(){
		for(uint i=0;i<vert;i++){
			delete[] prior[i];
		}
		delete[] length;
	}

	void initPriority(uint gn,pair<double,uint>** ginitsort,uint* outPL){
		vert=gn;
		prior=new pair<double,uint>*[vert];
		length=new uint[vert];
		for(uint i=0;i<vert;i++){
			//uint outdeg=g.getOutSize(i);
			uint outdeg=outPL[i+1]-outPL[i];
			prior[i]=new pair<double,uint>[outdeg];
			memcpy(prior[i],ginitsort[i],sizeof(pair<double,uint>)*outdeg);
			length[i]=outdeg;
			//if(i<5||i==485){
			//	for(uint j=0;j<outdeg;j++){
			//		cout<<"prior: "<<i<<" "<<prior[i][j].first<<" "<<prior[i][j].second<<endl;
			//	}
			//}
		}	
	}

	void resetPriority(pair<double,uint> *initsort_NID, uint resetNID, uint outdeg){
		memcpy(prior[resetNID],initsort_NID,sizeof(pair<double,uint>)*outdeg);
	}

	double findMin(uint NID){
		return prior[NID][0].first;
	}

	uint findMin_index(uint NID){
		return prior[NID][0].second;
	}

	int isbottom(uint NID, uint tmppos){
		uint tmplength=length[NID];
		if(2*tmppos>tmplength){
			return 0;
		}	
		else if((2*tmppos+1)>tmplength){
			return 1;
		}
		else{
			return 2;
		}
	}

	void increaseKey(uint NID, double newKey){//default, reset the key of the first element
		uint curtidx=1;
		prior[NID][0].first=newKey;
		while(isbottom(NID,curtidx)!=0){
			uint leftsonidx=curtidx*2;
			uint rightsonidx=curtidx*2+1;
			if(isbottom(NID,curtidx)==1){
				if(prior[NID][curtidx-1].first>prior[NID][leftsonidx-1].first){
					pair<double,uint> tmppair=prior[NID][curtidx-1];
					prior[NID][curtidx-1]=prior[NID][leftsonidx-1];
					prior[NID][leftsonidx-1]=tmppair;
				}
				curtidx=leftsonidx;
			}
			else{
				uint smallidx;
				if(prior[NID][leftsonidx-1].first<prior[NID][rightsonidx-1].first){
					smallidx=leftsonidx;
				}
				else{
					smallidx=rightsonidx;
				}
				if(prior[NID][curtidx-1].first>prior[NID][smallidx-1].first){
					pair<double,uint> tmppair=prior[NID][curtidx-1];
					prior[NID][curtidx-1]=prior[NID][smallidx-1];
					prior[NID][smallidx-1]=tmppair;
				}
				curtidx=smallidx;
			}
		}
	}	

};

#endif
