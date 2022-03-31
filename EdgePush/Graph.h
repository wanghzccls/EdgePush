#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
//#include <cstring>

using namespace std;

typedef unsigned int uint;

int mkpath(string s, mode_t mode=0755){
	size_t pre=0, pos;
	string dir;
	int mdret;
	if(s[s.size()-1]!='/'){
		s+='/';
	}
	while((pos=s.find_first_of('/',pre))!=string::npos){
		dir=s.substr(0,pos++);
		pre=pos;
		if(dir.size()==0) continue;
		if((mdret=::mkdir(dir.c_str(),mode)) && errno!=EEXIST){
			return mdret;
		}
	}
	return mdret;
}



class Graph{
public:
	uint n;	//number of nodes
	uint m;	//number of edges

	uint* outEL;
	uint* outPL;
	double* outWEL;//store the weight of each edge
	uint* inEL;
	uint* inPL;
	double* inWEL;
	
	double* outWeights;//the weight of each node
	double* inWeights;//the weight of each node
	double* sqrtoutWeights;//the sqrt of the degree of each node.	

	double sum_partrmax=0.0;//the common term in rmax(j):\sum_{j\in V} n_j^{2/3}\cdot d_j^{1/3}. 
	double totaldeg=0.0;

	double** rmax;
	pair<double,uint>** initsort;

	Graph(){
		
	}

	~Graph(){
		for(uint i=0;i<n;i++){
			delete[] rmax[i];
			delete[] initsort[i];
		}
		delete[] outEL;
		delete[] outPL;
		delete[] outWEL;
		delete[] outWeights;
		delete[] sqrtoutWeights;
		delete[] inEL;
		delete[] inPL;
		delete[] inWEL;
		delete[] inWeights;
	}


	void csrGraphChange(string filedir, string filelabel){
		cout<<"Change to csr format..."<<endl;
		original_inputGraph(filedir,filelabel);
		csr_convert(filedir,filelabel);
		//cout<<"n="<<n<<" m="<<m<<endl;
	}

	void inputGraph(string filedir, string filelabel,double eps){
		stringstream ss_dir,ss_attr,ss_outEL,ss_outPL,ss_outWEL,ss_inEL,ss_inPL,ss_inWEL;
	
		ss_attr<<filedir<<filelabel<<".attribute";
		ifstream in_attr;
		in_attr.open(ss_attr.str());
		if(!in_attr){
			csrGraphChange(filedir,filelabel);
			in_attr.open(ss_attr.str());
			if(!in_attr){
				cout<<"===Read graph ERROR!==="<<endl;
				return;
			}
		}
		cout<<"Read graph attributes..."<<endl;
		string tmp;
		in_attr>>tmp>>n;
		in_attr>>tmp>>m;
		cout<<"n="<<n<<" m="<<m<<endl;
	
		in_attr.close();

		cout<<"Read graph edges..."<<endl;
		ss_outEL<<filedir<<filelabel<<".outEdges";
		ss_outPL<<filedir<<filelabel<<".outPtr";
		ss_outWEL<<filedir<<filelabel<<".outWEdges";
		ss_inEL<<filedir<<filelabel<<".inEdges";
		ss_inPL<<filedir<<filelabel<<".inPtr";
		ss_inWEL<<filedir<<filelabel<<".inWEdges";
	
		outEL=new uint[m];
		outPL=new uint[n+1];
		outWEL=new double[m];
		outWeights=new double[n];
		sqrtoutWeights=new double[n]; 
		inEL=new uint[m];
		inPL=new uint[n+1];
		inWEL=new double[m];
		inWeights=new double[n];

		rmax=new double*[n];
		initsort=new pair<double,uint>*[n];

		ifstream outf(ss_outEL.str(),ios::in | ios::binary);
		outf.read((char *)&outEL[0],sizeof(outEL[0])*m);

		ifstream outpf(ss_outPL.str(),ios::in | ios::binary);
		outpf.read((char *)&outPL[0],sizeof(outPL[0])*(n+1));

		ifstream outef(ss_outWEL.str(),ios::in | ios::binary);
		outef.read((char *)&outWEL[0],sizeof(outWEL[0])*m);

		ifstream inf(ss_inEL.str(),ios::in | ios::binary);
		inf.read((char *)&inEL[0],sizeof(inEL[0])*m);
 
		ifstream inpf(ss_inPL.str(),ios::in | ios::binary);
		inpf.read((char *)&inPL[0],sizeof(inPL[0])*(n+1));
 
 		ifstream inef(ss_inWEL.str(),ios::in | ios::binary);
		inef.read((char *)&inWEL[0],sizeof(inWEL[0])*m);
		
		for(uint i=0;i<n;i++){
			uint outdeg=outPL[i+1]-outPL[i];
			double tmpoutw=0.0;
			for(uint j=0;j<outdeg;j++){
				tmpoutw+=outWEL[(outPL[i]+j)];
				sum_partrmax+=sqrt(outWEL[(outPL[i]+j)]);
			}
			
			totaldeg+=tmpoutw;
			outWeights[i]=tmpoutw;
			sqrtoutWeights[i]=sqrt(tmpoutw);

			uint indeg=inPL[i+1]-inPL[i];
			double tmpinw=0.0;
			for(uint j=0;j<indeg;j++){
				tmpinw+=inWEL[(inPL[i]+j)];
			}
			inWeights[i]=tmpinw;
		}
		sum_partrmax=totaldeg/sum_partrmax;

		
		vector<pair<double,uint> >tmpsort;
		for(uint i=0;i<n;i++){
			uint outdeg=outPL[i+1]-outPL[i];
			rmax[i]=new double[outdeg];
			initsort[i]=new pair<double,uint>[outdeg];
			
			for(uint j=0;j<outdeg;j++){
				double tmpedgew=outWEL[(outPL[i]+j)];
				tmpsort.push_back(make_pair(tmpedgew,j));
				rmax[i][j]=eps*sum_partrmax*sqrt(tmpedgew);
			}
			sort(tmpsort.begin(),tmpsort.end(),greater<pair<double,uint> >());
			
			for(int j=0;j<outdeg;j++){
				uint tmpindex=tmpsort[j].second;
				double tmpedgew=outWEL[(outPL[i]+tmpindex)];
				initsort[i][j]=make_pair(((rmax[i][tmpindex])/tmpedgew),tmpindex);
			}
			tmpsort.clear();	
		}

		outf.close();
		outpf.close();
		outef.close();
		inf.close();
		inpf.close();
		inef.close();

	}


	uint getInSize(uint vert){
		return (inPL[vert+1]-inPL[vert]);
	}
	uint getInVert(uint vert, uint pos){
		//return inAdjList[vert][pos];
		return inEL[(inPL[vert]+pos)];
	}
	double getInEdgeWeight(uint vert, uint pos){
		return inWEL[(inPL[vert]+pos)];
	}
	double getInVertWeight(uint vert){
		return inWeights[vert];
	}

	uint getOutSize(uint vert){
		return (outPL[vert+1]-outPL[vert]);
	}
	uint getOutVert(uint vert, uint pos){
		//return outAdjList[vert][pos];
		return outEL[(outPL[vert]+pos)];
	}
	double getOutEdgeWeight(uint vert, uint pos){
		return outWEL[(outPL[vert]+pos)];
	}
	double getOutVertWeight(uint vert){
		return outWeights[vert];
	}
	double getOutSqrtVertWeight(uint vert){
		return sqrtoutWeights[vert];
	}

private:	
	uint** inAdjList;
	uint** outAdjList;
	double** inAdjWList;
	double** outAdjWList;
	uint* indegree;
	uint* outdegree;
	
	void original_inputGraph(string filedir, string filelabel){
		m=0;
		string filename="dataset/"+filelabel+".txt";
		ifstream infile;
		infile.open(filename);
		if(!infile){
			cout<<"ERROR: unable to open source dataset: "<<filename<<endl;
			return;
		}

		cout<<"Read the original txt file..."<<endl;
		infile >> n;

		indegree=new uint[n];
		outdegree=new uint[n];
		for(uint i=0;i<n;i++){
			indegree[i]=0;
			outdegree[i]=0;
		}
		//read graph and get degree info
		uint from;
		uint to;
		double weight;
		while(infile>>from>>to>>weight){
			outdegree[from]++;
			indegree[to]++;
		}

		inAdjList=new uint*[n];
		outAdjList=new uint*[n];
		
		inAdjWList=new double*[n];
		outAdjWList=new double*[n];

		uint* pointer_in=new uint[n];
		uint* pointer_out=new uint[n];
		for(uint i=0;i<n;i++){
			inAdjList[i]=new uint[indegree[i]];
			outAdjList[i]=new uint[outdegree[i]];
			
			inAdjWList[i]=new double[indegree[i]];
			outAdjWList[i]=new double[outdegree[i]];

			pointer_in[i]=0;
			pointer_out[i]=0;
		}
		infile.clear();
		infile.seekg(0);

		clock_t t1=clock();
		infile >> n;
		//cout<<"Vertice num="<<n<<endl;
		while(infile>>from>>to>>weight){
			outAdjList[from][pointer_out[from]]=to;
			outAdjWList[from][pointer_out[from]]=weight;
			pointer_out[from]++;
			inAdjList[to][pointer_in[to]]=from;
			inAdjWList[to][pointer_in[to]]=weight;
			pointer_in[to]++;

			m++;
		}
		infile.close();
		clock_t t2=clock();
		//cout<<"Edge num="<<m<<endl;
		cout<<"reading in graph takes "<<(t2-t1)/(1.0*CLOCKS_PER_SEC)<<" s."<<endl;

		delete[] pointer_in;
		delete[] pointer_out;
	}

	uint gettxtInSize(uint vert){
		return indegree[vert];
	}
	uint gettxtInVert(uint vert, uint pos){
		return inAdjList[vert][pos];
	}
	double gettxtInVertWeight(uint vert, uint pos){
		return inAdjWList[vert][pos];
	}
	uint gettxtOutSize(uint vert){
		return outdegree[vert];
	}
	uint gettxtOutVert(uint vert, uint pos){
		return outAdjList[vert][pos];
	}
	double gettxtOutVertWeight(uint vert, uint pos){
		return outAdjWList[vert][pos];
	}

	void csr_convert(string filedir, string filelabel){
		cout<<"Write to csr file..."<<endl;
			
		uint *coutEL=new uint[m];
		uint *coutPL=new uint[n+1];
		double *coutWEL=new double[m];
		uint *cinEL=new uint[m];
		uint *cinPL=new uint[n+1];
		double *cinWEL=new double[m];
		coutPL[0]=0;
		uint outid=0;
		uint out_curnum=0;
  
		cinPL[0]=0;
		uint inid=0;
		uint in_curnum=0;
	
		for(uint i=0;i<n;i++){
			outid+=gettxtOutSize(i);
			coutPL[i+1]=outid;
			for(uint j=0;j<gettxtOutSize(i);j++){
				coutEL[out_curnum]=gettxtOutVert(i,j);
				coutWEL[out_curnum]=gettxtOutVertWeight(i,j);
				out_curnum+=1;
			}
			inid+=gettxtInSize(i);
			cinPL[i+1]=inid;
			for(uint j=0;j<gettxtInSize(i);j++){
				cinEL[in_curnum]=gettxtInVert(i,j);
				cinWEL[in_curnum]=gettxtInVertWeight(i,j);
				in_curnum+=1;
			}
		}

		stringstream ss_cdir,ss_cattr,ss_coutEL,ss_coutWEL,ss_coutPL,ss_cinEL,ss_cinPL,ss_cinWEL;
		ss_cdir<<filedir;
		mkpath(ss_cdir.str());    

		ss_cattr<<ss_cdir.str()<<filelabel<<".attribute";
		ofstream cout_attr;
		cout_attr.open(ss_cattr.str());
		if(!cout_attr){
			cout<<"ERROR: unable to open attribute file: "<<ss_cattr.str()<<endl;
			return;
		}

		cout_attr<<"n "<<n<<"\n";
		cout_attr<<"m "<<m<<"\n";
		cout_attr.close();
	
		ss_coutEL<<ss_cdir.str()<<filelabel<<".outEdges";
		ss_coutPL<<ss_cdir.str()<<filelabel<<".outPtr";
		ss_coutWEL<<ss_cdir.str()<<filelabel<<".outWEdges";
	
		ss_cinEL<<ss_cdir.str()<<filelabel<<".inEdges";
		ss_cinPL<<ss_cdir.str()<<filelabel<<".inPtr";
		ss_cinWEL<<ss_cdir.str()<<filelabel<<".inWEdges";
		
		ofstream foutEL(ss_coutEL.str(),ios::out | ios::binary);
		ofstream foutPL(ss_coutPL.str(),ios::out | ios::binary);
		ofstream foutWEL(ss_coutWEL.str(),ios::out | ios::binary);
		ofstream finEL(ss_cinEL.str(),ios::out | ios::binary);
		ofstream finPL(ss_cinPL.str(),ios::out | ios::binary);
		ofstream finWEL(ss_cinWEL.str(),ios::out | ios::binary);
	
		foutEL.write((char *)&coutEL[0],sizeof(outEL[0])*m);
		foutPL.write((char *)&coutPL[0],sizeof(outPL[0])*(n+1));
		foutWEL.write((char *)&coutWEL[0],sizeof(outWEL[0])*m);
		finEL.write((char *)&cinEL[0],sizeof(inEL[0])*m);
		finPL.write((char *)&cinPL[0],sizeof(inPL[0])*(n+1));
		finWEL.write((char *)&cinWEL[0],sizeof(inWEL[0])*m);
	
		foutEL.close();
		foutPL.close();
		foutWEL.close();
		finEL.close();
		finPL.close();
		finWEL.close();
	
		delete[] coutPL;
		delete[] coutEL;
		delete[] coutWEL;
		delete[] cinPL;
		delete[] cinEL;
		delete[] cinWEL;

		return;
	}
};


