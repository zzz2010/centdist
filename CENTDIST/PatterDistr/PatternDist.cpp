#include <cstring>
#include "stdafx.h"
#include "correlation.h"
#define MAXSEQLEN 10000
//#define NORMALIZE
#define SQ2PI 2.50662827
bool getThresh=false;
	float cutoff=0;
vector<unsigned long> rmPoslist;
vector<float> histdata;

void flushhistdata(double* CDhist,int length)
{
	histdata.clear();
	for(int i=0;i<length;i++)
		histdata.push_back(CDhist[i]);
}

void outputhistdata(PARAM * setting)
{
	string outfile=setting->inputFile.append(".hist");
	ofstream os(outfile.c_str());
	for(int i=0;i<histdata.size();i++)
		os<<histdata[i]<<endl;
	os.close();
}

/***********************************************/
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=0, i;
static char* syntax = "patterDistr -i inputPosFile  [-o outputDIR -findThresh  -ratio minSupportRatio -pdt specifiedThresh -rm rmfile -rmlen len -rt rmThresh -K order]";


//if(nargs < 3){ printf("%s\n",syntax); exit(0); }
p = new param_st();//(PARAM*)calloc(1,sizeof(PARAM));
p->seedlength= 5;
p->min_supp_ratio=0.05;
p->FDRthresh=1.e-7;
p->olThresh=0;
p->pdThresh=0.18;
p->outputDIR=string(".");
p->max_motif_length=26;
p->N_motif=2;
p->resolution=40;
p->startwSize=100;
p->PeakRange=4000;
p->weightFile="";

while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1){ p->inputFile=argv[++iarg]; }
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputDIR=argv[++iarg];
	else if(!strcmp(argv[iarg],"-findThresh") ) getThresh=true;
	else if(!strcmp(argv[iarg],"-rm") && iarg < nargs-1) p->weightFile=argv[++iarg];
	else if(!strcmp(argv[iarg],"-rt") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->olThresh);
	else if(!strcmp(argv[iarg],"-pdt") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&cutoff);
	else if(!strcmp(argv[iarg],"-ratio") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->min_supp_ratio);
	else if(!strcmp(argv[iarg],"-FDR") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->FDRthresh);
	else if(!strcmp(argv[iarg],"-maxlen") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->max_motif_length);
	else if(!strcmp(argv[iarg],"-K") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_motif);
	else if(!strcmp(argv[iarg],"-rmlen") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->seedlength);
	else if(!strcmp(argv[iarg],"-mbr") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->startwSize);
	else if(!strcmp(argv[iarg],"-rs") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->resolution);
	else if(!strcmp(argv[iarg],"-pr") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->PeakRange);
	else{
		//printf("ERROR: Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}

if(!p->inputFile[0]){ printf("ERROR: %s\n",syntax); exit(0); }
//
//
//if(p->N_motif<1){ printf("WARNING: Wrong number of output motifs. Using default N=20\n"); p->N_motif=20; }
//if(p->min_supp_ratio >1){ printf("WARNING: Minimum ratio should be <1. Using default=0.05\n"); p->min_supp_ratio=0.05; }
return(p);
}


int  split_string (char *string, char *items[], int num)
/***********************************************/
{
char *ch;
int i=0;

ch = strchr(string,'\n');
if(ch) *ch = '\0';
ch = string;
while(1){
	items[i] = ch;
	ch = strchr(ch,'\t');
	if(ch) *ch = '\0';
	else break;
	if(i>=num-1) break;
	ch++;
	i++;
}
return(i+1);
}

void LoadrmPoslist(PARAM * setting)
{

	FILE* fp = fopen(setting->weightFile.c_str(),"r");
	
	char buffer[300];
	char *items[10];

	int SEQLEN=0;
	int maxSeqNum=0;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		{
			float score;
			sscanf(items[2],"%f",&score);
			if(score<setting->olThresh)
				continue;
		}
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		if(pos>SEQLEN)
			SEQLEN=pos;
		rmPoslist.push_back((unsigned long)pos+seqid*MAXSEQLEN);
	}
	sort(rmPoslist.begin(),rmPoslist.end());
	//cout<<rmPoslist[rmPoslist.size()-1]<<endl;
		fclose(fp);
}

bool testInside(vector<unsigned long>& v, int extendlen,unsigned long searchVal,PARAM * setting)
{
	vector<unsigned long>::iterator low,up;
	low=lower_bound (v.begin(), v.end(), searchVal); //    
	
	if(low==v.end()||low!=v.begin()&&(unsigned long)*low!=searchVal)
		low--;
  up= upper_bound (v.begin(), v.end(), searchVal); //  
		if(up==v.end())
		up--;
		if((searchVal-(unsigned long)*low)<setting->seedlength)
			  return true;
		  if(((unsigned long)*up-searchVal)<setting->seedlength)
			  return true;
  return false;
}

void CenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}
		/***modified***/
		//else if(abs(lastpos-pos)<setting->resolution/2)
		//{
		//	//cout<<"a";
		//	//lastpos=pos;
		//	continue;
		//}
		//else
		//	lastpos=pos;

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;
//ap::real_1d_array orderNum;
//orderNum.setbounds(0,BINNUMBER-1);
//	for(int i=0;i<BINNUMBER;i++)
//		orderNum(i)=i;

lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
		CDHist[i]=0;
		
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE
			//if(count>maxSeqNum)
			//sumNorm=sumNorm/(maxSeqNum+1);
			//else
				//sumNorm=sumNorm/(positionlist.size());
			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;
				//cout<<CDHist[j]<<endl;
				//s+=CDHist[j];
			}
#endif
/**modified**/
			CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];

			//ap::real_1d_array CDrank;
			//CDrank.setcontent(0,BINNUMBER-1,CDHist);
			//double rscore=spearmanrankcorrelation(orderNum,CDrank,BINNUMBER);
		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-1;j++)
			{
				/**modified**/
				double temp=fabs(CDHist[j]-CDHist[j+1]);
				if(temp<0)
					temp=0;
				CDDiff[j]=temp;
			
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10)
					break;
			if(minbin<50)
			{
				setting->resolution=setting->resolution*2;
				BINNUMBER= ceil((double)SEQLEN/setting->resolution);
				BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));

				continue;
			}
			
	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/

/**********use 0-order************/
			for(int j=0;j<BINNUMBER-2;j++)
			{
				/**modified**/
				double temp=fabs(CDHist[j]-2*CDHist[j+1]+CDHist[j+2]);
				if(temp<0)
					temp=0;
				
				CDDiff2[j]=temp;
				ssr+=CDDiff2[j];
			}
/**********use 0-order************/
		
				for(int i=0;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				if(Z0<3)
				{
					Z1=0;
				//if(Z2>2*Z0||Z2>2*Z1)
					Z2=0;
				}
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
				//else
				//	up=!up;
				//if(up)
				//{
				//		step=step/2;
				//		for(int k=0;k<step;k++)
				//			ITER++;
				//}
				//else
				//{
				//	step=step/2;
				//		for(int k=0;k<step;k++)
				//			ITER--;
				//}
				

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			bestscore+=bestZ1score;
		if(setting->N_motif==3)
			bestscore+=bestZ1score+bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore/(setting->N_motif*setting->N_motif)<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void UDSmoothCenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
					CDDiff[i]=0;
				if(i<BINNUMBER-2)
					CDDiff2[i]=0;

			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];
for(int j=0;j<BINNUMBER-slidewin;j++)
{
	for(int kk=j;kk<j+slidewin;kk++)
		SmoothCDHist[j]+=CDHist[kk];
		
	SmoothCDHist[j]=SmoothCDHist[j]/slidewin;
	
}
for(int rr=BINNUMBER-slidewin;rr<BINNUMBER;rr++)
{
	SmoothCDHist[rr]=SmoothCDHist[BINNUMBER-slidewin-1];
}

		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin]);
				//cout<<temp<<endl;
				if(temp>0)
				CDDiff[j]=1;
				else
				CDDiff[j]=0;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;
			//if(minbin<10)
			//{
			//	setting->resolution=setting->resolution*2;
			//	BINNUMBER= ceil((double)SEQLEN/setting->resolution);
			//	BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
			//	slidewin=ceil((double)slidewin/2);

			//	continue;
			//}

	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/
ssr=0;
/**********use 0-order************/
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(SmoothCDHist[j]-2*SmoothCDHist[j+(int)slidewin/2]+SmoothCDHist[j+slidewin]);
				if(temp>0)//&&CDDiff[j+(int)slidewin/2]
				CDDiff2[j]=1;
				else
				CDDiff2[j]=0;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
			break;
/**********use 0-order************/
		//start i=1 , to avoid to few sample in the center!
				for(int i=1;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(binominalCDF(ratio2,cddiff , n2));//(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));invnormaldistribution
				double Z2=(binominalCDF(ratio3,cddiff2 , n3));//(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)
				//{
				//	Z1=0;
				////if(Z2>2*Z0||Z2>2*Z1)
				//	Z2=0;
				//}
/**modified**/
				//if(Z2>1000)
				//	Z2=0;
				//if(Z1>1000)
				//	Z1=0;
				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 
				//if(score<0)
				//	score=0;

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			bestscore=bestZ2score;
			//bestscore+=bestZ1score+bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}



void UDCenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
					CDDiff[i]=0;
				if(i<BINNUMBER-2)
					CDDiff2[i]=0;

			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];


		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDHist[j]-CDHist[j+slidewin]);
				//cout<<temp<<endl;
				if(temp>0)
				CDDiff[j]=1;
				else
				CDDiff[j]=0;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;
			//if(minbin<10)
			//{
			//	setting->resolution=setting->resolution*2;
			//	BINNUMBER= ceil((double)SEQLEN/setting->resolution);
			//	BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
			//	slidewin=ceil((double)slidewin/2);

			//	continue;
			//}

	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/
ssr=0;
/**********use 0-order************/
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDHist[j]-2*CDHist[j+(int)slidewin/2]+CDHist[j+slidewin]);
				if(temp>0)//&&CDDiff[j+(int)slidewin/2]
				CDDiff2[j]=1;
				else
				CDDiff2[j]=0;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
			break;
/**********use 0-order************/
		//start i=1 , to avoid to few sample in the center!
				for(int i=1;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=invnormaldistribution(binominalCDF(ratio2,cddiff , n2));//(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=invnormaldistribution(binominalCDF(ratio3,cddiff2 , n3));//(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)
				//{
				//	Z1=0;
				////if(Z2>2*Z0||Z2>2*Z1)
				//	Z2=0;
				//}
/**modified**/
				if(Z2>1000)
					Z2=0;
				if(Z1>1000)
					Z1=0;
				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 
				//if(score<0)
				//	score=0;

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			bestscore+=bestZ1score;
		if(setting->N_motif==3)
			bestscore+=bestZ1score+bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void SmoothCenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)200/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
					if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];
			double mean=0;
for(int j=0;j<BINNUMBER;j++)
{
	mean+=CDHist[j];
	
}
mean/=BINNUMBER;
double sd=0;
for(int j=0;j<BINNUMBER;j++)
{
	sd+=(CDHist[j]-mean)*(CDHist[j]-mean);
}
sd/=BINNUMBER;
double sqsd=SQ2PI*sqrt(sd);
for(int j=0;j<BINNUMBER-slidewin+1;j++)
{
	double gnorm=0;
	for(int jj=j;jj<j+slidewin;jj++)
	{
		double gauss=exp(-1.0*(j-jj)*(j-jj)/(2*sd))/(sqsd);
		SmoothCDHist[j]+=CDHist[jj]*gauss;
		gnorm+=gauss;
	}
	SmoothCDHist[j]/=gnorm;

	cout<<SmoothCDHist[j]<<"\t"<<CDHist[j]<<endl;
}
for(int j=0;j<slidewin-1;j++)
SmoothCDHist[BINNUMBER-j-1]=SmoothCDHist[BINNUMBER-slidewin];


		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin])/slidewin;
				//cout<<temp<<endl;
				if(temp<0)
					temp=0;
				CDDiff[j]=temp;
				ssr+=CDDiff[j];
				//cout<<CDDiff[j]<<endl;
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;
			if(minbin<10)
			{
				setting->resolution=setting->resolution*2;
				BINNUMBER= ceil((double)SEQLEN/setting->resolution);
				BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));

				continue;
			}
			
	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/
			ssr=0;

/**********use 0-order************/
			for(int j=0;j<BINNUMBER-slidewin-1;j++)
			{
				/**modified**/
				double temp=(CDDiff[j]-CDDiff[j+slidewin])/slidewin;//(SmoothCDHist[j]-2*SmoothCDHist[j+1]+SmoothCDHist[j+2]);
				if(temp<0)
					temp=0;
				CDDiff2[j]=temp;
				ssr+=CDDiff2[j];
			   //cout<<CDDiff2[j]<<endl;
			}
			if(ssr==0)
				break;
/**********use 0-order************/
		
				for(int i=0;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				if(Z0<10)
				{
					Z1=0;
				//if(Z2>2*Z0||Z2>2*Z1)
					Z2=0;
				}
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			bestscore+=bestZ1score;
		if(setting->N_motif==3)
			bestscore+=bestZ1score+bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore/(setting->N_motif*setting->N_motif)<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void PeakShapeScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			smooth(CDHist,SmoothCDHist,BINNUMBER);
//for(int j=0;j<BINNUMBER-slidewin;j++)
//{
//	for(int kk=j;kk<j+slidewin;kk++)
//		SmoothCDHist[j]+=CDHist[kk];
//		
//	SmoothCDHist[j]=SmoothCDHist[j]/slidewin;
//	//SmoothCDHist[j]=CDHist[j];
//	
//}
//for(int rr=BINNUMBER-slidewin;rr<BINNUMBER;rr++)
//{
//	SmoothCDHist[rr]=SmoothCDHist[BINNUMBER-slidewin-1];
//}

		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin]);
				//cout<<CDHist[j]<<endl;
		/*		if(temp<0)
					temp=0;*/
				CDDiff[j]=temp;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;


/**********use 0-order************/
			ssr=0;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDDiff[j]-fabs(CDDiff[j+slidewin]));//(SmoothCDHist[j]-2*SmoothCDHist[j+(int)slidewin/2]+SmoothCDHist[j+slidewin]);//

				if(temp<0)
					temp=0;
				CDDiff2[j]=temp;
				//cout<<temp<<endl;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
				break;
/**********use 0-order************/
			int i=400/setting->resolution;
				//for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
						//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=0;j<i-1;j++)
					{
						if(localmax< SmoothCDHist[j]&&j<2*slidewin)
						{
							localmax=SmoothCDHist[j];
							localmaxI=j;
						}
					}
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						if(CDDiff[j]>0&&j>=localmaxI)
						cddiff+=CDDiff[j];
						if(CDDiff2[j]>0&&j>=localmaxI)
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1&&CDDiff[j]>0)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
						//bgdiff2+=fabs(CDDiff[j])+fabs(CDDiff[j+slidewin]);
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1-localmaxI)/(BINNUMBER-slidewin-localmaxI);
				double ratio3=(double)(i+1-localmaxI)/(BINNUMBER-slidewin-localmaxI);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i));//SmoothCDHist[i+slidewin]);
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-slidewin));//CDDiff[i+slidewin]);//
				cddiff-=bias1/2;
				cddiff2-=bias2/2;

				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));


				score=Z0; //+ssr*ssr


					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	
				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif>1&&CDScore>9.6118)
			bestscore=1000+bestZ1score;
			//bestscore=bestZ1score;
		if(setting->N_motif>2&&CDScore>9.6118&&bestZ1score>9.6118)
			bestscore=2000+bestZ2score;
		//bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void WholeShapeScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)100/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			smooth(CDHist,SmoothCDHist,BINNUMBER);
 

		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin]);
				//cout<<CDHist[j]<<endl;
		/*		if(temp<0)
					temp=0;*/
				CDDiff[j]=temp;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;


/**********use 0-order************/
			ssr=0;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDDiff[j]-(CDDiff[j+slidewin]));//(SmoothCDHist[j]-2*SmoothCDHist[j+(int)slidewin/2]+SmoothCDHist[j+slidewin]);//

				if(temp<0)
					temp=0;
				CDDiff2[j]=temp;
				//cout<<temp<<endl;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
				break;
/**********use 0-order************/
			int i=400/setting->resolution;
				//for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
						//find peak location 
					int localmaxI=0;
					double localmax=0;
					int startpoint=0;
					double minacc=1000;
					for(int j=0;j<i-1;j++)
					{
						if(localmax< SmoothCDHist[j])
						{
							localmax=SmoothCDHist[j];
							localmaxI=j;
						}	
					}
					// find the transition start point
					//for(int j=localmaxI;j<i/2;j++)
					//{
					//	if(minacc>fabs(CDDiff[j]))
					//	{
					//		minacc=fabs(CDDiff2[j]);
					//		startpoint=j;
					//	}
					//}
					startpoint=i;
					bestwindowId=i=startpoint;
					int nc,pc;
					nc=pc=0;
					double opVal=0;
					double proVal=0;
					if(i>BINNUMBER-slidewin)
						i=BINNUMBER-slidewin;
					for(int j=0;j<bestwindowId+1;j++)
					{
						double temp=SmoothCDHist[j]-SmoothCDHist[j+1];
						if(temp*(j-localmaxI)<0)
						{
							nc++;
							opVal+=fabs(temp);
						}
						else
						{
							pc++;
							proVal+=fabs(temp);
						}

					}
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						if(CDDiff[j]>0&&j>=localmaxI)
						cddiff+=CDDiff[j];
						if(CDDiff2[j]>0&&j>=localmaxI)
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1&&CDDiff[j]>0)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
						//bgdiff2+=fabs(CDDiff[j])+fabs(CDDiff[j+slidewin]);
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=0.5;
				double ratio3=0.5;
				double n=bgcnter+cdcnter;
				//double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i));//SmoothCDHist[i+slidewin]);
				//double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-slidewin));//CDDiff[i+slidewin]);//
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double ratio4=(double)(i+1-localmaxI)/(BINNUMBER-slidewin-localmaxI);
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i));//SmoothCDHist[i+slidewin]);
				cddiff-=bias1;
			

				double n2=pc+nc;
				double n3=proVal+opVal;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));// binominalCDF(ratio,cdcnter,n);//
				//double Z1=(proVal-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));//binominalCDF(0.5,pc,pc+nc);//(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				//double Z2=(pc-n2*ratio2)/sqrt(n2*0.5*0.5);//binominalCDF(0.5,proVal,proVal+opVal);

				double n4=bgdiff+cddiff;
				double Z1=(cddiff-n4*ratio4)/sqrt(n4*ratio4*(1-ratio4));
				double Z2=(proVal-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));

				score=Z0; //+ssr*ssr


					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	
				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
	if(setting->N_motif>1)
			bestscore+=bestZ1score;
			//bestscore=bestZ1score;
		if(setting->N_motif>2)
			bestscore+=bestZ2score;
		//bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void FARSmoothCenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			smooth(CDHist,SmoothCDHist,BINNUMBER);
//for(int j=0;j<BINNUMBER-slidewin;j++)
//{
//	for(int kk=j;kk<j+slidewin;kk++)
//		SmoothCDHist[j]+=CDHist[kk];
//		
//	SmoothCDHist[j]=SmoothCDHist[j]/slidewin;
//	//SmoothCDHist[j]=CDHist[j];
//	
//}
//for(int rr=BINNUMBER-slidewin;rr<BINNUMBER;rr++)
//{
//	SmoothCDHist[rr]=SmoothCDHist[BINNUMBER-slidewin-1];
//}

		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin]);
				//cout<<CDHist[j]<<endl;
				if(temp<0)
					temp=0;
				CDDiff[j]=temp;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;


/**********use 0-order************/
			ssr=0;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDDiff[j]-fabs(CDDiff[j+slidewin]));//(SmoothCDHist[j]-2*SmoothCDHist[j+(int)slidewin/2]+SmoothCDHist[j+slidewin]);//

				if(temp<0)
					temp=0;
				CDDiff2[j]=temp;
				//cout<<temp<<endl;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
				break;
/**********use 0-order************/
		
				for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						if(CDDiff[j]>0)
						cddiff+=CDDiff[j];
						if(CDDiff2[j]>0)
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
						//bgdiff2+=fabs(CDDiff[j])+fabs(CDDiff[j+slidewin]);
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-slidewin);
				double ratio3=(double)(i+1)/(BINNUMBER-slidewin);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-slidewin));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-slidewin));
				cddiff-=bias1/2;
				cddiff2-=bias2/2;

				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));


				score=Z0; //+ssr*ssr


					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	
				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif>1&&CDScore>8.3452)
			bestscore=1000+bestZ1score;
			//bestscore=bestZ1score;
		if(setting->N_motif>2&&CDScore>8.3452&&bestZ1score>8.3452)
			bestscore=2000+bestZ2score;
		//bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void FARSmoothQuantiScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			smooth(CDHist,SmoothCDHist,BINNUMBER);
			//SmoothCDHist=CDHist;
			//CDHist=SmoothCDHist;
/**********use 0-order************/
		
				for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=0;j<i-1;j++)
					{
						if(localmax< SmoothCDHist[j]&&j<2*slidewin)
						{
							localmax=SmoothCDHist[j];
							localmaxI=j;
						}
					}
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
					}
					int localwin=i-localmaxI;
					double lwsize=(localwin*setting->resolution/2);
					double avgVel=0;
					double avgAcc=0;
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-localwin)
						{
							CDDiff[j]=(double)(SmoothCDHist[j]-SmoothCDHist[j+localwin])/lwsize;
							avgVel+=CDDiff[j];
						}
						if(j<BINNUMBER-2*localwin)
						{
						   CDDiff2[j]=(double)(SmoothCDHist[j]-SmoothCDHist[j+(int)localwin]+fabs(SmoothCDHist[j+2*localwin]-SmoothCDHist[j+(int)localwin]))/pow(lwsize,2);
						   avgAcc+=CDDiff2[j];
						}
					}
					avgVel/=(double)(BINNUMBER-localwin-i);
					avgAcc/=(double)(BINNUMBER-2*localwin-i);
					double sdVel=0;
					double sdAcc=0;
					for(int j=i+1;j<BINNUMBER-localwin;j++)
					{
						sdVel+=(CDDiff[j]-avgVel)*(CDDiff[j]-avgVel);
						sdAcc+=(CDDiff2[j]-avgAcc)*(CDDiff2[j]-avgAcc);
					}
					sdVel=sqrt(sdVel/(double)(BINNUMBER-localwin-i));
					sdAcc=sqrt(sdAcc/(double)(BINNUMBER-2*localwin-i));
				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(SmoothCDHist[localmaxI]-SmoothCDHist[i])/lwsize/sdVel;
				double Z2=(SmoothCDHist[localmaxI]-SmoothCDHist[i]+fabs(SmoothCDHist[i+localwin]-SmoothCDHist[i]))/pow(lwsize,2)/sdAcc;
/**modified**/
				//if(Z0<3)//(Z1>2*Z0||Z2>2*Z0)
				//{
				//	Z1=0;
				//	Z2=0;
				//}/*
				//if(Z2<0)
				//{
				//	Z2=0;
				//}*/
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;		
						if(setting->N_motif==3)
							bestwindowId=localwin;		
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			//bestscore+=bestZ1score+bestZ2score;
		bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void BasicQuantiScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			//smooth(CDHist,SmoothCDHist,BINNUMBER);
			//SmoothCDHist=CDHist;
			//CDHist=SmoothCDHist;
/**********use 0-order************/
		
				for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER-i-1);

					for(int j=0;j<BINNUMBER-1;j++)
					{
						double temp=CDHist[j]-avgVel;
						if(temp<0)
							CDDiff[j]=0;
						else
							CDDiff[j]=temp;

						if(j<BINNUMBER-2)
						CDDiff2[j]=fabs(temp);
						if(temp<0&&j<i+1)
							CDDiff2[j]=0;


					}
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}

					for(int j=i+1;j<BINNUMBER;j++)
					{
						
						if(j<BINNUMBER-1)
							bgdiff+=CDDiff[j];
						
						if(j<BINNUMBER-2)
						{
						   
						   bgdiff2+=CDDiff2[j];
						}
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)//(Z1>2*Z0||Z2>2*Z0)
				//{
				//	Z1=0;
				//	Z2=0;
				//}/*
				//if(Z2<0)
				//{
				//	Z2=0;
				//}*/
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;		
	
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			//bestscore+=bestZ1score+bestZ2score;
		bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}


void FixWindScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	int win1=50;
	int win2=200;
	int win3=500;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			//smooth(CDHist,SmoothCDHist,BINNUMBER);
			//SmoothCDHist=CDHist;
			//CDHist=SmoothCDHist;
/**********use 0-order************/
			getThresh=false;
		double Z0,Z1,Z2;
							double cdcnter=0;
					double bgcnter=0;	
					double score=0;
		bool rmNoise=0;
			int i=win1*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;

					//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=BINNUMBER/2;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER/2);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=BINNUMBER/2;j<BINNUMBER;j++)
						{
							double temp=CDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}
				double ratio=(double)(i+1)/(BINNUMBER/2+i+1);

				double n=bgcnter+cdcnter;
				Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==1)
					{
						score=Z0; //+ssr*ssr
						bestwindowId=i;
					}
				}
			rmNoise=1;
			 i=win1*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;
					//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=BINNUMBER/2;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER/2);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=BINNUMBER/2;j<BINNUMBER;j++)
						{
							double temp=CDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}  
				double ratio=(double)(i+1)/(BINNUMBER/2+i+1);

				double n=bgcnter+cdcnter;
				   Z1=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==2)
					{
						score=Z1; //+ssr*ssr
						bestwindowId=i;
					}
				}

				 i=win2*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;

					for(int j=BINNUMBER/2;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER/2);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=BINNUMBER/2;j<BINNUMBER;j++)
						{
							double temp=CDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}
				double ratio=(double)(i+1)/(BINNUMBER/2+i+1);


				double n=bgcnter+cdcnter;
				   Z2=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==3)
					{
						score=Z2; //+ssr*ssr
						bestwindowId=i;
					}
				}

					if(score>=bestscore&&cdcnter>50)
					{
								
	
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			//bestscore+=bestZ1score+bestZ2score;
		bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void UseredBGScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	int p=setting->inputFile.find_last_of("/")+1;
	string bgfile="/home/chipseq/public_html/Pomoda/userdata/ES/bg/";
	bgfile=bgfile.append(setting->inputFile.substr(p).c_str());
	FILE* bgfp = fopen(bgfile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
		vector<float> bgscorelist;
	vector<int> bgpositionlist;
	set<float> bgsortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	int win1=50;
	int win2=200;
	int win3=500;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}


SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* bgCDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;
				int bgmaxSeqNum=0;
/*********************************************************BG set*********************************************/
						while(fgets(buffer,300,bgfp)){
						int seqid=0;
						int pos=0;

						int n = split_string(&buffer[0],items,10);
						if(n<2)
							continue;
						sscanf(items[0],"%d",&seqid);
						sscanf(items[1],"%d",&pos);
						//cout<<seqid<<endl;
						//if(getThresh)
						if(bgmaxSeqNum<seqid)
							bgmaxSeqNum=seqid;
						
						if(rmPoslist.size()>0)
						{
							bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
							if(rmflag)
								continue;
						}
						//filter repeat
						if(lastseq!=seqid)
						{
							lastseq=seqid;
							lastpos=-1000;
				#ifdef NORMALIZE
							for(int i=0;i<sameSeqPos.size();i++)
							{
								int p=sameSeqPos[i];
								problist[p]=problist[p]/sumSameSeq;
							}

							sameSeqPos.clear();
							sumSameSeq=0;
				#endif
						}

						{
							float score;
							sscanf(items[2],"%f",&score);
							bgscorelist.push_back(score);
							bgsortedThresh.insert(score);
				#ifdef NORMALIZE
							sameSeqPos.push_back(scorelist.size()-1);
							double temp=exp(score*log10);
							sumSameSeq+=temp;
							problist.push_back(temp);
				#endif
						}

						if(pos>SEQLEN)
							SEQLEN=pos;
						bgpositionlist.push_back(pos);
					}
						fclose(bgfp);
/*********************************************************BG set*********************************************/

lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				bgCDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			double bgCNT=0;
			for(int i=0;i<bgpositionlist.size();i++)
			{
				int pos=bgpositionlist[i];
				if(bgscorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=maxSeqNum/bgmaxSeqNum; //important normalize
						bgCDHist[tt]+=step;
				}
			}
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			//smooth(CDHist,SmoothCDHist,BINNUMBER);
			//SmoothCDHist=CDHist;
			//CDHist=SmoothCDHist;
/**********use 0-order************/
			getThresh=false;
		double Z0,Z1,Z2;
							double cdcnter=0;
					double bgcnter=0;	
					double score=0;
		bool rmNoise=0;
			int i=win1*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;

					for(int j=0;j<BINNUMBER;j++)
					{
						bgcnter+=bgCDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=0;j<BINNUMBER;j++)
						{
							double temp=bgCDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}
				double ratio=(double)(i+1)/(BINNUMBER+i+1);

				double n=bgcnter+cdcnter;
				Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==1)
					{
						score=Z0; //+ssr*ssr
						bestwindowId=i;
					}
				}
			rmNoise=1;
			 i=win1*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;
					//find peak location 
					for(int j=0;j<BINNUMBER;j++)
					{
						bgcnter+=bgCDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=0;j<BINNUMBER;j++)
						{
							double temp=bgCDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}
				double ratio=(double)(i+1)/(BINNUMBER+i+1);

				double n=bgcnter+cdcnter;
				   Z1=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==2)
					{
						score=Z1; //+ssr*ssr
						bestwindowId=i;
					}
				}

				 i=win2*2/setting->resolution-1;
				{
					cdcnter=0;
					bgcnter=0;

				for(int j=0;j<BINNUMBER;j++)
					{
						bgcnter+=bgCDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER);
					for(int j=0;j<i+1;j++)
					{
						if(rmNoise)
						cdcnter+=CDHist[j]-avgVel;
						else
						cdcnter+=CDHist[j];

						
					}		  
					if(rmNoise)
					{
						bgcnter=0;
						for(int j=0;j<BINNUMBER;j++)
						{
							double temp=bgCDHist[j]-avgVel;
							if(temp>0)
								bgcnter+=temp;
						}
					}
				double ratio=(double)(i+1)/(BINNUMBER+i+1);


				double n=bgcnter+cdcnter;
				   Z2=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
					if(setting->N_motif==3)
					{
						score=Z2; //+ssr*ssr
						bestwindowId=i;
					}
				}

					if(score>=bestscore&&cdcnter>50)
					{
								
	
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			//bestscore+=bestZ1score+bestZ2score;
		bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}



void BasicPlusQuantiScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)80/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestS1,bestS2,bestS3,bestS4;
bestS1=bestS2=bestS3=bestS4=0;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
			for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)ceil(60.0/setting->resolution)];}
			smooth(CDHist,SmoothCDHist,BINNUMBER,slidewin);

					



			//SmoothCDHist=CDHist;
			//CDHist=SmoothCDHist;
/**********use 0-order************/
		
				for(int i=slidewin;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					//find peak location 
					int localmaxI=0;
					double localmax=0;
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
					}
					double avgVel=bgcnter/(double)(BINNUMBER-i-1);

					for(int j=0;j<BINNUMBER-1;j++)
					{
						double temp=CDHist[j]-avgVel;
						if(temp<0)
							CDDiff[j]=0;
						else
							CDDiff[j]=temp;

						if(j<BINNUMBER-2)
						CDDiff2[j]=fabs(temp);
						if(temp<0&&j<i+1)
							CDDiff2[j]=0;
					}
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}

					for(int j=i+1;j<BINNUMBER;j++)
					{
						
						if(j<BINNUMBER-1)
							bgdiff+=CDDiff[j];
						
						if(j<BINNUMBER-2)
						{
						   
						   bgdiff2+=CDDiff2[j];
						}
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-1);
				double ratio3=(double)(i+1)/(BINNUMBER-2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-1));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)//(Z1>2*Z0||Z2>2*Z0)
				//{
				//	Z1=0;
				//	Z2=0;
				//}/*
				//if(Z2<0)
				//{
				//	Z2=0;
				//}*/
/**modified**/

				score=Z0; //+ssr*ssr
			

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;		
	
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;


		/***********************************plus section*******************************/
								double scorePlUS=0;
							//plus section
							int pc=0;
							int nc=0;
							//find peak location 
							int localmaxI=0;
							double localmax=0;
							for(int j=0;j<bestwindowId-1;j++)
							{
								if(localmax< SmoothCDHist[j]&&j<2*slidewin)
								{
									localmax=SmoothCDHist[j];
									localmaxI=j;
								}
							}
							
							for(int j=0;j<bestwindowId+1;j++)
							{
								double temp=SmoothCDHist[j]-SmoothCDHist[j+1];
								if(temp*(j-localmaxI)<0)
									nc++;
								else
									pc++;

							}
							scorePlUS=(double)(pc)/(pc+nc);
							double sP2=(SmoothCDHist[localmaxI]-SmoothCDHist[localmaxI+10])/fabs(SmoothCDHist[localmaxI+10]-SmoothCDHist[localmaxI+20]+0.001);
							double sP3=0;
							double sP4=0;
							for(int j=localmaxI;j<bestwindowId+5;j++)
							{
								
								double temp=SmoothCDHist[j]-SmoothCDHist[j+1];
								sP3+=temp*fabs(temp);
								temp=SmoothCDHist[j]-SmoothCDHist[j+slidewin];
								sP4+=temp*fabs(temp);
							}

							bestS1=scorePlUS;
							bestS2=sP2;
							bestS3=sP3;
							bestS4=sP4;
		/***********************************plus section*******************************/
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
						

				if(!getThresh)
					break;
				ITER++;
				step--;
		}while(step>1);



		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
			bestscore=bestZ1score;
		if(setting->N_motif==3)
			//bestscore+=bestZ1score+bestZ2score;
		bestscore=bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			
			//cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<sP2<<"\t"<<scorePlUS<<endl; 
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<bestS1<<"\t"<<bestS2<<"\t"<<bestS3<<"\t"<<bestS4<<endl; 
}

double PWMScoreBIASZstat(map<float,float> scoredist,float split=-1)
{
	int i,j,k;
	float mu1,mu2,var1,var2,n1,n2;
	mu1=mu2=var1=var2=n1=n2=0;
	if(split==-1)
	{
		//half and half
		map<float,float>::iterator iter=scoredist.begin();
		for(i=0;i<scoredist.size()/2;i++)
			iter++;
		split=iter->first;
	}
	map<float,float>::iterator iter=scoredist.begin();
		for(i=0;i<scoredist.size();i++)
		{
			if(iter->first<=split)
			{
				n1++;
				mu1+=iter->second;
				var1+=iter->second*iter->second;

			}
			else
			{
				n2++;
				mu2+=iter->second;
				var2+=iter->second*iter->second;
			}
			iter++;
		}
		mu1=mu1/n1;
		mu2=mu2/n2;
		var1=var1/n1-mu1*mu1;
		var2=var2/n2-mu2*mu2;
		float score=(mu1-mu2)/(sqrt(var1/n1+var2/n2));
		return score;

}

double findCutoffFromPWMscoreDist(map<float,float> scoredist)
{
	int i;
	float split=-1;
	double cutoff=0;
	int step=scoredist.size()/2;

	float startpoint= scoredist.begin()->first;
		while(step<scoredist.size()-1)
		{
			int leftmv=0;
			map<float,float>::iterator iter=scoredist.begin();
			map<float,float>::iterator iter2=scoredist.begin();
			//set to middle point
			for(i=0;i<step;i++)
			{
					iter++;
					iter2++;
			}
			split=iter->first;
			//get min of right
			float mindis=1000;
			while(iter2!=scoredist.end())
			{
				mindis=min(iter2->second,mindis);
				iter2++;
			}
			float mindis2=1000;
			while(iter->first!=startpoint)
			{
				leftmv++;
				mindis2=min(iter->second,mindis2);
				if(mindis2<=mindis)
					break;
				iter--;
			}
			if(iter->first==startpoint)
			{
				cutoff=iter->first;
				break;
			}
			else
			{	
				iter++;
				startpoint=iter->first;
				step=scoredist.size()-(scoredist.size()-step+leftmv)/2;
				
			}
			

		}
		return cutoff;

}
map<float,float> mergeBin(map<float,float> scoredist,map<float,int> scorelist,int mincnt)
{
map<float,float> newmergedist;
map<float,int>::reverse_iterator iter=scorelist.rbegin();
int waitmergeCnt=0;
float waitmergeMass=0;
while(iter!=scorelist.rend())
{
//if(iter->second<mincnt)
{
waitmergeCnt+=iter->second;
waitmergeMass+=iter->second*scoredist[iter->first];
}
if(waitmergeCnt>mincnt)
{
newmergedist[iter->first]=waitmergeMass/waitmergeCnt;
waitmergeCnt=0;
waitmergeMass=0;
}
iter++;
}

return newmergedist;

}

void FARCenterDistribtionScore2(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	map<float,float> scoredist;
	map<float,int> scoreCount;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)200/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
		
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 

double *bgCDDiff=new double[BINNUMBER-1]; 
double *bgCDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double pwmbiasZstat=0;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		//assuming at most one occurrence per sequence
		//maxSeqNum=5000;
/*******************************consider the PWM Score-distance Distribution*******************/	
		if(sortedThresh.size()>10)
		{
			for(int i=0;i<sortedThresh.size();i++)
			{
				scoredist[*ITER]=0;
					scoreCount[*ITER]=0;
					ITER++;
				
			}
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				int tt=abs(pos-center);
				float pwmscore=scorelist[i];
				scoreCount[pwmscore]++;
				scoredist[pwmscore]+=tt;
			}
			ITER=sortedThresh.begin();
			for(int i=0;i<sortedThresh.size();i++)
			{
				scoredist[*ITER]=scoredist[*ITER]/scoreCount[*ITER];
				ITER++;
			}
			scoredist=mergeBin(scoredist,scoreCount,sortedThresh.size()/BINNUMBER);
			if(scoredist.size()>4)
			pwmbiasZstat=PWMScoreBIASZstat(scoredist);
			else
				pwmbiasZstat=1.96;
			ITER=sortedThresh.begin();
			if(step>maxSeqNum)
			{
				step=maxSeqNum;
				for(int i=0;i<sortedThresh.size()-maxSeqNum;i++)
					ITER++;
				
			}
		}
		else
			pwmbiasZstat=1.96;
/*******************************consider the PWM Score-distance Distribution*******************/

		//if(sortedThresh.size()>3000)
		//	for(int i=0;i<sortedThresh.size()-3000;i++)
		//		ITER++;
		//getThresh=false;
		//scoredist=mergeBin(scoredist,scoreCount,50);
		//if(pwmbiasZstat!=1.96)
		//pwmbiasZstat=PWMScoreBIASZstat(scoredist);
		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				{
				CDDiff[i]=0;
				bgCDDiff[i]=0;
				}
				if(i<BINNUMBER-2)
				{
				CDDiff2[i]=0;
				bgCDDiff2[i]=0;
				}
			}
			if(getThresh)
				cutoff=*ITER;
			
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(coverSeq);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
					for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)floor(60.0/setting->resolution)];}
			//CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];
					smooth(CDHist,SmoothCDHist,BINNUMBER,slidewin);
		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				
				/**modified**/
				double temp=(SmoothCDHist[j]-SmoothCDHist[j+slidewin])/slidewin;
				//cout<<temp<<endl;
				if(temp<0)
				{	
					bgCDDiff[j]=-temp;
				   CDDiff[j]=0;
				}
				else
				{
					CDDiff[j]=temp;
				   	bgCDDiff[j]=0;
				}
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;
			//if(minbin<10)
			//{
			//	//break;
			//	setting->resolution=setting->resolution*2;
			//	BINNUMBER= ceil((double)SEQLEN/setting->resolution);
			//	BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
			//	slidewin=ceil((double)slidewin/2);

			//	continue;
			//}
			
	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/
ssr=0;
/**********use 0-order************/
			for(int j=1;j<BINNUMBER-2*slidewin;j++)
			{
				/**modified**/
				double temp=(SmoothCDHist[j]-2*SmoothCDHist[j+(int)slidewin]+SmoothCDHist[j+2*slidewin])/((slidewin)*(slidewin));//(CDDiff[j]-CDDiff[j+slidewin])/(slidewin);//
				if(temp<0)
				{	
					bgCDDiff2[j]=-temp;
				    CDDiff2[j]=0;
				}
				else
				{
					CDDiff2[j]=temp;
				   	bgCDDiff2[j]=0;
				}
				//cout<<temp<<endl;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
			break;
/**********use 0-order************/
			//int i=960/setting->resolution;
				for(int i=0;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
						bgdiff+=bgCDDiff[j];
						bgdiff2+=bgCDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++) ///zzz modify
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						{
							bgdiff+=CDDiff[j];
							cddiff+=bgCDDiff[j];
						}
						if(j<BINNUMBER-2)
						{
							bgdiff2+=CDDiff2[j];
							cddiff2+=bgCDDiff2[j];
						}
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=0.5;
				double ratio3=0.5;
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-slidewin));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2*slidewin));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)
				//{
				//	Z1=0;
				////if(Z2>2*Z0||Z2>2*Z1) ||Z1<0||Z2<0
				//	Z2=0;
				//}
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif>1)
				//	score=Z0+Z1;
				//score=Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 
				//score=Z2;

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					flushhistdata(CDHist,BINNUMBER);
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
				else
					up=!up;
				if(up)
				{
					
						step=step/2;
						for(int k=0;k<step;k++)
							ITER++;
				}
				else
				{
					step=step/2;
						for(int k=0;k<step;k++)
							ITER--;
				}

				if(!getThresh)
					break;
				//ITER++;
				//step--;
				//break;
				
		}while(step>1);
		//double pvalue=
		if(setting->N_motif==2)
			//bestscore+=bestZ1score;
		bestscore+=bestZ1score;
		if(setting->N_motif==3)
			//bestscore=bestZ2score;
			bestscore+=bestZ1score+pwmbiasZstat;
			int window=(bestwindowId+1)*bestwin;
			
			//double pvalue=1-normaldistribution(min(bestZ1score,min(bestZ2score,CDScore)));
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
			outputhistdata(setting);
			//if(pwmbiasZstat>=0)
			//cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<pwmbiasZstat<<"\t"<<bestZ1score<<endl; 
			//else
			//cout<<setting->inputFile<<"\t"<<0<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<pwmbiasZstat<<"\t"<<bestZ1score<<endl; 
			//delete[] CDHist;
			//delete[] CDDiff;
			//delete[] CDDiff;
			//delete[] bgCDDiff;
			//delete[] bgCDDiff2;
			//delete[] SmoothCDHist;
}


void FARCenterDistribtionScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	const double log10=log(10.0);
	char buffer[300];
	char *items[10];
	vector<float> scorelist;
	vector<int> positionlist;
	set<float> sortedThresh;
	vector<int> sameSeqPos;
	int SEQLEN=0;
	int maxSeqNum=0;
	int lastpos=-1000;
	int lastseq=-1;
	double sumSameSeq=0;

/***debug**/
	int slidewin=(int)400/setting->resolution;//setting->N_motif;
//setting->resolution=20;

	vector<float> problist;
	while(fgets(buffer,300,fp)){
		int seqid=0;
		int pos=0;

		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		//if(getThresh)
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		
		if(rmPoslist.size()>0)
		{
			bool rmflag=testInside(rmPoslist,setting->seedlength,(unsigned long)pos+seqid*MAXSEQLEN,setting);
			if(rmflag)
				continue;
		}
		//filter repeat
		if(lastseq!=seqid)
		{
			lastseq=seqid;
			lastpos=-1000;
#ifdef NORMALIZE
			for(int i=0;i<sameSeqPos.size();i++)
			{
				int p=sameSeqPos[i];
				problist[p]=problist[p]/sumSameSeq;
			}

			sameSeqPos.clear();
			sumSameSeq=0;
#endif
		}

		{
			float score;
			sscanf(items[2],"%f",&score);
			scorelist.push_back(score);
			sortedThresh.insert(score);
#ifdef NORMALIZE
			sameSeqPos.push_back(scorelist.size()-1);
			double temp=exp(score*log10);
			sumSameSeq+=temp;
			problist.push_back(temp);
#endif
		}

		if(pos>SEQLEN)
			SEQLEN=pos;
		positionlist.push_back(pos);
	}
		fclose(fp);
		if(positionlist.size()<setting->min_supp_ratio*maxSeqNum)
		{
			cout<<setting->inputFile<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl; 
			return;
		}
SEQLEN=ceil(SEQLEN/20.0)*20;
int BINNUMBER= ceil((double)SEQLEN/setting->resolution);
int BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
int center=SEQLEN/2;
	double* CDHist=new double[BINNUMBER];
	double* SmoothCDHist=new double[BINNUMBER];
double *CDDiff=new double[BINNUMBER-1]; 
double *CDDiff2=new double[BINNUMBER-2]; 
double bestscore=-1000000000000;
int bestwin=1;
double bestZ1score=-1000000000000;
double bestZ2score=-1000000000000;
double bestCutoff=0;
bool up=true;
double lastbestscore=-1000000000000;
double sumNorm=0;
				double CDScore=-1000000000000;
				int bestbgcnt,bestcdcnt,bestwindowId;


lastpos=10000;
	int step=sortedThresh.size();
		set<float>::iterator ITER=sortedThresh.begin();
		//if(sortedThresh.size()>3000)
		//	for(int i=0;i<sortedThresh.size()-3000;i++)
		//		ITER++;

		do
		{
			for(int i=0;i<BINNUMBER;i++)
			{
				CDHist[i]=0;
				SmoothCDHist[i]=0;
				if(i<BINNUMBER-1)
				CDDiff[i]=0;
				if(i<BINNUMBER-2)
				CDDiff2[i]=0;
			}
			if(getThresh)
				cutoff=*ITER;
			int count=0;
			int coverSeq=0;
			for(int i=0;i<positionlist.size();i++)
			{
				int pos=positionlist[i];
				if(scorelist[i]>=cutoff){
					int tt=abs(pos-center);
						tt=tt/BINSIZE2;
						if(tt>=BINNUMBER)
						tt=BINNUMBER-1;
						double step=1;
#ifdef NORMALIZE
						step=scorelist[i];//exp(scorelist[i]*log10);
			
				sumNorm+=step;
				count++;
				if(lastpos>pos)
				{
					coverSeq++;	
				}

				lastpos=pos;
#endif

				CDHist[tt]+=step;
				}
			}
#ifdef NORMALIZE

			sumNorm=sumNorm/(count);
			//double s=0;
			for(int j=0;j<BINNUMBER;j++)
			{
				CDHist[j]/=sumNorm;

			}
#endif
/**modified**/
					for(int i=0;i<60/setting->resolution;i++)
			{CDHist[BINNUMBER-1-i]=CDHist[BINNUMBER-1-(int)floor(60.0/setting->resolution)];}
			//CDHist[BINNUMBER-1]=CDHist[BINNUMBER-2];

		double ssr=0;
			double minbin=1000;
			for(int j=0;j<BINNUMBER-slidewin;j++)
			{
				/**modified**/
				double temp=(CDHist[j]-CDHist[j+slidewin])/slidewin;
				//cout<<temp<<endl;
				if(temp<0)
					temp=0;
				CDDiff[j]=temp;
				ssr+=CDDiff[j];
				if(CDHist[j+1]<minbin)
					minbin=CDHist[j+1];
			}
			if(BINNUMBER<10||ssr==0)
					break;
			if(minbin<10)
			{
				setting->resolution=setting->resolution*2;
				BINNUMBER= ceil((double)SEQLEN/setting->resolution);
				BINSIZE2=ceil((double)SEQLEN/2/(BINNUMBER));
				slidewin=ceil((double)slidewin/2);

				continue;
			}
			
	/**********use absolute 1-order************/
			//for(int j=0;j<BINNUMBER-2;j++)
			//{
			//	CDDiff2[j]=fabs(CDDiff[j+1]-CDDiff[j]);
			//	ssr+=CDDiff2[j];
			//}
	/**********use absolute 1-order************/
ssr=0;
/**********use 0-order************/
			for(int j=0;j<BINNUMBER-2*slidewin;j++)
			{
				/**modified**/
				double temp=(CDDiff[j]-CDDiff[j+slidewin])/(slidewin);//(CDHist[j]-2*CDHist[j+(int)slidewin/2]+CDHist[j+slidewin]);
				if(temp<0)
					temp=0;
				CDDiff2[j]=temp;
				//cout<<temp<<endl;
				ssr+=CDDiff2[j];
			}
			if(ssr==0)
			break;
/**********use 0-order************/
			//int i=400/setting->resolution;
				for(int i=0;i<BINNUMBER/2;i++)
				{
					double cdcnter=0;
					double bgcnter=0;
					double cddiff=0;
						double cddiff2=0;
					double bgdiff=0;
						double bgdiff2=0;
					double score=0;
					for(int j=0;j<i+1;j++)
					{
						cdcnter+=CDHist[j];
						cddiff+=CDDiff[j];
						cddiff2+=CDDiff2[j];
					}
					for(int j=i+1;j<BINNUMBER;j++)
					{
						bgcnter+=CDHist[j];
						if(j<BINNUMBER-1)
						bgdiff+=CDDiff[j];
						if(j<BINNUMBER-2)
						bgdiff2+=CDDiff2[j];
					}

				  
				double ratio=(double)(i+1)/BINNUMBER;
				double ratio2=(double)(i+1)/(BINNUMBER-slidewin+1);
				double ratio3=(double)(i+1)/(BINNUMBER-2*slidewin+2);
				double n=bgcnter+cdcnter;
				double bias1=fabs(cdcnter/(i+1)-bgcnter/(BINNUMBER-i-slidewin));
				double bias2=fabs(cddiff/(i+1)-bgdiff/(BINNUMBER-i-2*slidewin));
				//cddiff-=bias1/2;
				//cddiff2-=bias2/2;
				double n2=bgdiff+cddiff;
				double n3=bgdiff2+cddiff2;
				double Z0=(cdcnter-n*ratio)/sqrt(n*ratio*(1-ratio));
				double Z1=(cddiff-n2*ratio2)/sqrt(n2*ratio2*(1-ratio2));
				double Z2=(cddiff2-n3*ratio3)/sqrt(n3*ratio3*(1-ratio3));
/**modified**/
				//if(Z0<3)
				//{
				//	Z1=0;
				////if(Z2>2*Z0||Z2>2*Z1) ||Z1<0||Z2<0
				//	Z2=0;
				//}
/**modified**/

				score=Z0; //+ssr*ssr
				//if(setting->N_motif==2)
				//	score=Z0+Z1;
				//if(setting->N_motif==3)
				//	score=Z0+Z1+Z2; 

					if(score>=bestscore&&cdcnter>50)
					{
						bestwindowId=i;				
						bestwin=setting->resolution;
						//cout<<(bestwindowId+1)*bestwin<<"\t"<<Z0<<"\t"<<Z1<<"\t"<<Z2<<endl;
						bestscore=score;
						CDScore=Z0;
						bestZ1score=Z1;
						bestZ2score=Z2;
						bestbgcnt=bgcnter;
						bestcdcnt=cdcnter;
						
					}
			
				
				}
				if(bestscore>lastbestscore)
				{
					lastbestscore=bestscore;
					bestCutoff=cutoff;	

				}
				else
					up=!up;
				if(up)
				{
						step=step/2;
						for(int k=0;k<step;k++)
							ITER++;
				}
				else
				{
					step=step/2;
						for(int k=0;k<step;k++)
							ITER--;
				}

				if(!getThresh)
					break;
				//ITER++;
				//step--;
				//break;
				
		}while(step>1);
		//double pvalue=
		if(setting->N_motif==2)
			bestscore+=bestZ1score;
		////bestscore=bestZ1score;
		if(setting->N_motif==3)
		//	//bestscore=bestZ2score;
			bestscore+=bestZ1score+bestZ2score;
			int window=(bestwindowId+1)*bestwin;
			//double pvalue=1-normaldistribution(min(bestZ1score,min(bestZ2score,CDScore)));
			cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<"\t"<<bestCutoff<<"\t"<<CDScore<<"\t"<<bestZ1score<<endl; 
}

void SimplePeakEnrichmentScore(PARAM * setting)
{
	FILE* fp = fopen(setting->inputFile.c_str(),"r");
	int binNum=(setting->PeakRange+setting->resolution-1)/setting->resolution;
	int* CDHist=new int[binNum];
	int maxSeqNum=0;
	for(int i=0;i<binNum;i++)
		CDHist[i]=0;
	char buffer[30];
	char *items[10];
	int center=setting->PeakRange/2;
	while(fgets(buffer,30,fp)){
		int seqid=0;
		int pos=0;
		int n = split_string(&buffer[0],items,10);
		if(n<2)
			continue;
		sscanf(items[0],"%d",&seqid);
		sscanf(items[1],"%d",&pos);
		if(maxSeqNum<seqid)
			maxSeqNum=seqid;
		CDHist[abs(pos-center)/(setting->resolution/2)]++;
	}
		fclose(fp);
		int windowsize=100;
		double CDScore;
		int bestbgcnt,bestcdcnt,bestwindowId;
		double bestscore=0;
	    for(int i=0;i<binNum/5;i++)
		{
			int j=0;
			double cdcnter=0;
			double bgcnter=0;
			
			double score=0;
			for(j=0;j<i+1;j++)
				cdcnter+=CDHist[j];
			for(j=i+1;j<binNum;j++)
				bgcnter+=CDHist[j];
			//only average on bg
			double avgvalue=(double)(bgcnter/(binNum-i-1));
			double ssr=0;
			for(j=i+1;j<binNum;j++)
					ssr+=(CDHist[j]-avgvalue)*(CDHist[j]-avgvalue);
			ssr=sqrt( (double)(ssr/(binNum-i-1)));
			bgcnter=ssr+avgvalue; //1.95996*

			score=(double)((cdcnter/(i+1))/bgcnter);
			bgcnter-=ssr;  //1.95996*
			bgcnter*=(binNum-i-1);
			windowsize=(i+1)*(setting->PeakRange/binNum);

			if(score>=bestscore&&cdcnter>maxSeqNum*setting->min_supp_ratio)
			{
				bestwindowId=i;		
				bestscore=score;
				CDScore=cdcnter;
				bestbgcnt=bgcnter;
				bestcdcnt=cdcnter;
			}
		
		}
		int window=bestwindowId*setting->resolution;
		cout<<setting->inputFile<<"\t"<<bestscore<<"\t"<<window<<endl;
}
int main(int argc, char* argv[])
{

	PARAM * setting=read_parameters (argc, argv);
	unsigned long int aa=1;
	//double data[10];
	//double* cc=new double[10];
	//for(int i=0;i<10;i++)
	//{
	//	data[i]=1+0.2*i;
	//}
	//smooth(data,cc,10);
	//for(int i=0;i<10;i++)
	//{
	//	cout<<cc[i]<<endl;
	//}
	
	//cout<<sizeof(aa)*8<<endl;
	if(setting->weightFile!="")
	{
		LoadrmPoslist(setting);
	}

	//SimplePeakEnrichmentScore(setting);
FARCenterDistribtionScore2(setting);
//UDSmoothCenterDistribtionScore(setting);
//FARSmoothCenterDistribtionScore(setting);
//UDSmoothCenterDistribtionScore(setting);
//FARSmoothQuantiScore(setting);
	 //BasicQuantiScore(setting);
	// BasicPlusQuantiScore(setting);
	//FixWindScore(setting);
	//UseredBGScore(setting);
	//FARSmoothCenterDistribtionScore(setting);
	//PeakShapeScore(setting);
	//WholeShapeScore(setting);
//cin>>aa;
	return 0;
	
	
}
