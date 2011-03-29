#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
int maxw;
typedef struct match_st{
	int seqid;
	float score;
	int pos;
} MATCH;

typedef struct matrix_st{
	char* mname;
	int mlen;
} MAT;
int partition(int* list,int left,int right,int pivotIndex){
	int pivotValue=list[pivotIndex];
	int temp = list[pivotIndex];
	list[pivotIndex]=list[right];
	list[right]=temp;
	int storeIndex=left;
	int i;
	for (i=left;i<right;i++){
		if (list[i]<pivotValue){
			temp=list[storeIndex];
			list[storeIndex]=list[i];
			list[i]=temp;
			storeIndex+=1;
		}
	}
	temp=list[storeIndex];
	list[storeIndex]=list[right];
	list[right]=temp;
	return storeIndex;
}

int findFirstK(int* list,int left,int right,int k){
	int pivotNewIndex;
	int i;
	if (right>left){
		pivotNewIndex=partition(list,left,right,(left+right)/2);
		if (pivotNewIndex>k){
			return findFirstK(list,left,pivotNewIndex-1,k);
		}
		else if (pivotNewIndex<k){
			return findFirstK(list,pivotNewIndex+1,right,k);
		}
		else return list[pivotNewIndex];
	}
}

int findFirstKiter(int* list,int left,int right,int k){
	int pivotNewIndex;
	int Left=left;
	int Right=right;
	int i;
	while (right>left){
		pivotNewIndex=partition(list,left,right,(left+right)/2);
		if (pivotNewIndex>k){
			right=pivotNewIndex-1;
		}
		else if (pivotNewIndex<k){
			left=pivotNewIndex+1;
		}else{
			break;
		}
	}
	return list[k];
}
/* A Quicksort for structures of type address. */
int sort_MATCH_pos(MATCH items[], int left, int right)
{
	register int i, j;
	//type of key to sort
	int x;
	MATCH temp;
	i = left; j = right;
	x = items[(left+right)/2].pos;
	do {
		while((items[i].pos<x) && (i < right)) i++;
		while((items[j].pos>x) && (j > left)) j--;
		if(i <= j) {
			temp = items[i];
			items[i] = items[j];
			items[j] = temp;
			i++; j--;
		}
	} while(i <= j);

	if(left < j) sort_MATCH_pos(items, left, j);
	if(i < right) sort_MATCH_pos(items, i, right);
}

/***********************************************/
void check (void *x)
/***********************************************/
{
	if (!x){
		fprintf(stderr,"Out of memory\n");
		exit(0);
	}
}

void str_replace(char * from, char * to, char * s,char *dest){
	char*t4=0;
	int len_dest=0;
	while(strstr(s,from)){
		t4=strstr(s,from);
		strncpy(dest+len_dest,s,t4-s);
		len_dest+=t4-s;
		dest[len_dest]=0;
		strcat(dest,to);
		len_dest+=strlen(to);
		t4+=strlen(from);
		s=t4;
	}
	if (!t4) strcpy(dest,s);
	else strcat(dest,t4);
}
int get_num_lines(char* f){
	char s[strlen(f)+100];
	sprintf(s,"wc -l \"%s\"",f);
	char s2[strlen(f)+100];
	str_replace("$","\\$",s,s2);
	FILE* fin=popen(s2,"r");
	int lines;
	fscanf(fin,"%d",&lines);
	fclose(fin);
	return lines;
}

int read_pos_file(char* fname,MATCH** Mp,int L,int mlen){
	int num_matches=get_num_lines(fname);
	int x,y;
	float z;
	check(*Mp = (MATCH*)malloc(num_matches*sizeof(MATCH)));
	int i=0;
	FILE* f=fopen(fname,"r");
	while (fscanf(f,"%d%d%f",&x,&y,&z)==3){
		(*Mp)[i].seqid=x;
		(*Mp)[i].score=z;
		(*Mp)[i++].pos=y-L/2+mlen;
	}
	fclose(f);
	return num_matches;
}


int read_matrix_file(char* fname,MAT** Mp){
	int num_lines=get_num_lines(fname);
	int x,y;
	check(*Mp = (MAT*)malloc(num_lines*sizeof(MAT)));
	int i=0;
	char  buffer[500];
	FILE* f=fopen(fname,"r");
	while (fscanf(f,"%s%d",&buffer,&y)==2){
		(*Mp)[i].mname=(char*)malloc(strlen(buffer)+1);
		strcpy((*Mp)[i].mname,buffer);
		(*Mp)[i].mlen=y;
		//fprintf(stderr,"%s\n%d\n",(*Mp)[i].mname,(*Mp)[i].mlen);
		i++;
	}
	fclose(f);
	return num_lines;
}

int count(MATCH* Mp,int num_matches,int*bgdist,int*fgdist,int N,int L,int w){
	for (int i=0; i<num_matches;i++){
		if (-w*3/2<=Mp[i].pos && Mp[i].pos<w*3/2){
			if (-w/2<=Mp[i].pos && Mp[i].pos<w/2){
				fgdist[Mp[i].seqid-1]+=1;
			}
			else{
				bgdist[N*(Mp[i].pos>0)+Mp[i].seqid-1]+=1;
			}
		}
	}
}
int computescore(int* fgcounts,int* bgcounts){
	int maxscore=0;
	int minscore=0;
	int cumsum=0;
	for (int i=maxw;i>=0;i--){
		cumsum+=fgcounts[i]-bgcounts[i];
		if (cumsum>maxscore){
			maxscore=cumsum;
		}
		if (cumsum<minscore){
			minscore=cumsum;
		}
	}
	return maxscore+minscore;
}
int countnum(int* fg,int fglen,int* counts){
	memset(counts,0,(maxw+1)*sizeof(int));
	for (int i=0;i<fglen;i++){
		counts[fg[i]]+=1;
	}
}
int* randindex;
int ri=0;
int maxri=0;
int genbgcounts(int* bgdist,int bgdistlen,int bglen,int* counts){
	memset(counts,0,(maxw+1)*sizeof(int));
	for (int j=0;j<bglen;j++){
		counts[bgdist[randindex[ri+j]]]+=1;
	}
}

int main(int argc,char** argv){
	if (argc==1){
		printf("arguments: numseq seqlen inputdir");
		exit(1);
	}
	int randtrials=101;
	int N=atoi(argv[1]);
	int L=atoi(argv[2]);
	char* inputdir=argv[3];
	maxw=(L-100)/3;
	FILE* f;
	MATCH* Mp;
	MAT* Matp;
	char matfile[strlen(inputdir)+50];
	sprintf(matfile,"%s/matlist.txt",inputdir);
	char resultsfile[strlen(inputdir)+50];
	sprintf(resultsfile,"%s/results.txt",inputdir);
	FILE* fresults=fopen(resultsfile,"w");

	int num_mat=read_matrix_file(matfile,&Matp);
	//printf("%d\t%d\n",Mp[num_matches-1].seqid,Mp[num_matches-1].pos);
	int bgdist[2*N];
	int fgdist[N];
	int bgcounts2[maxw+1];
	int bgcounts[maxw+1];
	
	int fgcounts[maxw+1];
	int bg[N];
	time_t start=clock();
	int i,w,j;
	int Matpi;
	//fprintf(stderr,"generating 1000*N random numbers\n");
	randindex=(int*)malloc(1000*N*sizeof(int));
	maxri=1000*N;
	for (i=0;i<maxri;i++){
		randindex[i]=rand()%(2*N);
	}
	//fprintf(stderr,"time elapsed=%lf\n",1.0*(clock()-start)/CLOCKS_PER_SEC);
	//fprintf(stderr,"counting\n");

	int sumscore=0;
	int sumscoresq=0;
	int minscore=999999999;
	double meanscore;
	double sdscore;
	double zscore;
	int pscore;
	int score=0;
	int score1,score2;
	int scores[randtrials];
	int st=0;
	int bestw;
	double wmean[1000];
	int wi;
	double wsd[1000];
	fprintf(fresults,"ID\tNAME\tENRICHMENT\tZSCORE\tRANDOMMEAN\tRANDOMERROR\tHALFWINDOW\n");
	for (Matpi=0;Matpi<num_mat;Matpi++){
		char posfile[strlen(inputdir)+strlen(Matp[Matpi].mname)+50];
		sprintf(posfile,"%s/%s.pos",inputdir,Matp[Matpi].mname);
		int mlen=Matp[Matpi].mlen;
		fprintf(fresults,"%d\t%s\t",Matpi,Matp[Matpi].mname);
		fprintf(stderr,"%d\t%s\t",Matpi,Matp[Matpi].mname);
		int num_matches=read_pos_file(posfile,&Mp,L,mlen); //create Mp
		double bestmeanscore=-99999;
		double bestzscore=-999999;
		int bestscore=-999999;
		double bestsdscore=1;
		int bestw=50;
		int score0;
		int score1;
		int score;
		wi=0;
		for(w=50;w<=maxw;w=ceil(w*1.25/50)*50){
			memset(fgdist,0,N*sizeof(int));
			memset(bgdist,0,2*N*sizeof(int));
			count(Mp,num_matches,bgdist,fgdist,N,L,w);
			countnum(fgdist,N,fgcounts);
			countnum(bgdist,N,bgcounts);
			countnum(bgdist+N,N,bgcounts2);
			if (bgcounts[0]==N){
				bgcounts[0]=N-1;
				bgcounts[1]=1;
			}
			if (bgcounts2[0]==N){
				bgcounts2[0]=N-1;
				bgcounts2[1]=1;
			}
			score0=computescore(fgcounts,bgcounts);
			score1=computescore(fgcounts,bgcounts2);
			score=score0<score1?score0:score1;
			sumscore=0;
			sumscoresq=0;
			for (i=0; i<randtrials;i++){
				ri=(ri+N)%maxri;
				genbgcounts(bgdist,2*N,N,fgcounts);
				score0=computescore(fgcounts,bgcounts);
				score1=computescore(fgcounts,bgcounts2);
				scores[i]=score0<score1?score0:score1;
				sumscore+=scores[i];
				sumscoresq+=scores[i]*scores[i];
			}
			meanscore=sumscore*1.0/randtrials;
			sdscore=sqrt((sumscoresq*1.0/randtrials-meanscore*meanscore)*randtrials/(randtrials-1));
			sdscore=sdscore>wsd[wi-1]?sdscore:wsd[wi-1];
			sdscore=sdscore>10?sdscore:10;
			wmean[wi]=meanscore;
			wsd[wi]=sdscore;
			wi++;
			meanscore=meanscore>0?meanscore:0;
			zscore=(score-meanscore)/sdscore;
			if (zscore>bestzscore) {
				bestzscore=zscore;
				bestmeanscore=meanscore;
				bestsdscore=sdscore;
				bestscore=score;
				bestw=w;
			}
		}
		fprintf(fresults,"%d\t%f\t%f\t%f\t%d",bestscore,bestzscore,bestmeanscore,bestsdscore,bestw/2);
		wi=0;
		/*for(w=50;w<=maxw;w+=50){
			fprintf(fresults,"\t%f\t%f",wmean[wi],wsd[wi]);
			printf("\t%f\t%f",wmean[wi],wsd[wi]);
			wi++;
		}*/
		fprintf(fresults,"\n");
		fprintf(stderr,"%d\t%f\t%f\t%f\t%d",bestscore,bestzscore,bestmeanscore,bestsdscore,bestw/2);
		fprintf(stderr,"\n");
		free(Mp);
	}
	fclose(fresults);
	//printf("time elapsed=%lf\n",1.0*(clock()-start)/CLOCKS_PER_SEC);
	return 0;
}
