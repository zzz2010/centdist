#pragma once
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <map>
#include "binomialdistr.h"
#include "blas.h"
#include "trinverse.h"
#include "cholesky.h"
#include "spdsolve.h"
#include "lbfgs.h"
#include "minlm.h"
#include <sstream>
#include <set>


#include <algorithm>

//#define FPMIN 1.0E-300
//#define EPS 1.0E-200
//#define MAXIT 100
//#define pi 3.141592653587932384626
//#define SQR(X) ((X)*(X))
//#define CUBE(X) ((X)*(X)*(X))
//#define MAX(X,Y) ( ((X)>(Y))? (X): (Y) )
//#define MIN(X,Y) ( ((X)>(Y))? (Y): (X) )
//#define MAXPWMLEN 50
//#define MAXSEQLEN 10000


using namespace std;

#define DEBUG 0
//#define Setting->seedlength 6
#define SWITCHLEN 2000
#define DNA_N_ONLY

#define VAL      long int
#define nextDouble (double)(rand()%1000)/1000

#ifndef isnan(x)
#define isnan(x) ((x) != (x))
#endif
#define MINSCORE -1000000


#ifndef NOTACGT(x)
#define NOTACGT(x) (x>3||x<0)
#endif

//#ifndef min
//#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
//#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
//#endif

	//int seedlength=6;
	//float min_supp_ratio=0.5;
	//float FDRthresh=1.e-7;
	//float olThresh=0.02;
	//float pdThresh=0.18;
	//long max_motif_length=26;
typedef struct param_st{
	//char inputFile[300];
	//char outputDIR[300];
	//char weightFile[300];
	string inputFile;
	string outputDIR;
	string weightFile;
	int N_motif;
	int seedlength;
	int resolution;
	int startwSize;
	int PeakRange;
	float min_supp_ratio;
	float FDRthresh;
	float olThresh;
	float pdThresh;
	long max_motif_length;
}PARAM;

static inline int acgt(char w)
{
	switch(w)
	{
	case 'a': return 0;
	case 'c': return 1;
	case 'g': return 2;
	case 't': return 3;
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'T': return 3;
	case 'R': return 4;
	case 'Y': return 5;
	case 'K': return 6;
	case 'M': return 7;
	case 'S': return 8;
	case 'W': return 9;
	case 'B': return 10;
	case 'D': return 11;
	case 'H': return 12;
	case 'V': return 13;
	case 'N': return 14;
	case 'n': return 14;
	case 'X': return 15;
		//RYKMSWBDHVN
	default: return -1;
	}

}

inline void smooth(double* data,double* c, int length,int span=5)
{

	int width=span-1+span%2; // force it to be odd
	if(width<3)
	{
		//do nothing...
		for(int j=0;j<length;j++)
			c[j]=data[j];
		return;
	}
    double* cbegin=new double[width-2];
	double* cend=new double[width-2];
	int head=(width-1)/2;
	for(int j=0;j<length;j++)
		 c[j]=0;
	for(int j=0;j<length-width+1;j++)
	{
		for(int kk=j;kk<j+width;kk++)
		c[j+head]+=data[kk];
		
		c[j+head]=c[j+head]/width;
	
	}
	for(int i=0;i<width-2;i++)
	{
		cbegin[i]=0;
		cend[i]=0;
	}
	for(int i=0;i<width-2;i++)
	{
		for(int j=i;j<width-2;j++)
		{
			cbegin[j]+=data[i];
			cend[j]+=data[length-i-1];
		}
		
	}
	for(int i=0;i<head;i++)
	{
		cbegin[i*2]=cbegin[i*2]/(1+i*2);
		c[i]=cbegin[i*2];

		cend[width-3-i*2]=cend[width-3-i*2]/(width-2-i*2);
		c[length-head+i]=cend[width-3-i*2];
	}
	
	delete[] cbegin;
	delete[] cend;

}

inline double binominal(double p,int k , int n)
{
	double prob=0;
	for(int i=0;i<k;i++)
	{
		prob+=log((double)n-i)-log((double)i+1);
	}
	prob+=k*log(p);
	prob+=(n-k)*log(1-p);
	return exp(prob);
}



 inline double getZscore(double val,double E, double var)
{
	//double temp=val-E;
	//temp=temp/sqrt(var);
	//return temp;
	return val/var;
}

inline long getReverseComplementHashing(long hash, int len)
{
	long ret=0;
	for(int i = 0; i < len; ++i)
	{
		ret<<=2;
		long temp=hash%4;
		hash>>=2;
		ret+=3-temp;
		
		
	}

	return ret;
}



inline long getHashing(const char* source, int start, int len)
{
		int i=0;
		long ret=0;
		
		ret = acgt( source[start]);
		if(ret<0)
			return -1;
		for (i = start + 1; i < start + len; i++)
		{ 
			ret<<=2;
			long temp=acgt(source[i]);
			if(temp<0||temp>3)
				return -1;
			ret+=temp;	

			
		}

			
			return ret;

}



inline char reverseC(char w)
{
	switch(w)
	{
	case 'a': return 't';
	case 'c': return 'g';
	case 'g': return 'c';
	case 't': return 'a';
	case 'A': return 'T';
	case 'C': return 'G';
	case 'G': return 'C';
	case 'T': return 'A';
	case 'N': return 'N';
	case 'n': return 'N';
	case 'R': return 'Y';
	case 'Y': return 'R';
	case 'W': return 'S';
	case 'S': return 'W';
	case 'V': return 'B';
	case 'B': return 'V';
	case 'K': return 'M';
	case 'M': return 'K';
	case 'H': return 'D';
	case 'D': return 'H';
	
	}
}

inline string reverseString(string Tag)
{
		string temp="";
		for(int i=0;i<Tag.size();i++)
		temp.push_back(reverseC(Tag[Tag.size()-1-i]));
		//temp.push_back(reverseC(Tag[i]));
		return temp;
}


inline string Int2String(int num)
{
	string s;
	stringstream ss(s);
	ss<<num;
	return ss.str();
}

inline VAL gcd(VAL  a, VAL  b)
{
     if  (b  ==   0 ) 
         return  a;
     else
         return  gcd(b, a % b);
}

inline VAL  C( int  a,  int  b)
{
    VAL  x = 1 , y = 1 ;
    VAL  r;
     int  i;
     for  (i = 0 ; i < b; i ++ )
     {
        x  *=  (a - i); 
        y  *=  (b - i);
        r  =  gcd(x, y);
        x  /=  r;
        y  /=  r;
    }
     return  x  /  y;
} 



//inline long double loggamma(long double x){
//	long double a1,a2,a3,a4,a5,f,xx,z;
//	a1=log(2*pi)/2;
//	a2=1.0/1680;
//	a3=1.0/1260;
//	a4=1.0/360;
//	a5=1.0/12;
//	f=1;
//	if (x<=0) return -1;
//	xx=x;
//	while (xx<7){
//		f*=xx;
//		xx+=1;
//	}
//	f=-log(f);
//	z=1.0/(xx*xx);
//	return (f+(xx-0.5)*log(xx)-xx+a1+(((-a2*z+a3)*z-a4)*z+a5)/xx);
//}
//
//inline long double betacf(int a,int b,long double x){
//	int m,m2;
//	long double qab,qap,qam,c,d,h,aa,delx;
//	qab=a+b;
//	qap=a+1.0;
//	qam=a-1.0;
//	c=1.0;
//	d=1.0-1.0*qab*x/qap;
//	if (fabs(d)<FPMIN) 
//		d = FPMIN;
//	d=1.0/d;
//	h=d;
//	m=1;
//	while(m<=MAXIT){
//		m2=2*m;
//		aa=1.0*m*(b-m)*x/((qam+m2)*(a+m2));
//		d=1.0+aa*d;
//		if (fabs(d)<FPMIN) 
//			d=FPMIN;
//		c=1.0+1.0*aa/c;
//		if (fabs(c)<FPMIN) 
//			c=FPMIN;
//		d=1.0/d;
//		h*=d*c;
//		aa=-1.0*(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
//		d=1.0+aa*d;
//		if (fabs(d) < FPMIN) 
//			d=FPMIN;
//		c=1.0+aa/c;
//		if (fabs(c) < FPMIN)
//			c=FPMIN;
//		d=1.0/d;
//		delx=d*c;
//		h*=delx;
//		if (fabs(delx-1.0)<EPS)
//			break;
//		m+=1;
//	}
//	return h;
//
//}
//
//inline long double betai(int a,int b,long double x){
//	long double bt;
//	if (a==0) 
//
//		return 1;
//	if (b<=0) 
//
//		return 0;
//	if (x<0.0 || x>1.0) 
//		
//		return -1;
//	if (x==0.0 || x==1.0) 
//		bt=0.0;
//	else
//		bt=exp(loggamma(a+b)-loggamma(a)-loggamma(b)+a*log(x)+b*log(1-x));
//
//	if (x<(a+1.0)/(a+b+2.0))
//		return bt*betacf(a,b,x)/a;
//	else
//		return 1.0-bt*betacf(b,a,1.0-x)/b;
//}


 inline long double binominalCDF(double p,int k , int n)
{
		long double prob=0;
		//for(int i=0;i<k+1;i++)
		//	prob+=binominal(p,i,n);
		//k=n-k;		
		prob=1.0-binomialcdistribution(k,n,p); //betai( k, n-k+1,p);
		return prob;
}

  inline long double binominalTail(double p,int k , int n)
{
		long double prob=0;
		//for(int i=0;i<k+1;i++)
		//	prob+=binominal(p,i,n);
		//k=n-k;		
		prob=binomialcdistribution(k,n,p); //betai( k, n-k+1,p);
		return prob;
}

  inline int invbinominalCDF(double y,double p, int n,int maxcount)
{
		long double prob=0;
		//for(int i=0;i<k+1;i++)
		//	prob+=binominal(p,i,n);
		//k=n-k;		
		double key=y;
		int last=maxcount;
		int first=0;
		   while (first <= last) {
			   int mid = (first + last) / 2;  // compute mid point.
			 double  prob=binomialcdistribution(mid,n,p); //betai( k, n-k+1,p);
			   if (key <  prob) 
				   first = mid + 1;  // repeat search in top half.
			   else if (key >  prob) 
				   last = mid - 1; // repeat search in bottom half.
			   else
				   return mid;     // found it. return position /////
		   }
		   return (first + 1);    // failed to find key

		
		
}

  inline bool KSTest(vector<VAL> cdpos,int seqlen)
{
	int i;
	for(i=0;i<cdpos.size();i++)
		cdpos[i]=cdpos[i]%seqlen;
	sort( cdpos.begin(),cdpos.end());
	//cout<<*(cdpos.end()-1)<<endl;
	double score=0;
	for(i=0;i<cdpos.size();i++)
	{
		double temp= fabs((double)i/cdpos.size()-cdpos[i]/(seqlen));
		if(temp>score)
			score=temp;
	}

		if(score>0.3288)
			return true;
		return false;
}

///fit exponential curve
   inline double fitParameter(double* data,int num,double & snratio,double & probE, double& noiseP)
   {
		 bool result;
    bool waserrors;
    bool referror;
    bool lin1error;
    bool lin2error;
    bool eqerror;
    bool converror;
    bool scerror;
    
    int nset;
    int n;
    int m;
    ap::real_1d_array x;
    ap::real_1d_array xe;
    ap::real_1d_array b;
    int i;
    int j;
    double v;
    ap::real_2d_array a;
    lmstate state;
    lmreport rep;

    waserrors = false;
    referror = false;
    lin1error = false;
    lin2error = false;
    eqerror = false;
    converror = false;
    scerror = false;
	    //
    // Reference problem.
    // RKind is a algorithm selector:
    // * 0 = FJ
    // * 1 = FGJ
    // * 2 = FGH
    //
	int rkind=2;

	x.setbounds(0, 2);
    n = 3;
    m = 3;

	            x(0) = 5.64e-02;
            x(1) = 0.883;
            x(2) = 0.407;
		if( rkind==0 )
        {
            minlmfj(n, m, x, 0.0, 0.0, 400, state);
        }
        if( rkind==1 )
        {
            minlmfgj(n, m, x, 0.0, 0.0, 400, state);
        }
        if( rkind==2 )
        {
            minlmfgh(n, x, 0.0, 0.0, 400, state);
        }

		 while(minlmiteration(state))
        {
			    if (state.x(0) < 0)
                    state.x(0) = 0.000000001;
                if (state.x(1) < 0)
                    state.x(1) = 0.000000001;
                if (state.x(2) < 0)
                    state.x(2) = 0.000000001;
                if (state.x(0) >1)
                    state.x(0) = 0.999999999;
                if ((state.x(2)/(1-state.x(0))) > 1)
                    state.x(2) = (1-state.x(0))-0.000000001;
            
            //
            // (a*exp(-b*x)+c-y)^2
            //
            state.f =0;
                for (int dd = 1; dd < num+1; dd++)
			    {
    			    double temp= state.x(0)*exp(-state.x(1)*dd)+state.x(2) -data[dd-1];
                    state.f += temp * temp;
			    }
            if( state.needfg||state.needfgh )
            {
				  state.g(0)= state.g(1)= state.g(2)=0;
                    for (int dd = 1; dd < num+1; dd++)
                    {
                        double temp = state.x(0) * exp(-state.x(1) * dd) + state.x(2) ;
                        double temp2 = temp - data[dd-1];
                        
                        temp=2*temp2;
                        state.g(0) += temp * (exp(-state.x(1) * dd));
                        state.g(1) += temp * (exp(-state.x(1) * dd)) * (-dd) * state.x(0);
                        state.g(2) += temp;
                    }
            }
            if( state.needfij )
            {
                state.fi(0) = state.x(0)-2;
                state.fi(1) = state.x(1);
                state.fi(2) = state.x(2)-state.x(0);
                state.j(0,0) = 1;
                state.j(0,1) = 0;
                state.j(0,2) = 0;
                state.j(1,0) = 0;
                state.j(1,1) = 1;
                state.j(1,2) = 0;
                state.j(2,0) = -1;
                state.j(2,1) = 0;
                state.j(2,2) = 1;
            }
            if( state.needfgh )
            {
					 for (int dd = 1; dd < num+1; dd++)
                    {
						double exp_bx2 = ap::sqr((exp(-state.x(1) * dd)));
                        double exp_bx = exp(-state.x(1) * dd);
                        double x2=dd*dd;
                        double T = state.x(0) * exp_bx + state.x(2) - data[dd-1];
                        state.h(0, 0) += 2 * exp_bx2;
                        state.h(0, 1) += 2 * (-state.x(0) * (-dd) * exp_bx2 + T * (-dd) * exp_bx);
                        state.h(0, 2) += 2*exp_bx;
                        state.h(1, 0) += 2 * (-state.x(0) * (-dd) * exp_bx2 + T * (-dd) * exp_bx);
                        state.h(1, 1) += 2 * (state.x(0) * state.x(0)*x2*exp_bx2+state.x(0) *x2*exp_bx);
                        state.h(1, 2) += 2*(-state.x(0) * (-dd) * exp_bx);
                        state.h(2, 0) += 2 * exp_bx;
                        state.h(2, 1) += 2 * (-state.x(0) * (-dd) * exp_bx);
                        state.h(2, 2) += 2;
                    }
            }
          /* scerror = scerror||!rkindvsstatecheck(rkind, state);*/
        }
        minlmresults(state, x, rep);
        referror = referror||rep.terminationtype<=0;
		snratio=state.x(0);
		probE=exp(-state.x(1));
		noiseP=state.x(2)/(1-state.x(0));
		if( referror )
			return -1;
		 
		return  state.f/num;
   
   }