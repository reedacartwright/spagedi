#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Xatools.h"

char wait_a_char(); // defined in Autocio.c

#define NR_END 1
#define FREE_ARG char* 
#ifndef ERRORFILE
	#define ERRORFILE "error.txt"
#endif
#ifndef MISSING
	#define MISSING -999.99F
#endif

#ifndef MISSVAL
	#define MISSVAL HUGE_VAL
#endif

#if defined(CURSES_HAVE_CURSES_H)
#	include <curses.h>
#elif defined(CURSES_HAVE_NCURSES_H)
#	include <ncurses.h>
#elif defined(CURSES_HAVE_NCURSES_NCURSES_H)
#	include <ncurses/ncurses.h>
#elif defined(CURSES_HAVE_NCURSES_CURSES_H)
#	include <ncurses/curses.h>
#endif

// write can be found in unistd and conflicts with the write in this file
#define write write_string

/*********************************************************************/

int *ivector(long int nl, long int nh)
/*allocate a int vector with subscript range[nl..nh]*/
{
	int *v;
	
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in ivector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_ivector(int *v, long int nl, long int nh)
/*free a float vector allocated with ivector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/*********************************************************************/

long *lvector(long int nl, long int nh)
/*allocate a int vector with subscript range[nl..nh]*/
{
	long *v;
	
	v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in lvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_lvector(long *v, long int nl, long int nh)
/*free a float vector allocated with ivector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/*********************************************************************/

float *vector(long int nl, long int nh)
/*allocate a float vector with subscript range[nl..nh]*/
{
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in vector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_vector(float *v, long int nl, long int nh)
/*free a float vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/********************************************************************/

double *dvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_dvector(double *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/********************************************************************/

char *cvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	char *v;

	v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_cvector(char *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/****************************************************************************/
int **imatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/*allocate pointers to rows*/
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in imatrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in imatrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_imatrix(int **m, long nrl,long nrh,long ncl, long nch)
/*free an int matrix allocated by imatrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

float **matrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/*allocate pointers to rows*/
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float *)));

	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_matrix(float **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

double **dmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/*allocate pointers to rows*/
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double *)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_dmatrix(double **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

char **cmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;
	
	/*allocate pointers to rows*/
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char *)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_cmatrix(char **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	
/****************************************************************************/

float ***fpmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float ***m;
	
	/*allocate pointers to rows*/
	m=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float **)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_fpmatrix(float ***m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

float ***f3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a float 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_f3tensor(float ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an float f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	


/****************************************************************************/

double ***d3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a float 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_d3tensor(double ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an float f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	

/****************************************************************************/

int ***i3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a int 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_i3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an int f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	

/****************************************************************************/

float bico(int n, int k)
/*returns the binomial coefficient Cn(haut)k(bas) = n!/(k!(n-k)!) as a floating-
point number*/
{
	float factln(int n);
	return (float)floor(0.5+exp(factln(n)-factln(k)-factln(n-k))); 
}

/****************************************************************************/

float factln(int n)
/*returns ln(n!)*/
{
	float gammln(float xx);
	static float a[101];
	if (n<0) {
		write(ERRORFILE,"\nNegative factorial in routine factln"); 
		exit(1);
	}
	if (n <=1) return 0.0f;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0f));
	else return gammln(n+1.0f);
}

/*******************************************************/

float gammln(float xx)
/*numerical recipes p.214, returns the value of ln(gamma(xx)) for xx>0*/
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677, 
		24.01409824083091,-1.231739572450155, 
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y; 
	return (float)(-tmp+log(2.5066282746310005*ser/x));
}

/*******************************************************/

float factrl(int n)
/*returns the value n! as a floating-point number*/
{
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;
	
	if(n<0) {
		write(ERRORFILE,"\nNegative factorial in routine factrl");
		exit(1);
	}
	if(n>32) return (float)exp(gammln(n+1.0f));
	
	while(ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}


/******************************************************************************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
		int j;
		long k;
		static long iy=0;
		static long iv[NTAB];
		float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j= (int) iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=(float) AM*iy) > RNMX) return (float) RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/********************************************************************************/

/* (C) Copr. 1986-92 Numerical Recipes Software . */
long timeseed(void)
/*output a seed for random generator based on system time by  multiplying sec*min*hours or give 1 if 
sec=min=hour=0*/
{
		struct tm *local;
		time_t t;
		long seed;

		seed=1;
		t = time(NULL);
		local = localtime(&t); 
		if(local->tm_sec) seed=local->tm_sec; 
		if(local->tm_min) seed *= local->tm_min; 
		if(local->tm_hour) seed *= local->tm_hour;
		return seed;
}		 



/****************************************************************************/
/****************************************************************************/

void write_string(char *filename,char *info)
/* print string info to display and writes string info to logfile*/
{
	FILE *fp;
	
	printf("%s",info);
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"%s",info);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		printf("error in writing %s : cannot open file",info); 
		wait_a_char();
		exit(3);
	}
}	  /*end of write*/

/****************************************************************************/

void write_tofile_only(char *filename,char *info)
/*  writes string info to logfile without display*/
{
	FILE *fp;
	
	/*printf("%s",info);  */
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"%s",info);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		printf("error in writing %s : cannot open file",info); 
		wait_a_char();
		exit(3);
	}
}	  /*end of write*/

/****************************************************************************/

double sqr(double x)
{
	return x*x;
}

/****************************************************************************/

void avestd(float data[], unsigned long n, float missing, float *ave, float *std)
/*given array data[1..n], returns its mean as ave and its standard deviation as std
value in array data equal to missing are not taken into account*/
{
	unsigned long j;
	float s,ep,var;
	unsigned long count=0;
	
	for (*ave=0.0F,j=1;j<=n;j++) if(data[j] != missing) {   /*skip over missing values*/
		*ave += data[j];
		count++;
	}
	*ave /= (float) count;
	var=ep=0.0F;
	for (j=1;j<=n;j++) if(data[j] != missing) {
		s=data[j]-(*ave);
		ep += s;
		var += s*s;
	}
	if(count>1) var =(var-ep*ep/(float) count)/(count-1.0F);
	else var = 0.0F;
	*std = (float) sqrt(var);
}  /*end of procedure avestd*/

/****************************************************************************/

void tptest(float data1[],float data2[],unsigned long n,float missing,float *t,float *prob)
/*given the paired arrays data1[1..n] and data2[1..n] this routine returns Student's t
for paired data set, and its significance as prob, small values of prob indicating a 
significant difference of means, numerical recipes p. 618
modified routine, values in data1 and data2 equal to missing are not taken into account
attention missing values should be missing both in data1 and data2 at the same positions
returns t=0 and prob=999 if not more than one non-missing value*/
{
	unsigned long j,count=0;
	float var1,var2,ave1,ave2,sd,df,cov=0.0F;
	
	avestd(data1,n,missing,&ave1,&var1);
	avestd(data2,n,missing,&ave2,&var2);
    var1 = (float)sqr(var1);	/*because avestd returns standard deviation and not variance*/
    var2 = (float)sqr(var2);
	for (j=1;j<=n;j++) if(data1[j] != missing) {  /*skip over missing values in data1 which must also correspond to missing values in data2*/
    	cov += (data1[j]-ave1)*(data2[j]-ave2);
    	count ++;
    }
    df = (float)(count-1);
    if(df>0) cov /= df;
    else {            /*if only one non missing value in data1*/
    	*t=0.0F;
    	*prob=999.0F;
    	write(ERRORFILE,"\n Error in tptest: only one non missing value in data1");
    	return;
    }
    sd=(float)sqrt((var1+var2-2.0*cov)/count);
    *t=(ave1-ave2)/sd; 
    *prob=0.0F;	/*routine betai not yet implemented*/
/*    *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));       */
   }	/*end procedure tptest*/	

/****************************************************************************/

void tutest(float data1[],unsigned long n1,float data2[],unsigned long n2,float *t,float *prob)
/*given the arrays data1[1..n1] and data2[1..n2] this routine returns Student's t
 with unequal variances and its significance as prob, small values of prob indicating a 
significant difference of means, numerical recipes p. 617*/
{
	unsigned long j,count=0;
	float var1,var2,ave1,ave2,sd,df,cov=0.0F;
	
	if((n1<=1) || (n2<=1)) {
		if((n1==0) || (n2==0)) 
			write(ERRORFILE,"\n Error in tutest: n1 or n2 equal to zero");
		else write(ERRORFILE,"\n Error in tutest: n1 or n2 equal to one");
		*t=0.0F;
		*prob=999.0F;
		return;
	}
	avestd(data1,n1,MISSING,&ave1,&var1);
	avestd(data2,n2,MISSING,&ave2,&var2);
    var1 = (float)sqr(var1);	/*because avestd returns standard deviation and not variance*/
    var2 = (float)sqr(var2);
    if(!  sqrt(var1/n1+var2/n2)) {
			write(ERRORFILE,"\n Error in tutest: var1 and var2 equal to zero");
		*t=0.0F;
		*prob=999.0F;
		return;
	}
    	
    *t=(float)((ave1-ave2)/sqrt(var1/n1+var2/n2));
    df = (float)(sqr(var1/n1+var2/n2)/(sqr(var1/n1)/(n1-1)+sqr(var2/n2)/(n2-1)));
    *prob=0.0F;	/*routine betai not yet implemented*/
/*    *prob=betai(0.5*df,0.5,df/(df+sqr(*t)));       */
}	/*end procedure tutest*/	



/****************************************************************************/

int arrondir(float x)
{
	int i=(int) x;
	
	if(x == fabs(x)) {
		if((x-i) < 0.5F) return i;
		else return i+1;
	} else {
		if((i-x) < 0.5F) return i;
		else return i-1;
	}
}



/*****************************************************************************/

void shellsort(unsigned long n, float a[])
/*sorts an array[1..n] into ascending numerical order by Shell's method
n is input, a is replaced on output by its sorted rearrangement
Numerical recipes in C p. 332
*/
{
	unsigned long i,j,inc;
	float v;
	inc=1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do {
		inc /= 3;
		for (i=inc+1;i<=n;i++) {
			v=a[i];
			j=i;
			while(a[j-inc] > v) {
				a[j]=a[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v;
		}
	}	while (inc > 1);
}

/*****************************************************************************/

void shellsortd(unsigned long n, double a[])
/*sorts an array[1..n] into ascending numerical order by Shell's method
n is input, a is replaced on output by its sorted rearrangement
Numerical recipes in C p. 332
*/
{
	unsigned long i,j,inc;
	double v;
	inc=1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do {
		inc /= 3;
		for (i=inc+1;i<=n;i++) {
			v=a[i];
			j=i;
			while(a[j-inc] > v) {
				a[j]=a[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v;
		}
	}	while (inc > 1);
}


/*****************************************************************************/

/*this procedure does not work, problem with first position*/
void shellsort_bis(unsigned long first,unsigned long last,float a[])
/*sorts an array[first..last] into ascending numerical order by Shell's method
n is input, a is replaced on output by its sorted rearrangement
Numerical recipes in C p. 332
*/
{
	unsigned long i,j,inc,n;
	float v;
	
	n=last-first+1;
	inc=1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do {
		inc /= 3;
		for (i=inc+first;i<=last;i++) {
			v=a[i];
			j=i;
			while(a[j-inc] > v) {
				a[j]=a[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v;
		}
	}	while (inc > 1);
}

/******************************************************************************/
/******************************************************************************/

void resample_init(int *index,int first,int last)
/*initializes index_vector by setting values respectively to first, first+1,first+2,...last*/
{
	int i;
	
	for(i=first;i<=last;i++) index[i]=i;
}

/**************************************************************************************/

void resample_shuffle(int *index, int first,int last, long *seed)
/*randomize the order of the values pointed to by index_vector by shuffling*/
{
	int i,n,posfirst,valuetemp,picki;
	int VERIF=0;
	long iok,icheck;
	char smess[200];
	
	n=last-first+1;
	posfirst=first;
	
	for(i=0;i<n-1;i++) {
		picki = posfirst + (int) ((float) (n-i) * ran1(seed));
		valuetemp = index[posfirst];
		index[posfirst]=index[picki];
		index[picki]=valuetemp;
		posfirst++;
	}
	
	if(VERIF) {
		for(i=first,iok=0;i<=last;i++) iok +=i;
		for(i=first,icheck=0;i<=last;i++) icheck +=index[i];
		if(icheck != iok) {
			write(ERRORFILE,"\n error after shuffling, see shuffled vector below:\ni\tvalue");
			for(i=first;i<=last;i++) {
				sprintf(smess,"\n%d\t%d",i,index[i]);
				write(ERRORFILE,smess);
				exit(1);
			}
		}  /*end if icheck different than iok*/
	}	/*end if VERIF*/ 	
}

/**************************************************************************************/

/*test of shuffling procedure

index=ivector(1,100);
store=ivector(1,100);
for(r=1;r<=100;r++) store[r]=0;
resample_init(index,1,100);
seed=timeseed();
seed = -seed;
sprintf(smess,"\nseed:%ld",seed);
write(ERRORFILE,smess);
for (l=0;l<100;l++) {
	resample_shuffle(index, 1,100, &seed);
	store[index[100]]++;
}
write(ERRORFILE,"\ni\tscore");
for (r=1;r<=100;r++) {
	sprintf(smess,"\n%d\t%d",r,store[r]);
	write(ERRORFILE,smess);
}

end of test of shuffling procedure*/

/**********************************************************************************/
struct resample_stat_type **resample_stat_type_matrix(long nrl, long nrh,long ncl, long nch)
/*allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	struct resample_stat_type **m;
	
	/*allocate pointers to rows*/
	m=(struct resample_stat_type **) malloc((size_t)((nrow+NR_END)*sizeof(struct resample_stat_type*)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in resample_stat_type_matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(struct resample_stat_type *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(struct resample_stat_type)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in resample_stat_type_matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}

/**********************************************************************************/

void free_resample_stat_type_matrix(struct resample_stat_type **m, long nrl,long nrh,long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/**********************************************************************************/


void resample_analyze(float *r_vector,unsigned long first,unsigned long last,
			struct resample_stat_type *r_stat,float observed)
/*sorts r_vector then take the 5% and 1% critical values*/
/*computes proportion of values higher or lower than observed value = P values*/
/*!!!does only work temporarily if r_vector is from 1 to n!!!*/
{
	unsigned long low95,high95,low99,high99;
	unsigned long n,i;
	double average=0.,var=0.;
	float ploworequal,phighorequal;
	
	r_stat->obs=observed;
	n=last-first+1;
	
	/*compute mean*/
	for(i=first;i<=last;i++) average+=r_vector[i];
	average/=n;
	r_stat->mean=(float)average;

	/*compute variance*/
	for(i=first;i<=last;i++) var+=(r_vector[i]-average)*(r_vector[i]-average);
	var/=(n-1);

	/*compute sd*/
	r_stat->sd=(float)sqrt(var);

	/*shellsort_bis(first,last,r_vector); */      /*this procedure doesn't work presently*/
/*does only work temporarily if r_vector is from 1 to n*/	
	shellsort(n,r_vector);

	/*count the # of different values*/
	r_stat->nd=1;
	for(i=first+1;i<=last;i++) if(r_vector[i]!=r_vector[i-1]) r_stat->nd++;


	low95=(int) ((float) n / 40.0F);	/*1/40=2.5%*/
	high95=(int) (39.0F*(float) n/40.0F)+1;
	if(first+high95 > last) high95=last-first;
	low99=(int) ((float) n / 200.0F);
	high99=(int) (199.0F*(float) n/200.0F)+1;
	if(first+high99 > last) high99=last-first;
    
    if(n>=40) {
	    r_stat->low95=r_vector[first+low95];
	    r_stat->high95=r_vector[first+high95];
	} else r_stat->low95=r_stat->high95=0.0F;	/*if insufficient number of samples*/
    if(n>=200) {
	    r_stat->low99=r_vector[first+low99];
	    r_stat->high99=r_vector[first+high99];
    } else r_stat->low99=r_stat->high99=0.0F;	/*if insufficient number of samples*/
    
	i=first;
    while((r_vector[i] <= observed) && (i<=last)) {
    	i++;
    }
    if(i>last) r_stat->plow = 1.0F;
    else if(i==first) r_stat->plow = 0.0F;
    else r_stat->plow = (float) (i-first+1)/(n+1);			/*see Box 18.5 p.816 Sokal & Rholf*/
	
	i=last;
    while((r_vector[i] >= observed) && (i>=first)) i--;
    if(i<first) r_stat->phigh = 1.0F;
    else if(i==last) r_stat->phigh = 0.0F;
    else r_stat->phigh = (float) (last-i+1)/(n+1);			/*see Box 18.5 p.816 Sokal & Rholf*/
	
	/*i=first;
    while((r_vector[i] < observed) && (i<=last)) i++;
    if(i>last) ploworequal = 1.0F;
    else if(i==first) ploworequal = 0.0F;
    else ploworequal = (float) (i-first+1)/(n+1);		
	i=last;
    while((r_vector[i] > observed) && (i>=first)) i--;
    if(i<first) phighorequal = 1.0F;
    else if(i==last) phighorequal = 0.0F;
    else phighorequal = (float) (last-i+1)/(n+1);
	if(ploworequal < phighorequal) r_stat->pbil=2.*ploworequal;
	else r_stat->pbil=2.*phighorequal;*/

	if(r_stat->plow < r_stat->phigh) r_stat->pbil=r_stat->plow + (1.0f-r_stat->phigh);
	else r_stat->pbil=r_stat->phigh + (1.0f-r_stat->plow);

	if(r_stat->pbil>1.0f) r_stat->pbil=1.0f;
}	 
/**********************************************************************************/

void resample_write(char *filename,float *r_vector,float actual_value, unsigned long first,unsigned long last,char *var_name) 
{
	FILE *fp;
	unsigned long i;
	
	if((fp=fopen(filename,"wt"))==NULL) {
		write(ERRORFILE,"\ncannot open file in resample_write");
		exit(1);
	}
	fprintf(fp,"Resampled statistics with samplesize : %ld",last-first+1);
	fprintf(fp,"\nActual value of %s:\t%8.4f",var_name,actual_value);
	fprintf(fp,"\n%s",var_name);
	for(i=first;i<=last;i++)	fprintf(fp,"\n%8.4f",r_vector[i]);
	fclose(fp);
}	/*end procedure resample_write*/

/**********************************************************************************/

void linear_fit(float x[],float y[],int ndata,float *a, float *b)
/*given a set of data points x[0..ndata],y[0..ndata] fit them to a straight line y=a+bx by 
minimizing chi2. Return only a, the intercept and b, the slope of the line from Numerical recipes p.665. */
{
		int i;
		float t,sxoss,sx=0.0F,sy=0.0F,st2=0.0F,ss;
		
		*b=0.0;
		for(i=0;i<ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss = (float)ndata;
		sxoss=sx/ss;
		for(i=0;i<ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
		*b /= st2;
		*a=(sy-sx*(*b))/ss;
}	    /* end of linear_fit */


