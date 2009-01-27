#pragma once
#ifndef XATOOLS_H
#define XATOOLS_H

int *ivector(long int nl, long int nh);
void free_ivector(int *v, long int nl, long int nh);
long *lvector(long int nl, long int nh);
void free_lvector(long *v, long int nl, long int nh);
float *vector(long int nl, long int nh);
void free_vector(float *v, long int nl, long int nh);
double *dvector(long int nl, long int nh);
void free_dvector(double *v, long int nl, long int nh);          
char *cvector(long int nl, long int nh);
void free_cvector(char *v, long int nl, long int nh);
int **imatrix(long nrl, long nrh,long ncl, long nch);
void free_imatrix(int **m, long nrl,long nrh,long ncl, long nch);
float **matrix(long nrl, long nrh,long ncl, long nch);
void free_matrix(float **m, long nrl,long nrh,long ncl, long nch);
double **dmatrix(long nrl, long nrh,long ncl, long nch);
void free_dmatrix(double **m, long nrl,long nrh,long ncl, long nch);
char **cmatrix(long nrl, long nrh,long ncl, long nch);
void free_cmatrix(char **m, long nrl,long nrh,long ncl, long nch);
float ***fpmatrix(long nrl, long nrh,long ncl, long nch);
void free_fpmatrix(float ***m, long nrl,long nrh,long ncl, long nch);
float ***f3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
void free_f3tensor(float ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
double ***d3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
void free_d3tensor(double ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
int ***i3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
void free_i3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh);

float bico(int n, int k);
float ran1(long *idum);
long timeseed(void);
float gammln(float xx);
float factrl(int n);
float factln(int n);

void write_string(char *filename,char *info);
void write_tofile_only(char *filename,char *info);
double sqr(double x);
void avestd(float data[], unsigned long n, float missing, float *ave, float *std);
void tptest(float data1[],float data2[],unsigned long n,float missing,float *t,float *prob);
void tutest(float data1[],unsigned long n1,float data2[],unsigned long n2,float *t,float *prob);
int arrondir(float x);
void shellsort(unsigned long n, float a[]);
void shellsortd(unsigned long n, double a[]);
void shellsort_bis(unsigned long first,unsigned long last, float a[]);
void resample_init(int *index,int first,int last);
void resample_shuffle(int *index, int first,int last, long *seed);

struct resample_stat_type {
	float obs;		/*observed value*/
	float mean;		/*average of resampling*/
	float sd;
	float low95;  	/*critical values at the 5%level*/
	float high95;
	float low99;    /*critical values at the 1%level*/
	float high99;
	float plow;    /*proportion of values lower than observed ones*/
	float phigh;   /*proportion of values higher than observed ones*/
	float pbil;		/*bilateral test*/ 
	int n;			/*# of resamplings*/
	int nd;			/*# of different values among all resamplings (with valid values)*/
	int mode;		/*mode of permutation (object permuted)*/
};
struct resample_stat_type **resample_stat_type_matrix(long nrl, long nrh,long ncl, long nch);
void free_resample_stat_type_matrix(struct resample_stat_type **m, long nrl,long nrh,long ncl, long nch);
void resample_analyze(float *r_vector,unsigned long first,unsigned long last,struct resample_stat_type *r_stat,float observed);
void resample_write(char *filename,float *r_vector,float actual_value, unsigned long first,unsigned long last,char *var_name);
void linear_fit(float x[],float y[],int ndata,float *a, float *b);

#endif

