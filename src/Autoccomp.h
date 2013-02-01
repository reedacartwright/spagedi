/************************************************************************* 
 * Copyright (c) 2002-2011 Olivier Hardy and Xavier Vekemans             *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *************************************************************************/

#pragma once
#ifndef AUTOCCOMP_H
#define AUTOCCOMP_H

#ifdef PACKAGE_STRING
#	define VERSION PACKAGE_STRING
#else
#	define  VERSION "SPAGeDi 1.3e-STABLE"
#endif

#define  VARTYPE 0  /*0 for locus with genotypes, 1 for quantitative variable*/
#define  NMAX 100000	/*max number of individuals*/
#define  NCOORDMAX 3    /*max number of spatial coordinates*/
#define  MMAX 10000    /*max number of loci or variables*/
#define  NDIGITMAX 3    /*max number of digits per allele*/
#define  PLOIDYMAX 8    /*max number of digits per allele*/
#define  MAXALLID 999 /*max number of alleles per locus*/
#define  MAXNOM 31  /*max number of characters read for individual names*/
#define  MAXINTERVALS 102	/*max number of classes of intervals*/
#define  NRESAMPLE 20001 	/*max number of resamplings*/
#define  SMAX 10000	/*max number of characters for strings*/
#define  SMAX2 1000000	/*max number of characters for strings*/
#define  ERRORTXT "error.txt"	/*name of file for error messages*/
#define  MISSVAL HUGE_VAL	/*value given for missing value*/
#define  MISSING 999.99F        /*missing value for avestd and ttest routines*/

#include <limits.h>
#include <stdlib.h>

#ifndef PATH_MAX
#	define PATH_MAX _MAX_PATH
#endif

#include "Xatools.h"

extern char errorfile[PATH_MAX];
#define  ERRORFILE errorfile /* name of variable that contains errorfile */

struct name{
	char n[MAXNOM];
};

/*in file main.c*/
void mainAnalysisBtwInd(int argc,int n,int ntot,double *xi,double *yi,double *zi,double **dij,int *sgi,int Nsg,int *skgi,int Nskg,
	int *catskg,int *cati,int Ncat,int *Nik,int nc,double *maxc,float dijmin,float dijmax,
	int m,int ndigit,int ploidy,int *ploidyi,int *Nallelel,int **Nallelekl,int **Nvalgenkl,int ***gilc,
	float ***Pkla,int **allelesizela,float **Masizekl,float **Vasizekl,float ***Mgdlaa,float givenF,
	struct name namei[],char namelocus[][MAXNOM],struct name namecat[],
	int TypeComp,int cat1,int cat2,int FreqRef,float **givenPla,int *Ngivenallelel,int JKest,
	int NS,int Stat[],int printdistmatrix,float sigmaest,float density,float dwidth,
	int Npermut[],int permutalleles,int writeresampdistr,int regdetails,int varcoef,int Rbtwloc,
	int permutdetails,char outputfile[]);
void mainAnalysisBtwPop(int argc,int StatType,int n,double *xp,double *yp,double *zp,double **dij,int *popi,int Npop,
	int *catp,int *cati,int Ncat,int nc,double *maxc,float dijmin,float dijmax,
	int m,int ndigit,int ploidy,int *ploidyi,int *Nallelel,int **Nvalgenkl,int ***gilc,
	float ***Pkla,int **allelesizela,float **Masizekl,float **Vasizekl,float ***Mgdlaa,
	struct name namepop[],char namelocus[][MAXNOM],struct name namecat[],
	int PWstat,int TypeComp,int cat1,int cat2,int FreqRef,int JKest,int NS,int Stat[],int printdistmatrix,
	int Npermut[],int permutalleles,int writeresampdistr,int regdetails,int varcoef,int Rbtwloc,int permutdetails,char outputfile[]);



/*in file autoccomp.c*/
void define_groups(int n,struct name namei[],double *xi,double *yi,double *zi,
		int Ncat,struct name namecat[],struct name namecati[],int *cati,int *Nik,
		int *Nsg,struct name namesg[],int *sgi,int *Nig,
		int *Nskg,struct name nameskg[],int *skgi,int *catskg,int *Niskg,
		int orig_order);
void reorder_data_by_category(int n,int Ncat,int *cati,int *cat1,int *cat2,
			int *ploidyi,int *sgi,int *skgi,double *xi,double *yi,double *zi,
			struct name *namei,struct name *namecati,int ***gilc);
void define_pop(int StatType,int n,int Npop,int Nsg,int Ncat,int Nskg,
		int *popi,int *sgi,int *cati,int *skgi,int *Nip,int *Nig,int *Nik,int *Niskg,
		int *catskg,int *catp,struct name *namepop,struct name *namesg,struct name *namecat,struct name *nameskg,
		double *xp,double *yp,double *zp,double *xi,double *yi,double *zi);
void checkdist(int n, int *nc, double *maxc, double *xi, double *yi, double *zi,double **Mdij,
			   int Nsg,int *sgi,int *cati,int StatType,int TypeComp,
			   double *mdc,double *mlndc, int *npc,float **indexpartic);
void compute_allele_freq(int n,int Ncat,int *cati,int m,
			int ndigit,int ploidy,int ***gilc,int *ploidyi,int *Nallelel,int **allelesizela,float ***Mgdlaa,
			int alleledist,float ***Pkla,int **Nallelekl,int **Nmissinggenotkl,int **Nincompletegenotkl,
			int **Nvalgenkl,float **Hekl,float **hTkl,float **uTkl,float **Dmkl,float **Dwmkl,float **Masizekl,float **Vasizekl);
void compute_pairwise_corr_F(int n,int ntot,int Ncat,int *cati,int m,int ndigit,int ploidy,
			float missdat,int ***gilc,int *Nallelel,int **allelesizela,float ***Distla1a2,
			float ***corrSlij[],int NS,int Stat[12],int FreqRef,float **givenPla,int *Ngivenallelel,
			int TypeComp,float givenF,int compute_inbreeding_coef_only,int JKl);
void compute_F_R_stat(int n,int Npop,int pop1,int pop2,int *popi,int m,int *Nallelel,
		int ploidy,int ***gilc,int **allelesizela,int NS,int Stat[],
		float **FstatSlr[],float ***corrSlij[],int Rstat_only,int JKest);
void compute_G_N_stat(int n,int Npop,int pop1,int pop2,int *popi,int m,int *Nallelel,
		int ploidy,int ***gilc,float ***Ppla,float ***Distla1a2,int NS,int Stat[],
		float **FstatSlr[],float ***corrSlij[],int JKest,int computeallelefreq);
void compute_corr_per_dist_class (int n,int m,int nc,double *maxc,int Ncat,int *cati,
		int StatType,int TypeComp,int FreqRef,int varcoef,double *xi,double *yi,double *zi,double **Mdij,
		int *sgi,float dijmin,float dijmax,float ***corrlij,float **corrlc,float ***Rll,int Rbtwloc,float ***V,float **R2pl,
		long *seed,int JKl);
void inter_locus_corr(int n,int m,float ***corrlij,float **Rll,float **V,float *R2pl,long *seed);
void estimate_sigma_2D_kinship (int n,int m,double *xi,double *yi,double *zi,double **Mdij,
		int *sgi,float ***corrlij,float **corrlc,int ploidy,int Stat,int JKl,float density,float dwidth);
void analyse_resampling(int m,int cinit,int nc,int Npermut,float ***corrlcs,
		struct resample_stat_type **r_statlc);

void compute_ANOCOVA_within_fam(int n,int m,double *xi,double *yi,int Nloc,
								 float missdat,float **qiv,float ***corrtv1v2);


void NestedANOVA(int a,int *bi,int **Nij,double ***Yijk,double SS[4],double MS[4],double s2[4]);

void permut_allelesizes
	(int m,int *Nallelel,int **allelesizela,int **allelesizelamix,long *seed);
void permut_genetic_distances
	(int m,int *Nallelel,float ***Mgdlaa,float ***Mgdlaamix,long *seed);
void permut_allelesizes_among_2pop_alleles
	(int m,int *Nallelel,float ***Ppla,int pop1,int pop2,int **allelesizela,int **allelesizelamix,long *seed);
void permut_genetic_distances_among_2pop_alleles
	(int m,int *Nallelel,float ***Ppla,int pop1,int pop2,float ***Mgdlaa,float ***Mgdlaamix,long *seed);
void permut_locations
	(int n,double *x,double *y,double *z,double **Mdij,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed);
void permut_locations_of_groups
	(int n,double *x,double *y,double *z,double **Mdij,int *groupi,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed);
void grumph
	(int n,double *x,double *y,double *z,double **Mdij,int *groupi,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed);
void permut_locations_within_cat
	(int n,int *cati,int Ncat,double *x,double *y,double *z,double **Mdij,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed);
void permut_locations_of_groups_within_cat
	(int n,int *cati,int Ncat,double *x,double *y,double *z,double **Mdij,int *groupi,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed);
void permut_genes_among_indiv
	(int n,int m,int ploidy,int *ploidyi,int ***gilc,int ***gilcmix,long *seed);
void permut_genes_among_indiv_within_pop
	(int n,int m,int ploidy,int *ploidyi,int *popi, int Npop,
	int ***gilc,int ***gilcmix,long *seed);
void permut_indiv_among_pop
	(int n,int *popi, int Npop,int *popimix,long *seed);
void permut_indiv_among_pop_within_categ
	(int n,int *popi, int Npop,int *cati,int Ncat,int *popimix,long *seed);

#endif

