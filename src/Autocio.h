/************************************************************************* 
 * Copyright (c) 2002-2009 Olivier Hardy and Xavier Vekemans             *
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
#ifndef AUTOCIO_H
#define AUTOCIO_H

#include "Xatools.h"

void commom_errors();
void get_input_output_file_names(int argc,char *argv[],char inputfile[],char outputfile[],char instrfile[]);
void import_data_file(char *inputfile);
void export_data_file(int StatType,char inputfile[],int n,
	double *xi,double *yi,double *zi,double *xp,double *yp,double *zp,double **Mdij,
	int *popi,int Npop,int m,int ndigit,int *Nallelel,int ***gilc,int **allelesizela,
	struct name namei[],struct name namepop[],char namelocus[][MAXNOM]);
void readbasicinfoF(char *inputfile,int *n,int *ncat,int *ncoord,int *mp,int *ndigit,int *ploidy,
		char catname[MAXNOM],char namecoord[] [MAXNOM],char namelocus[][MAXNOM],
		int *nc,double *maxc);
void readsecondinfoF(char *inputfile,int n,struct name namei[],int Ncat,struct name namecat[],
		struct name namecati[],int *cati,int ncoord,double *xi,double *yi,double *zi,
		int m,int ndigit,int ploidy,int ***gila,int *Nallelel,int **allelesizela,
		int *ploidyi,int Nip[]);
void displaybasicinfoF(int argc,char *inputfilename,char *outputfilename,
		int n,int Ncat,int *Nik,int ncoord,int m,int ndigit,
		int ploidy,int *ploidyi,int Nip[],char namelocus[] [MAXNOM],
		char namecoord[] [MAXNOM],struct name *namecat,
		int nc,double *maxc,int Nsg,int *Nig,int Nskg,int *Niskg);
void define_analysisF(int argc,char *instrfile,int n,int ploidy,int ndigit,int m,int Ncat,int Nsg,int Nskg,int ncoord,
	  int *StatType,int *NS,int Stat[12],int *TypeComp,int *cat1,int *cat2,
	  struct name *namecat,int *FreqRef,float *givenF,int *writeallelefreq,int *JKest,int *distmatrix,int Npermut[],
	  float *dijmin,float *dijmax,int *writeresampdistr,int *regdetails,int *varcoef,int *Rbtwloc,float *sigmaest,float *density,float *dwidth,
	  int *permutdetails,int *distm,char *inputfilename,char *distfilename,char *freqfilename,
	  int *definealleledist,char *alleledistfilename,int *export,long *seed);
void read_pairwise_spatial_distances(char *inputfile,char *distfile,int n,
		struct name *namei,double **Mdij);
void read_pairwise_genetic_distances(char *inputfile,char *distfile,int m,
		char namelocusp[] [MAXNOM],int *Nallelel,int **allelenamela,float ***Mgdlaa);
void read_allele_frequencies(char *freqfile,int m,char namelocus[][MAXNOM],int *Nallelel,
							 int **allelesizela, float **givenPla,int *Ngivenallelel);
void write_allele_freq(char *outputfilename,int n,int m,char namelocus[][MAXNOM],int ndigit,
		int Npop,int *Nip,struct name namepop[],int StatType,int NS,int Stat[],int **allelesizela,
		int **Nallelepl,float ***Ppla,float **Hepl,float **hTpl,float **uTpl,float **Dmpl,float **Dwmpl,float **Masizepl,float **Vasizepl,
		int **Nmissinggenotpl,int **Nincompletegenotpl,int **Nvalgenpl,int printallelefreq,int FreqRef,int *Ngivenallelel,float **givenPla);
void displaydist(int argc,char *outputfilename,int nc,double *maxc,double *mdc,double *mlndc,
				 int *npc,float **indexpartic);
void writeIndStatresults(char *outputfilename,int n,int Nsg,int m,char namelocus[][MAXNOM],int nc,double *maxc,
		int *npc,float **indexpartic,double *mdc,double *mlndc,float dijmin,float dijmax,float givenF,
		int TypeComp,int cat1,int cat2,struct name *namecat,int FreqRef,int NS,int Stat[12],float **corrSvc[10],float density,float dwidth,
		int JKest,int regdetails,int varcoef,int Rbtwloc,float ***RSll[12],float ***V[12],float **R2pl[12]);
void writePopStatresults(char *outputfilename,int Npop,int m,char namelocus[][MAXNOM],int nc,double *maxc,
		int *npc,float **indexpartic,double *mdc,double *mlndc,float dijmin,float dijmax,
		int StatType,int TypeComp,int cat1,int cat2,struct name *namecat,int NS,int Stat[12],float **corrSlc[],float **FstatSlr[],
		int JKest,int regdetails,int varcoef,int Rbtwloc,float ***RSll[12],float ***V[12],float **R2pl[12],int PWstat);
void writedistmatrices (char *outputfilename,int n,int m,float givenF,int TypeComp, int *cati,
		int printdistmatrix,double *xi,double *yi,double *zi,double **Mdij,int *sgi,
		int StatType,int NS,int Stat[12],
		float ***corrSlij[10],struct name namei[],char namelocus[][MAXNOM]);

void WriteIndPermutRes(char *outputfilename,int n,int ploidy,int Ncat,
		double *maxc,int m,char namelocus[][MAXNOM],int TypeComp,float givenF,
		int nc,int Nsg,int NS,int Stat[],int Npermut[],int permutalleles,
		int permutdetails,struct resample_stat_type **r_statSlc[10],long seedinit);
void WritePopPermutRes(char *outputfilename,int n,int ploidy,
		double *maxc,int m,char namelocus[][MAXNOM],int TypeComp,
		int nc,int NS,int Stat[],int Npermut[],int permutalleles,
		int permutdetails,struct resample_stat_type **r_statFSlr[],
		struct resample_stat_type **r_statSlc[],long seedinit,int PWstat);
void WriteAllelesPermutResForPopPair(char *outputfilename,int m,char namelocus[][MAXNOM],
		int p1,int p2,struct name *namepop,int NS,int Stat[],
		struct resample_stat_type **r_statSlc[10],int permutdetails,int comp);

char *nexttab(char *s, char *inputfile, int line);
char *nexttab_blank(char *s, char *inputfile, int line);
void readsfromfile(FILE *fp, char *s, char *inputfile, int *line);
void readsfromfile_no_end_of_file_check(FILE *fp, char *s, char *inputfile, int *line);
void readsfromstring(char *ins, char *outs,int nchar, char *inputfile, int line);
void readintfromstring(char *ins, int *outi, char *inputfile, int line);
void readintfromstring_bis(char *ins, int *outi, char *inputfile, int line);
void readfloatfromstring(char *ins, float *outf, char *inputfile, int line);
void readdoublefromstring(char *ins, double *outf, char *inputfile, int line);


//NON ANSI FUNCTION
char wait_a_char();

char *fgets_chomp(char * a, int b, FILE * c);

#endif

