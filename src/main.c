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
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Autoccomp.h"
#include "Autocio.h"
#include "Xatools.h"

long seed;


int main(int argc,char *argv[])
{
	char inputfile[PATH_MAX];			//name of datafile
	char outputfile[PATH_MAX];		//name of outputfile
	char instrfile[PATH_MAX];			//name of instruction file (for command line)
	char distfile[PATH_MAX];			//name of the file with spatial pairwise distances
	char alleledistfile[PATH_MAX];	//name of the file with distances btw alleles
	char freqfile[PATH_MAX];			//name of file with reference allele frequencies
	char smess[SMAX];
	int n,ntot;					//number of individuals
 	int Ncat;				//# of different categories of ind
	int *Nik;				//#	indiv in categ k, min and max
	char catname[MAXNOM];	//identifier of each category name
	int *cati;					//categ of ind i; categ of spatial group g
	int m;							//number of loci
	char namelocus[MMAX][MAXNOM];   //identifier of each locus (idem)
 	int ncoord;				//number of spatial coordinates
	char namecoord[NCOORDMAX+1] [MAXNOM];		//identifier of each variable (idem)
	int *npc;
	double *xi,*yi,*zi;			//position in x, y, z coordinates of individuals (i)
	double **Mdij;
	int distm;
	int Nsg;					//# of spatial groups
	int *sgi;					//spatial group of ind i
	int *Nig;		//# ind in group g; min and max values
	int Nskg;
	int *skgi;
	int *catskg;
	int *Niskg;
	int Npop,*popi,*Nip,*catp;
	double *xp,*yp,*zp;
	struct name *namei,*namecat,*namecati,*namesg,*nameskg,*namepop;
	int export;		//=1 if data to be exported to Genepop or Fstat format

	int ndigit;				//number of digits per allele
	int ploidy,*ploidyi,Niploidy[2+PLOIDYMAX];	//ploidy of the data for F; ploidy of ind i; #ind of ploidy p
	int maxnal=0;			//max # of alleles
	float missdat;			//value attributed for missing data
	int **Nmissinggenotpl,**Nincompletegenotpl;		//number of missing / incomplete genotypes per locus 
	float dijmin,dijmax,sigmaest,density,dwidth;			//minimal/max distance to take into account for spatial regression
	int nc;				//number of classes of distance intervals (does not include class=0 which is for Fis),class
	double maxc[MAXINTERVALS];//max value for each distance interval
	float givenF;

	int i,p,l,linit,k,newi,*newii,*oldii;	

	int ***gilc,a1,a2;			//genotype of ind i at locus l and chromosome a
	int *rl;
	int **allelesizela;	//allele lenght (number defined in the data file)
	float ***Mgdlaa;	//genetic distance between alleles at each locus
	float ***Ppla,**Hepl,**hTpl,**vTpl,**Dmpl,**Dwmpl,**Masizepl,**Vasizepl;		//matrix of expected heterozygosity, allele frequencies, mean of allele size, and variance of allele size for locus l and pop p
	int **Nallelepl,**Nvalgenpl;	//number of alleles per locus; tot # of defined gene per locus (=ploidy*#ind, except 0 values)
	int *Nallelel,Ngivenallelel[MMAX];
	int StatType;	//StatType=1 for pairs of ind, 2 for pairs of spatial groups, 3 for pairs of categorical groups
	int FreqRef;	//FreqRef=-1 when ref allele frequencies are given in data file, =0 when it reffers to the whole sample, =1 when it refers to a given category 
	float **givenPla = NULL;//reference allele frequencies given
	int TypeComp;	//TypeComp=0 when categories are not considered, =1 for within categories comparisons, =2 for between categories comparisons
	int cat1,cat2;	//category numbers for analysis within a given cat or between 2 given cats
	int NS,Stat[15];			//Stat[0 to NS-1] are the numbers corresponding to the chosen statistics; S is a counter for Stat; NS is the number of different statistics asked
	int JKest;				//JKest=1 if jackknife estimators are to be computed, =0 otherwise
	int PWstat;				//=1 if pairwise comparisons are to be made btw pop
	int Npermut[5];			//total number of permutations and counter

	int alleledist=0;		//=1 if genetic dist btw alleles to be defined
	int printallelefreq;	//=0 if allele freq are to be computed and written only for the whole sample, =1 if it has to be done for each spatial or categorical group
	int writeresampdistr=0;	// 1 if values obtained after each permutation are to be written
	int regdetails;			//!=0 if detailed account of regression analyses are to be given
	int varcoef;			//=1 if the actual variance of pairwise coef is to be computed (Ritland's method)
	int Rbtwloc;			//=1 if the mean inter-locus correlation coefficients for genetic distances must be computed
	int permutdetails;
	int printdistmatrix=0;	// 1 if genetic and spatial distances matrices are to be printed
	int estinbreeding=0,permutalleles=0;	//=0 if no estimate of inbreeding asked

	FILE *filep;

/*
argc=4;
strcpy(inputfile,"in.txt");
strcpy(outputfile,"out.txt");
strcpy(instrfile,"instruction.txt");
*/

	printf("\n%s - a program for Spatial Pattern Analysis of Genetic Diversity\n               Written by Olivier Hardy & Xavier Vekemans\n",VERSION);
	
	//get names of data and output files
	get_input_output_file_names(argc,argv,inputfile,outputfile,instrfile);

	//write common errors in error file
	commom_errors();

	readbasicinfoF(inputfile,&n,&Ncat,&ncoord,&m,&ndigit,&ploidy,catname,namecoord,namelocus,&nc,maxc);
	if(ndigit<=0) ploidy=1; //case of a dominant marker. ndigit is the value attributed to missing data
    // initializations
	ploidyi=ivector(0,n);
	cati=ivector(0,n);
	Nik=ivector(0,Ncat);
	xi=dvector(0,n);
	yi=dvector(0,n);
	zi=dvector(0,n);
	npc=ivector(0,abs(nc)+2);
	sgi=ivector(0,n);
	Nig=ivector(0,n);
	skgi=ivector(0,n);
	Niskg=ivector(0,n);
	catskg=ivector(0,n);
	Nallelel=ivector(0,m);

	gilc=i3tensor(0,n,0,m,0,ploidy-1);
	allelesizela=imatrix(0,m,0,2+(int)pow(10,ndigit));
	missdat=0.;
	if((namei=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(100);
	if((namecati=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(101);
	if((namecat=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(101);
	if((namesg=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(102);
	if((nameskg=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(103);

	readsecondinfoF(inputfile,n,namei,Ncat,namecat,namecati,cati,ncoord,xi,yi,zi,m,ndigit,ploidy,gilc,Nallelel,allelesizela,ploidyi,Niploidy);
	for(l=1;l<=m;l++) if(maxnal<Nallelel[l]) maxnal=Nallelel[l];   //define the max # of alleles
	define_groups(n,namei,xi,yi,zi,Ncat,namecat,namecati,cati,Nik,&Nsg,namesg,sgi,Nig,&Nskg,nameskg,skgi,catskg,Niskg,1);
	displaybasicinfoF(argc,inputfile,outputfile,n,Ncat,Nik,ncoord,m,ndigit,ploidy,ploidyi,Niploidy,namelocus,namecoord,namecat,nc,maxc,Nsg,Nig,Nskg,Niskg);	
	
	define_analysisF(argc,instrfile,n,ploidy,ndigit,m,Ncat,Nsg,Nskg,ncoord,&StatType,&NS,Stat,&TypeComp,&cat1,&cat2,namecat,&FreqRef,&givenF,&printallelefreq,&JKest,&printdistmatrix,Npermut,&dijmin,&dijmax,&writeresampdistr,&regdetails,&varcoef,&Rbtwloc,&sigmaest,&density,&dwidth,&permutdetails,&distm,inputfile,distfile,freqfile,&alleledist,alleledistfile,&export,&seed);
	
	
	
	if(alleledist){
		Mgdlaa=f3tensor(0,m,0,maxnal,0,maxnal);
		Mgdlaa[0][0][0]=1.;
		read_pairwise_genetic_distances(inputfile,alleledistfile,m,namelocus,Nallelel,allelesizela,Mgdlaa);
	}
	else{Mgdlaa=f3tensor(0,0,0,0,0,0); Mgdlaa[0][0][0]=0.;}
	
	if(TypeComp>=1){
		reorder_data_by_category(n,Ncat,cati,&cat1,&cat2,ploidyi,sgi,skgi,xi,yi,zi,namei,namecati,gilc);
		define_groups(n,namei,xi,yi,zi,Ncat,namecat,namecati,cati,Nik,&Nsg,namesg,sgi,Nig,&Nskg,nameskg,skgi,catskg,Niskg,0);
	}
	  
	if(StatType==3)
		Npop=Nsg;
	else if(StatType==4)
		Npop=Nskg;
	else // if(StatType<=2)
		Npop=Ncat; 

  	popi=ivector(0,n);
	if((namepop=(struct name*) malloc((n+1)*sizeof(struct name)))==NULL) exit(104);
	Nip=ivector(0,Npop);
	catp=ivector(0,Npop);
	xp=dvector(0,Npop);yp=dvector(0,Npop);zp=dvector(0,Npop);
	define_pop(StatType,n,Npop,Nsg,Ncat,Nskg,popi,sgi,cati,skgi,Nip,Nig,Nik,Niskg,catskg,catp,namepop,namesg,namecat,nameskg,xp,yp,zp,xi,yi,zi);

	//free memory
	free(namecati);
	free(namesg);
	free(nameskg);
	free_ivector(Nig,0,n);
	free_ivector(Niskg,0,n);

	//read pairwise distances
	if(distm==1){
		if(StatType==1) {
			Mdij=dmatrix(0,n,0,n);
			read_pairwise_spatial_distances(inputfile,distfile,n,namei,Mdij);
		} else /* if(StatType>=2)*/ {
			Mdij=dmatrix(0,Npop,0,Npop);
			read_pairwise_spatial_distances(inputfile,distfile,Npop,namepop,Mdij);
		}
		Mdij[0][0]=1.;
	} else {
		Mdij=dmatrix(0,0,0,0);
		Mdij[0][0]=0.;
	}

	//read reference allele frequencies
	if(FreqRef==-1){
		givenPla=matrix(1,m,0,2+(int)pow(10,ndigit));
		read_allele_frequencies(freqfile,m,namelocus,Nallelel,allelesizela,givenPla,Ngivenallelel);
	}

	PWstat=0;
	if(printdistmatrix) PWstat=1;
	if(ncoord || Mdij[0][0]==1.) PWstat=2;


	//export data if necessary
	if(export) export_data_file(StatType,inputfile,n,xi,yi,zi,xp,yp,zp,Mdij,popi,Npop,m,ndigit,Nallelel,gilc,allelesizela,namei,namepop,namelocus);

	//init for allele freq analyses
	Ppla=f3tensor(-1,Npop,0,m,0,maxnal);
	Nmissinggenotpl=imatrix(0,Npop,0,m);
	Nincompletegenotpl=imatrix(0,Npop,0,m);
	Nallelepl=imatrix(0,Npop,0,m);
	Nvalgenpl=imatrix(0,Npop,0,m);
	Hepl=matrix(0,Npop,0,m);
	hTpl=matrix(0,Npop,0,m);
	vTpl=matrix(0,Npop,0,m);
	Dmpl=matrix(0,Npop,0,m);
	Dwmpl=matrix(0,Npop,0,m);
	Masizepl=matrix(0,Npop,0,m);
	Vasizepl=matrix(0,Npop,0,m);
	rl=ivector(0,m);

	printf("\n\n\nThe program is now doing the analyses\nTo stop it before the end, press 'Ctrl' + 'c'");
	printf("\n\nComputing allele frequencies. Please, wait.");
	compute_allele_freq(n,Npop,popi,m,ndigit,ploidy,gilc,ploidyi,Nallelel,allelesizela,Mgdlaa,alleledist,Ppla,Nallelepl,Nmissinggenotpl,Nincompletegenotpl,Nvalgenpl,Hepl,hTpl,vTpl,Dmpl,Dwmpl,Masizepl,Vasizepl);
	write_allele_freq(outputfile,n,m,namelocus,ndigit,Npop,Nip,namepop,StatType,NS,Stat,allelesizela,Nallelepl,Ppla,Hepl,hTpl,vTpl,Dmpl,Dwmpl,Masizepl,Vasizepl,Nmissinggenotpl,Nincompletegenotpl,Nvalgenpl,printallelefreq,FreqRef,Ngivenallelel,givenPla);
	
	ntot=n;	//restrict n (#indiv) to the indiv used
	if(cat1 && !cat2) { n=Nik[cat1]; for(p=1;p<=Npop;p++) if(catp[p]>1) Npop=p-1; }
	if(cat1 && cat2) { n=Nik[cat1]+Nik[cat2];	for(p=1;p<=Npop;p++) if(catp[p]>2) Npop=p-1; }

	if(StatType==1) mainAnalysisBtwInd(argc,n,ntot,xi,yi,zi,Mdij,sgi,Nsg,skgi,Nskg,catskg,cati,Ncat,Nik,nc,maxc,dijmin,dijmax,m,ndigit,ploidy,ploidyi,Nallelel,Nallelepl,Nvalgenpl,gilc,Ppla,allelesizela,Masizepl,Vasizepl,Mgdlaa,givenF,namei,namelocus,namecat,TypeComp,cat1,cat2,FreqRef,givenPla,Ngivenallelel,JKest,NS,Stat,printdistmatrix,sigmaest,density,dwidth,Npermut,permutalleles,writeresampdistr,regdetails,varcoef,Rbtwloc,permutdetails,outputfile);
	if(StatType>=2) mainAnalysisBtwPop(argc,StatType,n,xp,yp,zp,Mdij,popi,Npop,catp,cati,Ncat,nc,maxc,dijmin,dijmax,m,ndigit,ploidy,ploidyi,Nallelel,Nvalgenpl,gilc,Ppla,allelesizela,Masizepl,Vasizepl,Mgdlaa,namepop,namelocus,namecat,PWstat,TypeComp,cat1,cat2,FreqRef,JKest,NS,Stat,printdistmatrix,Npermut,permutalleles,writeresampdistr,regdetails,varcoef,Rbtwloc,permutdetails,outputfile);

	sprintf(smess,"\n\n ");
	write_tofile_only(outputfile,smess);
	printf("\n\nThe program has finished normally. \nLet's now open the file %c%s%c with a worksheet program\n... and digest SPAGeDi's results. \nPress any key to exit this program.\n",'"',outputfile,'"');
	if(argc<4) wait_a_char();			
	
	return 0;
}	//end of main



void mainAnalysisBtwPop(int argc,int StatType,int n,double *xp,double *yp,double *zp,double **Mdij,int *popi,int Npop,
	int *catp,int *cati,int Ncat,int nc,double *maxc,float dijmin,float dijmax,
	int m,int ndigit,int ploidy,int *ploidyi,int *Nallelel,int **Nvalgenkl,int ***gilc,
	float ***Ppla,int **allelesizela,float **Masizekl,float **Vasizekl,float ***Mgdlaa,
	struct name namepop[],char namelocus[][MAXNOM],struct name namecat[],
	int PWstat,int TypeComp,int cat1,int cat2,int FreqRef,int JKest,int NS,int Stat[],int printdistmatrix,
	int Npermut[],int permutalleles,int writeresampdistr,int regdetails,int varcoef,int Rbtwloc,int permutdetails,char outputfile[])
{
	int c,l,S,r,p1,p2;
	int linit,estinbreeding,maxnal;
	
	double *mdc,*mlndc;
	int *npc;
	float **indexpartic;

	float ***corrSlij[12],**corrSlc[12],**FstatSlr[12],***RSll[12],***V[12],**R2pl[12];

	float ***corrSlcp[12],***FstatSlrp[12];
	int Np,p,**allelesizelamix,***gilcmix,*popimix;
	float ***Mgdlaamix;
	double *xpmix,*ypmix,*zpmix,**Mdijmix;
	struct resample_stat_type **r_statSlc[12],**r_statFSlr[12]; //results of resampling for each loci and each dist class
	int comp,permutbypairs;
	char resampfile[50];			//name of file for resampling distrib
	float missdat=0.F;
	long seedinit;
	int F_Rstat,G_Nstat,pairwGst_Nst,pairwGij_Nij;  //=0 or 1 according to computations to do


	//init for dist classes analysis
	indexpartic=matrix(0,abs(nc)+2,0,1);
	npc=ivector(0,abs(nc)+2);
	mdc=dvector(0,abs(nc)+2);
	mlndc=dvector(0,abs(nc)+2);

	printf("\nComputing distance intervals. Please, wait.");
	checkdist(Npop,&nc,maxc,xp,yp,zp,Mdij,1,popi,catp,StatType,TypeComp,mdc,mlndc,npc,indexpartic);
	if(dijmax<=dijmin) dijmax=(float)maxc[nc];
	if(PWstat>1) displaydist(argc,outputfile,nc,maxc,mdc,mlndc,npc,indexpartic);

	//init for statistics computations
	linit=0;
	if(JKest) linit=-m;
	if(m==1) linit=1;
	for(S=1;S<=NS;S++){
		corrSlij[S]=f3tensor(linit,m,0,Npop,0,Npop);
		if(JKest) corrSlc[S]=matrix(-m,m+4,-20,nc+2);
		else corrSlc[S]=matrix(0,m+4,-20,nc+2);
		if(Rbtwloc) {RSll[S]=f3tensor(0,1,0,m,0,m); V[S]=f3tensor(0,1,0,5,0,2); R2pl[S]=matrix(0,1,-m,m);}
		FstatSlr[S]=matrix(0,m+2,0,4);
		for(l=0;l<=m+2;l++)for(r=0;r<=4;r++) FstatSlr[S][l][r]=(float)MISSVAL;
	}
	if(m==1) linit=1;
	else linit=0;

	if(ploidy>1) estinbreeding=1;

	F_Rstat=G_Nstat=pairwGst_Nst=pairwGij_Nij=0;
	for(S=1;S<=NS;S++){
		if(Stat[S]<=8) F_Rstat=1;
		if(Stat[S]>=9) G_Nstat=1;
		if(Stat[S]==9 || Stat[S]==11) pairwGst_Nst=1;
		if(Stat[S]==10 || Stat[S]==12) pairwGij_Nij=1;
	}


	//F-stat for pairs of spatial groups
	if(PWstat){
		printf("\nComputing pairwise differentiation between populations. Please, wait.");
		if(F_Rstat){
			for(p1=1;p1<Npop;p1++)for(p2=p1+1;p2<=Npop;p2++)if(TypeComp==0 || (TypeComp==1 && catp[p1]==catp[p2]) || (TypeComp==2 && catp[p1]!=catp[p2])) compute_F_R_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,allelesizela,NS,Stat,FstatSlr,corrSlij,0,JKest);
		}
		if(G_Nstat) compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,JKest,0);
		
		for(S=1;S<=NS;S++)compute_corr_per_dist_class (Npop,m,nc,maxc,Ncat,catp,StatType,TypeComp,FreqRef,varcoef,xp,yp,zp,Mdij,popi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],Rbtwloc,V[S],R2pl[S],&seed,JKest);
	}
	//global F-stat
	printf("\nComputing global differentiation among populations. Please, wait.");
	if(F_Rstat) compute_F_R_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,allelesizela,NS,Stat,FstatSlr,corrSlij,0,JKest);
	if(G_Nstat && !PWstat) compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,JKest,0);
	writePopStatresults(outputfile,Npop,m,namelocus,nc,maxc,npc,indexpartic,mdc,mlndc,dijmin,dijmax,StatType,TypeComp,cat1,cat2,namecat,NS,Stat,corrSlc,FstatSlr,JKest,regdetails,varcoef,Rbtwloc,RSll,V,R2pl,PWstat);

	if(printdistmatrix && !Npermut[0]) writedistmatrices(outputfile,Npop,m,0.,TypeComp,catp,printdistmatrix,xp,yp,zp,Mdij,popi,StatType,NS,Stat,corrSlij,namepop,namelocus);

	//PERMUTATIONS
	if(Npermut[0]){
		seedinit=seed;
		printf("\nComputing permutation tests. Please, wait.");
		for(S=1;S<=NS;S++){
			corrSlcp[S]=f3tensor(0,m,-2,nc+1,0,Npermut[0]);
			FstatSlrp[S]=f3tensor(0,m,1,4,0,Npermut[0]);
			r_statSlc[S]=resample_stat_type_matrix(0,m,-2,nc+1);
			r_statFSlr[S]=resample_stat_type_matrix(0,m,1,4);
			for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++) r_statSlc[S][l][c].mode=0;
			for(l=0;l<=m;l++)for(r=1;r<=4;r++) r_statFSlr[S][l][r].n=r_statFSlr[S][l][r].nd=0;
			for(l=0;l<=m;l++)for(r=1;r<=4;r++) {
				r_statFSlr[S][l][r].mode=(int)MISSVAL;
				r_statFSlr[S][l][r].obs=r_statFSlr[S][l][r].mean=r_statFSlr[S][l][r].sd=r_statFSlr[S][l][r].low95=r_statFSlr[S][l][r].high95=r_statFSlr[S][l][r].plow=r_statFSlr[S][l][r].phigh=r_statFSlr[S][l][r].pbil=(float)MISSVAL;
			}
		}
		for(S=1;S<=NS;S++) 
		for(S=1;S<=NS;S++) 

		//put obs val into p=0
		for(S=1;S<=NS;S++)for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++) corrSlcp[S][l][c][0]=corrSlc[S][l][c];
		for(S=1;S<=NS;S++)for(l=0;l<=m;l++) for(r=1;r<=4;r++) FstatSlrp[S][l][r][0]=FstatSlr[S][l][r];

		//LOCATIONS PERMUTATIONS
		if(PWstat)if((Np=Npermut[1])){
			printf("\nPermutations of spatial locations (%i)\n",Np);
			xpmix=dvector(0,n);
			ypmix=dvector(0,n);
			zpmix=dvector(0,n);
			if(Mdij[0][0]==1.) {Mdijmix=dmatrix(0,n,0,n); Mdijmix[0][0]=Mdij[0][0];}
			else Mdijmix=Mdij;
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(TypeComp!=1) permut_locations(Npop,xp,yp,zp,Mdij,xpmix,ypmix,zpmix,Mdijmix,&seed);
				if(TypeComp==1) permut_locations_within_cat(Npop,catp,Ncat,xp,yp,zp,Mdij,xpmix,ypmix,zpmix,Mdijmix,&seed);
				for(S=1;S<=NS;S++){
					compute_corr_per_dist_class (Npop,m,nc,maxc,Ncat,catp,StatType,TypeComp,FreqRef,0,xpmix,ypmix,zpmix,Mdijmix,popi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
					for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++){if(c==0)c=1; corrSlcp[S][l][c][p]=corrSlc[S][l][c];}
				}
			}
			free_dvector(xpmix,0,n);
			free_dvector(ypmix,0,n);
			free_dvector(zpmix,0,n);
			if(Mdij[0][0]==1.) free_dmatrix(Mdijmix,0,n,0,n);
			if(TypeComp!=1) for(S=1;S<=NS;S++)for(c=-2;c<=nc+1;c++){if(c==0)c=1; r_statSlc[S][0][c].mode=60;}
			if(TypeComp==1) for(S=1;S<=NS;S++)for(c=-2;c<=nc+1;c++){if(c==0)c=1; r_statSlc[S][0][c].mode=61;}
		}
		for(S=1;S<=NS;S++) analyse_resampling(m,-2,nc+1,Npermut[1],corrSlcp[S],r_statSlc[S]);

		//INDIVIDUALS PERMUTATIONS
		if((Np=Npermut[2])){
			popimix=ivector(0,n);
			printf("\nPermutations of individuals (%i)\n",Np);
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(TypeComp!=1) permut_indiv_among_pop(n,popi,Npop,popimix,&seed);
				if(TypeComp==1) permut_indiv_among_pop_within_categ(n,popi,Npop,cati,Ncat,popimix,&seed);
				if(F_Rstat) compute_F_R_stat(n,Npop,0,0,popimix,m,Nallelel,ploidy,gilc,allelesizela,NS,Stat,FstatSlr,corrSlij,0,0);
				if(G_Nstat) compute_G_N_stat(n,Npop,0,0,popimix,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,0,1);
				for(l=0;l<=m;l++)for(S=1;S<=NS;S++){FstatSlrp[S][l][3][p]=FstatSlr[S][l][3];FstatSlrp[S][l][4][p]=FstatSlr[S][l][4];}//Fst, Rho
			}
			free_ivector(popimix,0,n);
			if(TypeComp!=1)for(S=1;S<=NS;S++){r_statFSlr[S][0][3].mode=32; if(Stat[S]==3 || Stat[S]==4) r_statFSlr[S][0][4].mode=32;}
			if(TypeComp==1)for(S=1;S<=NS;S++){r_statFSlr[S][0][3].mode=33; if(Stat[S]==3 || Stat[S]==4) r_statFSlr[S][0][4].mode=33;}
		}
		for(S=1;S<=NS;S++) {analyse_resampling(m,3,3,Npermut[2],FstatSlrp[S],r_statFSlr[S]); if(Stat[S]==3 || Stat[S]==4) analyse_resampling(m,4,4,Npermut[2],FstatSlrp[S],r_statFSlr[S]);}

		//GENES PERMUTATIONS
		if((Np=Npermut[3])){ //make permutations for inbreeding coef
			gilcmix=i3tensor(0,n,0,m,0,ploidy-1);
			printf("\nPermutations of genes (%i)\n",Np);
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(TypeComp!=1) permut_genes_among_indiv(n,m,ploidy,ploidyi,gilc,gilcmix,&seed);
				if(TypeComp==1) /*permut genes among ind within cat*/permut_genes_among_indiv_within_pop(n,m,ploidy,ploidyi,cati,Ncat,gilc,gilcmix,&seed);
				if(F_Rstat) compute_F_R_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilcmix,allelesizela,NS,Stat,FstatSlr,corrSlij,0,0);
				for(S=1;S<=NS;S++)if(Stat[S]<=8)for(l=0;l<=m;l++) FstatSlrp[S][l][1][p]=FstatSlr[S][l][1];//Fit
				permut_genes_among_indiv_within_pop(n,m,ploidy,ploidyi,popi,Npop,gilc,gilcmix,&seed);
				if(F_Rstat) compute_F_R_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilcmix,allelesizela,NS,Stat,FstatSlr,corrSlij,0,0);
				for(S=1;S<=NS;S++)if(Stat[S]<=8)for(l=0;l<=m;l++) FstatSlrp[S][l][2][p]=FstatSlr[S][l][2];//Fis
			}
			free_i3tensor(gilcmix,0,n,0,m,0,ploidy-1);
			for(S=1;S<=NS;S++)if(Stat[S]<=8){
				if(TypeComp!=1) r_statFSlr[S][0][1].mode=20;
				if(TypeComp==1) r_statFSlr[S][0][1].mode=21;
				r_statFSlr[S][0][2].mode=22;
			}
		} //end of if estinbreeding
		for(S=1;S<=NS;S++)if(Stat[S]<=8) analyse_resampling(m,1,2,Npermut[3],FstatSlrp[S],r_statFSlr[S]);

		//WRITE LOC,POP,IND & GENE PERMUT TESTS
		if(Npermut[1] || Npermut[2] || Npermut[3]) WritePopPermutRes(outputfile,n,ploidy,maxc,m,namelocus,TypeComp,nc,NS,Stat,Npermut,0,permutdetails,r_statFSlr,r_statSlc,seedinit,PWstat);

		//ALLELE PERMUTATIONS
		if((Np=abs(Npermut[4]))){
			if(Npermut[4]<0) {Npermut[4]=abs(Npermut[4]); permutbypairs=1;}
			else permutbypairs=0;
			printf("\nPermutations of alleles (%i)\n",Np);
			seed=seedinit;
			maxnal=0;		//define the max # of alleles
			for(l=1;l<=m;l++) if(maxnal<Nallelel[l]) maxnal=Nallelel[l];
			allelesizelamix=imatrix(0,m,0,maxnal);
			Mgdlaamix=f3tensor(0,m,0,maxnal,0,maxnal);

			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(F_Rstat) permut_allelesizes(m,Nallelel,allelesizela,allelesizelamix,&seed);
				if(G_Nstat) permut_genetic_distances(m,Nallelel,Mgdlaa,Mgdlaamix,&seed);

				if(G_Nstat) compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaamix,NS,Stat,FstatSlr,corrSlij,0,0);
				if(F_Rstat) compute_F_R_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,allelesizelamix,NS,Stat,FstatSlr,corrSlij,1,0);
				for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
					for(l=0;l<=m;l++) for(r=1;r<=3;r++) FstatSlrp[S][l][r][p]=FstatSlr[S][l][r];
				}

				if(PWstat){
					if(F_Rstat) for(p1=1;p1<Npop;p1++)for(p2=p1+1;p2<=Npop;p2++)if(TypeComp==0 || (TypeComp==1 && catp[p1]==catp[p2]) || (TypeComp==2 && catp[p1]!=catp[p2])){
						compute_F_R_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,allelesizelamix,NS,Stat,FstatSlr,corrSlij,1,0);
					}
					for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
						compute_corr_per_dist_class (Npop,m,nc,maxc,Ncat,catp,StatType,TypeComp,FreqRef,0,xp,yp,zp,Mdij,popi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
						for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++) corrSlcp[S][l][c][p]=corrSlc[S][l][c];
					}
				}
			}

			for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
				for(c=-2;c<=nc+1;c++){if(Stat[S]<=8) r_statSlc[S][0][c].mode=10;  else r_statSlc[S][0][c].mode=11;}
				for(r=1;r<=3;r++){if(Stat[S]<=8) r_statFSlr[S][0][r].mode=10; else r_statFSlr[S][0][r].mode=11;}
				analyse_resampling(m,-2,nc+1,Npermut[4],corrSlcp[S],r_statSlc[S]);
				analyse_resampling(m,1,3,Npermut[4],FstatSlrp[S],r_statFSlr[S]);
			}
			WritePopPermutRes(outputfile,n,ploidy,maxc,m,namelocus,TypeComp,nc,NS,Stat,Npermut,1,permutdetails,r_statFSlr,r_statSlc,seedinit,PWstat);

			if(permutbypairs){
				if(npc[nc+1]==0) npc[nc+1]=(Npop*(Npop-1))/2;
				printf("\nPermutation of alleles for each pair of populations (%i pairs)\n",npc[nc+1]);
				comp=0;
				for(p1=1;p1<Npop;p1++)for(p2=p1+1;p2<=Npop;p2++)if(TypeComp==0 || (TypeComp==1 && catp[p1]==catp[p2]) || (TypeComp==2 && catp[p1]!=catp[p2])){
					for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){if(Stat[S]<=8) r_statSlc[S][0][0].mode=10; else  r_statSlc[S][0][0].mode=11;}
					comp++;
					if(npc[nc+1]>=100){if((comp%(npc[nc+1]/10))==0 || (comp<=(npc[nc+1]/10) && (comp%(npc[nc+1]/100))==0) ) printf(" %i",comp);}
					else if(npc[nc+1]>=10) {if((comp%(npc[nc+1]/10))==0) printf(" %i",comp);}
					else printf(" %i",comp);

 					if(F_Rstat) compute_F_R_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,allelesizela,NS,Stat,FstatSlr,corrSlij,1,0); //compute obs val
					if(pairwGst_Nst) compute_G_N_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,0,1);
					for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11)for(l=linit;l<=m;l++) corrSlcp[S][l][0][0]=corrSlij[S][l][p1][p2];	//transfer obs val into corrSlcp[S][l][0][0]
					if(pairwGij_Nij){
						compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,0,0);
						for(S=1;S<=NS;S++)if(Stat[S]==12) for(l=linit;l<=m;l++) corrSlcp[S][l][0][0]=corrSlij[S][l][p1][p2];
					}
					for(p=1;p<=Np;p++){
						if(F_Rstat){
							permut_allelesizes_among_2pop_alleles(m,Nallelel,Ppla,p1,p2,allelesizela,allelesizelamix,&seed);   
							compute_F_R_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,allelesizelamix,NS,Stat,FstatSlr,corrSlij,1,0);
							for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8) for(l=linit;l<=m;l++) corrSlcp[S][l][0][p]=corrSlij[S][l][p1][p2];
						}
						if(pairwGij_Nij){ 
							permut_genetic_distances(m,Nallelel,Mgdlaa,Mgdlaamix,&seed);
 							compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaamix,NS,Stat,FstatSlr,corrSlij,0,0);
							for(S=1;S<=NS;S++)if(Stat[S]==12) for(l=linit;l<=m;l++) corrSlcp[S][l][0][p]=corrSlij[S][l][p1][p2];
						}
						if(pairwGst_Nst){
							permut_genetic_distances_among_2pop_alleles(m,Nallelel,Ppla,p1,p2,Mgdlaa,Mgdlaamix,&seed);
							compute_G_N_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaamix,NS,Stat,FstatSlr,corrSlij,0,1);
							for(S=1;S<=NS;S++)if(Stat[S]==11) for(l=linit;l<=m;l++) corrSlcp[S][l][0][p]=corrSlij[S][l][p1][p2];
						}
					}
 					for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12) analyse_resampling(m,0,0,Npermut[4],corrSlcp[S],r_statSlc[S]);
					WriteAllelesPermutResForPopPair(outputfile,m,namelocus,p1,p2,namepop,NS,Stat,r_statSlc,permutdetails,comp);
				}
			}
 			free_imatrix(allelesizelamix,0,m,0,maxnal);
			free_f3tensor(Mgdlaamix,0,m,0,maxnal,0,maxnal);
		}

		//free memory
		for(S=1;S<=NS;S++){
			free_f3tensor(corrSlcp[S],0,m,-2,nc+1,0,Npermut[0]);
			free_f3tensor(FstatSlrp[S],0,m,1,4,0,Npermut[0]);
			free_resample_stat_type_matrix(r_statSlc[S],0,m,-2,nc+1);
			free_resample_stat_type_matrix(r_statFSlr[S],0,m,1,9);
		}

		if(printdistmatrix){
			if(G_Nstat) compute_G_N_stat(n,Npop,0,0,popi,m,Nallelel,ploidy,gilc,Ppla,Mgdlaa,NS,Stat,FstatSlr,corrSlij,JKest,0);
			if(F_Rstat) for(p1=1;p1<Npop;p1++)for(p2=p1+1;p2<=Npop;p2++)if(TypeComp==0 || (TypeComp==1 && catp[p1]==catp[p2]) || (TypeComp==2 && catp[p1]!=catp[p2])){
				compute_F_R_stat(n,2,p1,p2,popi,m,Nallelel,ploidy,gilc,allelesizela,NS,Stat,FstatSlr,corrSlij,0,JKest);
			}
			for(S=1;S<=NS;S++)compute_corr_per_dist_class (Npop,m,nc,maxc,Ncat,catp,StatType,TypeComp,FreqRef,varcoef,xp,yp,zp,Mdij,popi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
			writedistmatrices(outputfile,Npop,m,0.,TypeComp,catp,printdistmatrix,xp,yp,zp,Mdij,popi,StatType,NS,Stat,corrSlij,namepop,namelocus);
		}

	}//end of PERMUTATIONS



	linit=0;
	if(JKest) linit=-m;
	if(m==1) linit=1;
	for(S=1;S<=NS;S++){
		free_f3tensor(corrSlij[S],linit,m,0,Npop,0,Npop);
//		free_matrix(corrSlc[S],0,m+4,-20,nc+2);
//		free_matrix(FstatSlr[S],0,m+2,0,4);
		if(Rbtwloc) {free_f3tensor(RSll[S],0,1,0,m,0,m); free_f3tensor(V[S],0,1,0,5,0,2);}
	}

}


void mainAnalysisBtwInd(int argc,int n,int ntot,double *xi,double *yi,double *zi,double **Mdij,int *sgi,int Nsg,int *skgi,int Nskg,
	int *catskg,int *cati,int Ncat,int *Nik,int nc,double *maxc,float dijmin,float dijmax,
	int m,int ndigit,int ploidy,int *ploidyi,int *Nallelel,int **Nallelekl,int **Nvalgenkl,int ***gilc,
	float ***Ppla,int **allelesizela,float **Masizekl,float **Vasizekl,float ***Mgdlaa,float givenF,
	struct name namei[],char namelocus[][MAXNOM],struct name namecat[],
	int TypeComp,int cat1,int cat2,int FreqRef,float **givenPla,int *Ngivenallelel,int JKest,int NS,int Stat[],int printdistmatrix,float sigmaest,float density,float dwidth,
	int Npermut[],int permutalleles,int writeresampdistr,int regdetails,int varcoef,int Rbtwloc,int permutdetails,char outputfile[])
{
	int i,c,l,S,ni,nf;
	int linit,estinbreeding,maxnal;
	
	double *mdc,*mlndc;
	int *npc;
	float **indexpartic;

	float ***corrSlij[12],**corrSlc[12],***RSll[12],***V[12],**R2pl[12];

	float ***corrSlcp[12];
	int Np,p,Statbis[2],**allelesizelamix,***gilcmix,*sgimix;
	double *xmix,*ymix,*zmix,*xpmix,*ypmix,*zpmix,*xg,*yg,*zg,**Mdijmix;
	float ***Mgdlaamix;
	struct resample_stat_type **r_statSlc[12]; //results of resampling for each loci and each dist class

	char resampfile[50];			//name of file for resampling distrib
	float missdat=0.F;
	long seedinit;

	//add an intra-group distance class
	if(Nsg>1 && nc>=0) {
		for(c=nc;c>=1;c--) maxc[c+1]=maxc[c];
		maxc[1]=-1.;
		nc++;
	}
	//init for dist classes analysis
	indexpartic=matrix(0,abs(nc)+2,0,1);
	npc=ivector(0,abs(nc)+2);
	mdc=dvector(0,abs(nc)+2);
	mlndc=dvector(0,abs(nc)+2);

	//dist computations
	printf("\nComputing distance intervals. Please, wait.");
	checkdist(n,&nc,maxc,xi,yi,zi,Mdij,Nsg,sgi,cati,1,TypeComp,mdc,mlndc,npc,indexpartic);
	if(dijmax<=dijmin) dijmax=(float)maxc[nc];
	displaydist(argc,outputfile,nc,maxc,mdc,mlndc,npc,indexpartic);

	//init for statistics computations
	linit=0;
	if(JKest) linit=-m;
	if(m==1) linit=1;
	if(TypeComp!=1){
		if(NS*(m-linit+1)*n*n*4 > 250000000) {printf("\n\nSPAGeDi will attempt to allocate ca. %i Mo of RAM and may fail,\nin which case you should use a computer with more memory,\nselect just one statistic, or reduce your data set\nPress RETURN to go on",(int)(NS*(m-linit+1)*n*n*4./1000000.)); wait_a_char();}  
		for(S=1;S<=NS;S++) corrSlij[S]=f3tensor(linit,m,0,n,0,n);
	}
	if(TypeComp==1) for(S=1;S<=NS;S++){//create set of small matrices for within cat i,j pairs
		corrSlij[S]=fpmatrix(linit,m,0,n);
		for(l=linit;l<=m;l++) corrSlij[S][l][0]=vector(0,n);
		ni=1;nf=Nik[1];
		for(l=linit;l<=m;l++) corrSlij[S][l][1]=vector(ni,nf);
		for(i=2;i<=n;i++){
			if(cati[i]!=cati[i-1]){ni+=Nik[cati[i-1]];nf+=Nik[cati[i]];}
			for(l=linit;l<=m;l++) corrSlij[S][l][i]=vector(ni,nf);
		}
	}
	if(JKest) for(S=1;S<=NS;S++) corrSlc[S]=matrix(-m,m+4,-28,abs(nc)+2);
	else for(S=1;S<=NS;S++) corrSlc[S]=matrix(0,m+4,-28,abs(nc)+2);
	if(Rbtwloc)for(S=1;S<=NS;S++) {RSll[S]=f3tensor(0,1,0,m,0,m); V[S]=f3tensor(0,1,0,5,0,2); R2pl[S]=matrix(0,1,-m,m);}
	estinbreeding=0;
	for(S=1;S<=NS;S++) if(ploidy>1 && (Stat[S]==1 || Stat[S]==2 || Stat[S]==5)) estinbreeding=1; 
	
	//stat computations
	printf("\nComputing pairwise statistics between individuals. Please, wait.");
	compute_pairwise_corr_F(n,ntot,Ncat,cati,m,ndigit,ploidy,missdat,gilc,Nallelel,allelesizela,Mgdlaa,corrSlij,NS,Stat,FreqRef,givenPla,Ngivenallelel,TypeComp,givenF,0,JKest);
	for(S=1;S<=NS;S++) compute_corr_per_dist_class (n,m,nc,maxc,Ncat,cati,1,TypeComp,FreqRef,varcoef,xi,yi,zi,Mdij,sgi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],Rbtwloc,V[S],R2pl[S],&seed,JKest);
	if(sigmaest)for(S=1;S<=NS;S++){
		if(Stat[S]==1 || Stat[S]==2 || Stat[S]==11 || Stat[S]==4) estimate_sigma_2D_kinship (n,m,xi,yi,zi,Mdij,sgi,corrSlij[S],corrSlc[S],ploidy,Stat[S],JKest,density,dwidth);
	}
	writeIndStatresults(outputfile,n,Nsg,m,namelocus,nc,maxc,npc,indexpartic,mdc,mlndc,dijmin,dijmax,givenF,TypeComp,cat1,cat2,namecat,FreqRef,NS,Stat,corrSlc,density,dwidth,JKest,regdetails,varcoef,Rbtwloc,RSll,V,R2pl);
	if(printdistmatrix && !Npermut[0]) writedistmatrices(outputfile,n,m,givenF,TypeComp,cati,printdistmatrix,xi,yi,zi,Mdij,sgi,1,NS,Stat,corrSlij,namei,namelocus);

	//PERMUTATIONS
	if(Npermut[0]){
		printf("\nComputing permutation tests. Please, wait.");
		seedinit=seed;
		for(S=1;S<=NS;S++) corrSlcp[S]=f3tensor(0,m,-2,nc+1,0,Npermut[0]);
		for(S=1;S<=NS;S++) r_statSlc[S]=resample_stat_type_matrix(0,m,-2,nc+1);
		for(S=1;S<=NS;S++) for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++) r_statSlc[S][l][c].mode=0;


		//put obs val into p=0
		for(S=1;S<=NS;S++)for(l=0;l<=m;l++)for(c=-2;c<=nc;c++) corrSlcp[S][l][c][0]=corrSlc[S][l][c];

		//LOCATIONS PERMUTATIONS
		if((Np=Npermut[1])){
			printf("\nPermutations of spatial locations (%i)\n",Np);
			xmix=dvector(0,n);ymix=dvector(0,n);zmix=dvector(0,n);
			if(Mdij[0][0]==1.) {Mdijmix=dmatrix(0,n,0,n); Mdijmix[0][0]=Mdij[0][0];}
			else Mdijmix=Mdij;
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(Nsg<=1){
					if(TypeComp==1)	permut_locations_within_cat(n,cati,Ncat,xi,yi,zi,Mdij,xmix,ymix,zmix,Mdijmix,&seed);
					else permut_locations(n,xi,yi,zi,Mdij,xmix,ymix,zmix,Mdijmix,&seed);
				}
				if(Nsg>1){
					if(TypeComp==1) permut_locations_of_groups_within_cat (n,cati,Ncat,xi,yi,zi,Mdij,skgi,xmix,ymix,zmix,Mdijmix,&seed);
					else permut_locations_of_groups(n,xi,yi,zi,Mdij,sgi,xmix,ymix,zmix,Mdijmix,&seed);
				}
				for(S=1;S<=NS;S++){
					compute_corr_per_dist_class (n,m,nc,maxc,Ncat,cati,1,TypeComp,FreqRef,0,xmix,ymix,zmix,Mdijmix,sgi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
					for(l=0;l<=m;l++)for(c=-2;c<=nc;c++){if(c==0 && Nsg<=1)c=1;if(c==0 && Nsg>1)c=2; corrSlcp[S][l][c][p]=corrSlc[S][l][c];}
				}
			}

			free_dvector(xmix,0,n);free_dvector(ymix,0,n);free_dvector(zmix,0,n);
			if(Mdij[0][0]==1.) free_dmatrix(Mdijmix,0,n,0,n);
		
			if(Nsg<=1 && TypeComp!=1) for(c=-2;c<=nc;c++){if(c==0) c=1;for(S=1;S<=NS;S++) r_statSlc[S][0][c].mode=40;}
			if(Nsg<=1 && TypeComp==1) for(c=-2;c<=nc;c++){if(c==0) c=1;for(S=1;S<=NS;S++) r_statSlc[S][0][c].mode=41;}
			if(Nsg>1 && TypeComp!=1) for(c=-2;c<=nc;c++){if(c==0) c=2;for(S=1;S<=NS;S++) r_statSlc[S][0][c].mode=50;}
			if(Nsg>1 && TypeComp==1) for(c=-2;c<=nc;c++){if(c==0) c=2;for(S=1;S<=NS;S++) r_statSlc[S][0][c].mode=51;}

		}
		for(S=1;S<=NS;S++) analyse_resampling(m,-2,-1,Npermut[1],corrSlcp[S],r_statSlc[S]);
		if(Nsg<=1) for(S=1;S<=NS;S++) analyse_resampling(m,1,1,Npermut[1],corrSlcp[S],r_statSlc[S]);
		for(S=1;S<=NS;S++) analyse_resampling(m,2,nc,Npermut[1],corrSlcp[S],r_statSlc[S]);

		//INDIVIDUALS PERMUTATIONS
		if((Np=Npermut[2])){//permut of ind(positions) among spatial groups
			sgimix=ivector(0,n);
			printf("\nPermutations of individuals (%i)\n",Np);
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(TypeComp==1) permut_indiv_among_pop_within_categ(n,sgi,Nsg,cati,Ncat,sgimix,&seed);
				else permut_indiv_among_pop(n,sgi,Nsg,sgimix,&seed);
				for(S=1;S<=NS;S++){
					compute_corr_per_dist_class (n,m,nc,maxc,Ncat,cati,1,TypeComp,FreqRef,0,xi,yi,zi,Mdij,sgimix,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
					for(l=0;l<=m;l++) corrSlcp[S][l][1][p]=corrSlc[S][l][1];
				}
			}
			free_ivector(sgimix,0,n);
			if(TypeComp!=1) for(S=1;S<=NS;S++) r_statSlc[S][0][1].mode=30;
			else for(S=1;S<=NS;S++) r_statSlc[S][0][1].mode=31;
		}
		if(Nsg>1) for(S=1;S<=NS;S++) analyse_resampling(m,1,1,Npermut[2],corrSlcp[S],r_statSlc[S]);
	
		//GENES PERMUTATIONS
		if((Np=Npermut[3])){ //make permutations for inbreeding coef
			gilcmix=i3tensor(0,ntot,0,m,0,ploidy-1);
			for(i=1;i<=ntot;i++)for(l=1;l<=m;l++)for(c=0;c<ploidy;c++) gilcmix[i][l][c]=gilc[i][l][c];
			printf("\nPermutations of genes (%i)\n",Np);
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				if(TypeComp==1) permut_genes_among_indiv_within_pop(n,m,ploidy,ploidyi,cati,Ncat,gilc,gilcmix,&seed);
				else permut_genes_among_indiv(n,m,ploidy,ploidyi,gilc,gilcmix,&seed);
				compute_pairwise_corr_F(n,ntot,Ncat,cati,m,ndigit,ploidy,missdat,gilcmix,Nallelel,allelesizela,Mgdlaa,corrSlij,NS,Stat,FreqRef,givenPla,Ngivenallelel,TypeComp,givenF,1,0);
				for(S=1;S<=NS;S++){
					if(Stat[S]<=2 || Stat[S]==5){
						compute_corr_per_dist_class (n,m,-1,maxc,Ncat,cati,1,TypeComp,FreqRef,0,xi,yi,zi,Mdij,sgi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
						for(l=0;l<=m;l++) corrSlcp[S][l][0][p]=corrSlc[S][l][0];
					}
					else for(l=0;l<=m;l++) corrSlcp[S][l][0][p]=(float)MISSVAL;
				}
			}
			free_i3tensor(gilcmix,0,ntot,0,m,0,ploidy-1);
			if(TypeComp!=1) for(S=1;S<=NS;S++) r_statSlc[S][0][0].mode=20;
			if(TypeComp==1) for(S=1;S<=NS;S++) r_statSlc[S][0][0].mode=21;

		} //end of if estinbreeding
		for(S=1;S<=NS;S++) analyse_resampling(m,0,0,Npermut[3],corrSlcp[S],r_statSlc[S]);
		
		//WRITE RESULTS OF LOC, IND & GENE PERMUT TESTS
		if(Npermut[1] || Npermut[2] || Npermut[3]) WriteIndPermutRes(outputfile,n,ploidy,Ncat,maxc,m,namelocus,TypeComp,givenF,nc,Nsg,NS,Stat,Npermut,0,permutdetails,r_statSlc,seedinit);

		//ALLELE PERMUTATIONS
		if((Np=Npermut[4])){ 
			seed=seedinit;
			printf("\nPermutations of allele sizes (%i)\n",Np);
			Statbis[1]=5;
			maxnal=0;//define the max # of alleles
			for(l=1;l<=m;l++) if(maxnal<Nallelel[l]) maxnal=Nallelel[l];
			allelesizelamix=imatrix(0,m,0,maxnal);
			Mgdlaamix=f3tensor(0,m,0,maxnal,0,maxnal);
			for(p=1;p<=Np;p++){
				if(Np>=100){if((p%(Np/10))==0 || (p<=(Np/10) && (p%(Np/100))==0) ) printf(" %i",p);}
				else if((p%(Np/10))==0) printf(" %i",p);
				permut_allelesizes(m,Nallelel,allelesizela,allelesizelamix,&seed);
				permut_genetic_distances(m,Nallelel,Mgdlaa,Mgdlaamix,&seed);
				compute_pairwise_corr_F(n,ntot,Ncat,cati,m,ndigit,ploidy,missdat,gilc,Nallelel,allelesizelamix,Mgdlaamix,corrSlij,NS,Stat,FreqRef,givenPla,Ngivenallelel,TypeComp,givenF,0,0);
				for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==14){
					compute_corr_per_dist_class (n,m,nc,maxc,Ncat,cati,1,TypeComp,FreqRef,0,xi,yi,zi,Mdij,sgi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
					for(l=0;l<=m;l++)for(c=-2;c<=nc+1;c++) corrSlcp[S][l][c][p]=corrSlc[S][l][c];
				}
			}
			free_imatrix(allelesizelamix,0,m,0,maxnal);
			free_f3tensor(Mgdlaamix,0,m,0,maxnal,0,maxnal);
			
			for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==14){
				if(Stat[S]==5) for(c=-2;c<=nc+1;c++) r_statSlc[S][0][c].mode=10;
				if(Stat[S]==14) for(c=-2;c<=nc+1;c++) r_statSlc[S][0][c].mode=11;
				analyse_resampling(m,-2,nc,Npermut[4],corrSlcp[S],r_statSlc[S]);
			}
			WriteIndPermutRes(outputfile,n,ploidy,Ncat,maxc,m,namelocus,TypeComp,givenF,nc,Nsg,NS,Stat,Npermut,1,permutdetails,r_statSlc,seedinit);
		}

		//free memory
		for(S=1;S<=NS;S++) free_f3tensor(corrSlcp[S],0,m,-2,nc+1,0,Npermut[0]);
		for(S=1;S<=NS;S++) free_resample_stat_type_matrix(r_statSlc[S],0,m,-2,nc+1);

		if(printdistmatrix){
			compute_pairwise_corr_F(n,ntot,Ncat,cati,m,ndigit,ploidy,missdat,gilc,Nallelel,allelesizela,Mgdlaa,corrSlij,NS,Stat,FreqRef,givenPla,Ngivenallelel,TypeComp,givenF,0,0);
			for(S=1;S<=NS;S++) compute_corr_per_dist_class (n,m,nc,maxc,Ncat,cati,1,TypeComp,FreqRef,varcoef,xi,yi,zi,Mdij,sgi,dijmin,dijmax,corrSlij[S],corrSlc[S],RSll[S],0,V[S],R2pl[S],&seed,0);
			writedistmatrices(outputfile,n,m,givenF,TypeComp,cati,printdistmatrix,xi,yi,zi,Mdij,sgi,1,NS,Stat,corrSlij,namei,namelocus);
		}
	}//end of PERMUTATIONS


	linit=0;
	if(JKest) linit=-m;
	if(m==1) linit=1;
	if(TypeComp!=1) for(S=1;S<=NS;S++) free_f3tensor(corrSlij[S],linit,m,0,n,0,n);
//	for(S=1;S<=NS;S++) free_matrix(corrSlc[S],0,m+4,-20,abs(nc)+2);
	if(Rbtwloc)for(S=1;S<=NS;S++) {free_f3tensor(RSll[S],0,1,0,m,0,m); free_f3tensor(V[S],0,1,0,5,0,2);}
}




