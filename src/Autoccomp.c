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
#include "Xatools.h"
#include "Autoccomp.h"
#include "Autocio.h"

#define write write_string

void define_groups(int n,struct name namei[],double *xi,double *yi,double *zi,
		int Ncat,struct name namecat[],struct name namecati[],int *cati,int *Nik,
		int *Nsg,struct name namesg[],int *sgi,int *Nig,
		int *Nskg,struct name nameskg[],int *skgi,int *catskg,int *Niskg,
		int orig_order)
{
	int i,j,g,k,nsg,nskg;
	int ncatobs,newcat;
	int *newgroup;
	char smess[SMAX];

	/*define and count the different spatial groups and attribute a spatial group number to each indiv*/
	if(orig_order){
		sgi[1]=1;
 		nsg=1;
		Nig[1]=1;
		strncpy(namesg[1].n,namei[1].n,MAXNOM-1);
		for(i=2;i<=n;i++){
			if(xi[i]!=xi[i-1] || yi[i]!=yi[i-1] || zi[i]!=zi[i-1]){
				nsg++;
				sgi[i]=nsg;
				Nig[nsg]=1;
				sprintf(smess,namei[i].n,nsg);
				strncpy(namesg[nsg].n,smess,MAXNOM-1);
			}
			else{
				sgi[i]=nsg;
				Nig[nsg]++;
			}
		}
		*Nsg=nsg;
		if(*Nsg==n) *Nsg=0;
	}

	if(!orig_order){
		newgroup=ivector(0,n);
		for(g=0;g<=n;g++) newgroup[g]=0;

		newgroup[sgi[1]]=1;
		nsg=1;
		strncpy(namesg[1].n,namei[1].n,MAXNOM-1);
		for(i=2;i<=n;i++){
			if(newgroup[sgi[i]]==0){
				nsg++;
				newgroup[sgi[i]]=nsg;
				strncpy(namesg[nsg].n,namei[i].n,MAXNOM-1);
			}
		}
		*Nsg=nsg;

		for(g=0;g<=n;g++)Nig[g]++;
		for(i=1;i<=n;i++){
			sgi[i]=newgroup[sgi[i]];
			Nig[sgi[i]]++;
		}

		free_ivector(newgroup,0,n);

		if(*Nsg==n) *Nsg=0;
	}					 


	/*define and count the # of categories and attribute a category number to each indiv*/
	if(Ncat){
		ncatobs=1;
		cati[1]=1;
		strncpy(namecat[1].n,namecati[1].n,MAXNOM-1);
		for(i=2;i<=n;i++){
			newcat=1;
			for(k=1;k<=ncatobs;k++){
				if(strncmp(namecati[i].n,namecat[k].n,MAXNOM-1)==0) {cati[i]=k;newcat=0;}
			}
			if(newcat){
				ncatobs++;
				strncpy(namecat[ncatobs].n,namecati[i].n,MAXNOM-1);
				cati[i]=ncatobs;
			}
		}
		if(Ncat!=ncatobs && orig_order){
			printf("\n\nWARNING: %i categories observed for %i expected\nCheck the different categories found in the file error.txt\nPress any key to stop the program now.",ncatobs,Ncat);
			sprintf(smess,"\nWARNING: %i categories observed for %i expected",ncatobs,Ncat);
			write(ERRORFILE,smess);
			sprintf(smess,"\nThe following category names were found:");
			write(ERRORFILE,smess);
			for(k=1;k<=ncatobs;k++) {
				sprintf(smess,"\n%s",namecat[k].n);
				write(ERRORFILE,smess);
			}
			wait_a_char();
			exit(2);
		}

		Nik[0]=n;
		for(k=1;k<=Ncat;k++) Nik[k]=0;
		for(i=1;i<=n;i++) Nik[cati[i]]++;

	}/*end of if(Ncat)*/
	else for(i=1;i<=n;i++) cati[i]=0;


	/*count and define the spatio-categorical groups*/
	if(Ncat>1 && *Nsg>1){
		nskg=1;
		skgi[1]=1;
		Niskg[1]=1;
		catskg[1]=cati[1];
		sprintf(smess,"%s-%s",namesg[1].n,namecat[1].n);
		strncpy(nameskg[1].n,smess,MAXNOM-1);

		for(i=2;i<=n;i++){
			skgi[i]=0;
			if(sgi[i]==sgi[i-1]){
				j=i-1;
				do{
					if(cati[i]==cati[j]){
						skgi[i]=skgi[j];
						Niskg[skgi[i]]++;
					}
					j--;
				}while(sgi[j]==sgi[i] && j>0 && skgi[i]==0);
			}
			if(skgi[i]==0){
				nskg++;
				skgi[i]=nskg;
				catskg[nskg]=cati[i];
				Niskg[nskg]=1;
				sprintf(smess,"%s-%s",namesg[sgi[i]].n,namecat[cati[i]].n);
				strncpy(nameskg[nskg].n,smess,MAXNOM-1);
			}
		}
		*Nskg=nskg;

	}
	else *Nskg=0;


}/*end of define_groups*/

/*******************************************************************************/

/*function that reorder the data (i.e. individuals i=1->n) begining	with 
category 1,..., or by the category(ies) selected for within 1 cat or between 2 cat
analyses*/
void reorder_data_by_category(int n,int Ncat,int *cati,int *cat1,int *cat2,
			int *ploidyi,int *sgi,int *skgi,double *xi,double *yi,double *zi,
			struct name *namei,struct name *namecati,int ***gilc) 
{
	int *newii,*oldii,newi;
	int i,k;
	
	/*reorder individuals by categories*/
	newii=ivector(0,n);
	oldii=ivector(0,n);
	newi=0;
	/*put cat1 and cat2 in first place if necessary*/
	if(*cat1)for(i=1;i<=n;i++)if(cati[i]==*cat1) {newii[i]=(++newi);	oldii[newi]=i;}
	if(*cat2)for(i=1;i<=n;i++)if(cati[i]==*cat2) {newii[i]=(++newi);	oldii[newi]=i;}
	for(k=1;k<=Ncat;k++)if(k!=*cat1 && k!=*cat2)for(i=1;i<=n;i++)if(cati[i]==k) {newii[i]=(++newi);	oldii[newi]=i;}
	if(*cat1) *cat1=1;
	if(*cat2) *cat2=2;
	for(i=1;i<=n;i++){ /*loop over newi*/
		ploidyi[0]=ploidyi[i];ploidyi[i]=ploidyi[oldii[i]];ploidyi[oldii[i]]=ploidyi[0];
		cati[0]=cati[i];cati[i]=cati[oldii[i]];cati[oldii[i]]=cati[0];
		sgi[0]=sgi[i];sgi[i]=sgi[oldii[i]];sgi[oldii[i]]=sgi[0];
		skgi[0]=skgi[i];skgi[i]=skgi[oldii[i]];skgi[oldii[i]]=skgi[0];
		xi[0]=xi[i];xi[i]=xi[oldii[i]];xi[oldii[i]]=xi[0];
		yi[0]=yi[i];yi[i]=yi[oldii[i]];yi[oldii[i]]=yi[0];
		zi[0]=zi[i];zi[i]=zi[oldii[i]];zi[oldii[i]]=zi[0];
		strcpy(namei[0].n,namei[i].n);strcpy(namei[i].n,namei[oldii[i]].n);strcpy(namei[oldii[i]].n,namei[0].n);
		strcpy(namecati[0].n,namecati[i].n);strcpy(namecati[i].n,namecati[oldii[i]].n);strcpy(namecati[oldii[i]].n,namecati[0].n);
		gilc[0]=gilc[i];gilc[i]=gilc[oldii[i]];gilc[oldii[i]]=gilc[0];
		
		oldii[newii[i]]=oldii[i];
		newii[oldii[i]]=newii[i];
	}
	free_ivector(newii,0,n);
	free_ivector(oldii,0,n);
	/*limit the data to specified cat*/
}					

/*******************************************************************************/
/*function that defines "populations" according to StatType (i.e. pop are either
spatial groups, categorical groups, spatio-categorical groups. 
For StatType=1 (analysis at ind level), "populations" are defined as categories*/

void define_pop(int StatType,int n,int Npop,int Nsg,int Ncat,int Nskg,
		int *popi,int *sgi,int *cati,int *skgi,int *Nip,int *Nig,int *Nik,int *Niskg,
		int *catskg,int *catp,struct name *namepop,struct name *namesg,struct name *namecat,struct name *nameskg,
		double *xp,double *yp,double *zp,double *xi,double *yi,double *zi)
{
	int i,p;


	if(StatType<=2) for(i=1;i<=n;i++) popi[i]=cati[i];
	if(StatType==3) for(i=1;i<=n;i++) popi[i]=sgi[i];
	if(StatType==4) for(i=1;i<=n;i++) popi[i]=skgi[i];


	if(StatType<=2) for(p=1;p<=Npop;p++) strncpy(namepop[p].n,namecat[p].n,MAXNOM-1); 
	if(StatType==3) for(p=1;p<=Npop;p++) strncpy(namepop[p].n,namesg[p].n,MAXNOM-1); 
	if(StatType==4) for(p=1;p<=Npop;p++) strncpy(namepop[p].n,nameskg[p].n,MAXNOM-1); 

	if(StatType<=2) for(p=1;p<=Npop;p++) Nip[p]=Nik[p];  
	if(StatType==3) for(p=1;p<=Npop;p++) Nip[p]=Nig[p];  
	if(StatType==4) for(p=1;p<=Npop;p++) {Nip[p]=Niskg[p]; catp[p]=catskg[p]; }
	Nip[0]=0;
	for(p=1;p<=Npop;p++) Nip[0]+=Nip[p];
	if(StatType==1) Nip[0]=n;
	/*compute (average) spatial coordinates of each pop*/
	for(p=0;p<=Npop;p++) xp[p]=yp[p]=zp[p]=0.;
	for(i=1;i<=n;i++) {xp[popi[i]]+=xi[i];yp[popi[i]]+=yi[i];zp[popi[i]]+=zi[i];}
	for(p=0;p<=Npop;p++) if(Nip[p]) {xp[p]/=Nip[p];yp[p]/=Nip[p];zp[p]/=Nip[p];}

}

/*******************************************************************************/

/*Function that checks if the last class of distance interval includes 
the largest possible distance among individuals (if not it adds
one class and add_class=1); it computes the number of pairs of individuals for
each class of distance intervals (npTc); and it computes the mean observed
distance between pairs of individual for each class (mdTc)*/  
void checkdist(int n, int *nc, double *maxc, double *xi, double *yi, double *zi,double **Mdij,
			   int Nsg,int *sgi,int *cati,int StatType,int TypeComp,
			   double *mdc,double *mlndc, int *npc,float **indexpartic)
{
	int i,j,c;	/*current individual (i,j), class of distance interval (c)*/
	int add_classSup=0;
	double dij,lndij,dmax;	/*distance btw i & j, ln(dij),maximal distance between individuals*/
	double sumd[MAXINTERVALS+3],sumlnd[MAXINTERVALS+3];/*sum of distances (or ln(dist)) within class*/
	int np,nplnc[MAXINTERVALS+3];	   
	int **participic;
	int p;
	double *distp;
	
	/*TypeComp=0 -> comp btw all ind; TypeComp=1 -> comp btw ind of same category; TypeComp=2 -> comp btw ind of different categories*/ 
	if((*nc)<0){ /*define -nc class with equal # pairs of individuals per class*/
		(*nc)=-(*nc);
		distp=dvector(0,(n*(n-1))/2);
		p=1;
		for(i=1;i<n;i++)for(j=i+1;j<=n;j++){
			if(TypeComp==1)if(cati[i]!=cati[j]) continue;
			if(TypeComp==2)if(cati[i]==cati[j]) continue;
			if(StatType==1 && Mdij[0][0]==0)if(sgi[i]==sgi[j]) continue; /*don't take into account intra-group class*/
			if(Mdij[0][0]==1.){
				if(Mdij[i][j]!=(float)MISSVAL) dij=Mdij[i][j];
				else continue;
			}
			else if(Mdij[0][0]==-1.) dij=acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
			else dij=sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
			
			distp[p++]=dij;
		}
		p--;
		shellsortd(p,distp);
		if(Nsg>1 && StatType==1){
			maxc[1]=-1.;  /*add an intra-group class*/
			(*nc)++; 
			for(c=2;c<=(*nc);c++) maxc[c]=distp[(int)(p*(c-1)/((*nc)-1))];
		}
		else for(c=1;c<=(*nc);c++) maxc[c]=distp[(int)(p*c/(*nc))];
		maxc[(*nc)]*=1.0001;
		free_dvector(distp,0,(n*(n-1))/2);
	}

	/*initialization*/
	for(c=0;c<=(*nc)+1;c++){
		sumd[c]=sumlnd[c]=0.0F;
		npc[c]=nplnc[c]=0;
	}
	dmax=0.0F;
 	participic=imatrix(0,n,0,(*nc)+1);
	for(i=0;i<=n;i++)for(c=0;c<=(*nc)+1;c++)participic[i][c]=0;

	/*consider all pairs of individuals i & j*/
	for(i=1;i<n;i++){
		for(j=(i+1);j<=n;j++){
			if(TypeComp==1)if(cati[i]!=cati[j]) continue;
			if(TypeComp==2)if(cati[i]==cati[j]) continue;

			if(Mdij[0][0]==1.){
				if(Mdij[i][j]!=(float)MISSVAL) dij=Mdij[i][j];
				else continue;
			}
			else if(Mdij[0][0]==-1.) dij=acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
			else dij=sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );

			/*check if a new class interval must be added*/
			if((dij>maxc[(*nc)])&&(add_classSup==0) /*&&(*nc>0)*/ ){
				(*nc)++;
				add_classSup=1;
			}
			/*check the maximal distance class*/
			if(dij>dmax){
				dmax=dij*1.000001;
				if(add_classSup) maxc[(*nc)]=dmax;/*redefine the added class interval*/
				/* if(*nc==0) maxc[1]=dmax;		*/
			}
			/*define the distance interval containing dij*/
			if(StatType==1 && (sgi[i]==sgi[j]) && Mdij[0][0]==0) c=1;	/*intra-group class*/
			else{
				c=1; 
				while(dij>maxc[c] && c<=(*nc)) c++;
			}
			/*add one pair to this class and the distance to the sum*/
			npc[c]++;
			sumd[c]+=dij;
			if(dij>0.){
				lndij=log(dij);
				nplnc[c]++;
				sumlnd[c]+=lndij;
			}
			participic[i][c]++;
			participic[j][c]++;
		}/*end loop j*/
	}/*end loop i*/

	npc[(*nc)+1]=0;
	for(c=1;c<=(*nc);c++){
		if(npc[c]) mdc[c]=sumd[c]/(double)npc[c];
		else mdc[c]=(float)MISSVAL;
		if(nplnc[c]) mlndc[c]=sumlnd[c]/(double)nplnc[c];
		else mlndc[c]=(float)MISSVAL;
		npc[(*nc)+1]+=npc[c];
	}
	mdc[0]=mlndc[0]=(float)MISSVAL;


	/*count the % of individuals participating to each class and the CV of the participating rate*/
	for(c=1;c<=(*nc);c++)if(npc[c]){
		indexpartic[c][0]=indexpartic[c][1]=0.;
		for(i=1;i<=n;i++){
			if(participic[i][c]) indexpartic[c][0]++;
			indexpartic[c][1]+=(float)(sqr(participic[i][c]-(2.*npc[c]/(double)n))/(double)(n-1)); /*var of # de fois qu'on partic*/
		}
		indexpartic[c][0]/=n;
		indexpartic[c][1]=(float)(sqrt(indexpartic[c][1])/(2.*npc[c]/(double)n)); /*coef var du # de fois que chaque ind partic*/
	}
	indexpartic[0][0]=1.;
	indexpartic[0][1]=0.;
	
	free_imatrix(participic,0,n,0,(*nc)+1);	  

}	/*end procedure checkdist*/




/****************************************************************************/


/*Function computing allele frequencies*/ 
void compute_allele_freq(int n,int Ncat,int *cati,int m,
			int ndigit,int ploidy,int ***gilc,int *ploidyi,int *Nallelel,int **allelesizela,float ***Mgdlaa,
			int alleledist,float ***Pkla,int **Nallelekl,int **Nmissinggenotkl,int **Nincompletegenotkl,
			int **Nvalgenkl,float **Hekl,float **hTkl,float **vTkl,float **Dmkl,float **Dwmkl,float **Masizekl,float **Vasizekl)
{
	int i,l,a,a1,a2,ncomp,k,k1,k2,npairs,nallelepairs,nloci;/*counter for individuals (i,j), locus (l), group (g), allele (a), category*/ 
	float Navalid,nvalidpop,sumweight;		/*number of individuals with valid data*/
	float Sasize,SSasize;	/*for estimators of allele size coef*/
	float **Jkl;
	int *Nik;
	int maxnal;


	//Pkla[k][l][a]: freq of allele a at locus l in pop k, over all pop when k=0 (weighted by # of defined genes per pop) and k=-1 (arithmetic average over pop). When a=0, it gives the number of defined genes.

	/*count # ind per cat*/
	Nik=ivector(0,Ncat);
	for(k=0;k<=Ncat;k++) Nik[k]=0;
	for(i=1;i<=n;i++){
		Nik[0]++;
		// TODO: where should the else statement attach?  I'm resolving it like the compiler would
		if(Ncat) {
			if(cati[i]<=Ncat)
				Nik[cati[i]]++;
			else
				exit(133);
		}
	}
	maxnal=0;
	for(l=1;l<=m;l++) if(Nallelel[l]>maxnal) maxnal=Nallelel[l];
	Jkl=matrix(0,Ncat,0,m);
	
	/*determine allele frequencies (allele '0' correspond to a missing data)*/
	for(k=0;k<=Ncat;k++) for(l=0;l<=m;l++) Nvalgenkl[k][l]=Nmissinggenotkl[k][l]=Nincompletegenotkl[k][l]=0;
	for(k=-1;k<=Ncat;k++) for(l=0;l<=m;l++) for(a=0;a<=Nallelel[l];a++) Pkla[k][l][a]=0.;

	for(i=1;i<=n;i++){
		for(l=1;l<=m;l++){
			Navalid=(float)ploidyi[i];
			for(a=0;a<ploidyi[i];a++) if(gilc[i][l][a]==0) Navalid--;
			if(Navalid){
				for(a=0;a<ploidyi[i];a++)if(gilc[i][l][a]){
					/*count # of valid alleles*/
					Nvalgenkl[0][l]++;
					if(Ncat) Nvalgenkl[cati[i]][l]++;
					/*add ind allele freq*/
					Pkla[0][l][gilc[i][l][a]]+=1.0f/Navalid;
					if(Ncat>1) Pkla[cati[i]][l][gilc[i][l][a]]+=1.0f/Navalid;
				}
			}
			else{	/*count number of missing genotypes*/
				Pkla[0][l][0]++;
				Nmissinggenotkl[0][l]++;
				if(Ncat>1) Pkla[cati[i]][l][0]++; 
			}
			if(Navalid>0 && Navalid<ploidyi[i]) Nincompletegenotkl[0][l]++;
			if(Ncat){
				if(Navalid==0) Nmissinggenotkl[cati[i]][l]++;
				if(Navalid>0 && Navalid<ploidyi[i]) Nincompletegenotkl[cati[i]][l]++;
			}
		}
	}/*end of loop i*/

	for(l=1;l<=m;l++) for(a=1;a<=Nallelel[l];a++) for(k=0;k<=Ncat;k++) Pkla[k][l][a]/=(Nik[k]-Pkla[k][l][0]);

	for(l=1;l<=m;l++) for(k=0;k<=Ncat;k++) Pkla[k][l][0]=(float)Nvalgenkl[k][l];

	//compute mean allele freq over pop without weighting by the #ind/pop (put in k=-1)
	for(l=1;l<=m;l++){
		nvalidpop=0.;
		for(k=1;k<=Ncat;k++)if(Nvalgenkl[k][l]){
			nvalidpop++;
			for(a=1;a<=Nallelel[l];a++) Pkla[-1][l][a]+=Pkla[k][l][a];
		}
		for(a=1;a<=Nallelel[l];a++) Pkla[-1][l][a]/=nvalidpop;
	}


	/*determine the number of alleles per locus*/
	for(l=0;l<=m;l++) for(k=0;k<=Ncat;k++) Nallelekl[k][l]=0;
	for(l=1;l<=m;l++) for(a=1;a<=Nallelel[l];a++)for(k=0;k<=Ncat;k++)if(Pkla[k][l][a]) Nallelekl[k][l]++;
		
	/*determine the expected heterozygosity (gene diversity) per locus Hekl*/
	for(k=0;k<=Ncat;k++) Hekl[k][0]=0.;
	for(l=1;l<=m;l++)for(k=0;k<=Ncat;k++){ 
		Jkl[k][l]=0.;
		for(a=1;a<=Nallelel[l];a++) Jkl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a];/*homozygosity*/
		Hekl[k][l]=1.0f-Jkl[k][l];
		Hekl[k][l]*=(float)(Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.0f));//corrected heterozygosity		 
		Hekl[k][0]+=Hekl[k][l]/m;
	}

	//determine diversity of unordered and ordered alleles per locus hTkl, vTkl
	for(l=1;l<=m;l++){
		//within each pop
		for(k=1;k<=Ncat;k++){
			if(Nvalgenkl[k][l]>1){
				hTkl[k][l]=vTkl[k][l]=Dmkl[k][l]=Dwmkl[k][l]=sumweight=0.0f;
				nallelepairs=0;
				for(a1=1;a1<=Nallelel[l];a1++) if(Pkla[k][l][a1]) for(a2=1;a2<=Nallelel[l];a2++){
					if(a1!=a2) hTkl[k][l]+=Pkla[k][l][a1]*Pkla[k][l][a2];
					if(alleledist){
						vTkl[k][l]+=Mgdlaa[l][a1][a2]*Pkla[k][l][a1]*Pkla[k][l][a2];
						if(a1!=a2 && Pkla[k][l][a1] && Pkla[k][l][a2]){
							nallelepairs++;
							Dmkl[k][l]+=Mgdlaa[l][a1][a2];
							Dwmkl[k][l]+=Mgdlaa[l][a1][a2]*Pkla[k][l][a1]*Pkla[k][l][a2];
							sumweight+=Pkla[k][l][a1]*Pkla[k][l][a2];
						}
					}
				}
				hTkl[k][l]*=(float)(Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.));//corrected diversity		 
				vTkl[k][l]*=(float)(Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.));//corrected diversity
				if(nallelepairs) Dmkl[k][l]/=nallelepairs;
				else Dmkl[k][l]=(float)MISSVAL;
				if(sumweight) Dwmkl[k][l]/=sumweight;
				else Dwmkl[k][l]=(float)MISSVAL;
			} 
			else hTkl[k][l]=vTkl[k][l]=Dmkl[k][l]=Dwmkl[k][l]=(float)MISSVAL;
		}
		
		//total
		hTkl[0][l]=vTkl[0][l]=Dmkl[0][l]=Dwmkl[0][l]=0.0f;
		npairs=0;
		for(k1=1;k1<Ncat;k1++)for(k2=k1+1;k2<=Ncat;k2++){
			if(Nvalgenkl[k1][l]>1 && Nvalgenkl[k2][l]>1){
				npairs++;
				for(a1=1;a1<=Nallelel[l];a1++) if(Pkla[k1][l][a1]) for(a2=1;a2<=Nallelel[l];a2++) if(Pkla[k2][l][a2]){
					if(a1!=a2) hTkl[0][l]+=Pkla[k1][l][a1]*Pkla[k2][l][a2];
					if(alleledist) vTkl[0][l]+=Mgdlaa[l][a1][a2]*Pkla[k1][l][a1]*Pkla[k2][l][a2];
				}
			}
		}
		hTkl[0][l]/=npairs;	
		vTkl[0][l]/=npairs;
		
		nallelepairs=0;
		sumweight=0.0f;
		if(alleledist) for(a1=1;a1<=Nallelel[l];a1++) for(a2=1;a2<=Nallelel[l];a2++){
			if(a1!=a2 && Pkla[-1][l][a1] && Pkla[-1][l][a2]){
				nallelepairs++;
				Dmkl[0][l]+=Mgdlaa[l][a1][a2];
				Dwmkl[0][l]+=Mgdlaa[l][a1][a2]*Pkla[-1][l][a1]*Pkla[-1][l][a2];
				sumweight+=Pkla[-1][l][a1]*Pkla[-1][l][a2];
			}
		}
		if(nallelepairs) Dmkl[0][l]/=nallelepairs;
		else Dmkl[0][l]=(float)MISSVAL;
		if(sumweight) Dwmkl[0][l]/=sumweight;
		else Dwmkl[0][l]=(float)MISSVAL;

	} //end loop over locus

	for(k=0;k<=Ncat;k++){
		nloci=0;
		hTkl[k][0]=vTkl[k][0]=Dmkl[k][0]=Dwmkl[k][0]=0.0f;	
		for(l=1;l<=m;l++)if(Nvalgenkl[k][l]>1){
			nloci++;
			hTkl[k][0]+=hTkl[k][l];
			vTkl[k][0]+=vTkl[k][l];
			Dmkl[k][0]+=Dmkl[k][l];
			Dwmkl[k][0]+=Dwmkl[k][l];
		}
		hTkl[k][0]/=nloci;
		vTkl[k][0]/=nloci;
		Dmkl[k][0]/=nloci;
		Dwmkl[k][0]/=nloci;
	}


	/*determine the mean and variance of allele size*/
	for(l=1;l<=m;l++){
		for(k=0;k<=Ncat;k++){
			Sasize=SSasize=0.;
			for(a=1;a<=Nallelel[l];a++){
				Sasize+=Pkla[k][l][a]*allelesizela[l][a];
				SSasize+=Pkla[k][l][a]*allelesizela[l][a]*allelesizela[l][a];
			}
			Masizekl[k][l]=Sasize;
			Vasizekl[k][l]=(float)((SSasize-Sasize*Sasize)*((double)Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.)));
		}
	}

	free_matrix(Jkl,0,Ncat,0,m);
	free_ivector(Nik,0,Ncat);


}/*end procedure compute_allele_freq*/


/****************************************************************************/

/*Function computing the pairwise kinship (or related coef =stat) for each pair of 
individuals i,j (j>i) and each locus (corrSlij[stat][l][i][j]), 
as well as for all loci in corrSlij[stat][0][][], and all loci but locus l in corrSlij[stat][-l][][]. 
The inbreeding coefficient for ind i is in corrSlij[][][i][i]. 

If no correlation can be computed (no variation, missing value,...), value =(float)MISSVAL.*/ 

void compute_pairwise_corr_F(int n,int ntot,int Ncat,int *cati,int m,int ndigit,int ploidy,
			float missdat,int ***gilc,int *N1allelel,int **allelesizela,float ***Distla1a2,
			float ***corrSlij[],int NS,int Stat[12],int FreqRef,float **givenPla,int *Ngivenallelel,
			int TypeComp,float givenF,int compute_inbreeding_coef_only,int JKl)
{
	int i,j,k,l,linit,g,ci,cj,a,a1,a2,maxa,S,c;//counter for individuals (i,j), locus (l), gene (g), allele (a), max value for a at locus l, type of statistic (S=1 for Loiselle, 2 for Ritland original, 3 for Wright coef of relationship, 4 for Rousset) 
	int newSumNall; //sums of the number of alleles per locus
	float Pnew;	//corrected (Ritland original estimator) allele frequency
	int newNloci;	//number of valid loci excluded those monomorphic
	int Ncomp,NcompRitl,Ncompasc[MMAX],NcompNij[MMAX];		//comp=1 if there is no missing dat, otherwize=0
	double Ngivenallele,Nvalidallele;
	float ***Pkla,**Masizekl,**Vasizekl,**Plak[MMAX];
	int **Nallelekl,**Nvalgenkl,*Nivalidk,Ncvalid,ploidyi,Nallelel[MMAX];
	double Sasize,SSasize;
	float SoP,SoPt,Ritl=(float)MISSVAL,Ritlt,RitlN,RitlD;
	float Rous=(float)MISSVAL,**LoisD=NULL,Lois=(float)MISSVAL,Loist,LoisDt,LoisN,LoisNt,**DivN=NULL,Nij=(float)MISSVAL,SNij,SNijcomp;	//homozygosity per locus, denominator in Loiselle eq
	float **h2kl,**RdomD,Rdom=(float)MISSVAL,Rdomt,RdomDt,RdomN,RdomNt;
	float **Hokl,Hw,Hwt,Hb,Hbt,Hot;   //for estimator of Rousset, 1999
	int **Nvalidwpairkl;
	float ascor=(float)MISSVAL,asc,Sasc,SNcomp;	//for estimators of allele size coef
	float Pa,Pb,Pc,Pd,Sab,Sac,Sad,Sbc,Sbd,Sca,Scb,Scd,Sda,Sdb,rLynch=(float)MISSVAL,dLynch=(float)MISSVAL,rLynch1,rLynch2,dLynch1,dLynch2,WrLynch,WdLynch,WrLynch1,WrLynch2,WdLynch1,WdLynch2,SrLynch1,SrLynch2,SWrLynch,SWrLynch1,SWrLynch2,SWdLynch,SdLynch1,SdLynch2,SWdLynch1,SWdLynch2,*SrLynch1l,*SrLynch2l,*SWrLynch1l,*SWrLynch2l,*SdLynch1l,*SdLynch2l,*SWdLynch1l,*SWdLynch2l;
	float So[MMAX],Sij,rli,Srli[MMAX],SSrli[MMAX],Nrli[MMAX],Var,MLrli,SumW; 
	float rQueller=(float)MISSVAL,rQueller1,rQueller2,WrQueller,WrQueller1,WrQueller2,SrQueller1,SrQueller2,SWrQueller1,SWrQueller2,*SrQueller1l,*SrQueller2l,*SWrQueller1l,*SWrQueller2l;
	int Elois,Eritl,Erous,Erela,Easc,ENij,Erlynch,Edlynch,Erqueller,Edwang,Erwang,Ekinshipdom,Erelatdom,Erli;   //define if the statistic is asked (0=no), value>0 = value of S
	float **Plai[MMAX],***varkla,pi,qi,*Spik,*SSpik,cov,covt,var,vart,Mori=(float)MISSVAL;
	int **nvalidkl;
	int FRef;  //=0 if ref pop = all sample, =#categories if ref pop=cat, =-1 if ref pop = given allele freq
	double **Sokl; //for coef of Li 1993
	double **a2kl,**a3kl,**a4kl,**a22kl,**ukl,*Uk,**Wkl,*WWang,**bkl,**ckl,**dkl,**ekl,**fkl,**gkl,**Vkl,**N1kl,**N2kl,**N3kl,**N4kl,**N5kl,**N6kl,**N7kl,**N8kl,*P1,*P2,*P3,phiWang,rWang=(float)MISSVAL,dWang=(float)MISSVAL;
	double a2temp,a3temp,a4temp,a22temp,btemp,ctemp,dtemp,etemp,ftemp,gtemp,Vtemp,N1temp,N2temp,N3temp,N4temp,N5temp,N6temp,N7temp,N8temp;
	double WWangJK,P1JK,P2JK,P3JK,a2JK,a3JK,a4JK,a22JK,bJK,cJK,dJK,eJK,fJK,gJK,VJK,N1JK,N2JK,N3JK,N4JK,N5JK,N6JK,N7JK,N8JK;
	int Swang,Nvalidloci;

	//define if the statistic is asked (0=no), value>0 = value of S
	Elois=Eritl=Erous=Erela=Easc=Erlynch=Edlynch=Erqueller=Erwang=Edwang=Ekinshipdom=Erelatdom=Erli=ENij=0;  
	for(S=1;S<=NS;S++){
		if(Stat[S]==1) Elois=S;		//estimator of kinship coef of Loiselle et al.,1995
		if(Stat[S]==2) Eritl=S;		//estimator of kinship coef of Ritland,1996
		if(Stat[S]==3) Erela=S;		//estimator of relationship coef 
		if(Stat[S]==4) Erous=S;		//estimator of Rousset,2000
		if(Stat[S]==5) Easc=S;		//estimator of allele size correlation
		if(Stat[S]==6) Erlynch=S;	//estimator of relationship coef following Lynch&Ritland 1999
		if(Stat[S]==7) Edlynch=S;	//estimator of fraternity following Lynch&Ritland 1999
		if(Stat[S]==8) Erqueller=S;	//estimator of relationship coef following Lynch&Ritland 1999
		if(Stat[S]==9) Erwang=S;	//estimator of relationship coef following Wang 2002
		if(Stat[S]==10) Edwang=S;	//estimator of fraternity coef following Wang 2002
		if(Stat[S]==11) Ekinshipdom=S;	//estimator of kinship for dominant marker
		if(Stat[S]==12) Erelatdom=S;	//estimator of relationship coef for dominant marker
		if(Stat[S]==13) Erli=S;		//estimator of relationship coef for dominant marker
		if(Stat[S]==14) ENij=S;		//estimator of Nij (kinship equivalent taking into account the phylogenetic distance between alleles)
	}

	if(m==1) linit=1;
	else linit=0;

	//determine # alleles/locus	 (FreqRef=-1 when reference allele frequencies used)
	if(FreqRef==-1) for(l=1;l<=m;l++) Nallelel[l]=Ngivenallelel[l];
	else for(l=1;l<=m;l++) Nallelel[l]=N1allelel[l];

	maxa=0;
	for(l=1;l<=m;l++)if(Nallelel[l]>maxa) maxa=Nallelel[l]; 

	if(FreqRef==1) FRef=Ncat;	   //in case of local ref freq
	else FRef=0;				   //in case of global (or given) ref freq

	//compute allele frequencies
	Pkla=f3tensor(0,Ncat,0,m,0,maxa);
	Nvalgenkl=imatrix(0,Ncat,0,m);
	Nivalidk=ivector(0,Ncat);
	if(FreqRef!=-1){
		for(k=0;k<=FRef;k++) for(l=1;l<=m;l++)for(a=0;a<=Nallelel[l];a++) Pkla[k][l][a]=0.;
		for(k=0;k<=FRef;k++) for(l=1;l<=m;l++) Nvalgenkl[k][l]=0;
		for(l=1;l<=m;l++){
			for(k=0;k<=FRef;k++) Nivalidk[k]=0;
			for(i=1;i<=ntot;i++){
				ploidyi=Ncvalid=0;
				for(c=ploidy-1;c>=0;c--) if(gilc[i][l][c]){
					if(ploidyi==0) ploidyi=c+1;
					Ncvalid++;
				}
				if(Ncvalid) for(c=0;c<ploidyi;c++) if(gilc[i][l][c]){
					Nvalgenkl[0][l]++;
					if(FRef) Nvalgenkl[cati[i]][l]++;
					Pkla[0][l][gilc[i][l][c]]+=1.0f/Ncvalid;
					if(FRef) Pkla[cati[i]][l][gilc[i][l][c]]+=1.0f/Ncvalid;
				}
				if(Ncvalid) Nivalidk[0]++;
				if(FRef && Ncvalid) Nivalidk[cati[i]]++; 
			}
			for(a=1;a<=Nallelel[l];a++) for(k=0;k<=FRef;k++) Pkla[k][l][a]/=Nivalidk[k];
		}
	}
	if(FreqRef==-1) for(l=1;l<=m;l++)for(a=0;a<=Nallelel[l];a++) Pkla[0][l][a]=givenPla[l][a]; 
	
	//compute tot # alleles
	Nallelekl=imatrix(0,Ncat,0,m);
	for(l=1;l<=m;l++)for(k=0;k<=FRef;k++){
		Nallelekl[k][l]=0;
		for(a=1;a<=Nallelel[l];a++) if(Pkla[k][l][a]) Nallelekl[k][l]++;
	}

	/*compute mean and variance of allele lenght*/
	Masizekl=matrix(0,Ncat,0,m);
	Vasizekl=matrix(0,Ncat,0,m);
	for(l=1;l<=m;l++)for(k=0;k<=FRef;k++){
		Sasize=SSasize=0.;
		for(a=1;a<=Nallelel[l];a++){
			Sasize+=Pkla[k][l][a]*allelesizela[l][a];
			SSasize+=Pkla[k][l][a]*allelesizela[l][a]*allelesizela[l][a];
		}
		Masizekl[k][l]=(float)Sasize;
		if(FreqRef!=-1) Vasizekl[k][l]=(float)((SSasize-Sasize*Sasize)*((double)Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.)));
		else Vasizekl[k][l]=(float)(SSasize-Sasize*Sasize);
	}


	
	/*KINSHIP ESTIMATORS (estimator Loiselle, estimator Ritland, allele size correlation)*/
	/*compute inbreeding for each indiv*/
	if(Elois || Eritl || Easc || ENij){
		/*compute denominators for kinship estimator of Loiselle*/
		if(Elois){
			LoisD=matrix(0,FRef,0,m);
			for(k=0;k<=FRef;k++)for(l=0;l<=m;l++) LoisD[k][l]=0.;
			for(k=0;k<=FRef;k++)for(l=1;l<=m;l++)for(a=1;a<=Nallelel[l];a++) LoisD[k][l]+=Pkla[k][l][a]*(1.0f-Pkla[k][l][a]);
		}
		//compute phylodiversity for Nij estimator
		if(ENij){
			DivN=matrix(0,FRef,0,m);
			for(k=0;k<=FRef;k++)for(l=0;l<=m;l++) DivN[k][l]=0.0f;
			for(k=0;k<=FRef;k++)for(l=1;l<=m;l++){
				for(a1=1;a1<=Nallelel[l];a1++)for(a2=1;a2<=Nallelel[l];a2++){
					DivN[k][l]+=Pkla[k][l][a1]*Pkla[k][l][a2]*Distla1a2[l][a1][a2];
				}
				if(FreqRef!=-1 && Nvalgenkl[k][l]>1) DivN[k][l]*=Nvalgenkl[k][l]/(Nvalgenkl[k][l]-1.0f);
			}
		}


		/*compute weight per locus*/
		if(Elois) for(l=1;l<=m;l++) corrSlij[Elois][l][0][0]=LoisD[0][l];
		if(Eritl) for(l=1;l<=m;l++) corrSlij[Eritl][l][0][0]=Nallelekl[0][l]-1.0f;
		if(Easc) for(l=1;l<=m;l++) corrSlij[Easc][l][0][0]=Vasizekl[0][l];
		if(ENij) for(l=1;l<=m;l++) corrSlij[ENij][l][0][0]=DivN[0][l];
		if(ploidy>1){  /*compute inbreeding coef from kinship between genes within indiv*/
			for(i=1;i<=n;i++){
				SoPt=LoisDt=LoisNt=Sasc=SNcomp=SNij=SNijcomp=0.;
				newNloci=newSumNall=0;

				if(FreqRef>0) k=cati[i];
				else k=0;

				for(l=1;l<=m;l++){
					NcompRitl=Ncomp=0;
					LoisN=SoP=ascor=0.;
					Nij=0.;

					for(ci=0;ci<(ploidy-1);ci++) for(cj=(ci+1);cj<ploidy;cj++){ 
						if(gilc[i][l][ci]&&gilc[i][l][cj]){	/*if no missing value*/
							NcompRitl++; Ncomp++; 
							/*make the allele comparisons between i & j*/
							if(gilc[i][l][ci]==gilc[i][l][cj]){ 
								if(Eritl){

									SoP+=1.0f/Pkla[k][l][gilc[i][l][ci]];	  

					/*				if(FreqRef==-1) SoP+=1./Pkla[0][l][gilc[i][l][ci]];
									else{ //compute the corrected (elimination of the compared ind) allele freq (for Ritland's estimator)
										Ngivenallele=Nvalidallele=0.;
										for(g=0;g<ploidy;g++)if(gilc[i][l][g]){
											Nvalidallele++;
											if(gilc[i][l][g]==gilc[i][l][ci]) Ngivenallele++;
										}
										Pnew=(Pkla[k][l][gilc[i][l][ci]]*Nvalgenkl[k][l]-Ngivenallele)/(Nvalgenkl[k][l]-Nvalidallele);
										if((int)(Pkla[k][l][gilc[i][l][ci]]*Nvalgenkl[k][l])>(int)Ngivenallele) SoP+=1./Pnew;
										else NcompRitl--;
									}
							*/	}
								LoisN+=1.;
							}//end of if(gilc[i][l][ci]==gilc[j][l][cj])
							
							if(Elois)LoisN+=((1.0f-LoisD[k][l])-Pkla[k][l][gilc[i][l][ci]]-Pkla[k][l][gilc[i][l][cj]]);
							
							if(Easc && Vasizekl[k][l]){
								asc=(allelesizela[l][gilc[i][l][ci]]-Masizekl[k][l])*(allelesizela[l][gilc[i][l][cj]]-Masizekl[k][l])/Vasizekl[k][l];
								if(FreqRef!=-1) asc+=1.0f/(Nvalgenkl[k][l]-1.0f);
								ascor+=asc;
								Sasc+=asc*Vasizekl[k][l];
							}

							if(ENij){
								Nij+=(DivN[k][l]-Distla1a2[l][gilc[i][l][ci]][gilc[i][l][cj]])/DivN[k][l];
				//				Nij+=Distla1a2[l][gilc[i][l][ci]][gilc[i][l][cj]];
								SNij+=Nij*DivN[k][l];
							}
						}//end of if missing value
					}//end of loops ci and cj
					
					if(Elois && Ncomp && LoisD[k][l]){ 
						LoisN/=Ncomp;
						Lois=(LoisN/LoisD[k][l]);
						if(FreqRef!=-1) Lois+=1.0f/(Nvalgenkl[k][l]-1.0f); //sampling bias correction
						LoisNt+=Lois*LoisD[k][l];
						LoisDt+=LoisD[k][l];
					}
					else Lois=(float)MISSVAL;
					if(Elois) corrSlij[Elois][l][i][i]=Lois;//kinship of Loiselle
			
					if(Eritl && NcompRitl && (Nallelekl[k][l]-1)){ 
						SoP/=NcompRitl;
						Ritl=(SoP-1.0f)/(Nallelekl[k][l]-1.0f);
						SoPt+=SoP;
						newNloci++;
						newSumNall+=Nallelekl[k][l];
					}
					else Ritl=(float)MISSVAL;
 					if(Eritl) corrSlij[Eritl][l][i][i]=Ritl;//kinship of Ritland

					if(Easc && Ncomp && Vasizekl[k][l]) {
						ascor/=Ncomp; 
						SNcomp+=Ncomp*Vasizekl[k][l];
						Ncompasc[l]=(int)(Ncomp*Vasizekl[k][l]);
					}
					else ascor=(float)MISSVAL;
					if(Easc) corrSlij[Easc][l][i][i]=ascor;//allele size correlation

					if(ENij && Ncomp && DivN[k][l]) {
						Nij/=Ncomp; 
						SNijcomp+=Ncomp*DivN[k][l];
						NcompNij[l]=(int)(Ncomp*DivN[k][l]);
					}
					else Nij=(float)MISSVAL;
					if(ENij) corrSlij[ENij][l][i][i]=Nij;//inbreeding based on allele distances
				
				}/*end of loop l*/

				if(m>1){
					if(Elois && LoisDt) Loist=LoisNt/LoisDt;
					else Loist=(float)MISSVAL;
					if(Elois) corrSlij[Elois][0][i][i]=Loist;//kinship of Loiselle

					if(Eritl && newNloci) Ritlt=(SoPt-newNloci)/(newSumNall-newNloci);
					else Ritlt=(float)MISSVAL;
					if(Eritl) corrSlij[Eritl][0][i][i]=Ritlt;//kinship of Ritland

					if(Easc && SNcomp) ascor=Sasc/SNcomp; 
					else ascor=(float)MISSVAL;
 					if(Easc) corrSlij[Easc][0][i][i]=ascor;//value used for allele size similarity (identical alleles included)

					if(ENij && SNijcomp) Nij=SNij/SNijcomp; 
					else Nij=(float)MISSVAL;
 					if(ENij) corrSlij[ENij][0][i][i]=Nij;
				}
				//Compute multilocus estimator for all loci-1 (to compute jackknife estimators in other routines)
				if(JKl && m>1)for(l=1;l<=m;l++){
					if(Elois){
						if(corrSlij[Elois][0][i][i]!=(float)MISSVAL && corrSlij[Elois][l][i][i]!=(float)MISSVAL && (LoisDt-LoisD[k][l])){
							LoisN=corrSlij[Elois][0][i][i]*LoisDt-corrSlij[Elois][l][i][i]*LoisD[k][l];
							Lois=LoisN/(LoisDt-LoisD[k][l]);
						}
						if(corrSlij[Elois][l][i][i]==(float)MISSVAL) Lois=corrSlij[Elois][0][i][i];
						if(LoisDt==LoisD[k][l]) Lois=(float)MISSVAL;
						corrSlij[Elois][-l][i][i]=Lois;
					}

					if(Eritl){
						if(corrSlij[Eritl][0][i][i]!=(float)MISSVAL && corrSlij[Eritl][l][i][i]!=(float)MISSVAL && (newNloci>1)){
							RitlN=corrSlij[Eritl][0][i][i]*(newSumNall-newNloci)-corrSlij[Eritl][l][i][i]*(Nallelekl[k][l]-1.0f);
							RitlD=(newSumNall-newNloci)-(Nallelekl[k][l]-1.0f);
							Ritl=RitlN/RitlD;
						}
						if(corrSlij[Eritl][l][i][i]==(float)MISSVAL) Ritl=corrSlij[Eritl][0][i][i];
						if((newSumNall-newNloci)==(Nallelekl[k][l]-1.0f)) Ritl=(float)MISSVAL;
						corrSlij[Eritl][-l][i][i]=Ritl;
					}

					if(Easc){
						if(corrSlij[Easc][0][i][i]!=(float)MISSVAL && corrSlij[Easc][l][i][i]!=(float)MISSVAL && (SNcomp-Ncompasc[l])){
							ascor=(corrSlij[Easc][0][i][i]*SNcomp-corrSlij[Easc][l][i][i]*Ncompasc[l])/(SNcomp-Ncompasc[l]); 
						}
						if(corrSlij[Easc][l][i][i]==(float)MISSVAL) ascor=corrSlij[Easc][0][i][i];
						if(SNcomp==Ncompasc[l]) ascor=(float)MISSVAL;
 						corrSlij[Easc][-l][i][i]=ascor;
					}

					if(ENij){
						if(corrSlij[ENij][0][i][i]!=(float)MISSVAL && corrSlij[ENij][l][i][i]!=(float)MISSVAL && (SNijcomp-NcompNij[l])){
							Nij=(corrSlij[ENij][0][i][i]*SNijcomp-corrSlij[ENij][l][i][i]*NcompNij[l])/(SNijcomp-NcompNij[l]); 
						}
						if(corrSlij[ENij][l][i][i]==(float)MISSVAL) Nij=corrSlij[ENij][0][i][i];
						if(SNijcomp==NcompNij[l]) Nij=(float)MISSVAL;
 						corrSlij[ENij][-l][i][i]=Nij;
					}
				}//end of JKl
			}//end of loop i
		
		}//end of if(ploidy>1)
		else for(S=1;S<=NS;S++) for(i=1;i<=n;i++) for(l=linit;l<=m;l++) corrSlij[S][l][i][i]=(float)MISSVAL;

		//compute kinship for each pair of different indiv 
		if(compute_inbreeding_coef_only==0) for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];	//define ref allele freq
				else k=0;

				SoPt=LoisDt=LoisNt=Sasc=SNcomp=SNij=SNijcomp=0.;
				newNloci=newSumNall=0;
				for(l=1;l<=m;l++){
					Ncomp=NcompRitl=0;
					LoisN=SoP=ascor=0.;
					Nij=0.;
					
					for(ci=0;ci<ploidy;ci++) for(cj=0;cj<ploidy;cj++){ 
						if(gilc[i][l][ci]&&gilc[j][l][cj]){	//if no missing value
							Ncomp++; NcompRitl++;
							//make the allele comparisons between i & j
							if(gilc[i][l][ci]==gilc[j][l][cj]){ 
								LoisN+=1.0f;
								if(Eritl){

									SoP+=1.0f/Pkla[0][l][gilc[i][l][ci]];		 

			/*						if(FreqRef==-1) SoP+=1./Pkla[0][l][gilc[i][l][ci]];
									else{ //compute the corrected (elimination of the compared ind) allele freq
										Ngivenallele=Nvalidallele=0.;
										for(g=0;g<ploidy;g++){
											if(gilc[i][l][g]) Nvalidallele++;
											if(gilc[i][l][g]==gilc[i][l][ci]) Ngivenallele++;
											if(gilc[j][l][g]) Nvalidallele++;
											if(gilc[j][l][g]==gilc[i][l][cj]) Ngivenallele++;
										}
										Pnew=(Pkla[k][l][gilc[i][l][ci]]*Nvalgenkl[k][l]-Ngivenallele)/(Nvalgenkl[k][l]-Nvalidallele);
										if((int)(Pkla[k][l][gilc[i][l][ci]]*Nvalgenkl[k][l])>(int)Ngivenallele) SoP+=1./Pnew;

										else NcompRitl--;
									}
							*/	}
							}//end of if(gilc[i][l][ci]==gilc[j][l][cj])
							if(Elois)LoisN+=((1.0f-LoisD[k][l])-Pkla[k][l][gilc[i][l][ci]]-Pkla[k][l][gilc[j][l][cj]]);

							if(Easc && Vasizekl[k][l]){
								asc=(allelesizela[l][gilc[i][l][ci]]-Masizekl[k][l])*(allelesizela[l][gilc[j][l][cj]]-Masizekl[k][l])/Vasizekl[k][l];
								if(FreqRef!=-1) asc+=1.0f/(Nvalgenkl[k][l]-1.0f);
								ascor+=asc;
								Sasc+=asc*Vasizekl[k][l];;
							}

							if(ENij){
								Nij+=(DivN[k][l]-Distla1a2[l][gilc[i][l][ci]][gilc[j][l][cj]])/DivN[k][l];
				//				Nij+=Distla1a2[l][gilc[i][l][ci]][gilc[j][l][cj]];
								SNij+=Nij*DivN[k][l];
							}
						}//end of missing value
					}//end of loops ci and cj
					
					if(Elois){
						if(Ncomp && LoisD[k][l]){ 
							LoisN/=Ncomp;
							Lois=(LoisN/LoisD[k][l]);
							if(FreqRef!=-1) Lois+=1.0f/(Nvalgenkl[k][l]-1.0f);
							LoisNt+=Lois*LoisD[k][l];
							LoisDt+=LoisD[k][l];
   							corrSlij[Elois][l][i][j]=Lois;//kinship of Loiselle
   							corrSlij[Elois][l][j][i]=LoisD[k][l];//weight
						}
						else{
   							corrSlij[Elois][l][i][j]=(float)MISSVAL;//kinship of Loiselle
   							corrSlij[Elois][l][j][i]=0.;//weight
						}
					}

					if(Eritl){
						if(NcompRitl && (Nallelekl[k][l]-1)){ 
							SoP/=NcompRitl;
							Ritl=(SoP-1.0f)/(Nallelekl[k][l]-1.0f);
							SoPt+=SoP;
							newNloci++;
							newSumNall+=Nallelekl[k][l];
   							corrSlij[Eritl][l][i][j]=Ritl;//kinship of Ritland
   							corrSlij[Eritl][l][j][i]=(Nallelekl[k][l]-1.0f);//weight
						}
						else{
   							corrSlij[Eritl][l][i][j]=(float)MISSVAL;//kinship of Ritland
   							corrSlij[Eritl][l][j][i]=0.;//weight
						}
					}
					
					if(Easc){
						if(Ncomp && Vasizekl[k][l]) {
							ascor/=Ncomp; 
							SNcomp+=Ncomp*Vasizekl[k][l];
   							corrSlij[Easc][l][i][j]=ascor;
   							corrSlij[Easc][l][j][i]=Ncomp*Vasizekl[k][l];//weight
						}
						else{
   							corrSlij[Easc][l][i][j]=(float)MISSVAL;
   							corrSlij[Easc][l][j][i]=0.;//weight
						}
					}

					if(ENij){
						if(Ncomp && DivN[k][l]) {
							Nij/=Ncomp; 
							SNijcomp+=Ncomp*DivN[k][l];
   							corrSlij[ENij][l][i][j]=Nij;
   							corrSlij[ENij][l][j][i]=Ncomp*DivN[k][l];//weight
						}
						else{
   							corrSlij[ENij][l][i][j]=(float)MISSVAL;
   							corrSlij[ENij][l][j][i]=0.;//weight
						}
					}
				}//end of loop l

				if(m>1){
					if(Elois){
						if(LoisDt){
							Loist=LoisNt/LoisDt;
							corrSlij[Elois][0][i][j]=Loist;//kinship of Loiselle
							corrSlij[Elois][0][j][i]=LoisDt;
						}
						else{
							corrSlij[Elois][0][i][j]=(float)MISSVAL;
							corrSlij[Elois][0][j][i]=0.;
						}
					}
			
					if(Eritl){
						if(newNloci){
							Ritlt=(SoPt-newNloci)/(newSumNall-newNloci);
							corrSlij[Eritl][0][i][j]=Ritlt;//kinship of Ritland
							corrSlij[Eritl][0][j][i]=(float)(newSumNall-newNloci);
						}
						else{
							corrSlij[Eritl][0][i][j]=(float)MISSVAL;
							corrSlij[Eritl][0][j][i]=0.0f;
						}
					}

					if(Easc){
						if(SNcomp){
							ascor=Sasc/SNcomp; 
							corrSlij[Easc][0][i][j]=ascor;//value used for allele size similarity (identical alleles included)
							corrSlij[Easc][0][j][i]=SNcomp;
						}
						else{
							corrSlij[Easc][0][i][j]=(float)MISSVAL;
							corrSlij[Easc][0][j][i]=0.;
						}
					}

					if(ENij){
						if(SNijcomp){
							Nij=SNij/SNijcomp; 
							corrSlij[Easc][0][i][j]=Nij;
							corrSlij[Easc][0][j][i]=SNijcomp;
						}
						else{
							corrSlij[Easc][0][i][j]=(float)MISSVAL;
							corrSlij[Easc][0][j][i]=0.;
						}
					}

				}//end of if(m>1)
				
				//Compute multilocus estimator for all loci-1 (to compute jackknife estimators in other routines)
		/**/	if(JKl && m>1)for(l=1;l<=m;l++){

					for(S=1;S<=NS;S++)if(S==Elois || S==Eritl || S==Easc || S==ENij){
						if(corrSlij[S][0][i][j]!=(float)MISSVAL && corrSlij[S][l][i][j]!=(float)MISSVAL && (corrSlij[S][0][j][i]-corrSlij[S][l][j][i])){
							LoisN=corrSlij[S][0][i][j]*corrSlij[S][0][j][i]-corrSlij[S][l][i][j]*corrSlij[S][l][j][i];
							Lois=LoisN/(corrSlij[S][0][j][i]-corrSlij[S][l][j][i]);
						}
						if(corrSlij[S][l][i][j]==(float)MISSVAL) Lois=corrSlij[S][0][i][j];
						if((corrSlij[S][0][j][i]-corrSlij[S][l][j][i])==0) Lois=(float)MISSVAL;
						corrSlij[S][-l][i][j]=Lois;
					}
	
		/*			
					if(Elois){
						if(corrSlij[Elois][0][i][j]!=(float)MISSVAL && corrSlij[Elois][l][i][j]!=(float)MISSVAL && (LoisDt-LoisD[k][l])){
							LoisN=corrSlij[Elois][0][i][j]*LoisDt-corrSlij[Elois][l][i][j]*LoisD[k][l];
							Lois=LoisN/(LoisDt-LoisD[k][l]);
						}
						if(corrSlij[Elois][l][i][j]==(float)MISSVAL) Lois=corrSlij[Elois][0][i][j];
						if((LoisDt-LoisD[k][l])==0) Lois=(float)MISSVAL;
						corrSlij[Elois][-l][i][j]=Lois;
					}

					if(Eritl){
						if(corrSlij[Eritl][0][i][j]!=(float)MISSVAL && corrSlij[Eritl][l][i][j]!=(float)MISSVAL && (newNloci>1)){
							RitlN=corrSlij[Eritl][0][i][j]*(newSumNall-newNloci)-corrSlij[Eritl][l][i][j]*(Nallelekl[k][l]-1.);
							RitlD=(newSumNall-newNloci)-(Nallelekl[k][l]-1.);
							Ritl=RitlN/RitlD;
						}
						if(corrSlij[Eritl][l][i][j]==(float)MISSVAL) Ritl=corrSlij[Eritl][0][i][j];
						if((newSumNall-newNloci)==(Nallelekl[k][l]-1.)) Ritl=(float)MISSVAL;
						corrSlij[Eritl][-l][i][j]=Ritl;
					}

					if(Easc){
						if(corrSlij[Easc][0][i][j]!=(float)MISSVAL && corrSlij[Easc][l][i][j]!=(float)MISSVAL && (SNcomp-Ncompasc[l])){
							ascor=(corrSlij[Easc][0][i][j]*SNcomp-corrSlij[Easc][l][i][j]*Ncompasc[l])/(SNcomp-Ncompasc[l]); 
						}
						if(corrSlij[Easc][l][i][j]==(float)MISSVAL) ascor=corrSlij[Easc][0][i][j];
						if(SNcomp==Ncompasc[l]) ascor=(float)MISSVAL;
 						corrSlij[Easc][-l][i][j]=ascor;
					}
		*/   	 
				}//end of JKl		   
		  	 
			}//end of loop j
		}//end of loop i
		if(Elois) free_matrix(LoisD,0,FRef,0,m);
		if(ENij) free_matrix(DivN,0,FRef,0,m);
		if(compute_inbreeding_coef_only) return;
	}//end of if(Elois || Eritl || Easc)



	//KINSHIP / RELATIONSHIP ESTIMATORS FOR DOMINANT MARKERS (estimator of Hardy)
	//compute inbreeding for each indiv
	if(Ekinshipdom || Erelatdom){
		//compute heritability of dominant markers
		h2kl=matrix(0,FRef,0,m);
		for(k=0;k<=FRef;k++)for(l=1;l<=m;l++){
			if(Pkla[k][l][1] && (1.0f-Pkla[k][l][1])) h2kl[k][l]=(float)(2.0/(1.0+givenF)*(sqrt(givenF*givenF+4.0*(1.0-Pkla[k][l][2])*(1.0-givenF))+givenF)/(sqrt(givenF*givenF+4.*(1.0-Pkla[k][l][2])*(1.0-givenF))+2.0-givenF));
			else h2kl[k][l]=0.;
		}
		
		//compute denominators for kinship estimator for dominant marker
		RdomD=matrix(0,FRef,0,m);
		for(k=0;k<=FRef;k++)for(l=1;l<=m;l++) RdomD[k][l]=Pkla[k][l][2]*(1.0f-Pkla[k][l][2])*h2kl[k][l];

		//compute weight per locus
		if(Ekinshipdom) for(l=1;l<=m;l++) corrSlij[Ekinshipdom][l][0][0]=RdomD[0][l];
		if(Erelatdom) for(l=1;l<=m;l++) corrSlij[Erelatdom][l][0][0]=RdomD[0][l];

		//compute kinship/relationship for each pair of different indiv
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];	//define ref allele freq
				else k=0;

				RdomDt=RdomNt=0.;
				for(l=1;l<=m;l++){
					if(RdomD[k][l] && gilc[i][l][0] && gilc[j][l][0]){	//if no missing value
						//make the allele comparisons between i & j
						if(gilc[i][l][0]==1 && gilc[j][l][0]==1) RdomN=Pkla[k][l][2]*Pkla[k][l][2];
						else if(gilc[i][l][0]==2 && gilc[j][l][0]==2) RdomN=Pkla[k][l][1]*Pkla[k][l][1];
						else RdomN=-1.0f*Pkla[k][l][1]*Pkla[k][l][2];

						Rdom=(RdomN/RdomD[k][l]);
						if(FreqRef!=-1) Rdom+=1.0f/(Nvalgenkl[k][l]-1.0f);
						RdomNt+=Rdom*RdomD[k][l];
						RdomDt+=RdomD[k][l];
   						if(Erelatdom) corrSlij[Erelatdom][l][i][j]=Rdom;//relationship
   						if(Erelatdom) corrSlij[Erelatdom][l][j][i]=RdomD[k][l];//weight
   						if(Ekinshipdom) corrSlij[Ekinshipdom][l][i][j]=Rdom*(1.0f+givenF)/2.0f;//kinship
   						if(Ekinshipdom) corrSlij[Ekinshipdom][l][j][i]=RdomD[k][l];//weight
					}
					else{
   						if(Erelatdom) corrSlij[Erelatdom][l][i][j]=(float)MISSVAL;//relationship
   						if(Erelatdom) corrSlij[Erelatdom][l][j][i]=0.;//weight
   						if(Ekinshipdom) corrSlij[Ekinshipdom][l][i][j]=(float)MISSVAL;//kinship
   						if(Ekinshipdom) corrSlij[Ekinshipdom][l][j][i]=0.;//weight
					}
				}//end of loop l

				if(m>1){
					if(RdomDt){
						Rdomt=RdomNt/RdomDt;
   						if(Erelatdom) corrSlij[Erelatdom][0][i][j]=Rdomt;//relationship
   						if(Erelatdom) corrSlij[Erelatdom][0][j][i]=RdomDt;//weight
						if(Ekinshipdom) corrSlij[Ekinshipdom][0][i][j]=Rdomt*(1.0f+givenF)/2.0f;//kinship
						if(Ekinshipdom) corrSlij[Ekinshipdom][0][j][i]=RdomDt;
					}
					else{
   						if(Erelatdom) corrSlij[Erelatdom][0][i][j]=(float)MISSVAL;//relationship
   						if(Erelatdom) corrSlij[Erelatdom][0][j][i]=0.;//weight
   						if(Ekinshipdom) corrSlij[Ekinshipdom][0][i][j]=(float)MISSVAL;//kinship
   						if(Ekinshipdom) corrSlij[Ekinshipdom][0][j][i]=0.;//weight
					}

				}//end of if(m>1)
				
				//Compute multilocus estimator for all loci-1 (to compute jackknife estimators in other routines)
		/**/	if(JKl && m>1)for(l=1;l<=m;l++){

					for(S=1;S<=NS;S++)if(S==Ekinshipdom || S==Erelatdom){
						if(corrSlij[S][0][i][j]!=(float)MISSVAL && corrSlij[S][l][i][j]!=(float)MISSVAL && (corrSlij[S][0][j][i]-corrSlij[S][l][j][i])){
							RdomN=corrSlij[S][0][i][j]*corrSlij[S][0][j][i]-corrSlij[S][l][i][j]*corrSlij[S][l][j][i];
							Rdom=RdomN/(corrSlij[S][0][j][i]-corrSlij[S][l][j][i]);
						}
						if(corrSlij[S][l][i][j]==(float)MISSVAL) Rdom=corrSlij[S][0][i][j];
						if((corrSlij[S][0][j][i]-corrSlij[S][l][j][i])==0) Rdom=(float)MISSVAL;
						corrSlij[S][-l][i][j]=Rdom;
					}
	
				}//end of JKl		   
		  	 
			}//end of loop j
		}//end of loop i
		free_matrix(h2kl,0,FRef,0,m);
		free_matrix(RdomD,0,FRef,0,m);	  	
	}//end of if(Ekinshipdom)


	//ESTIMATOR OF ROUSSET
	if(Erous && ploidy>1){
		//compute observed (=ref) heterozygosity within ind
		Hokl=matrix(0,FRef,0,m);
		Nvalidwpairkl=imatrix(0,FRef,0,m);
		for(k=0;k<=FRef;k++) for(l=0;l<=m;l++) {
			Hokl[k][l]=0.0f;
			Nvalidwpairkl[k][l]=0;
		}

		for(i=1;i<=n;i++){
			newNloci=0;
			if(FreqRef>0) k=cati[i];
			else k=0;

			for(l=1;l<=m;l++){
				Ncomp=0;
				Hw=0.;
				for(ci=0;ci<(ploidy-1);ci++) for(cj=(ci+1);cj<ploidy;cj++){ 
					if(gilc[i][l][ci]&&gilc[i][l][cj]){	//if no missing value
						Ncomp++; 
						if(gilc[i][l][ci]!=gilc[i][l][cj]) Hw+=1.;
					}//end of if missing value
				}//end of loops ci and cj
				if(Ncomp && (Nallelekl[k][l]-1)) newNloci++; 
				if(Erous && Ncomp && (Nallelekl[k][l]-1)){ 
					Hw/=Ncomp;
					Hokl[k][l]+=Hw;
					Nvalidwpairkl[k][l]++;
				}
				else Hw=(float)MISSVAL;
				corrSlij[Erous][l][i][i]=Hw;
			}//end of loop l
		}//end of loop i
		for(k=0;k<=FRef;k++) for(l=1;l<=m;l++){	//mean obs heterozygosity per locus and ref pop
			if(Nvalidwpairkl[k][l]) Hokl[k][l]/=Nvalidwpairkl[k][l];
			else Hokl[k][l]=(float)MISSVAL;
		}
		//weight per locus
		for(l=1;l<=m;l++) corrSlij[Erous][l][0][0]=Hokl[0][l];


		//compute Rousset's estimator for each pair of different indiv 
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];	//define ref pop
				else k=0;

				Hbt=Hwt=Hot=SNcomp=0.;
				newNloci=0;
				for(l=1;l<=m;l++){
					Ncomp=0;
					Hb=0.;
					
					for(ci=0;ci<ploidy;ci++) for(cj=0;cj<ploidy;cj++){ 
						if(gilc[i][l][ci]&&gilc[j][l][cj]){	//if no missing value
							Ncomp++;
							if(gilc[i][l][ci]!=gilc[j][l][cj]) Hb++;
						}//end of missing value
					}//end of loops ci and cj
					
					if(Ncomp && (Nallelekl[k][l]-1)) newNloci++;

					if(Ncomp && Hokl[k][l]!=(float)MISSVAL && corrSlij[Erous][l][i][i]!=(float)MISSVAL && corrSlij[Erous][l][j][j]!=(float)MISSVAL){ 
						Hb/=Ncomp;
						Hw=(corrSlij[Erous][l][i][i]+corrSlij[Erous][l][j][j])/2.0f;
					/*	Rous=(Hb-Hw)/Hokl[k][l];	*/				// estimator "a*" of Rousset
					/*	Rous=(Hb-Hokl[k][l])/Hokl[k][l];*/			//estimator of Hardy
						Rous=(Hb-0.5f*Hw-0.5f*Hokl[k][l])/Hokl[k][l]; //estimator "a" of Rousset
						Hwt+=Hw;
						Hbt+=Hb;
						Hot+=Hokl[k][l];

						corrSlij[Erous][l][i][j]=Rous;//Rousset's estimator
						corrSlij[Erous][l][j][i]=Hokl[k][l]; //weight
					}
					else{
						corrSlij[Erous][l][i][j]=(float)MISSVAL;
						corrSlij[Erous][l][j][i]=0.;
					}
				}//end of loop l

				if(m>1){
					if(newNloci && Hot){
					/*	Rous=(Hbt-Hwt)/Hot; */ 	   //estimator "a*" of Rousset
					/*	Rous=(Hbt-Hot)/Hot;	*/	   //estimator of Hardy
						Rous=(Hbt-0.5f*Hwt-0.5f*Hot)/Hot; //estimator "a" of Rousset

						corrSlij[Erous][0][i][j]=Rous;//Rousset's estimator
						corrSlij[Erous][0][j][i]=Hot; //weight
					}
					else{
						corrSlij[Erous][0][i][j]=(float)MISSVAL;
						corrSlij[Erous][0][j][i]=0.;
					}
				}//end of if(m>1)

				//Compute multilocus estimator for all loci-1 (to compute jackknife estimators in other routines)
		/**/	if(JKl && m>1)for(l=1;l<=m;l++){
					if(corrSlij[Erous][0][i][j]!=(float)MISSVAL && corrSlij[Erous][l][i][j]!=(float)MISSVAL && (Hot-Hokl[k][l])){
						Rous=(corrSlij[Erous][0][i][j]*Hot-corrSlij[Erous][l][i][j]*Hokl[k][l])/(Hot-Hokl[k][l]);
					}
					if(corrSlij[Erous][l][i][j]==(float)MISSVAL) Rous=corrSlij[Erous][0][i][j];
					if(Hot==Hokl[k][l]) Rous=(float)MISSVAL;
					corrSlij[Erous][-l][i][j]=Rous;
				}
		
			}//end of loop j
		}//end of loop i

		free_matrix(Hokl,0,FRef,0,m);
		free_imatrix(Nvalidwpairkl,0,FRef,0,m);

	}//end of ESTIMATOR OF ROUSSET
	

	
	//COEF RELATIONSHIP (MORAN'S I on indiv allele freq)
	if(Erela){   
		//compute allele freq for each ind and variances in allele freq
		for(l=1;l<=m;l++) Plai[l]=matrix(0,Nallelel[l],0,n);
		for(l=1;l<=m;l++) Plak[l]=matrix(0,Nallelel[l],0,FRef);
		varkla=f3tensor(0,FRef,0,m,0,maxa);
		Spik=vector(0,FRef);
		SSpik=vector(0,FRef);
		nvalidkl=imatrix(0,FRef,0,m);
		for(l=1;l<=m;l++) for(a=1;a<=Nallelel[l];a++){
			for(k=0;k<=FRef;k++) {
				Spik[k]=SSpik[k]=0.0f;
				nvalidkl[k][l]=0;
			}
			k=0;
			for(i=1;i<=n;i++){
				if(FRef) k=cati[i];
				pi=qi=0.;
				for(g=0;g<ploidy;g++){
					if(gilc[i][l][g]==a) pi++;
					else if(gilc[i][l][g]) qi++;
				}
				if(pi+qi){
					pi/=(pi+qi);
					Plai[l][a][i]=pi;
					Spik[k]+=pi;
					SSpik[k]+=pi*pi;
					nvalidkl[k][l]++;
				}
				else Plai[l][a][i]=(float)MISSVAL;
			}
			for(k=0;k<=FRef;k++){
				if(nvalidkl[k][l]) varkla[k][l][a]=SSpik[k]/nvalidkl[k][l]-Spik[k]*Spik[k]/(nvalidkl[k][l]*nvalidkl[k][l]);
				else varkla[k][l][a]=(float)MISSVAL;
				Plak[l][a][k]=Spik[k]/nvalidkl[k][l];
			}
		}
		for(k=0;k<=FRef;k++)for(l=1;l<=m;l++) {
			varkla[k][l][0]=0.;
			for(a=1;a<=Nallelel[l];a++) varkla[k][l][0]+=varkla[k][l][a];
		}
		//weight per locus
		for(l=1;l<=m;l++) corrSlij[Erela][l][0][0]=varkla[0][l][0];

		//compute correl matrix for each pair of ind
		for(i=1;i<n;i++){
			for(j=i+1;j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef>0 && cati[i]==cati[j]) k=cati[i];
				else k=0;

				covt=vart=0.;
				for(l=1;l<=m;l++){
					if(Plai[l][1][i]!=(float)MISSVAL && Plai[l][1][j]!=(float)MISSVAL){
						cov=var=0.;
						for(a=1;a<=Nallelel[l];a++){
					//		cov+=(Plai[l][a][i]-Pkla[k][l][a])*(Plai[l][a][j]-Pkla[k][l][a]);
							cov+=(Plai[l][a][i]-Plak[l][a][k])*(Plai[l][a][j]-Plak[l][a][k]);
							if((nvalidkl[k][l]-1.)>0) cov+=varkla[k][l][a]/(nvalidkl[k][l]-1.0f); //bias correction
							var+=varkla[k][l][a];	 
						}
						covt+=cov;
						vart+=var;

						corrSlij[Erela][l][i][j]=cov/var;
						corrSlij[Erela][l][j][i]=var;
					} 
					else{
						corrSlij[Erela][l][i][j]=(float)MISSVAL;
						corrSlij[Erela][l][j][i]=0.;
					}
				}
				if(m>1){
					if(vart){
						corrSlij[Erela][0][i][j]=covt/vart;
						corrSlij[Erela][0][j][i]=vart;
					}
					else{
						corrSlij[Erela][0][i][j]=(float)MISSVAL;
						corrSlij[Erela][0][j][i]=0.;
					}
				}
				/*Compute multilocus estimator for all loci-1 (to compute jackknife estimators in other routines)*/
		/**/	if(JKl && m>1)for(l=1;l<=m;l++){
					if(corrSlij[Erela][0][i][j]!=(float)MISSVAL && corrSlij[Erela][l][i][j]!=(float)MISSVAL && (vart-varkla[k][l][0])){
						Mori=(corrSlij[Erela][0][i][j]*vart-corrSlij[Erela][l][i][j]*varkla[k][l][0])/(vart-varkla[k][l][0]);
					}
					if(corrSlij[Erela][l][i][j]==(float)MISSVAL) Mori=corrSlij[Erela][0][i][j]; 
					if(vart==varkla[k][l][0]) Mori=(float)MISSVAL;
					corrSlij[Erela][-l][i][j]=Mori;
				}
		
			}/*end of loop j*/
		}/*end of loop i*/
		for(l=1;l<=m;l++) free_matrix(Plai[l],0,Nallelel[l],0,n);
		for(l=1;l<=m;l++) free_matrix(Plak[l],0,Nallelel[l],0,FRef);
		free_f3tensor(varkla,0,FRef,0,m,0,maxa);
		free_vector(Spik,0,FRef);
		free_vector(SSpik,0,FRef);
		free_imatrix(nvalidkl,0,FRef,0,m);

	}/*COEF RELATIONSHIP (MORAN'S I on indiv allele freq)*/





	/*COEF RELATIONSHIP AND FRATERNITY (Lynch & Ritland 1999)*/
	if((Erlynch || Edlynch) && ploidy==2){
		SrLynch1l=vector(0,m);
		SWrLynch1l=vector(0,m);
		SrLynch2l=vector(0,m);
		SWrLynch2l=vector(0,m);
		SdLynch1l=vector(0,m);
		SWdLynch1l=vector(0,m);
		SdLynch2l=vector(0,m);
		SWdLynch2l=vector(0,m);

		/*weight per locus*/
		if(Erlynch) for(l=1;l<=m;l++) corrSlij[Erlynch][l][0][0]=(float)MISSVAL;
		if(Edlynch) for(l=1;l<=m;l++) corrSlij[Edlynch][l][0][0]=(float)MISSVAL;

		/*compute kinship for each pair of different indiv */
		k=0;
		
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];
				
				SrLynch1=SrLynch2=SWrLynch1=SWrLynch2=SdLynch1=SdLynch2=SWdLynch1=SWdLynch2=0.;
				for(l=1;l<=m;l++){
					Sac=Sca=Sad=Sda=Sbc=Scb=Sbd=Sdb=Sab=Scd=0.;
					if(gilc[i][l][0]&&gilc[i][l][1]&&gilc[j][l][0]&&gilc[j][l][1]){	/*if no missing value*/
						/*make the allele comparisons between i & j*/
						if(gilc[i][l][0]==gilc[j][l][0]) Sac=Sca=1.;
						if(gilc[i][l][0]==gilc[j][l][1]) Sad=Sda=1.;
						if(gilc[i][l][1]==gilc[j][l][0]) Sbc=Scb=1.;
						if(gilc[i][l][1]==gilc[j][l][1]) Sbd=Sdb=1.;
						if(gilc[i][l][0]==gilc[i][l][1]) Sab=1.;
						if(gilc[j][l][0]==gilc[j][l][1]) Scd=1.;
						
						Pa=Pkla[k][l][gilc[i][l][0]];
						Pb=Pkla[k][l][gilc[i][l][1]];
						Pc=Pkla[k][l][gilc[j][l][0]];
						Pd=Pkla[k][l][gilc[j][l][1]];

						if(Erlynch){
							WrLynch1=((1+Sab)*(Pa+Pb)-4*Pa*Pb)/(2*Pa*Pb);
							if(WrLynch1!=0.) rLynch1=(Pa*(Sbc+Sbd)+Pb*(Sac+Sad)-4*Pa*Pb)/((1+Sab)*(Pa+Pb)-4*Pa*Pb);
							else rLynch1=(float)MISSVAL;
							WrLynch2=((1+Scd)*(Pc+Pd)-4*Pc*Pd)/(2*Pc*Pd);
							if(WrLynch2!=0.) rLynch2=(Pc*(Sda+Sdb)+Pd*(Sca+Scb)-4*Pc*Pd)/((1+Scd)*(Pc+Pd)-4*Pc*Pd);
							else rLynch2=(float)MISSVAL;
							WrLynch=WrLynch1+WrLynch2;
							rLynch=(rLynch1+rLynch2)/2.0f;  							 
						/*	rLynch=(rLynch1*WrLynch1 + rLynch2*WrLynch2)/WrLynch;	 */
							if(rLynch1==(float)MISSVAL && rLynch2!=(float)MISSVAL) rLynch=rLynch2;
							if(rLynch2==(float)MISSVAL && rLynch1!=(float)MISSVAL) rLynch=rLynch1;

							SrLynch1+=(Pa*(Sbc+Sbd)+Pb*(Sac+Sad)-4*Pa*Pb)/(2*Pa*Pb);
							SWrLynch1+=WrLynch1;
							SrLynch2+=(Pc*(Sda+Sdb)+Pd*(Sca+Scb)-4*Pc*Pd)/(2*Pc*Pd);
							SWrLynch2+=WrLynch2;
			/**/			if(JKl){
								SrLynch1l[l]=(Pa*(Sbc+Sbd)+Pb*(Sac+Sad)-4*Pa*Pb)/(2*Pa*Pb);
								SWrLynch1l[l]=((1+Sab)*(Pa+Pb)-4*Pa*Pb)/(2*Pa*Pb);
								SrLynch2l[l]=(Pc*(Sda+Sdb)+Pd*(Sca+Scb)-4*Pc*Pd)/(2*Pc*Pd);
								SWrLynch2l[l]=((1+Scd)*(Pc+Pd)-4*Pc*Pd)/(2*Pc*Pd);
							}
			  
							corrSlij[Erlynch][l][i][j]=rLynch;
							corrSlij[Erlynch][l][j][i]=WrLynch;
						}
						if(Edlynch){
							WdLynch1=((1+Sab)*(1-Pa-Pb)+2*Pa*Pb)/(2*Pa*Pb);
							if(WdLynch1!=0.) dLynch1=(2.0f*Pa*Pb-Pa*(Sbc+Sbd)-Pb*(Sac+Sad)+Sac*Sbd+Sad*Sbc)/((1+Sab)*(1-Pa-Pb)+2*Pa*Pb);
							else dLynch1=(float)MISSVAL;
							WdLynch2=((1+Scd)*(1-Pc-Pd)+2*Pc*Pd)/(2*Pc*Pd);
							if(WdLynch2!=0.) dLynch2=(2.0f*Pc*Pd-Pc*(Sda+Sdb)-Pd*(Sca+Scb)+Sca*Sdb+Scb*Sda)/((1+Scd)*(1-Pc-Pd)+2*Pc*Pd);
							else dLynch2=(float)MISSVAL;
							WdLynch=WdLynch1+WdLynch2;
							dLynch=(dLynch1+dLynch2)/2.0f;   								 
						/*	dLynch=(dLynch1*WdLynch1 + dLynch2*WdLynch2)/WdLynch;		 */
							if(dLynch1==(float)MISSVAL && dLynch2!=(float)MISSVAL) dLynch=dLynch2;
							if(dLynch2==(float)MISSVAL && dLynch1!=(float)MISSVAL) dLynch=dLynch1;
						
							SdLynch1+=(2.0f*Pa*Pb-Pa*(Sbc+Sbd)-Pb*(Sac+Sad)+Sac*Sbd+Sad*Sbc)/(2*Pa*Pb);
							SWdLynch1+=WdLynch1;
							SdLynch2+=(2.0f*Pc*Pd-Pc*(Sda+Sdb)-Pd*(Sca+Scb)+Sca*Sdb+Scb*Sda)/(2*Pc*Pd);
							SWdLynch2+=WdLynch2;
				/**/		if(JKl){
								SdLynch1l[l]=(2.0f*Pa*Pb-Pa*(Sbc+Sbd)-Pb*(Sac+Sad)+Sac*Sbd+Sad*Sbc)/(2*Pa*Pb);
								SWdLynch1l[l]=((1+Sab)*(1-Pa-Pb)+2*Pa*Pb)/(2*Pa*Pb);
								SdLynch2l[l]=(2.0f*Pc*Pd-Pc*(Sda+Sdb)-Pd*(Sca+Scb)+Sca*Sdb+Scb*Sda)/(2*Pc*Pd);
								SWdLynch2l[l]=((1+Scd)*(1-Pc-Pd)+2*Pc*Pd)/(2*Pc*Pd);
							}
				 
							corrSlij[Edlynch][l][i][j]=dLynch;
							corrSlij[Edlynch][l][j][i]=WdLynch;
						}


					}/*end of missing value*/
					else{
						if(Erlynch) { corrSlij[Erlynch][l][i][j]=(float)MISSVAL;	corrSlij[Erlynch][l][j][i]=0.;}
						if(Edlynch) { corrSlij[Edlynch][l][i][j]=(float)MISSVAL;	corrSlij[Edlynch][l][j][i]=0.;}
					}
				}/*end of loop l*/

				if(m>1){
					if(Erlynch){
						SWrLynch=SWrLynch1 + SWrLynch2;
						if(/*SWrLynch1 + SWrLynch2*/SWrLynch1 && SWrLynch2){
							rLynch=(SrLynch1/SWrLynch1 + SrLynch2/SWrLynch2)/2.0f;	 
					/*		rLynch=(SrLynch1 + SrLynch2)/(SWrLynch1 + SWrLynch2);	 */
						}
						if(SWrLynch1 && !SWrLynch2) rLynch=SrLynch1/SWrLynch1;
						if(!SWrLynch1 && SWrLynch2) rLynch=SrLynch2/SWrLynch2;
						if(!SWrLynch1 && !SWrLynch2) rLynch=(float)MISSVAL;
						corrSlij[Erlynch][0][i][j]=rLynch;
						corrSlij[Erlynch][0][j][i]=SWrLynch;
					}

					if(Edlynch){
						SWdLynch=SWdLynch1 + SWdLynch2;
						if(/*SWdLynch1+SWdLynch2*/ SWdLynch1 && SWdLynch2){
							dLynch=(SdLynch1/SWdLynch1 + SdLynch2/SWdLynch2)/2.0f;	  
					/*		dLynch=(SdLynch1 + SdLynch2)/(SWdLynch1 + SWdLynch2);	  */
						}
						if(SWdLynch1 && !SWdLynch2) dLynch=SdLynch1/SWdLynch1;
						if(!SWdLynch1 && SWdLynch2) dLynch=SdLynch2/SWdLynch2;
						if(!SWdLynch1 && !SWdLynch2) dLynch=(float)MISSVAL;
						corrSlij[Edlynch][0][i][j]=dLynch;
						corrSlij[Edlynch][0][j][i]=SWdLynch;
					}
				}/*end of if(m>1)*/

		/**/	if(JKl && m>1)for(l=1;l<=m;l++){
					if((SWrLynch1-SWrLynch1l[l]) && (SWrLynch2-SWrLynch2l[l])) rLynch=((SrLynch1-SrLynch1l[l])/(SWrLynch1-SWrLynch1l[l]) + (SrLynch2-SrLynch2l[l])/(SWrLynch2-SWrLynch2l[l]))/2.0f;
					if(Erlynch)if(corrSlij[Erlynch][l][i][j]==(float)MISSVAL) rLynch=corrSlij[Erlynch][0][i][j];
					if((SWrLynch1==SWrLynch1l[l]) || (SWrLynch2==SWrLynch2l[l])) rLynch=(float)MISSVAL;

					if((SWdLynch1-SWdLynch1l[l]) && (SWdLynch2-SWdLynch2l[l])) dLynch=((SdLynch1-SdLynch1l[l])/(SWdLynch1-SWdLynch1l[l]) + (SdLynch2-SdLynch2l[l])/(SWdLynch2-SWdLynch2l[l]))/2.0f;
					if(Edlynch)if(corrSlij[Edlynch][l][i][j]==(float)MISSVAL) dLynch=corrSlij[Edlynch][0][i][j];
					if((SWdLynch1==SWdLynch1l[l]) || (SWdLynch2==SWdLynch2l[l])) dLynch=(float)MISSVAL;

					if(Erlynch) corrSlij[Erlynch][-l][i][j]=rLynch;
					if(Edlynch) corrSlij[Edlynch][-l][i][j]=dLynch;
				}
		 
			}/*end of loop j*/
		}/*end of loop i*/
		free_vector(SrLynch1l,0,m);
		free_vector(SWrLynch1l,0,m);
		free_vector(SrLynch2l,0,m);
		free_vector(SWrLynch2l,0,m);
		free_vector(SdLynch1l,0,m);
		free_vector(SWdLynch1l,0,m);
		free_vector(SdLynch2l,0,m);
		free_vector(SWdLynch2l,0,m);

	}/*end of COEF RELATIONSHIP AND DOMINANCE (Lynch & Ritland 1999)*/


	/*COEF RELATIONSHIP (Queller & Goodnight 1989)*/
	if(Erqueller && ploidy==2){
		SrQueller1l=vector(0,m);
		SWrQueller1l=vector(0,m);
		SrQueller2l=vector(0,m);
		SWrQueller2l=vector(0,m);

		/*weight per locus*/
		if(Erqueller) for(l=1;l<=m;l++) corrSlij[Erqueller][l][0][0]=(float)MISSVAL;

		/*compute coef for each pair of different indiv */
		k=0;
		
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];
				
				SrQueller1=SrQueller2=SWrQueller1=SWrQueller2=0.;
				for(l=1;l<=m;l++){
					Sac=Sca=Sad=Sda=Sbc=Scb=Sbd=Sdb=Sab=Scd=0.;
					if(gilc[i][l][0]&&gilc[i][l][1]&&gilc[j][l][0]&&gilc[j][l][1]){	/*if no missing value*/
						/*make the allele comparisons between i & j*/
						if(gilc[i][l][0]==gilc[j][l][0]) Sac=Sca=1.;
						if(gilc[i][l][0]==gilc[j][l][1]) Sad=Sda=1.;
						if(gilc[i][l][1]==gilc[j][l][0]) Sbc=Scb=1.;
						if(gilc[i][l][1]==gilc[j][l][1]) Sbd=Sdb=1.;
						if(gilc[i][l][0]==gilc[i][l][1]) Sab=1.;
						if(gilc[j][l][0]==gilc[j][l][1]) Scd=1.;
						
						Pa=Pkla[k][l][gilc[i][l][0]];
						Pb=Pkla[k][l][gilc[i][l][1]];
						Pc=Pkla[k][l][gilc[j][l][0]];
						Pd=Pkla[k][l][gilc[j][l][1]];

						WrQueller1=1+Sab-Pa-Pb;
						if(Nallelekl[k][l]==2 && Sab==0.) WrQueller1=0.;
						if(WrQueller1!=0.) rQueller1=(0.5f*(Sbc+Sbd+Sac+Sad)-Pa-Pb)/WrQueller1;
						else rQueller1=(float)MISSVAL;

						WrQueller2=1+Scd-Pc-Pd;
						if(Nallelekl[k][l]==2 && Scd==0.) WrQueller2=0.;
						if(WrQueller2!=0.) rQueller2=(0.5f*(Sda+Sdb+Sca+Scb)-Pc-Pd)/WrQueller2;
						else rQueller2=(float)MISSVAL;

						WrQueller=WrQueller1+WrQueller2;
						rQueller=(rQueller1+rQueller2)/2.0f;								  
					/*	rQueller=(rQueller1*WrQueller1 + rQueller2*WrQueller2)/WrQueller; */

						if(rQueller1==(float)MISSVAL && rQueller2!=(float)MISSVAL) rQueller=rQueller2;
						if(rQueller2==(float)MISSVAL && rQueller1!=(float)MISSVAL) rQueller=rQueller1;

						SrQueller1+=(0.5f*(Sbc+Sbd+Sac+Sad)-Pa-Pb);
						SWrQueller1+=WrQueller1;
						SrQueller2+=(0.5f*(Sda+Sdb+Sca+Scb)-Pc-Pd);
						SWrQueller2+=WrQueller2;
			/**/		if(JKl){
							SrQueller1l[l]=(0.5f*(Sbc+Sbd+Sac+Sad)-Pa-Pb);
							SWrQueller1l[l]=(1+Sab-Pa-Pb);
							SrQueller2l[l]=(0.5f*(Sda+Sdb+Sca+Scb)-Pc-Pd);
							SWrQueller2l[l]=(1+Scd-Pc-Pd);
						}
			
						corrSlij[Erqueller][l][i][j]=rQueller;
						corrSlij[Erqueller][l][j][i]=WrQueller;
					}/*end of missing value*/
					else{
						corrSlij[Erqueller][l][i][j]=(float)MISSVAL;
						corrSlij[Erqueller][l][j][i]=0.;
					}/*end of missing value*/
				}/*end of loop l*/

				if(m>1){
					WrQueller=SWrQueller1 + SWrQueller2;

					if(/*SWrQueller1 + SWrQueller2*/ SWrQueller1 && SWrQueller2){
						rQueller=(SrQueller1/SWrQueller1 + SrQueller2/SWrQueller2)/2.0f;	 
				/*		rQueller=(SrQueller1 + SrQueller2)/(SWrQueller1 + SWrQueller2);	 */
					}
					if(SWrQueller1 && !SWrQueller2) rQueller=SrQueller1/SWrQueller1;	 
					if(!SWrQueller1 && SWrQueller2) rQueller=SrQueller2/SWrQueller2;	 
					if(!SWrQueller1 && !SWrQueller2) rQueller=(float)MISSVAL;	 

					corrSlij[Erqueller][0][i][j]=rQueller;
					corrSlij[Erqueller][0][j][i]=WrQueller;
				}/*end of if(m>1)*/

				if(JKl && m>1)for(l=1;l<=m;l++){
					if((SWrQueller1-SWrQueller1l[l]) && (SWrQueller2-SWrQueller2l[l])) rQueller=((SrQueller1-SrQueller1l[l])/(SWrQueller1-SWrQueller1l[l]) + (SrQueller2-SrQueller2l[l])/(SWrQueller2-SWrQueller2l[l]))/2.0f;
					if(corrSlij[Erqueller][l][i][j]==(float)MISSVAL) rQueller=corrSlij[Erqueller][0][i][j];
					if((SWrQueller1==SWrQueller1l[l]) || (SWrQueller2==SWrQueller2l[l])) rQueller=(float)MISSVAL;

					corrSlij[Erqueller][-l][i][j]=rQueller;
				}
	/*	*/
			}/*end of loop j*/
		}/*end of loop i*/
		free_vector(SrQueller1l,0,m);
		free_vector(SWrQueller1l,0,m);
		free_vector(SrQueller2l,0,m);
		free_vector(SWrQueller2l,0,m);

	}/*end of COEF RELATIONSHIP (Queller & Goodnight 1989)*/




	//COEF RELATIONSHIP (Li et al. 1993 + sampling bias correction of Wang(2002))
	if(Erli && ploidy==2){

		a2kl=dmatrix(0,FRef,0,m);
		a3kl=dmatrix(0,FRef,0,m);
		Sokl=dmatrix(0,FRef,0,m);

		//compute coef for each pair of different indiv 
		for(l=1;l<=m;l++) Srli[l]=SSrli[l]=Nrli[l]=0.;
		
		for(k=0;k<=FRef;k++) for(l=1;l<=m;l++){
			a2kl[k][l]=a3kl[k][l]=Sokl[k][l]=0.;
			for(a=1;a<=Nallelel[l];a++) {
				a2kl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a];
				a3kl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a]*Pkla[k][l][a];
			}
			//sample size corrected a values (following Wang (2002)
			if(FreqRef!=-1){
				a2kl[k][l]=(Nvalgenkl[k][l]*a2kl[k][l]-1.)/(Nvalgenkl[k][l]-1.);
				a3kl[k][l]=(Nvalgenkl[k][l]*Nvalgenkl[k][l]*a3kl[k][l]-3.*(Nvalgenkl[k][l]-1.)*a2kl[k][l]-1.)/((Nvalgenkl[k][l]-1.)*(Nvalgenkl[k][l]-2.));
			}
			Sokl[k][l]= 2.*a2kl[k][l]-a3kl[k][l];
		}

	//	for(l=1;l<=m;l++) for(a=1;a<=Nallelel[l];a++) So[l]+= Pkla[0][l][a]*Pkla[0][l][a]*(2.-Pkla[0][l][a]);
		
		k=0;
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];

	//			if(k) for(l=1;l<=m;l++) for(a=1;a<=Nallelel[l];a++) So[l]+= Pkla[k][l][a]*Pkla[k][l][a]*(2.-Pkla[k][l][a]);
				
				for(l=1;l<=m;l++){
					if(gilc[i][l][0]&&gilc[i][l][1]&&gilc[j][l][0]&&gilc[j][l][1]){	/*if no missing value*/
						/*make the allele comparisons between i & j*/
						if((gilc[i][l][0]==gilc[j][l][0] && gilc[i][l][1]==gilc[j][l][1])||(gilc[i][l][0]==gilc[j][l][1] && gilc[i][l][1]==gilc[j][l][0])) Sij=1.;
						else if((gilc[i][l][0]==gilc[i][l][1] && (gilc[i][l][0]==gilc[j][l][0] || gilc[i][l][0]==gilc[j][l][1])) || (gilc[j][l][0]==gilc[j][l][1] && (gilc[j][l][0]==gilc[i][l][0] || gilc[j][l][0]==gilc[i][l][1]))) Sij=0.75;
						else if(gilc[i][l][0]==gilc[j][l][0] || gilc[i][l][0]==gilc[j][l][1] || gilc[i][l][1]==gilc[j][l][0] || gilc[i][l][1]==gilc[j][l][1]) Sij=0.5;
						else Sij=0.;

						rli=(float)((Sij-Sokl[k][l])/(1.0-Sokl[k][l]));
			
						corrSlij[Erli][l][i][j]=rli;

						Nrli[l]++;
						Srli[l]+=rli;
						SSrli[l]+=rli*rli;

					}/*end of missing value*/
					else{
						corrSlij[Erli][l][i][j]=(float)MISSVAL;
						corrSlij[Erli][l][j][i]=0.;
					}/*end of missing value*/
				}/*end of loop l*/

			}/*end of loop j*/
		}/*end of loop i*/

		if(m>1){
			//empirically derived weights
			for(l=1;l<=m;l++) {
				if(Nrli[l]>1){
					Var=(SSrli[l] - Srli[l]*Srli[l]/Nrli[l])/(Nrli[l]-1.0f);
					if(Var) corrSlij[Erli][l][0][0]=1.0f/Var;
					else corrSlij[Erli][l][0][0]=0.;
				}
				else corrSlij[Erli][l][0][0]=0.;
			}
			
			for(i=1;i<n;i++){
				for(j=(i+1);j<=n;j++){
					if(TypeComp==1)if(cati[i]!=cati[j]) continue;
					if(TypeComp==2)if(cati[i]==cati[j]) continue;
					if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];

					//multilocus estimates
					MLrli=SumW=0.;
					for(l=1;l<=m;l++)if(corrSlij[Erli][l][i][j]!=(float)MISSVAL){
						corrSlij[Erli][l][j][i]=corrSlij[Erli][l][0][0];
						MLrli+=corrSlij[Erli][l][0][0]*corrSlij[Erli][l][i][j];
						SumW+=corrSlij[Erli][l][0][0];
					}
					corrSlij[Erli][0][j][i]=SumW;
					if(SumW) corrSlij[Erli][0][i][j]=MLrli/SumW;
					else corrSlij[Erli][0][i][j]=(float)MISSVAL;

					//jackknife
					if(JKl)for(l=1;l<=m;l++){
						if(SumW==corrSlij[Erli][l][0][0]) corrSlij[Erli][-l][i][j]=(float)MISSVAL;
						else if(corrSlij[Erli][l][i][j]!=(float)MISSVAL) corrSlij[Erli][-l][i][j]=(MLrli-corrSlij[Erli][l][i][j]*corrSlij[Erli][l][0][0])/(SumW-corrSlij[Erli][l][0][0]);
						else corrSlij[Erli][-l][i][j]=corrSlij[Erli][0][i][j]; 
					}
				}//loop j
			}//loop i
		}//end of if(m>1)

		free_dmatrix(a2kl,0,FRef,0,m);
		free_dmatrix(a3kl,0,FRef,0,m);
		free_dmatrix(Sokl,0,FRef,0,m);


	}//end of COEF RELATIONSHIP (Li 1993)



	/*COEF RELATIONSHIP AND FRATERNITY (Wang 2002)*/
	if((Erwang || Edwang) && ploidy==2){
		a2kl=dmatrix(0,FRef,0,m);
		a3kl=dmatrix(0,FRef,0,m);
		a4kl=dmatrix(0,FRef,0,m);
		a22kl=dmatrix(0,FRef,0,m);
		ukl=dmatrix(0,FRef,0,m);
		Wkl=dmatrix(0,FRef,0,m);
		bkl=dmatrix(0,FRef,0,m);
		ckl=dmatrix(0,FRef,0,m);
		dkl=dmatrix(0,FRef,0,m);
		ekl=dmatrix(0,FRef,0,m);
		fkl=dmatrix(0,FRef,0,m);
		gkl=dmatrix(0,FRef,0,m);
		Vkl=dmatrix(0,FRef,0,m);
		N1kl=dmatrix(0,FRef,0,m);
		N2kl=dmatrix(0,FRef,0,m);
		N3kl=dmatrix(0,FRef,0,m);
		N4kl=dmatrix(0,FRef,0,m);
		N5kl=dmatrix(0,FRef,0,m);
		N6kl=dmatrix(0,FRef,0,m);
		N7kl=dmatrix(0,FRef,0,m);
		N8kl=dmatrix(0,FRef,0,m);
		Uk=dvector(0,m);
		WWang=dvector(0,m);
		P1=dvector(0,m);
		P2=dvector(0,m);
		P3=dvector(0,m);

		//Compute constant values according to allele freq (a-g coef + others)
		for(k=0;k<=FRef;k++){
			//single locus a coef
			for(l=1;l<=m;l++){
				a2kl[k][l]=a3kl[k][l]=a4kl[k][l]=0.;
				for(a=1;a<=Nallelel[l];a++) {
					a2kl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a];
					a3kl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a]*Pkla[k][l][a];
					a4kl[k][l]+=Pkla[k][l][a]*Pkla[k][l][a]*Pkla[k][l][a]*Pkla[k][l][a];
				}
				//sample size corrected a values
				if(FreqRef!=-1){
					a2kl[k][l]=(Nvalgenkl[k][l]*a2kl[k][l]-1.)/(Nvalgenkl[k][l]-1.);
					a3kl[k][l]=(Nvalgenkl[k][l]*Nvalgenkl[k][l]*a3kl[k][l]-3.*(Nvalgenkl[k][l]-1.)*a2kl[k][l]-1.)/((Nvalgenkl[k][l]-1.)*(Nvalgenkl[k][l]-2.));
					a4kl[k][l]=(Nvalgenkl[k][l]*Nvalgenkl[k][l]*Nvalgenkl[k][l]*a4kl[k][l]-6.*(Nvalgenkl[k][l]-1.)*(Nvalgenkl[k][l]-2.)*a3kl[k][l]-7.*(Nvalgenkl[k][l]-1.)*a2kl[k][l]-1.)/(Nvalgenkl[k][l]*Nvalgenkl[k][l]*Nvalgenkl[k][l]-6.*Nvalgenkl[k][l]*Nvalgenkl[k][l]+11.*Nvalgenkl[k][l]-6.);
				}

				a22kl[k][l]=a2kl[k][l]*a2kl[k][l];
				ukl[k][l]=2.*a2kl[k][l]-a3kl[k][l];
			}

			//weights per locus
			Uk[k]=0.;
			for(l=1;l<=m;l++) Uk[k]+=1./ukl[k][l];
			for(l=1;l<=m;l++) Wkl[k][l]=1./(Uk[k]*ukl[k][l]);

			//multilocus a coef
 			a2kl[k][0]=a3kl[k][0]=a4kl[k][0]=a22kl[k][0]=0.;
			for(l=1;l<=m;l++){
				a2kl[k][0]+=Wkl[k][l]*a2kl[k][l];
				a3kl[k][0]+=Wkl[k][l]*a3kl[k][l];
				a4kl[k][0]+=Wkl[k][l]*a4kl[k][l];
				a22kl[k][0]+=Wkl[k][l]*a22kl[k][l];
			}

			//multilocus and singlelocus b-g and other coef
			for(l=0;l<=m;l++){
				bkl[k][l]=2.*a22kl[k][l]-a4kl[k][l];
				ckl[k][l]=a2kl[k][l]-2.*a22kl[k][l]+a4kl[k][l];
				dkl[k][l]=4.*(a3kl[k][l]-a4kl[k][l]);
				ekl[k][l]=2.*(a2kl[k][l]-3.*a3kl[k][l]+2.*a4kl[k][l]);
				fkl[k][l]=4.*(a2kl[k][l]-a22kl[k][l]-2.*a3kl[k][l]+2.*a4kl[k][l]);
				gkl[k][l]=1.0f-7.*a2kl[k][l]+4.*a22kl[k][l]+10.*a3kl[k][l]-8.*a4kl[k][l];

				Vkl[k][l]=((1.0f-bkl[k][l])*(1.0f-bkl[k][l]))*(ekl[k][l]*ekl[k][l]*fkl[k][l]+dkl[k][l]*gkl[k][l]*gkl[k][l])-(1.0f-bkl[k][l])*(ekl[k][l]*fkl[k][l]-dkl[k][l]*gkl[k][l])*(ekl[k][l]*fkl[k][l]-dkl[k][l]*gkl[k][l])+2.*ckl[k][l]*dkl[k][l]*fkl[k][l]*(1.0f-bkl[k][l])*(gkl[k][l]+ekl[k][l])+ckl[k][l]*ckl[k][l]*dkl[k][l]*fkl[k][l]*(dkl[k][l]+fkl[k][l]);
				N1kl[k][l]=dkl[k][l]*fkl[k][l]*((ekl[k][l]+gkl[k][l])*(1.0f-bkl[k][l])+ckl[k][l]*(dkl[k][l]+fkl[k][l]));
				N2kl[k][l]=dkl[k][l]*(1.0f-bkl[k][l])*(gkl[k][l]*(1.0f-bkl[k][l]-dkl[k][l])+fkl[k][l]*(ckl[k][l]+ekl[k][l]));
				N3kl[k][l]=fkl[k][l]*(1.0f-bkl[k][l])*(ekl[k][l]*(1.0f-bkl[k][l]-fkl[k][l])+dkl[k][l]*(ckl[k][l]+gkl[k][l]));
				N4kl[k][l]=ckl[k][l]*dkl[k][l]*fkl[k][l]*(ekl[k][l]+gkl[k][l]);
				N5kl[k][l]=(1.0f-bkl[k][l])*(fkl[k][l]*ekl[k][l]*ekl[k][l]+dkl[k][l]*gkl[k][l]*gkl[k][l])-(ekl[k][l]*fkl[k][l]-dkl[k][l]*gkl[k][l])*(ekl[k][l]*fkl[k][l]-dkl[k][l]*gkl[k][l]);
				N6kl[k][l]=ckl[k][l]*(dkl[k][l]*gkl[k][l]-ekl[k][l]*fkl[k][l]);
				N7kl[k][l]=ckl[k][l]*ckl[k][l]*dkl[k][l]*fkl[k][l];
				N8kl[k][l]=ckl[k][l]*(1.0f-bkl[k][l]);
			}


		}/*end loop k*/


		/*weight per locus*/
		if(Erwang) for(l=1;l<=m;l++) corrSlij[Erwang][l][0][0]=(float)Wkl[0][l];
		if(Edwang) for(l=1;l<=m;l++) corrSlij[Edwang][l][0][0]=(float)Wkl[0][l];

		/*compute kinship for each pair of different indiv */
		k=0;
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];
				
				WWang[0]=a2temp=a3temp=a4temp=a22temp=P1[0]=P2[0]=P3[0]=0.;
				Nvalidloci=0;
				for(l=1;l<=m;l++){
					Swang=0;
					P1[l]=P2[l]=P3[l]=0.0;
					if(gilc[i][l][0]&&gilc[i][l][1]&&gilc[j][l][0]&&gilc[j][l][1]){	/*if no missing value*/
						Nvalidloci++;
						/*make the allele comparisons between i & j*/
						if(gilc[i][l][0]==gilc[j][l][0] || gilc[i][l][0]==gilc[j][l][1]) Swang++;
						if(gilc[i][l][1]==gilc[j][l][0] || gilc[i][l][1]==gilc[j][l][1]) Swang++;
						if(gilc[j][l][0]==gilc[i][l][0] || gilc[j][l][0]==gilc[i][l][1]) Swang++;
						if(gilc[j][l][1]==gilc[i][l][0] || gilc[j][l][1]==gilc[i][l][1]) Swang++;
						if(Swang==4) {P1[l]=1.;	P1[0]+=Wkl[k][l];}
						if(Swang==3) {P2[l]=1.; P2[0]+=Wkl[k][l];}
 						if(Swang==2) {P3[l]=1.; P3[0]+=Wkl[k][l];}

						WWang[0]+=Wkl[k][l];
						a2temp+=Wkl[k][l]*a2kl[k][l];
						a3temp+=Wkl[k][l]*a3kl[k][l];
						a4temp+=Wkl[k][l]*a4kl[k][l];
						a22temp+=Wkl[k][l]*a22kl[k][l];

						
						phiWang=(N1kl[k][l]*(P1[l]-1.)+N2kl[k][l]*P3[l]+N3kl[k][l]*P2[l])/Vkl[k][l];
						dWang=(N4kl[k][l]*(P1[l]+1.0f-2.*bkl[k][l])+N5kl[k][l]*(P1[l]-bkl[k][l])+N6kl[k][l]*(dkl[k][l]*P3[l]-fkl[k][l]*P2[l])-N7kl[k][l]*(P3[l]+P2[l]-dkl[k][l]-fkl[k][l])-N8kl[k][l]*(dkl[k][l]*gkl[k][l]*P3[l]+ekl[k][l]*fkl[k][l]*P2[l]))/Vkl[k][l];
						rWang=phiWang/2.+dWang;
						WWang[l]=Wkl[k][l];
					}
					else {rWang=dWang=(float)MISSVAL; WWang[l]=0.;}

					if(Erwang){corrSlij[Erwang][l][i][j]=(float)rWang;	corrSlij[Erwang][l][j][i]=(float)WWang[l];}
					if(Edwang){corrSlij[Edwang][l][i][j]=(float)dWang;	corrSlij[Edwang][l][j][i]=(float)WWang[l];}

				}/*end of loop l*/

				if(m>1){
					if(Nvalidloci){
						if(Nvalidloci==m){
							phiWang=(N1kl[k][0]*(P1[0]-1.)+N2kl[k][0]*P3[0]+N3kl[k][0]*P2[0])/Vkl[k][0];
							dWang=(N4kl[k][0]*(P1[0]+1.0f-2.*bkl[k][0])+N5kl[k][0]*(P1[0]-bkl[k][0])+N6kl[k][0]*(dkl[k][0]*P3[0]-fkl[k][0]*P2[0])-N7kl[k][0]*(P3[0]+P2[0]-dkl[k][0]-fkl[k][0])-N8kl[k][0]*(dkl[k][0]*gkl[k][0]*P3[0]+ekl[k][0]*fkl[k][0]*P2[0]))/Vkl[k][0];
							rWang=phiWang/2.+dWang;
						}
						else{
							P1[0]/=WWang[0];
							P2[0]/=WWang[0];
							P3[0]/=WWang[0];
							a2temp/=WWang[0];
							a3temp/=WWang[0];
							a4temp/=WWang[0];
							a22temp/=WWang[0];

							btemp=2.*a22temp-a4temp;
							ctemp=a2temp-2.*a22temp+a4temp;
							dtemp=4.*(a3temp-a4temp);
							etemp=2.*(a2temp-3.*a3temp+2.*a4temp);
							ftemp=4.*(a2temp-a22temp-2.*a3temp+2.*a4temp);
							gtemp=1.0f-7.*a2temp+4.*a22temp+10.*a3temp-8.*a4temp;

							Vtemp=((1.0f-btemp)*(1.0f-btemp))*(etemp*etemp*ftemp+dtemp*gtemp*gtemp)-(1.0f-btemp)*(etemp*ftemp-dtemp*gtemp)*(etemp*ftemp-dtemp*gtemp)+2.*ctemp*dtemp*ftemp*(1.0f-btemp)*(gtemp+etemp)+ctemp*ctemp*dtemp*ftemp*(dtemp+ftemp);
							N1temp=dtemp*ftemp*((etemp+gtemp)*(1.0f-btemp)+ctemp*(dtemp+ftemp));
							N2temp=dtemp*(1.0f-btemp)*(gtemp*(1.0f-btemp-dtemp)+ftemp*(ctemp+etemp));
							N3temp=ftemp*(1.0f-btemp)*(etemp*(1.0f-btemp-ftemp)+dtemp*(ctemp+gtemp));
							N4temp=ctemp*dtemp*ftemp*(etemp+gtemp);
							N5temp=(1.0f-btemp)*(ftemp*etemp*etemp+dtemp*gtemp*gtemp)-(etemp*ftemp-dtemp*gtemp)*(etemp*ftemp-dtemp*gtemp);
							N6temp=ctemp*(dtemp*gtemp-etemp*ftemp);
							N7temp=ctemp*ctemp*dtemp*ftemp;
							N8temp=ctemp*(1.0f-btemp);

							phiWang=(N1temp*(P1[0]-1.)+N2temp*P3[0]+N3temp*P2[0])/Vtemp;
							dWang=(N4temp*(P1[0]+1.0f-2.*btemp)+N5temp*(P1[0]-btemp)+N6temp*(dtemp*P3[0]-ftemp*P2[0])-N7temp*(P3[0]+P2[0]-dtemp-ftemp)-N8temp*(dtemp*gtemp*P3[0]+etemp*ftemp*P2[0]))/Vtemp;
							rWang=phiWang/2.+dWang;

						}/*end of else*/

					}
					else rWang=dWang=(float)MISSVAL;

					if(Erwang){corrSlij[Erwang][0][i][j]=(float)rWang;	corrSlij[Erwang][0][j][i]=(float)WWang[0];}
					if(Edwang){corrSlij[Edwang][0][i][j]=(float)dWang;	corrSlij[Edwang][0][j][i]=(float)WWang[0];}

				}/*end of if(m>1)*/

		/**/	if(JKl && m>1)for(l=1;l<=m;l++){  //l is the locus removed
					if(WWang[0]>Wkl[k][l]){
						if(WWang[l]){

							WWangJK=WWang[0]-Wkl[k][l];
							P1JK=(P1[0]*WWang[0]-Wkl[k][l]*P1[l])/WWangJK;
							P2JK=(P2[0]*WWang[0]-Wkl[k][l]*P2[l])/WWangJK;
							P3JK=(P3[0]*WWang[0]-Wkl[k][l]*P3[l])/WWangJK;
							a2JK=(a2temp*WWang[0]-Wkl[k][l]*a2kl[k][l])/WWangJK;
							a3JK=(a3temp*WWang[0]-Wkl[k][l]*a3kl[k][l])/WWangJK;
							a4JK=(a4temp*WWang[0]-Wkl[k][l]*a4kl[k][l])/WWangJK;
							a22JK=(a22temp*WWang[0]-Wkl[k][l]*a22kl[k][l])/WWangJK;

							bJK=2.*a22JK-a4JK;
							cJK=a2JK-2.*a22JK+a4JK;
							dJK=4.*(a3JK-a4JK);
							eJK=2.*(a2JK-3.*a3JK+2.*a4JK);
							fJK=4.*(a2JK-a22JK-2.*a3JK+2.*a4JK);
							gJK=1.0f-7.*a2JK+4.*a22JK+10.*a3JK-8.*a4JK;

							VJK=((1.0f-bJK)*(1.0f-bJK))*(eJK*eJK*fJK+dJK*gJK*gJK)-(1.0f-bJK)*(eJK*fJK-dJK*gJK)*(eJK*fJK-dJK*gJK)+2.*cJK*dJK*fJK*(1.0f-bJK)*(gJK+eJK)+cJK*cJK*dJK*fJK*(dJK+fJK);
							N1JK=dJK*fJK*((eJK+gJK)*(1.0f-bJK)+cJK*(dJK+fJK));
							N2JK=dJK*(1.0f-bJK)*(gJK*(1.0f-bJK-dJK)+fJK*(cJK+eJK));
							N3JK=fJK*(1.0f-bJK)*(eJK*(1.0f-bJK-fJK)+dJK*(cJK+gJK));
							N4JK=cJK*dJK*fJK*(eJK+gJK);
							N5JK=(1.0f-bJK)*(fJK*eJK*eJK+dJK*gJK*gJK)-(eJK*fJK-dJK*gJK)*(eJK*fJK-dJK*gJK);
							N6JK=cJK*(dJK*gJK-eJK*fJK);
							N7JK=cJK*cJK*dJK*fJK;
							N8JK=cJK*(1.0f-bJK);

							phiWang=(N1JK*(P1JK-1.)+N2JK*P3JK+N3JK*P2JK)/VJK;
							dWang=(N4JK*(P1JK+1.0f-2.*bJK)+N5JK*(P1JK-bJK)+N6JK*(dJK*P3JK-fJK*P2JK)-N7JK*(P3JK+P2JK-dJK-fJK)-N8JK*(dJK*gJK*P3JK+eJK*fJK*P2JK))/VJK;
							rWang=phiWang/2.+dWang;

						}
						else{  //when l was not informative so that removing it do not change the estimate
							if(Erwang)rWang=corrSlij[Erwang][0][i][j];
							if(Edwang)dWang=corrSlij[Edwang][0][i][j];
						}
					}
					else{	//when l was the sole informative locus
						rWang=(float)MISSVAL;
						dWang=(float)MISSVAL;
					}

					if(Erwang) corrSlij[Erwang][-l][i][j]=(float)rWang;
					if(Edwang) corrSlij[Edwang][-l][i][j]=(float)dWang;
				}
		 
			}/*end of loop j*/
		}/*end of loop i*/

		free_dmatrix(a2kl,0,FRef,0,m);
		free_dmatrix(a3kl,0,FRef,0,m);
		free_dmatrix(a4kl,0,FRef,0,m);
		free_dmatrix(a22kl,0,FRef,0,m);
		free_dmatrix(ukl,0,FRef,0,m);
		free_dmatrix(Wkl,0,FRef,0,m);
		free_dmatrix(bkl,0,FRef,0,m);
		free_dmatrix(ckl,0,FRef,0,m);
		free_dmatrix(dkl,0,FRef,0,m);
		free_dmatrix(ekl,0,FRef,0,m);
		free_dmatrix(fkl,0,FRef,0,m);
		free_dmatrix(gkl,0,FRef,0,m);
		free_dmatrix(Vkl,0,FRef,0,m);
		free_dmatrix(N1kl,0,FRef,0,m);
		free_dmatrix(N2kl,0,FRef,0,m);
		free_dmatrix(N3kl,0,FRef,0,m);
		free_dmatrix(N4kl,0,FRef,0,m);
		free_dmatrix(N5kl,0,FRef,0,m);
		free_dmatrix(N6kl,0,FRef,0,m);
		free_dmatrix(N7kl,0,FRef,0,m);
		free_dmatrix(N8kl,0,FRef,0,m);
		free_dvector(Uk,0,m);
		free_dvector(WWang,0,m);
		free_dvector(P1,0,m);
		free_dvector(P2,0,m);
		free_dvector(P3,0,m);

	}/*end of COEF RELATIONSHIP AND DOMINANCE (Wang 2002)*/

	
	
	
	free_f3tensor(Pkla,0,Ncat,0,m,0,maxa);
	free_imatrix(Nvalgenkl,0,Ncat,0,m);
	free_imatrix(Nallelekl,0,Ncat,0,m);
	free_ivector(Nivalidk,0,Ncat);
	free_matrix(Masizekl,0,Ncat,0,m);
	free_matrix(Vasizekl,0,Ncat,0,m);






	//Multilocus estimates (all estimators) with empirically derived weights
/*
	if(m>1)for(S=1;S<=NS;S++){
		//empirically derived weights
		for(l=1;l<=m;l++) Srli[l]=SSrli[l]=Nrli[l]=0.;
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];

				for(l=1;l<=m;l++) {
					rli=corrSlij[S][l][i][j];
					if(rli!=(float)MISSVAL){
						Nrli[l]++;
						Srli[l]+=rli;
						SSrli[l]+=rli*rli;
					}
				}
			}//loop j
		}//loop i

		for(l=1;l<=m;l++) {
			if(Nrli[l]>1){
				Var=(SSrli[l] - Srli[l]*Srli[l]/Nrli[l])/(Nrli[l]-1.);
				if(Var) corrSlij[S][l][0][0]=1./Var;
				else corrSlij[S][l][0][0]=0.;
			}
			else corrSlij[S][l][0][0]=0.;
		}
		
		for(i=1;i<n;i++){
			for(j=(i+1);j<=n;j++){
				if(TypeComp==1)if(cati[i]!=cati[j]) continue;
				if(TypeComp==2)if(cati[i]==cati[j]) continue;
				if(FreqRef==1 && cati[i]==cati[j]) k=cati[i];

				//multilocus estimates
				MLrli=SumW=0.;
				for(l=1;l<=m;l++)if(corrSlij[S][l][i][j]!=(float)MISSVAL){
					corrSlij[S][l][j][i]=corrSlij[S][l][0][0];
					MLrli+=corrSlij[S][l][0][0]*corrSlij[S][l][i][j];
					SumW+=corrSlij[S][l][0][0];
				}
				corrSlij[S][0][j][i]=SumW;
				if(SumW) corrSlij[S][0][i][j]=MLrli/SumW;
				else corrSlij[S][0][i][j]=(float)MISSVAL;

				//jackknife
				if(JKl)for(l=1;l<=m;l++){
					if(SumW==corrSlij[S][l][0][0]) corrSlij[S][-l][i][j]=(float)MISSVAL;
					else if(corrSlij[S][l][i][j]!=(float)MISSVAL) corrSlij[S][-l][i][j]=(MLrli-corrSlij[S][l][i][j]*corrSlij[S][l][0][0])/(SumW-corrSlij[S][l][0][0]);
					else corrSlij[S][-l][i][j]=corrSlij[S][0][i][j]; 
				}
			}//loop j
		}//loop i
	}//end of if(m>1)

*/	//end of Multilocus estimates (all estimators) with empirically derived weights


}/*end procedure compute_pairwise_corr_F*/



/****************************************************************************/


/*compute F-stat or R-stat between Npop>2 populations (Npop=# pop, pop1=pop2=0) or just
2 populations (Npop=2, pop1 & pop2 >0).
n :					# individuals (total) 
popi[i] :			n° of the population to which i belongs (popi >0).
m :					# of loci
Nallelel[l] :		# of alleles at locus l
gilc[i][l][a] :		allele of ind i at locus l for chromosome a
Allelesizela[l][a]: size of allele a at locus l (for microsatellites)

Stat :				defines the statistic to be computed and returned:
	Stat=1 : computes Fit, Fis, Fst and returns Fst
	Stat=2 : computes Fit, Fis, Fst and returns Fst/(1-Fst)
	Stat=3 : computes Fit, Fis, Fst and returns Rho
	Stat=4 : computes Fit, Fis, Fst and returns Rho/(1-Rho)
	Stat=5 : computes Rit, Ris, Rst and returns Rst
	Stat=6 : computes Rit, Ris, Rst and returns Rst/(1-Rst)
FstatSlr[S][l][r] values of Fit (r=1), Fis (r=2), Fst (r=3), Rho (r=4), 
						  Rit (r=1), Ris (r=2), Rst (r=3) at locus l (l=0 for multilocus)				

*/

void compute_F_R_stat(int n,int Npop,int pop1,int pop2,int *popi,int m,int *Nallelel,
		int ploidy,int ***gilc,int **allelesizela,int NS,int Stat[],
		float **FstatSlr[],float ***corrSlij[],int Rstat_only,int JKest)
{	
	int i,p,l,c,a,a1,a2,r,S,linit;	//counters for indv(i=1 to n), pop (p=1 to Npop), loci(l=1 to m), chromosomes (c=0 to ploidyi[i]-1), alleles (a=1 to Nallelel[l], statistic (r=1 to 4)
	int *Nip,Nipmax;	//# ind in pop i; max # over all pop
	int **Ncpi;			//#	chromosomes	with valid gene in ind i from pop  p
	int *nip, nc;		//ind n° in pop p (nip=1-Nip[p]; chromosome n° (nc=1-Ncpi[c][i])
	int *ploidyi;
	int **Nlpa[MMAX];		//# of allele a at locus l in pop p (totals for p=0 and a=0)
	double ***Gpic;		//value of the gene (indicator variable or allele size or breeding value) on chromosome a from ind i of pop p
	double **NumFlr, **DenFlr, **NumRlr, **DenRlr;	//numerator and denominator for F-stat / R-stat
	double SS[4],MS[4],s2[4];
	float **Flr,**Rlr;		//Flr =Fit (r=1), =Fis (r=2), =Fst (r=3), =Rho (r=4); Rlr =Rit (r=1), Ris (r=2), =Rst (r=3)
	float *Neil,*dm2l;	//Nei and DelaMuSquare distances
	int nJK;
	double S1,SS1,S2,SS2,D1,V1,V2,J1[MMAX],J2[MMAX],J12[MMAX],x1,x2,SJ1,SJ2,SJ12;//variables for Nei and Golstein distances
	int nvl = 0;

	//check whether Fstat (i.e. based on indicator variables) are to be computed
	if(Rstat_only==0){
		Rstat_only=1;
		for(S=1;S<=NS;S++)if(Stat[S]<=4 || Stat[S]==7) Rstat_only=0;
	}


	//check individual ploidy level
	ploidyi=ivector(0,n);
	for(i=1;i<=n;i++){
		ploidyi[i]=0;
		for(l=1;l<=m;l++){
			for(c=ploidy-1;c>=0;c--) if(gilc[i][l][c]) break;
			if(ploidyi[i]<(c+1)) ploidyi[i]=c+1;
		}
	}


	//compute the # of ind in each pop
	Nip=ivector(0,Npop);

	for(p=0;p<=Npop;p++) Nip[p]=0;
	if(pop1==0 && pop2==0){	  //for all pop
		for(i=1;i<=n;i++){
			if(popi[i]<=Npop) Nip[popi[i]]++;
			else {printf("\n\nERROR in ""compute_F-R_stat"": popi[%i]>Npop: %i>%i",i,popi[i],Npop); wait_a_char(); exit(33);} 
		}
		Nipmax=0;
		for(p=1;p<=Npop;p++) if(Nip[p]>Nipmax) Nipmax=Nip[p];
	}
	else if(Npop==2 && pop1 && pop2){	//for 2 pop defined as pop1 and pop2
		for(i=1;i<=n;i++){
			if(popi[i]==pop1) Nip[1]++;
			if(popi[i]==pop2) Nip[2]++;
		}
		if(Nip[1]>Nip[2]) Nipmax=Nip[1];
		else Nipmax=Nip[2];
	}
	else {printf("\n\nERROR in ""compute_F-R_stat"": pop1=%i, pop2=%i and Npop=%i",pop1,pop2,Npop); wait_a_char(); exit(35);}

	nip=ivector(0,n);
	NumFlr=dmatrix(0,m,0,4);
	DenFlr=dmatrix(0,m,0,4);
	NumRlr=dmatrix(0,m,0,4);
	DenRlr=dmatrix(0,m,0,4);
	Flr=matrix(-m,m+2,0,4);
	Rlr=matrix(-m,m+2,0,4);
	Neil=vector(-m,m+2);
	dm2l=vector(-m,m+2);
	Gpic=d3tensor(0,Npop,0,Nipmax,0,ploidy);	//Gpic[p][i][c] = value of the gene (indicator variable or allele size or breeding value) on chromosome a from ind i of pop p
	Ncpi=imatrix(0,Npop,0,Nipmax);			//Ncpi[p][i] = # of valid genes (chromosomes) in ind i from pop p
	for(l=0;l<=m;l++) Nlpa[l]=imatrix(0,Npop,0,Nallelel[l]);	 //counter for allele freq

	SJ1=SJ2=SJ12=0.;
	//Define the values for the Nested ANOVA	
	for(l=1;l<=m;l++){
		for(r=0;r<=4;r++) NumFlr[l][r]=DenFlr[l][r]=NumRlr[l][r]=DenRlr[l][r]=0.;
		for(p=0;p<=Npop;p++) for(a=0;a<=Nallelel[l];a++) Nlpa[l][p][a]=0;//allele freq counter

		for(a=0;a<=Nallelel[l];a++){  //loop over alleles (when a=0, R-stat are computed)
			if(Rstat_only) if(a>0) break;
			for(p=0;p<=Npop;p++) nip[p]=0;
			for(i=1;i<=n;i++){
				if(pop1==0 && pop2==0) p=popi[i];
				else{
					if(popi[i]==pop1) p=1;
					else if(popi[i]==pop2) p=2;
					else continue;
				}

				nip[p]++;
				nc=0;
				Gpic[p][nip[p]][0]=0.; //reset mean individual allele a freq
				for(c=0;c<ploidyi[i];c++){	 //loop over chromosomes
					if(gilc[i][l][c]) nc++;
					else continue;

					if(a==0){  //Rstat computed
						Gpic[p][nip[p]][nc]=allelesizela[l][gilc[i][l][c]];
						Nlpa[l][p][gilc[i][l][c]]++;  //for allele freq
					}
					if(a>0){  //Fstat computed, with Gpic as indicator variable for allele a
						if(gilc[i][l][c]==a) Gpic[p][nip[p]][nc]=1.;
						else Gpic[p][nip[p]][nc]=0.;

						Gpic[p][nip[p]][0]+=Gpic[p][nip[p]][nc]; //sum of allele freq to compute mean ind allele a freq
					}
				}
				if(nc){
					Ncpi[p][nip[p]]=nc;
					Gpic[p][nip[p]][0]/=nc;	
				}
				else nip[p]--;	//remove one ind as there is no data for it
			} //end of loop i

			NestedANOVA(Npop,nip,Ncpi,Gpic,SS,MS,s2);

			if(a==0){
				if((s2[2]+s2[3])!=(float)MISSVAL){			   
					NumRlr[l][2]=(s2[2]);			   //Ris
					DenRlr[l][2]=(s2[2]+s2[3]);
				}
				else NumRlr[l][2]=DenRlr[l][2]=(float)MISSVAL;

				if((s2[1]+s2[2]+s2[3])!=(float)MISSVAL){
					NumRlr[l][1]=(s2[1]+s2[2]);			//Rit
					DenRlr[l][1]=(s2[1]+s2[2]+s2[3]);
					NumRlr[l][3]=(s2[1]);				//Rst
					DenRlr[l][3]=(s2[1]+s2[2]+s2[3]);
				}
				else if((s2[1]+s2[2])!=(float)MISSVAL){  //Rst in case of haploids
					NumRlr[l][3]+=(s2[1]);			//Rst
					DenRlr[l][3]+=(s2[1]+s2[2]);
				}
				else NumRlr[l][1]=DenRlr[l][1]=NumRlr[l][3]=DenRlr[l][3]=(float)MISSVAL;
			}
			if(a>0){
				if((s2[2]+s2[3])!=(float)MISSVAL){
					NumFlr[l][2]+=(s2[2]);			   //Fis
					DenFlr[l][2]+=(s2[2]+s2[3]);

							//unweighted allele F-stat
							/*	if((s2[2]+s2[3])){	  
									Fbis[2]+=(s2[2]/(s2[2]+s2[3]));				
									naFbis[2]+=1.;
								}	*/
				}
				if((s2[1]+s2[2]+s2[3])!=(float)MISSVAL){
					NumFlr[l][1]+=(s2[1]+s2[2]);	   //Fit
					DenFlr[l][1]+=(s2[1]+s2[2]+s2[3]);
					NumFlr[l][3]+=(s2[1]);			   //Fst
					DenFlr[l][3]+=(s2[1]+s2[2]+s2[3]);

						 //unweighted allele F-stat
						/*	if((s2[1]+s2[2]+s2[3])){	  
									Fbis[1]+=((s2[1]+s2[2])/(s2[1]+s2[2]+s2[3]));  
									naFbis[1]+=1.;
									Fbis[3]+=(s2[1]/(s2[1]+s2[2]+s2[3]));
									naFbis[3]+=1.;
								} */
				}
				else if((s2[1]+s2[2])!=(float)MISSVAL){  //Fst in case of haploids
					NumFlr[l][3]+=(s2[1]);			//Fst
					DenFlr[l][3]+=(s2[1]+s2[2]);
				}

				//Rho stat
				for(p=1;p<=Npop;p++)for(i=1;i<=Nip[p];i++){
					if(Gpic[p][i][0]!=(float)MISSVAL){
						Gpic[p][i][1]=Gpic[p][i][0];
						Ncpi[p][i]=1;
					}
					else Ncpi[p][i]=0;
				}
				NestedANOVA(Npop,nip,Ncpi,Gpic,SS,MS,s2); 
				if((s2[1]+s2[2])!=(float)MISSVAL){
					NumFlr[l][4]+=(s2[1]);		   //Rho
					DenFlr[l][4]+=(s2[1]+s2[2]);
				}
			}
		}//end loop over alleles
		for(r=1;r<=4;r++){
			if(DenFlr[l][r] && DenFlr[l][r]!=(float)MISSVAL){
				Flr[l][r]=(float)(NumFlr[l][r]/DenFlr[l][r]);
			}
			else Flr[l][r]=(float)MISSVAL;

			if(DenRlr[l][r] && DenRlr[l][r]!=(float)MISSVAL){
				Rlr[l][r]=(float)(NumRlr[l][r]/DenRlr[l][r]);
			}
			else Rlr[l][r]=(float)MISSVAL;

		}

		//compute Nei and Goldstein dist
		if(pop1 && pop2){
			for(a=1;a<=Nallelel[l];a++) {Nlpa[l][1][0]+=Nlpa[l][1][a]; Nlpa[l][2][0]+=Nlpa[l][2][a];}
			if(Nlpa[l][1][0]>1 && Nlpa[l][2][0]>1){
				//Nei's unbiased standard distance
				J1[l]=J2[l]=J12[l]=0;
				for(a=1;a<=Nallelel[l];a++){
					x1=1.*Nlpa[l][1][a]/Nlpa[l][1][0];
					x2=1.*Nlpa[l][2][a]/Nlpa[l][2][0];
					J1[l]+=x1*x1;
					J2[l]+=x2*x2;
					J12[l]+=x1*x2;
				}
				J1[l]=(Nlpa[l][1][0]*J1[l]-1.)/(Nlpa[l][1][0]-1.);
				J2[l]=(Nlpa[l][2][0]*J2[l]-1.)/(Nlpa[l][2][0]-1.);
				Neil[l]=(float)-log(J12[l]/sqrt(J1[l]*J2[l]));
				SJ1+=J1[l];
				SJ2+=J2[l];
				SJ12+=J12[l];

				//Goldstein's distance
				S1=S2=SS1=SS2=0.;
				for(a=1;a<=Nallelel[l];a++){
					S1+=allelesizela[l][a]*Nlpa[l][1][a];
					S2+=allelesizela[l][a]*Nlpa[l][2][a];
					SS1+=allelesizela[l][a]*allelesizela[l][a]*Nlpa[l][1][a];
					SS2+=allelesizela[l][a]*allelesizela[l][a]*Nlpa[l][2][a];
				}
				V1=Nlpa[l][1][0]/(Nlpa[l][1][0]-1.)*( (SS1/Nlpa[l][1][0])-(S1/Nlpa[l][1][0])*(S1/Nlpa[l][1][0]) );
				V2=Nlpa[l][2][0]/(Nlpa[l][2][0]-1.)*( (SS2/Nlpa[l][2][0])-(S2/Nlpa[l][2][0])*(S2/Nlpa[l][2][0]) );
				D1=0.;
				for(a1=1;a1<=Nallelel[l];a1++) for(a2=1;a2<=Nallelel[l];a2++) D1+=(allelesizela[l][a1]-allelesizela[l][a2])*(allelesizela[l][a1]-allelesizela[l][a2])*(1.*Nlpa[l][1][a1]/Nlpa[l][1][0])*(1.*Nlpa[l][2][a2]/Nlpa[l][2][0]);
				dm2l[l]=(float)(D1-V1-V2);

			}
			else Neil[l]=dm2l[l]=(float)MISSVAL;
		}

	}//end loop over loci

	//multilocus estimates
	for(r=1;r<=4;r++) NumFlr[0][r]=DenFlr[0][r]=NumRlr[0][r]=DenRlr[0][r]=0.;
	for(l=1;l<=m;l++)for(r=1;r<=4;r++){
		if(DenFlr[l][r] && DenFlr[l][r]!=(float)MISSVAL){
			NumFlr[0][r]+=NumFlr[l][r];
			DenFlr[0][r]+=DenFlr[l][r];
		}
		if(DenRlr[l][r] && DenRlr[l][r]!=(float)MISSVAL){
			NumRlr[0][r]+=NumRlr[l][r];
			DenRlr[0][r]+=DenRlr[l][r];
		}
	}
	for(r=1;r<=4;r++){
		if(DenFlr[0][r]) Flr[0][r]=(float)(NumFlr[0][r]/DenFlr[0][r]);	//multilocus average
		else Flr[0][r]=(float)MISSVAL;
		if(DenRlr[0][r]) Rlr[0][r]=(float)(NumRlr[0][r]/DenRlr[0][r]);	//multilocus average
		else Rlr[0][r]=(float)MISSVAL;
	}
	if(pop1 && pop2){  //Nei and Goldstein dist
		Neil[0]=(float)-log(SJ12/sqrt(SJ1*SJ2));
		nvl=0;
		dm2l[0]=0.0f;
		for(l=1;l<=m;l++) if(dm2l[l]!=(float)MISSVAL){ nvl++; dm2l[0]+=dm2l[l];}
		if(nvl) dm2l[0]/=nvl;
		else dm2l[0]=(float)MISSVAL;
	}

	//Jackknife over loci for global analysis
	if(JKest && pop1==0 && pop2==0)for(r=1;r<=4;r++){
		if(Flr[0][r]!=(float)MISSVAL){//Fstat
			nJK=0;
			for(l=1;l<=m;l++){ //compute stat for all loci but locus l
				if(Flr[l][r]!=(float)MISSVAL){
					Flr[-l][r]=(float)((NumFlr[0][r]-NumFlr[l][r])/(DenFlr[0][r]-DenFlr[l][r]));
					nJK++;
				}
				else Flr[-l][r]=(float)MISSVAL;
			}
			if(nJK>1){
				Flr[m+1][r]=Flr[m+2][r]=Rlr[m+1][r]=Rlr[m+2][r]=0.;
				for(l=1;l<=m;l++)if(Flr[-l][r]!=(float)MISSVAL){ //compute pseudovalues
					Flr[-l][r]=nJK*Flr[0][r]-(nJK-1)*Flr[-l][r];
					Flr[m+1][r]+=Flr[-l][r];   
				}
				Flr[m+1][r]/=nJK;   //mean JK estimate
				for(l=1;l<=m;l++)if(Flr[-l][r]!=(float)MISSVAL){ //compute variance of pseudoval
					Flr[m+2][r]+=(Flr[-l][r]-Flr[m+1][r])*(Flr[-l][r]-Flr[m+1][r]);
				}
				Flr[m+2][r]=(float)sqrt((Flr[m+2][r]/(nJK-1))/nJK); //standard error JK estimate
			}
			else Flr[m+1][r]=Flr[m+2][r]=(float)MISSVAL;
		}
		else Flr[m+1][r]=Flr[m+2][r]=(float)MISSVAL;

		if(Rlr[0][r]!=(float)MISSVAL){  //Rsat
			nJK=0;
			for(l=1;l<=m;l++){ //compute stat for all loci but locus l
				if(Rlr[l][r]!=(float)MISSVAL){
					Rlr[-l][r]=(float)((NumRlr[0][r]-NumRlr[l][r])/(DenRlr[0][r]-DenRlr[l][r]));
					nJK++;
				}
				else Rlr[-l][r]=(float)MISSVAL;
			}
			if(nJK>1){
				for(l=1;l<=m;l++)if(Rlr[-l][r]!=(float)MISSVAL){ //compute pseudovalues
					Rlr[-l][r]=nJK*Rlr[0][r]-(nJK-1)*Rlr[-l][r];
					Rlr[m+1][r]+=Rlr[-l][r];   
				}
				Rlr[m+1][r]/=nJK;   //mean JK estimate
				for(l=1;l<=m;l++)if(Rlr[-l][r]!=(float)MISSVAL){ //compute variance of pseudoval
					Rlr[m+2][r]+=(Rlr[-l][r]-Rlr[m+1][r])*(Rlr[-l][r]-Rlr[m+1][r]);
				}
				Rlr[m+2][r]=(float)sqrt((Rlr[m+2][r]/(nJK-1))/nJK); //standard error JK estimate
			}
			else Rlr[m+1][r]=Rlr[m+2][r]=(float)MISSVAL;
		}
		else Rlr[m+1][r]=Rlr[m+2][r]=(float)MISSVAL;
	}

	//Multilocus estimates for all loci but one, for pairwise analysis
	if(JKest && pop1 && pop2){
		for(r=1;r<=4;r++)for(l=1;l<=m;l++){ //compute stat for all loci but locus l
			if(Flr[l][r]!=(float)MISSVAL) Flr[-l][r]=(float)((NumFlr[0][r]-NumFlr[l][r])/(DenFlr[0][r]-DenFlr[l][r]));
			else Flr[-l][r]=Flr[0][r];	
			if(Rlr[l][r]!=(float)MISSVAL) Rlr[-l][r]=(float)((NumRlr[0][r]-NumRlr[l][r])/(DenRlr[0][r]-DenRlr[l][r]));
			else Rlr[-l][r]=Rlr[0][r];	
		}

		for(l=1;l<=m;l++){
			if(Neil[l]!=(float)MISSVAL) Neil[-l]=(float)-log( (SJ12-J12[l])/sqrt((SJ1-J1[l])*(SJ2-J2[l])) );
			else Neil[-l]=Neil[0];

			if(dm2l[l]!=(float)MISSVAL) dm2l[-l]=(dm2l[0]*nvl-dm2l[l])/(nvl-1.0f);
			else dm2l[-l]=dm2l[0];
		}
	}

	//for global analysis
	if(pop1==0 && pop2==0)for(r=1;r<=3;r++)for(l=0;l<=m+2;l++)for(S=1;S<=NS;S++)if(Stat[S]<=4 || Stat[S]==7) FstatSlr[S][l][r]=Flr[l][r];
	if(pop1==0 && pop2==0) for(l=0;l<=m+2;l++)for(S=1;S<=NS;S++)if(Stat[S]==3 || Stat[S]==4) FstatSlr[S][l][4]=Flr[l][4];
	if(pop1==0 && pop2==0)for(r=1;r<=3;r++)for(l=0;l<=m+2;l++)for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8) FstatSlr[S][l][r]=Rlr[l][r];
	//for pairwise analysis
	if(JKest)linit=-m;
	else linit=0;
	if(m==1) linit=1;
	if(pop1 && pop2)for(l=linit;l<=m;l++)for(S=1;S<=NS;S++){//for pairwise comparisons
		if(Stat[S]==1) corrSlij[S][l][pop1][pop2]=Flr[l][3];				//Fst
		if(Stat[S]==2) corrSlij[S][l][pop1][pop2]=Flr[l][3]/(1.0f-Flr[l][3]);	//Fst/(1-Fst)
		if(Stat[S]==3) corrSlij[S][l][pop1][pop2]=Flr[l][4];				//Rho
		if(Stat[S]==4) corrSlij[S][l][pop1][pop2]=Flr[l][4]/(1.0f-Flr[l][4]);	//Rho/(1-Rho)
		if(Stat[S]==5) corrSlij[S][l][pop1][pop2]=Rlr[l][3];				//Rst
		if(Stat[S]==6) corrSlij[S][l][pop1][pop2]=Rlr[l][3]/(1.0f-Rlr[l][3]);	//Rst/(1-Rst)
		if(Stat[S]==7) corrSlij[S][l][pop1][pop2]=Neil[l];	//Nei's standard dist
		if(Stat[S]==8) corrSlij[S][l][pop1][pop2]=dm2l[l];	//Goldstein distance

		//weights of per locus estimates
		if(l>=0){
			if(DenFlr[l][3]==(float)MISSVAL) DenFlr[l][3]=0.;
			if(DenFlr[l][4]==(float)MISSVAL) DenFlr[l][4]=0.;
			if(DenRlr[l][3]==(float)MISSVAL) DenRlr[l][3]=0.;
			if(Stat[S]==1) corrSlij[S][l][pop2][pop1]=(float)DenFlr[l][3];	//Fst
			if(Stat[S]==2) corrSlij[S][l][pop2][pop1]=(float)DenFlr[l][3];	//Fst/(1-Fst)
			if(Stat[S]==3) corrSlij[S][l][pop2][pop1]=(float)DenFlr[l][4];	//Rho
			if(Stat[S]==4) corrSlij[S][l][pop2][pop1]=(float)DenFlr[l][4];	//Rho/(1-Rho)
			if(Stat[S]==5) corrSlij[S][l][pop2][pop1]=(float)DenRlr[l][3];	//Rst
			if(Stat[S]==6) corrSlij[S][l][pop2][pop1]=(float)DenRlr[l][3];	//Rst/(1-Rst)
		}
	}


	free_ivector(ploidyi,0,n);
	free_ivector(nip,0,n);
	free_matrix(Flr,-m,m+2,0,4);
	free_matrix(Rlr,-m,m+2,0,4);
	free_vector(Neil,-m,m+2);
	free_vector(dm2l,-m,m+2);
	free_dmatrix(NumRlr,0,m,0,4);
	free_dmatrix(DenRlr,0,m,0,4);
	free_dmatrix(NumFlr,0,m,0,4);
	free_dmatrix(DenFlr,0,m,0,4);
	free_ivector(Nip,0,Npop);
	free_d3tensor(Gpic,0,Npop,0,Nipmax,0,ploidy);
	free_imatrix(Ncpi,0,Npop,0,Nipmax);
	for(l=0;l<=m;l++) free_imatrix(Nlpa[l],0,Npop,0,Nallelel[l]);

}


/****************************************************************************/

	
/*compute G-stat or N-stat between populations (Npop=# pop, pop1=pop2=0) .
n :					# individuals (total) 
popi[i] :			n° of the population to which i belongs (popi >0).
m :					# of loci
Nallelel[l] :		# of alleles at locus l
gilc[i][l][a] :		allele of ind i at locus l for chromosome a
Allelesizela[l][a]: size of allele a at locus l (for microsatellites)

Stat :				defines the statistic to be computed and returned:
	Stat=9 : computes Gstat and pairwise Gst
	Stat=10 : computes Gstat and pairwise Gij
	Stat=11 : computes Nstat and pairwise Nst
	Stat=12 : computes Nstat and pairwise Nij

G-Nstatlr[l][r] values of Git (r=1), Gis (r=2), Gst (r=3),  
						  Nit (r=5), Nis (r=6), Nst (r=7) at locus l (l=0 for multilocus)				

*/

void compute_G_N_stat(int n,int Npop,int pop1,int pop2,int *popi,int m,int *Nallelel,
		int ploidy,int ***gilc,float ***Ppla,float ***Distla1a2,int NS,int Stat[],
		float **FstatSlr[],float ***corrSlij[],int JKest,int computeallelefreq)
{	
	int i,p,p1,p2,l,c,a,a1,a2,r,S,linit,Ncvalid;	//counters for indv(i=1 to n), pop (p=1 to Npop), loci(l=1 to m), chromosomes (c=0 to ploidyi[i]-1), alleles (a=1 to Nallelel[l], statistic (r=1 to 7)
	int **Ncpl = NULL, *Nip = NULL;	//#	chromosomes	with valid gene in pop  p at locus l
	int **Nlpa[MMAX];	//# of allele a at locus l in pop p (totals for p=0 and a=0)
	double **Plpa[MMAX];		//freq of allele a at locus l in pop p 
	double **NumGlr, **DenGlr, **NumNlr, **DenNlr;	//numerator and denominator for G-stat / N-stat
	float **Glr;		//=Git (r=1), =Gis (r=2), =Gst (r=3) 
	float **Nlr;		//=Nit (r=1), Nis (r=2), =Nst (r=3)
	double ***DivNlp1p2,***DivGlp1p2,**DivNrefp1p2,**DivGrefp1p2; //diversity estimates	per pair of pop
	double DivNwl[MMAX],DivGwl[MMAX],DivNtl[MMAX],DivGtl[MMAX];//diversity within pop and total
	double ***Nlp1p2,***Glp1p2,***Nstlp1p2,***Gstlp1p2;		//Nij, Gij, and pairwise Nst and Gst values per pop pair
	int Nvalpop,Nvalloc;
	int nJK;


	
	//Compute allele freq	
	for(l=0;l<=m;l++) Plpa[l]=dmatrix(0,Npop,0,Nallelel[l]);
	if(computeallelefreq){
		Ncpl=imatrix(0,Npop,0,m);			//Ncpl[p][l] = # of valid genes (chromosomes) in pop p at locus l
		Nip=ivector(0,Npop);				//#ind with valid genes per pop p
		for(l=0;l<=m;l++) Nlpa[l]=imatrix(0,Npop,0,Nallelel[l]);	 //counter for allele freq

		for(l=1;l<=m;l++){
			for(p=0;p<=Npop;p++) for(a=0;a<=Nallelel[l];a++) Nlpa[l][p][a]=0;
			for(p=0;p<=Npop;p++) Ncpl[p][l]=Nip[p]=0;

			for(i=1;i<=n;i++){
				if(pop1==0 && pop2==0) p=popi[i];
				else{
					if(popi[i]==pop1) p=1;
					else if(popi[i]==pop2) p=2;
					else continue;
				}
				Ncvalid=0;
				for(c=0;c<ploidy;c++)if(gilc[i][l][c]) Ncvalid++;
				if(Ncvalid){
					Nip[p]++;
					for(c=0;c<ploidy;c++)if(gilc[i][l][c]){	 //loop over chromosomes
						Ncpl[p][l]++;
						// TODO: Is this a bug?
						Nlpa[l][p][gilc[i][l][c]]+=(int)(1./Ncvalid);
					}
				}
			} //end of loop i

			for(a=1;a<=Nallelel[l];a++)for(p=1;p<=Npop;p++){
				if(Nip[p]) Plpa[l][p][a]=(double)Nlpa[l][p][a]/(double)Nip[p];
				else Plpa[l][p][a]=(float)MISSVAL;
			}
			for(p=1;p<=Npop;p++) Plpa[l][p][0]=Ncpl[p][l]; 	//put # of chromosomes with defined genes in Plpa[l][p][0]

		}//end loop over loci
	}
	else for(l=1;l<=m;l++)for(p=1;p<=Npop;p++)for(a=0;a<=Nallelel[l];a++) Plpa[l][p][a]=Ppla[p][l][a];


	//Compute diversity measures
	DivNlp1p2=d3tensor(0,m,0,Npop,0,Npop);
	DivGlp1p2=d3tensor(0,m,0,Npop,0,Npop);
	DivNrefp1p2=dmatrix(0,Npop,0,Npop);
	DivGrefp1p2=dmatrix(0,Npop,0,Npop);

	for(l=1;l<=m;l++){
		//diversity btw each pair of pop, including within each pop
		for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			
			if(Plpa[l][p1][0]>1 && Plpa[l][p2][0]>1){
				DivNlp1p2[l][p1][p2]=DivGlp1p2[l][p1][p2]=0.;
				for(a1=1;a1<=Nallelel[l];a1++)if(Plpa[l][p1][a1])for(a2=1;a2<=Nallelel[l];a2++){
					if(Distla1a2[0][0][0]) DivNlp1p2[l][p1][p2]+=Distla1a2[l][a1][a2]*Plpa[l][p1][a1]*Plpa[l][p2][a2];
					if(a1!=a2) DivGlp1p2[l][p1][p2]+=1.*Plpa[l][p1][a1]*Plpa[l][p2][a2];
				}
				if(p1==p2){
					DivNlp1p2[l][p1][p1]*=Plpa[l][p1][0]/(Plpa[l][p1][0]-1.);
					DivGlp1p2[l][p1][p1]*=Plpa[l][p1][0]/(Plpa[l][p1][0]-1.);
				}
			}   
			else DivNlp1p2[l][p1][p2]=DivGlp1p2[l][p1][p2]=(float)MISSVAL;
		}
		//mean diversity per pop
		DivNwl[l]=DivGwl[l]=0.;
		Nvalpop=0;
		for(p=1;p<=Npop;p++)if(DivNlp1p2[l][p][p]!=(float)MISSVAL){
			DivNwl[l]+=DivNlp1p2[l][p][p];
			DivGwl[l]+=DivGlp1p2[l][p][p];
			Nvalpop++;
		}
		if(Nvalpop){DivNwl[l]/=Nvalpop;	DivGwl[l]/=Nvalpop;}
		else DivNwl[l]=DivGwl[l]=(float)MISSVAL;
		//total diversity
		DivNtl[l]=DivGtl[l]=0.;
		Nvalpop=0;
		for(p1=1;p1<Npop;p1++)for(p2=p1+1;p2<=Npop;p2++)if(DivNlp1p2[l][p1][p2]!=(float)MISSVAL){
			DivNtl[l]+=DivNlp1p2[l][p1][p2];
			DivGtl[l]+=DivGlp1p2[l][p1][p2];
			Nvalpop++;
		}
		if(Nvalpop){DivNtl[l]/=Nvalpop; DivGtl[l]/=Nvalpop;}
		else DivNtl[l]=DivGtl[l]=(float)MISSVAL;

	}
	//multilocus estimators
	//for each pair of pop
	for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
		Nvalloc=0;
		DivNlp1p2[0][p1][p2]=DivGlp1p2[0][p1][p2]=DivNrefp1p2[p1][p2]=DivGrefp1p2[p1][p2]=0.;
		for(l=1;l<=m;l++) if(DivNlp1p2[l][p1][p2]!=(float)MISSVAL){
			DivNlp1p2[0][p1][p2]+=DivNlp1p2[l][p1][p2];
			DivNrefp1p2[p1][p2]+=DivNtl[l];
			DivGlp1p2[0][p1][p2]+=DivGlp1p2[l][p1][p2];
			DivGrefp1p2[p1][p2]+=DivGtl[l];
			Nvalloc++;
		}
/*		DivNlp1p2[0][p1][p2]/=Nvalloc;
		DivGlp1p2[0][p1][p2]/=Nvalloc;		*/
	}
	//for within pop
	Nvalloc=0;
	DivNwl[0]=DivGwl[0]=0.;
	for(l=1;l<=m;l++) if(DivNwl[l]!=(float)MISSVAL){
		DivNwl[0]+=DivNwl[l];
		DivGwl[0]+=DivGwl[l];
		Nvalloc++;
	}
/*	DivNwl[0]/=Nvalloc;
	DivGwl[0]/=Nvalloc;				*/
	// for total diversity
	Nvalloc=0;
	DivNtl[0]=DivGtl[0]=0.;
	for(l=1;l<=m;l++) if(DivNtl[l]!=(float)MISSVAL){
		DivNtl[0]+=DivNtl[l];
		DivGtl[0]+=DivGtl[l];
		Nvalloc++;
	}
/*	DivNtl[0]/=Nvalloc;
	DivGtl[0]/=Nvalloc;				*/


	//Compute differentiations
	NumGlr=dmatrix(0,m,0,4);
	DenGlr=dmatrix(0,m,0,4);
	NumNlr=dmatrix(0,m,0,4);
	DenNlr=dmatrix(0,m,0,4);
	Glr=matrix(-m,m+2,0,4);
	Nlr=matrix(-m,m+2,0,4);
	for(l=0;l<=m+2;l++)for(r=1;r<=4;r++) Glr[l][r]=Nlr[l][r]=(float)MISSVAL;
	Glp1p2=d3tensor(-m,m,0,Npop,0,Npop);
	Nlp1p2=d3tensor(-m,m,0,Npop,0,Npop);
	Gstlp1p2=d3tensor(-m,m,0,Npop,0,Npop);
	Nstlp1p2=d3tensor(-m,m,0,Npop,0,Npop);

	for(l=0;l<=m;l++)if(DivNtl[l]!=(float)MISSVAL){
		for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			//pairwise Gij
			if(DivGlp1p2[l][p1][p2]!=(float)MISSVAL){
				if(l==0) Glp1p2[l][p1][p2]=1.0f-DivGlp1p2[l][p1][p2]/DivGrefp1p2[p1][p2];
				else Glp1p2[l][p1][p2]=1.0f-DivGlp1p2[l][p1][p2]/DivGtl[l];
			}
			else Glp1p2[l][p1][p2]=(float)MISSVAL;
			//pairwise Nij
			if(DivNlp1p2[l][p1][p2]!=(float)MISSVAL){
				if(l==0) Nlp1p2[l][p1][p2]=1.0f-DivNlp1p2[l][p1][p2]/DivNrefp1p2[p1][p2];
				else Nlp1p2[l][p1][p2]=1.0f-DivNlp1p2[l][p1][p2]/DivNtl[l];
			}
			else Nlp1p2[l][p1][p2]=(float)MISSVAL;

			//pairwise Gst
			if(DivGlp1p2[l][p1][p1]!=(float)MISSVAL && DivGlp1p2[l][p2][p2]!=(float)MISSVAL && DivGlp1p2[l][p1][p2]!=(float)MISSVAL && DivGlp1p2[l][p1][p2]){
				Gstlp1p2[l][p1][p2]=1.0f-(0.5*DivGlp1p2[l][p1][p1]+0.5*DivGlp1p2[l][p2][p2])/DivGlp1p2[l][p1][p2];
			}
			else Gstlp1p2[l][p1][p2]=(float)MISSVAL;
			//pairwise Nst
			if(DivNlp1p2[l][p1][p1]!=(float)MISSVAL && DivNlp1p2[l][p2][p2]!=(float)MISSVAL && DivNlp1p2[l][p1][p2]!=(float)MISSVAL && DivNlp1p2[l][p1][p2]){
				Nstlp1p2[l][p1][p2]=1.0f-(0.5*DivNlp1p2[l][p1][p1]+0.5*DivNlp1p2[l][p2][p2])/DivNlp1p2[l][p1][p2];
			}
			else Nstlp1p2[l][p1][p2]=(float)MISSVAL;

		}
		//global Nst and Gst
		if(DivNwl[l]!=(float)MISSVAL) Nlr[l][3]=(float)(1.0-DivNwl[l]/DivNtl[l]);
		else Nlr[l][3]=(float)MISSVAL;

		if(DivGwl[l]!=(float)MISSVAL) Glr[l][3]=(float)(1.0-DivGwl[l]/DivGtl[l]);
		else Glr[l][3]=(float)MISSVAL;

		NumGlr[l][3]=DivGtl[l]-DivGwl[l];
		DenGlr[l][3]=DivGtl[l];
		NumNlr[l][3]=DivNtl[l]-DivNwl[l];
		DenNlr[l][3]=DivNtl[l];
	}


	//Jackknife over loci for global analysis
	if(JKest && pop1==0 && pop2==0)for(r=1;r<=3;r++){
		//Gstat
		if(Glr[0][r]!=(float)MISSVAL){
			nJK=0;
			for(l=1;l<=m;l++){ //compute stat for all loci but locus l
				if(Glr[l][r]!=(float)MISSVAL){
					Glr[-l][r]=(float)((NumGlr[0][r]-NumGlr[l][r])/(DenGlr[0][r]-DenGlr[l][r]));
					nJK++;
				}
				else Glr[-l][r]=(float)MISSVAL;
			}
			if(nJK>1){
				Glr[m+1][r]=Glr[m+2][r]=Nlr[m+1][r]=Nlr[m+2][r]=0.;
				for(l=1;l<=m;l++)if(Glr[-l][r]!=(float)MISSVAL){ //compute pseudovalues
					Glr[-l][r]=nJK*Glr[0][r]-(nJK-1)*Glr[-l][r];
					Glr[m+1][r]+=Glr[-l][r];   
				}
				Glr[m+1][r]/=nJK;   //mean JK estimate
				for(l=1;l<=m;l++)if(Glr[-l][r]!=(float)MISSVAL){ //compute variance of pseudoval
					Glr[m+2][r]+=(Glr[-l][r]-Glr[m+1][r])*(Glr[-l][r]-Glr[m+1][r]);
				}
				Glr[m+2][r]=(float)sqrt((Glr[m+2][r]/(nJK-1))/nJK); //standard error JK estimate
			}
			else Glr[m+1][r]=Glr[m+2][r]=(float)MISSVAL;
		}
		else Glr[m+1][r]=Glr[m+2][r]=(float)MISSVAL;
		//Nstat
		if(Nlr[0][r]!=(float)MISSVAL){
			nJK=0;
			for(l=1;l<=m;l++){ //compute stat for all loci but locus l
				if(Nlr[l][r]!=(float)MISSVAL){
					Nlr[-l][r]=(float)((NumNlr[0][r]-NumNlr[l][r])/(DenNlr[0][r]-DenNlr[l][r]));
					nJK++;
				}
				else Nlr[-l][r]=(float)MISSVAL;
			}
			if(nJK>1){
				Nlr[m+1][r]=Nlr[m+2][r]=0.;
				for(l=1;l<=m;l++)if(Nlr[-l][r]!=(float)MISSVAL){ //compute pseudovalues
					Nlr[-l][r]=nJK*Nlr[0][r]-(nJK-1)*Nlr[-l][r];
					Nlr[m+1][r]+=Nlr[-l][r];   
				}
				Nlr[m+1][r]/=nJK;
				for(l=1;l<=m;l++)if(Glr[-l][r]!=(float)MISSVAL){ //compute variance of pseudoval
					Nlr[m+2][r]+=(Nlr[-l][r]-Nlr[m+1][r])*(Nlr[-l][r]-Nlr[m+1][r]);
				}
				Nlr[m+2][r]=(float)sqrt((Nlr[m+2][r]/(nJK-1))/nJK);
			}
			else Nlr[m+1][r]=Nlr[m+2][r]=(float)MISSVAL;
		}
		else Nlr[m+1][r]=Nlr[m+2][r]=(float)MISSVAL;
	}

	//Multilocus estimates for all loci but one, for pairwise Gst, Nst analysis
	if(JKest && pop1 && pop2){
		for(r=1;r<=3;r++)for(l=1;l<=m;l++){ //compute stat for all loci but locus l
			if(Glr[l][r]!=(float)MISSVAL) Glr[-l][r]=(float)((NumGlr[0][r]-NumGlr[l][r])/(DenGlr[0][r]-DenGlr[l][r]));
			else Glr[-l][r]=Glr[0][r];	
			if(Nlr[l][r]!=(float)MISSVAL) Nlr[-l][r]=(float)((NumNlr[0][r]-NumNlr[l][r])/(DenNlr[0][r]-DenNlr[l][r]));
			else Nlr[-l][r]=Nlr[0][r];	
		}

	}

	//Multilocus estimates for all loci but one, for pairwise Gij, Nij, Gst, Nst analyses
	if(JKest && pop1==0 && pop2==0){
		for(l=1;l<=m;l++){
			for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
				if(DivNlp1p2[l][p1][p2]!=(float)MISSVAL) Nlp1p2[-l][p1][p2]=1.0f-(DivNlp1p2[0][p1][p2]-DivNlp1p2[l][p1][p2])/(DivNrefp1p2[p1][p2]-DivNtl[l]);
				else Nlp1p2[-l][p1][p2]=Nlp1p2[0][p1][p2];
				if(DivGlp1p2[l][p1][p2]!=(float)MISSVAL) Glp1p2[-l][p1][p2]=1.0f-(DivGlp1p2[0][p1][p2]-DivGlp1p2[l][p1][p2])/(DivGrefp1p2[p1][p2]-DivGtl[l]);
				else Glp1p2[-l][p1][p2]=Glp1p2[0][p1][p2];

				if(DivNlp1p2[l][p1][p2]!=(float)MISSVAL) Nstlp1p2[-l][p1][p2]=1.0f-(0.5*DivNlp1p2[0][p1][p1]+0.5*DivNlp1p2[0][p2][p2]-0.5*DivNlp1p2[l][p1][p1]-0.5*DivNlp1p2[l][p2][p2])/(DivNlp1p2[0][p1][p2]-DivNlp1p2[l][p1][p2]);
				else Nstlp1p2[-l][p1][p2]=Nstlp1p2[0][p1][p2];
				if(DivGlp1p2[l][p1][p2]!=(float)MISSVAL) Gstlp1p2[-l][p1][p2]=1.0f-(0.5*DivGlp1p2[0][p1][p1]+0.5*DivGlp1p2[0][p2][p2]-0.5*DivGlp1p2[l][p1][p1]-0.5*DivGlp1p2[l][p2][p2])/(DivGlp1p2[0][p1][p2]-DivGlp1p2[l][p1][p2]);
				else Gstlp1p2[-l][p1][p2]=Gstlp1p2[0][p1][p2];
			}
		}

	}



	//for global analysis
	if(pop1==0 && pop2==0)for(r=1;r<=3;r++)for(l=0;l<=m+2;l++)for(S=1;S<=NS;S++)if(Stat[S]==9 || Stat[S]==10) FstatSlr[S][l][r]=Glr[l][r];
	if(pop1==0 && pop2==0)for(r=1;r<=3;r++)for(l=0;l<=m+2;l++)for(S=1;S<=NS;S++)if(Stat[S]==11 || Stat[S]==12) FstatSlr[S][l][r]=Nlr[l][r];

	//for pairwise Gst or Nst for one pair of pop
	if(JKest)linit=-m;
	else linit=0;
	if(m==1) linit=1;
	if(pop1 && pop2)for(l=linit;l<=m;l++)for(S=1;S<=NS;S++){
		if(Stat[S]==9) corrSlij[S][l][pop1][pop2]=Glr[l][3];				//Gst
		if(Stat[S]==11) corrSlij[S][l][pop1][pop2]=Nlr[l][3];				//Nst

		//weights of per locus estimates
		if(l>=0){
			if(DenGlr[l][3]==(float)MISSVAL) DenGlr[l][3]=0.;
			if(DenNlr[l][3]==(float)MISSVAL) DenNlr[l][3]=0.;
			if(Stat[S]==9) corrSlij[S][l][pop2][pop1]=(float)DenGlr[l][3];	
			if(Stat[S]==11) corrSlij[S][l][pop2][pop1]=(float)DenNlr[l][3];	
		}
	}

	//for pairwise Gij, Gst, Nij, Nst
	if(pop1==0 && pop2==0)for(S=1;S<=NS;S++){
		if(Stat[S]==9)for(l=linit;l<=m;l++)for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			corrSlij[S][l][p1][p2]=(float)Gstlp1p2[l][p1][p2];
		}

		if(Stat[S]==10)for(l=linit;l<=m;l++)for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			corrSlij[S][l][p1][p2]=(float)Glp1p2[l][p1][p2];
	/*		if(p1!=p2 && l==0) corrSlij[S][l][p2][p1]=DivGrefp1p2[p1][p2];
			if(p1!=p2 && l!=0) corrSlij[S][l][p2][p1]=DenGlr[l][3];			   */
		}

		if(Stat[S]==11)for(l=linit;l<=m;l++)for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			corrSlij[S][l][p1][p2]=(float)Nstlp1p2[l][p1][p2];
		}

		if(Stat[S]==12)for(l=linit;l<=m;l++)for(p1=1;p1<=Npop;p1++)for(p2=p1;p2<=Npop;p2++){
			corrSlij[S][l][p1][p2]=(float)Nlp1p2[l][p1][p2];
	/*		if(p1!=p2 && l==0) corrSlij[S][l][p2][p1]=DivNrefp1p2[p1][p2];
			if(p1!=p2 && l!=0) corrSlij[S][l][p2][p1]=DenNlr[l][3];			   */
		}
	}



	for(l=0;l<=m;l++) free_dmatrix(Plpa[l],0,Npop,0,Nallelel[l]);
	if(computeallelefreq){
		free_imatrix(Ncpl,0,Npop,0,m);
		free_ivector(Nip,0,Npop);
		for(l=0;l<=m;l++) free_imatrix(Nlpa[l],0,Npop,0,Nallelel[l]);
	}
	free_matrix(Glr,-m,m+2,0,4);
	free_matrix(Nlr,-m,m+2,0,4);
	free_dmatrix(NumGlr,0,m,0,4);
	free_dmatrix(DenGlr,0,m,0,4);
	free_dmatrix(NumNlr,0,m,0,4);
	free_dmatrix(DenNlr,0,m,0,4);
	free_d3tensor(DivNlp1p2,0,m,0,Npop,0,Npop);
	free_d3tensor(DivGlp1p2,0,m,0,Npop,0,Npop);
	free_dmatrix(DivNrefp1p2,0,Npop,0,Npop);
	free_dmatrix(DivGrefp1p2,0,Npop,0,Npop);
	free_d3tensor(Nlp1p2,-m,m,0,Npop,0,Npop);
	free_d3tensor(Glp1p2,-m,m,0,Npop,0,Npop);
	free_d3tensor(Nstlp1p2,-m,m,0,Npop,0,Npop);
	free_d3tensor(Gstlp1p2,-m,m,0,Npop,0,Npop);

}



/****************************************************************************/




void compute_corr_per_dist_class (int n,int m,int nc,double *maxc,int Ncat,int *cati,
		int StatType,int TypeComp,int FreqRef,int varcoef,double *xi,double *yi,double *zi,double **Mdij,
		int *sgi,float dijmin,float dijmax,float ***corrlij,float **corrlc,float ***Rll,int Rbtwloc,float ***V,float **R2pl,
		long *seed,int JKl)
{
	int i,j,l,l2,linit,k,K; 		/*counter for individuals (i,j), locus (l),categ*/ 
	int c,inbreedingonly;			/*class of distance interval (c)*/ 
	double **statlc;
	double dij,logd;
	double *SPvdlin,*SPvdlog,*Sdlin,*Sdlog,*SSdlin,*SSdlog,*Svlin,*Svlog,*SSvlin,*SSvlog,val;
	double Sw,Sw2,Swcorr,Swcorr2,**SSwcorrcl,**Sdevcl;
	int *npairslin,*npairslog,*np,**npbiscl;
	double Mv,Vv,Md,Vd,SSDvlin,SSDdlin,SPDvdlin,SSDvlog,SSDdlog,SPDvdlog,a,b,r2;
	int nJK; /*#valid loci for jackknifying*/
	double mJK,vJK,*JK;
	float ***rescorrlij;

	if(nc==-1){
		inbreedingonly=1;
		nc=0;
	}
	else inbreedingonly=0;

	np=ivector(0,nc+1);
	npbiscl=imatrix(0,nc+1,0,m);
	SSwcorrcl=dmatrix(0,nc+1,0,m);
	Sdevcl=dmatrix(0,nc+1,0,m);

	if(TypeComp==1) K=Ncat;
	else K=0;
	SPvdlin=dvector(0,K);
	SPvdlog=dvector(0,K);
	Sdlin=dvector(0,K);
	Sdlog=dvector(0,K);
	SSdlin=dvector(0,K);
	SSdlog=dvector(0,K);
	Svlin=dvector(0,K);
	SSvlin=dvector(0,K);
	Svlog=dvector(0,K);
	SSvlog=dvector(0,K);
	npairslin=ivector(0,K);
	npairslog=ivector(0,K);

	linit=0;
	if(JKl) linit=-m;
	if(m==1) linit=1;
   	statlc=dmatrix(linit,m+2,-20,nc+2);

	for(c=0;c<=nc+1;c++)for(l=0;l<=m;l++){
		npbiscl[c][l]=0;
		SSwcorrcl[c][l]=Sdevcl[c][l]=0.;
	}

	for(l=linit;l<=m;l++){
		/*initialize*/
		for(c=0;c<=nc+1;c++){
			np[c]=0;
			statlc[l][c]=0.;
		}
		for(k=0;k<=K;k++) {
			npairslin[k]=npairslog[k]=0;
			Sdlin[k]=Sdlog[k]=Svlin[k]=Svlog[k]=SSdlin[k]=SSdlog[k]=SSvlin[k]=SSvlog[k]=SPvdlin[k]=SPvdlog[k]=0.;
		}
		k=0;
		for(i=1;i<=n;i++){
			if(corrlij[l][i][i]!=(float)MISSVAL){	/*if there is no missing value...*/
				np[0]++;				/*add a pair of ind/pop within class 0*/
				statlc[l][0]+=corrlij[l][i][i];
			}
			if(!inbreedingonly)for(j=i+1;j<=n;j++) {
				if(TypeComp==1)if(cati[i]!=cati[j]) {corrlij[l][i][j]=(float)MISSVAL;continue;}
				if(TypeComp==2)if(cati[i]==cati[j]) {corrlij[l][i][j]=(float)MISSVAL;continue;}

				val=corrlij[l][i][j];
			/*if(l<0){
					val=(corrlij[0][i][j]*corrlij[0][j][i]-corrlij[-l][i][j]*corrlij[-l][j][i])/(corrlij[0][j][i]-corrlij[-l][j][i]);
					if(corrlij[0][i][j]==(float)MISSVAL) val=(float)MISSVAL;
					if(corrlij[-l][i][j]==(float)MISSVAL) val=corrlij[0][i][j];
				}
			   */
				if(val==(float)MISSVAL) continue;

				if(Mdij[0][0]==1.){
					if(Mdij[i][j]!=(float)MISSVAL) dij=Mdij[i][j];
					else continue;
				}
				else if(Mdij[0][0]==-1.) dij=acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
				else dij=sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
			
 				if(StatType==1 && Mdij[0][0]==0 && (sgi[i]==sgi[j])) {c=1; dij=-1.;}	/*intra-group class*/
				else{c=1; while(dij>maxc[c] && c<=nc) c++;}

				np[c]++;				/*add a pair of allele within class c*/
				statlc[l][c]+=val;

				/*compute actual var of relatedness within each dist class*/
				if(varcoef && ( (l==0 && m>=2) || (m>=3 && l>0) ) ){
					/*compute total weight over loci*/
					Sw=0.;
					for(l2=1;l2<=m;l2++) if(l2!=l && corrlij[l2][i][j]!=(float)MISSVAL) Sw+=corrlij[l2][j][i];
					Sw2=Swcorr=Swcorr2=0.;
					for(l2=1;l2<=m;l2++) if(l2!=l && corrlij[l2][i][j]!=(float)MISSVAL){
						Swcorr+=(corrlij[l2][j][i]/Sw)*corrlij[l2][i][j];
						Swcorr2+=(corrlij[l2][j][i]/Sw)*corrlij[l2][i][j]*(corrlij[l2][j][i]/Sw)*corrlij[l2][i][j];
						Sw2+=(corrlij[l2][j][i]/Sw)*(corrlij[l2][j][i]/Sw);
					}
					if(Sw2<1.){	 /*estimated var for l=0 and jackknife values when locus l is removed*/
						npbiscl[c][l]++;
						SSwcorrcl[c][l]+=Swcorr;
						Sdevcl[c][l]+=(Swcorr*Swcorr-Swcorr2)/(1.0f-Sw2);
						npbiscl[nc+1][l]++;
						SSwcorrcl[nc+1][l]+=Swcorr;
						Sdevcl[nc+1][l]+=(Swcorr*Swcorr-Swcorr2)/(1.0f-Sw2);
					}
					else{
						
					}
					
/*					Sw=0.;
					for(l2=1;l2<=m;l2++) if(l2!=l && corrlij[l2][i][j]!=(float)MISSVAL) Sw+=corrlij[l2][0][0];
					Sw2=Swcorr=Swcorr2=0.;
					for(l2=1;l2<=m;l2++) if(l2!=l && corrlij[l2][i][j]!=(float)MISSVAL){
						Swcorr+=(corrlij[l2][0][0]/Sw)*corrlij[l2][i][j];
						Swcorr2+=(corrlij[l2][0][0]/Sw)*corrlij[l2][i][j]*(corrlij[l2][0][0]/Sw)*corrlij[l2][i][j];
						Sw2+=(corrlij[l2][0][0]/Sw)*(corrlij[l2][0][0]/Sw);
					}
				
					if(Sw2<1. && JKl){	 //estimated var for l=0 and jackknife values when locus l is removed
						npbiscl[c][l]++;
						SSwcorrcl[c][l]+=Swcorr;
						Sdevcl[c][l]+=(Swcorr*Swcorr-Swcorr2)/(1.0f-Sw2);
						npbiscl[nc+1][l]++;
						SSwcorrcl[nc+1][l]+=Swcorr;
						Sdevcl[nc+1][l]+=(Swcorr*Swcorr-Swcorr2)/(1.0f-Sw2);
					}
*/ 
				}
				
				/*compute correlation/regression btw pairwise autocorrelation coef & dist*/
				if(dij>-1. && dij>=dijmin && dij<=dijmax){
					if(TypeComp==1) k=cati[i];
					/*lin reg*/
					npairslin[k]++;
					Svlin[k]+=val;
					SSvlin[k]+=val*val;
					Sdlin[k]+=dij;
					SSdlin[k]+=dij*dij;
					SPvdlin[k]+=dij*val;
					/*log reg*/
					if(dij){
						logd=log(dij);
						npairslog[k]++;
						Svlog[k]+=val;
						SSvlog[k]+=val*val;
						Sdlog[k]+=logd;
						SSdlog[k]+=logd*logd;
						SPvdlog[k]+=logd*val;
					}
				}
			} //end loop j
		} //end loop i


		/*compute coef for each dist class and average coef over all dist classes*/
		//for inbreeding coef
		if(np[0]) statlc[l][0]=statlc[l][0]/np[0];
		else statlc[l][0]=(float)MISSVAL;
		//for other dist classes
		if(!inbreedingonly){
			if(nc>0){
				statlc[l][nc+1]=0.0;
				np[nc+1]=0;
				for(c=1;c<=nc;c++){
					if(np[c]){
						statlc[l][nc+1]+=statlc[l][c];
						np[nc+1]+=np[c];

						statlc[l][c]=statlc[l][c]/np[c];
					}
					else statlc[l][c]=(float)MISSVAL;
				}
				if(np[nc+1]) statlc[l][nc+1]=statlc[l][nc+1]/np[nc+1];
				else statlc[l][nc+1]=(float)MISSVAL;
			}
			else{
				if(np[1]) statlc[l][1]=statlc[l][1]/np[1];
				else statlc[l][1]=(float)MISSVAL;
			}
		}


		/*compute correlation/regression btw pairwise coef & dist*/
		if(TypeComp==1)for(k=1;k<=K;k++){
			npairslin[0]+=npairslin[k]; npairslog[0]+=npairslog[k]; Svlin[0]+=Svlin[k]; SSvlin[0]+=SSvlin[k]; Svlog[0]+=Svlog[k]; SSvlog[0]+=SSvlog[k]; Sdlin[0]+=Sdlin[k]; SSdlin[0]+=SSdlin[k]; Sdlog[0]+=Sdlog[k]; SSdlog[0]+=SSdlog[k]; SPvdlin[0]+=SPvdlin[k]; SPvdlog[0]+=SPvdlog[k];
		}	
		if(nc>0){
			/*reg with lin d*/
			statlc[l][-17]=npairslin[0];		   /*# pairs*/
			if(npairslin[0]){
				Mv=Svlin[0]/npairslin[0];
				SSDvlin=SSvlin[0]-Svlin[0]*Svlin[0]/npairslin[0];
				if(SSDvlin<0.000000001) SSDvlin=0.;/*when val is constant, Vv might however take some value due to error of round*/
				Vv=SSDvlin/(npairslin[0]-1.);

				Md=Sdlin[0]/npairslin[0];
				SSDdlin=SSdlin[0]-Sdlin[0]*Sdlin[0]/npairslin[0];
				if(SSDdlin<0.000000001) SSDdlin=0.;/*when val is constant, Vd might however take some value due to error of round*/
				Vd=SSDdlin/(npairslin[0]-1.);

				statlc[l][-16]=Md;	
				statlc[l][-15]=Vd;	
				statlc[l][-14]=Mv;	
				statlc[l][-13]=Vv;	

				if(SSDdlin && SSDvlin){
					SPDvdlin=SPvdlin[0]-Svlin[0]*Sdlin[0]/npairslin[0];
					r2=(SPDvdlin*SPDvdlin)/(SSDdlin*SSDvlin);
					b=SPDvdlin/SSDdlin;
					a=Mv-b*Md;

					statlc[l][-2]=statlc[l][-20]=b;	
					statlc[l][-19]=a;	
					statlc[l][-18]=r2;	
				}
				else statlc[l][-2]=statlc[l][-20]=statlc[l][-19]=statlc[l][-18]=(float)MISSVAL;

			}
			else statlc[l][-2]=statlc[l][-20]=statlc[l][-19]=statlc[l][-18]=statlc[l][-16]=statlc[l][-15]=statlc[l][-14]=statlc[l][-13]=(float)MISSVAL;	

			/*reg with log d*/
			statlc[l][-7]=npairslog[0];		   /*# pairs*/
			if(npairslog[0]){
				Mv=Svlog[0]/npairslog[0];
				SSDvlog=SSvlog[0]-Svlog[0]*Svlog[0]/npairslog[0];
				if(SSDvlog<0.000000001) SSDvlog=0.;/*when val is constant, Vv might however take some value due to error of round*/
				Vv=SSDvlog/(npairslog[0]-1.);

				Md=Sdlog[0]/npairslog[0];
				SSDdlog=SSdlog[0]-Sdlog[0]*Sdlog[0]/npairslog[0];
				if(SSDdlog<0.000000001) SSDdlog=0.;/*when val is constant, Vd might however take some value due to error of round*/
				Vd=SSDdlog/(npairslog[0]-1.);

				statlc[l][-6]=Md;	
				statlc[l][-5]=Vd;	
				statlc[l][-4]=Mv;	
				statlc[l][-3]=Vv;	

				if(SSDdlog && SSDvlog){
					SPDvdlog=SPvdlog[0]-Svlog[0]*Sdlog[0]/npairslog[0];
					r2=(SPDvdlog*SPDvdlog)/(SSDdlog*SSDvlog);
					b=SPDvdlog/SSDdlog;
					a=Mv-b*Md;

					statlc[l][-1]=statlc[l][-10]=b;	
					statlc[l][-9]=a;	
					statlc[l][-8]=r2;	
				}
				else statlc[l][-1]=statlc[l][-10]=statlc[l][-9]=statlc[l][-8]=(float)MISSVAL;

			}
			else statlc[l][-1]=statlc[l][-10]=statlc[l][-9]=statlc[l][-8]=statlc[l][-6]=statlc[l][-5]=statlc[l][-4]=statlc[l][-3]=(float)MISSVAL;	

			if(TypeComp==1 && FreqRef==1){
  				statlc[l][-17]=npairslin[0];
				statlc[l][-7]=npairslog[0];
				if(npairslin[0] || npairslog[0]){
					SSDvlin=SSDdlin=SPDvdlin=SSDvlog=SSDdlog=SPDvdlog=0.;
					for(k=1;k<=Ncat;k++){
						if(npairslin[k]){
							SSDvlin+=SSvlin[k]-Svlin[k]*Svlin[k]/npairslin[k];
							SSDdlin+=SSdlin[k]-Sdlin[k]*Sdlin[k]/npairslin[k];
							SPDvdlin+=SPvdlin[k]-Svlin[k]*Sdlin[k]/npairslin[k];
						}
						if(npairslog[k]){
							SSDvlog+=SSvlog[k]-Svlog[k]*Svlog[k]/npairslog[k];
							SSDdlog+=SSdlog[k]-Sdlog[k]*Sdlog[k]/npairslog[k];
							SPDvdlog+=SPvdlog[k]-Svlog[k]*Sdlog[k]/npairslog[k];
						}
					}
					/*reg with lin dist*/
					if(npairslin[0]){
						Mv=Svlin[0]/npairslin[0];
						Md=Sdlin[0]/npairslin[0];
						Vv=SSDvlin/(npairslin[0]-Ncat);
						Vd=SSDdlin/(npairslin[0]-Ncat);
					}
					else Mv=Md=Vv=Vd=(float)MISSVAL;
					if(SSDvlin<0.000000001) SSDvlin=0.;
					if(SSDdlin<0.000000001) SSDdlin=0.;
					if(SSDdlin && SSDvlin){
						r2=(SPDvdlin*SPDvdlin)/(SSDdlin*SSDvlin);
						b=SPDvdlin/SSDdlin;
						a=Mv-b*Md;
					}
					else r2=a=b=(float)MISSVAL;

					statlc[l][-2]=statlc[l][-20]=b;	
					statlc[l][-19]=a;	
					statlc[l][-18]=r2;
					statlc[l][-16]=Md;	
					statlc[l][-15]=Vd;	
					statlc[l][-14]=Mv;	
					statlc[l][-13]=Vv;	

					/*reg with log dist*/
					if(npairslog[0]){
						Mv=Svlog[0]/npairslog[0];
						Md=Sdlog[0]/npairslog[0];
						Vv=SSDvlog/(npairslin[0]-Ncat);
						Vd=SSDdlog/(npairslin[0]-Ncat);
					}
					else Mv=Md=Vv=Vd=(float)MISSVAL;
					if(SSDvlog<0.000000001) SSDvlog=0.;
					if(SSDdlog<0.000000001) SSDdlog=0.;
					if(SSDdlog && SSDvlog){
						r2=(SPDvdlog*SPDvdlog)/(SSDdlog*SSDvlog);
						b=SPDvdlog/SSDdlog;
						a=Mv-b*Md;
					}
					else r2=a=b=(float)MISSVAL;

					statlc[l][-1]=statlc[l][-10]=b;	
					statlc[l][-9]=a;	
					statlc[l][-8]=r2;
					statlc[l][-6]=Md;	
					statlc[l][-5]=Vd;	
					statlc[l][-4]=Mv;	
					statlc[l][-3]=Vv;	
			
				} /*end of if(npairslin[0] || npairslog[0])*/ 
				else statlc[l][-2]=statlc[l][-20]=statlc[l][-19]=statlc[l][-18]=statlc[l][-16]=statlc[l][-15]=statlc[l][-14]=statlc[l][-13]=statlc[l][-1]=statlc[l][-10]=statlc[l][-9]=statlc[l][-8]=statlc[l][-6]=statlc[l][-5]=statlc[l][-4]=statlc[l][-3]=(float)MISSVAL;	

			} /*end of if(TypeComp=1)*/
			statlc[l][-12]=statlc[l][-11]=(float)MISSVAL; /*memory not used*/
		}/*end of if(nc>0)*/
		else for(c=-20;c<=-1;c++) statlc[l][c]=(float)MISSVAL;	

	}/*end of l counter*/

	/*compute actual var of relatedness and standard error by jacknife*/
	if(varcoef && m>=2){
		for(c=0;c<=nc+1;c++){
			if(npbiscl[c][0]) corrlc[m+3][c]=(float)((Sdevcl[c][0]/npbiscl[c][0])-(SSwcorrcl[c][0]/npbiscl[c][0])*(SSwcorrcl[c][0]/npbiscl[c][0]));
			else corrlc[m+3][c]=(float)MISSVAL;
	
			if(m>=3){/*jackknife*/
				JK=dvector(0,m);
				if(npbiscl[c][0]){
					nJK=0;
					for(l=1;l<=m;l++) if(npbiscl[c][l]) nJK++;
					
					if(nJK>=2){
						mJK=vJK=0;
						for(l=1;l<=m;l++)if(npbiscl[c][l]){
							JK[l]=(Sdevcl[c][l]/npbiscl[c][l])-(SSwcorrcl[c][l]/npbiscl[c][l])*(SSwcorrcl[c][l]/npbiscl[c][l]);
							JK[l]=nJK*corrlc[m+3][c]-(nJK-1)*JK[l];/*pseudoval*/
							mJK+=JK[l]/nJK;
						}
						for(l=1;l<=m;l++)if(npbiscl[c][l]) vJK+=(JK[l]-mJK)*(JK[l]-mJK)/(nJK-1);
						corrlc[m+4][c]=(float)sqrt(vJK/nJK);
					}
					else corrlc[m+4][c]=(float)MISSVAL;
				}
				else corrlc[m+4][c]=(float)MISSVAL;
				free_dvector(JK,0,m);
			}
		}
	}

	/*transfer results*/
	for(l=1;l<=m;l++)for(c=-20;c<=nc+1;c++) corrlc[l][c]=(float)statlc[l][c];
	if(m>1)for(c=-20;c<=nc+1;c++) corrlc[0][c]=(float)statlc[0][c]; 
	if(JKl) for(l=1;l<=m;l++)for(c=-20;c<=nc+1;c++) corrlc[-l][c]=(float)statlc[-l][c];

	/*Jackknife over loci*/
	if(JKl && m>1){
		for(c=-20;c<=nc+1;c++){
			if(statlc[0][c]!=(float)MISSVAL){
				nJK=0;
				for(l=1;l<=m;l++) if(statlc[-l][c]!=(float)MISSVAL) nJK++;
				if(nJK>1){
					statlc[m+1][c]=0.;
					for(l=1;l<=m;l++){
						if(statlc[-l][c]!=(float)MISSVAL){ /*compute pseudovalues*/
							statlc[-l][c]=nJK*statlc[0][c]-(nJK-1)*statlc[-l][c];
							statlc[m+1][c]+=statlc[-l][c];
						}
					}
					statlc[m+1][c]/=nJK;  /*mean JK estimate*/
					statlc[m+2][c]=0.;
					for(l=1;l<=m;l++)if(statlc[-l][c]!=(float)MISSVAL){ /*compute variance JK estimate*/
						statlc[m+2][c]+=(statlc[-l][c]-statlc[m+1][c])*(statlc[-l][c]-statlc[m+1][c]);
					}
					statlc[m+2][c]=sqrt((statlc[m+2][c]/(nJK-1))/nJK); /*standard error JK estimate*/
				}
				else statlc[m+1][c]=statlc[m+2][c]=(float)MISSVAL;
			}
			else statlc[m+1][c]=statlc[m+2][c]=(float)MISSVAL;

			corrlc[m+1][c]=(float)statlc[m+1][c];
			corrlc[m+2][c]=(float)statlc[m+2][c];
		}
	}


	//compute inter-locus correlation coefficients
	if(Rbtwloc==1 && m>1){
		inter_locus_corr(n,m,corrlij,Rll[0],V[0],R2pl[0],seed);

		//inter-locus corr of residuals
		rescorrlij=f3tensor(0,m,0,n,0,n);
		for(i=1;i<n;i++)for(j=i+1;j<=n;j++){
			if(TypeComp==1)if(cati[i]!=cati[j]) continue;
			if(TypeComp==2)if(cati[i]==cati[j]) continue;

			//find dist class
			if(Mdij[0][0]==1.){
				if(Mdij[i][j]!=(float)MISSVAL) dij=Mdij[i][j];
				else continue;
			}
			else if(Mdij[0][0]==-1.) dij=acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
			else dij=sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
		
 			if(StatType==1 && Mdij[0][0]==0 && (sgi[i]==sgi[j])) {c=1; dij=-1.;}	/*intra-group class*/
			else{c=1; while(dij>maxc[c] && c<=nc) c++;}

			//get residuals
			for(l=1;l<=m;l++){
				if(corrlij[l][i][j]!=(float)MISSVAL)	rescorrlij[l][i][j]=corrlij[l][i][j]-corrlc[l][c];
	//			if(corrlij[l][i][j]!=(float)MISSVAL)	rescorrlij[l][i][j]=corrlij[l][i][j]-corrlc[0][c];
				else rescorrlij[l][i][j]=(float)MISSVAL;
				rescorrlij[l][j][i]=corrlij[l][j][i];//locus specific weights
			}
		}
		inter_locus_corr(n,m,rescorrlij,Rll[1],V[1],R2pl[1],seed);
		free_f3tensor(rescorrlij,0,m,0,n,0,n);
	}




	free_dmatrix(statlc,linit,m+2,-20,nc+2);
	free_ivector(np,0,nc+1);
	free_imatrix(npbiscl,0,nc+1,0,m);
	free_dmatrix(SSwcorrcl,0,nc+1,0,m);
	free_dmatrix(Sdevcl,0,nc+1,0,m);
	free_dvector(SPvdlin,0,K);
	free_dvector(SPvdlog,0,K);
	free_dvector(Sdlin,0,K);
	free_dvector(Sdlog,0,K);
	free_dvector(SSdlin,0,K);
	free_dvector(SSdlog,0,K);
	free_dvector(Svlin,0,K);
	free_dvector(SSvlin,0,K);
	free_dvector(Svlog,0,K);
	free_dvector(SSvlog,0,K);
	free_ivector(npairslin,0,K);
	free_ivector(npairslog,0,K);

	if(m>1) corrlc[0][-11]=0.; //contains sigma estim if asked 
	else corrlc[1][-11]=0.;


}	/*end procedure compute_corr_per_dist_class*/

/*************************************************************************************/
void inter_locus_corr(int n,int m,float ***corrlij,float **Rll,float **V,float *R2pl,long *seed)
{
	int i,j,l,l1,l2,Ncomp,Nlcomp[MMAX],nJK;
	double SSl1,SSl2,Sl1,Sl2,SP,Vl1,Vl2,Cov,Cor;
	double *Sp,***Sl,*S,*T,*SCEl,*SCEp,*SCEt,*SCEe,*SCEr,*CMl,*CMp,*CMe,*CMr,*Ve,*Vl,*Vp1,*Vp2,*Rl1,*Rl2;

	int *shufl,ljk,nrand,Nrand,NRand;
	double mean,mean1,mean2,sw1,sw2;

	
	for(l=0;l<=m;l++) {
		Nlcomp[l]=0;
		Rll[0][l]=0.;
	}

	for(l1=1;l1<m;l1++){
		for(l2=l1+1;l2<=m;l2++){
			SP=Sl1=Sl2=SSl1=SSl2=0.;
			Ncomp=0;
			for(i=1;i<n;i++)for(j=i+1;j<=n;j++)if(corrlij[l1][i][j]!=(float)MISSVAL && corrlij[l2][i][j]!=(float)MISSVAL){
				Ncomp++;
				Sl1+=corrlij[l1][i][j];
				Sl2+=corrlij[l2][i][j];
				SSl1+=corrlij[l1][i][j]*corrlij[l1][i][j];
				SSl2+=corrlij[l2][i][j]*corrlij[l2][i][j];
				SP+=corrlij[l1][i][j]*corrlij[l2][i][j];
			}
			if(Ncomp){
				Vl1=SSl1/Ncomp-(Sl1*Sl1/(Ncomp*Ncomp));
				Vl2=SSl2/Ncomp-(Sl2*Sl2/(Ncomp*Ncomp));

				Cov=SP/Ncomp-(Sl1*Sl2/(Ncomp*Ncomp));
				if(Vl1>0 && Vl2>0){
					Nlcomp[0]++;
					Rll[l1][l2]=(float)(Cov/(sqrt(Vl1)*sqrt(Vl2)));		//correlation between loci 1 & 2				
					Rll[0][0]+=Rll[l1][l2] /* *Rll[l1][l2] */ ;	//mean correlation between loci
					//jackknife average
					for(l=1;l<=m;l++){
						if(l==l1 || l==l2) continue;
						else{
							Nlcomp[l]++;
							Rll[0][l]+=Rll[l1][l2]; //mean correlation between loci excepted locus l
						}
					}
				}
				else Rll[l1][l2]=(float)MISSVAL;
			}
			else Rll[l1][l2]=(float)MISSVAL;
		}
	}
	for(l=0;l<=m;l++){
		if(Nlcomp[l]) Rll[0][l]/=Nlcomp[l];
		else Rll[0][l]=(float)MISSVAL;
	}
	//jackknife
	nJK=0;
	for(l=1;l<=m;l++) if(Rll[0][l]!=(float)MISSVAL) nJK++;
	if(nJK>1){
		Rll[1][0]=0.;
		for(l=1;l<=m;l++){
			if(Rll[0][l]!=(float)MISSVAL){ /*compute pseudovalues*/
				Rll[0][l]=nJK*Rll[0][0]-(nJK-1)*Rll[0][l];
				Rll[1][0]+=Rll[0][l];
			}
		}
		Rll[1][0]/=nJK;  /*mean JK estimate*/
		Rll[2][0]=0.;
		for(l=1;l<=m;l++)if(Rll[0][l]!=(float)MISSVAL){ /*compute variance JK estimate*/
			Rll[2][0]+=(Rll[0][l]-Rll[1][0])*(Rll[0][l]-Rll[1][0]);
		}
		Rll[2][0]=(float)sqrt((Rll[2][0]/(nJK-1))/nJK); /*standard error JK estimate*/
	}
	else Rll[1][0]=Rll[2][0]=(float)MISSVAL;





	//Inter-locus correlation based on 2 pools of loci
	shufl=ivector(0,m);
	for(ljk=1;ljk<=m;ljk++){
		for(l=1;l<ljk;l++) shufl[l]=l;
		for(l=ljk+1;l<=m;l++) shufl[l-1]=l;
		
		NRand=0;
		R2pl[-ljk]=0.0f;
		if(m<20) Nrand=10;
		else Nrand=1;
		for(nrand=1;nrand<=Nrand;nrand++){
			resample_shuffle(shufl,1,m-1,seed);
			SP=Sl1=Sl2=SSl1=SSl2=0.;
			Ncomp=0;
			for(i=1;i<n;i++)for(j=i+1;j<=n;j++){
				//average over 1st pool of loci
				mean1=sw1=0.;
				for(l=1;l<=(int)((m-1)/2);l++){
					l1=shufl[l];
					if(corrlij[l1][i][j]!=(float)MISSVAL){
						mean1+=corrlij[l1][i][j];
						sw1+=corrlij[l1][j][i];
					}
				}
				if(sw1) mean1/=sw1;
				else mean1=(float)MISSVAL;
				//average over 2nd pool of loci
				mean2=sw2=0.;
				for(l=1+(int)((m-1)/2);l<=(m-1);l++){
					l2=shufl[l];
					if(corrlij[l2][i][j]!=(float)MISSVAL){
						mean2+=corrlij[l2][i][j];
						sw2+=corrlij[l2][j][i];
					}
				}
				if(sw2) mean2/=sw2;
				else mean2=(float)MISSVAL;

				if(mean1!=(float)MISSVAL && mean2!=(float)MISSVAL){
					Ncomp++;
					Sl1+=mean1;
					Sl2+=mean2;
					SSl1+=mean1*mean1;
					SSl2+=mean2*mean2;
					SP+=mean1*mean2;
				}
			}
			if(Ncomp){
				Vl1=SSl1/Ncomp-(Sl1*Sl1/(Ncomp*Ncomp));
				Vl2=SSl2/Ncomp-(Sl2*Sl2/(Ncomp*Ncomp));

				Cov=SP/Ncomp-(Sl1*Sl2/(Ncomp*Ncomp));
				if(Vl1>0 && Vl2>0){
					Nlcomp[0]++;
					Cor=Cov/(sqrt(Vl1)*sqrt(Vl2));		//correlation				
				}
				else Cor=(float)MISSVAL;
			}
			else Cor=(float)MISSVAL;

			if(Cor!=(float)MISSVAL){NRand++; R2pl[-ljk]+=(float)Cor;}
		}
		if(NRand) R2pl[-ljk]/=NRand;
		else R2pl[-ljk]=(float)MISSVAL;

	}
	//mean over all
	Ncomp=0;
	R2pl[0]=0.;
	for(l=1;l<=m;l++) if(R2pl[-l]!=(float)MISSVAL){ Ncomp++; R2pl[0]+=R2pl[-l];}
	R2pl[0]/=Ncomp;
	//jackknife
	nJK=0;
	for(l=1;l<=m;l++) if(R2pl[-l]!=(float)MISSVAL) nJK++;
	if(nJK>1){
		mean=0.;
		for(l=1;l<=m;l++){
			if(R2pl[-l]!=(float)MISSVAL){ /*compute pseudovalues*/
				R2pl[-l]=nJK*R2pl[0]-(nJK-1)*R2pl[-l];
				mean+=R2pl[-l];
			}
		}
		R2pl[1]=(float)mean/nJK;  /*mean JK estimate*/
		R2pl[2]=0.;
		for(l=1;l<=m;l++)if(R2pl[-l]!=(float)MISSVAL){ /*compute variance JK estimate*/
			R2pl[2]+=(R2pl[-l]-(float)mean)*(R2pl[-l]-(float)mean);
		}
		R2pl[2]=(float)sqrt((R2pl[2]/(nJK-1))/nJK); /*standard error JK estimate*/
	}
	else R2pl[1]=R2pl[2]=(float)MISSVAL;

	free_ivector(shufl,0,m);





		






	//ANOVA2 à 1obs
/**/	Sl=d3tensor(0,n,0,n,-m,m);
	Sp=dvector(-m,m);
	S=dvector(-m,m);
	T=dvector(-m,m);
	SCEl=dvector(-m,m);
	SCEp=dvector(-m,m);
	SCEt=dvector(-m,m);
	SCEe=dvector(-m,m);
	SCEr=dvector(-m,m);
	CMl=dvector(-m,m);
	CMp=dvector(-m,m);
	CMe=dvector(-m,m);
	CMr=dvector(-m,m);
	Ve=dvector(-m,m);
	Vl=dvector(-m,m);
	Vp1=dvector(-m,m);
	Vp2=dvector(-m,m);
	Rl1=dvector(-m,m);
	Rl2=dvector(-m,m);


	//sum over all pairs within locus
	for(l=1;l<=m;l++){
		Sp[l]=0.;
		for(i=1;i<n;i++)for(j=i+1;j<=n;j++)if(corrlij[l][i][j]!=(float)MISSVAL) Sp[l]+=corrlij[l][i][j];
	}
	Sp[0]=0.;
	for(l=1;l<=m;l++) Sp[0]+=Sp[l];
	for(l=1;l<=m;l++) Sp[-l]=Sp[0]-Sp[l]; //for jackknife
	//sum over all loci within pairs
	for(i=1;i<n;i++)for(j=i+1;j<=n;j++){
		Sl[i][j][0]=0.;
		for(l=1;l<=m;l++)if(corrlij[l][i][j]!=(float)MISSVAL) Sl[i][j][0]+=corrlij[l][i][j];
		for(l=1;l<=m;l++) Sl[i][j][-l]=Sl[i][j][0]-corrlij[l][i][j]; //for jackknife
	}
	//overall sum
	S[0]=0.;
	for(l=1;l<=m;l++) S[0]+=Sp[l];
	for(l=1;l<=m;l++) S[-l]=S[0]-Sp[l]; //for jackknife
	//SSG	
	for(l=1;l<=m;l++){
		T[l]=0.;
		for(i=1;i<n;i++)for(j=i+1;j<=n;j++) if(corrlij[l][i][j]!=(float)MISSVAL) T[l]+=corrlij[l][i][j]*corrlij[l][i][j];
	}
	T[0]=0.;
	for(l=1;l<=m;l++) T[0]+=T[l];
	for(l=1;l<=m;l++) T[-l]=T[0]-T[l]; //for jackknife

	//SCEl
	SCEl[0]=0.;
	for(l=1;l<=m;l++)if(Sp[l]!=(float)MISSVAL) SCEl[0]+=Sp[l]*Sp[l];
	for(l=1;l<=m;l++){if(Sp[l]!=(float)MISSVAL) SCEl[-l]=SCEl[0]-Sp[l]*Sp[l]; else SCEl[-l]=SCEl[0];} //for jackknife
	SCEl[0]=SCEl[0]/(n*(n-1.)/2.) - S[0]*S[0]/(m*n*(n-1.)/2.);
	for(l=1;l<=m;l++) SCEl[-l]=SCEl[-l]/(n*(n-1.)/2.) - S[-l]*S[-l]/((m-1)*n*(n-1.)/2.);
	//SCEp
	for(l=0;l<=m;l++) SCEp[-l]=0.;
	for(i=1;i<n;i++)for(j=i+1;j<=n;j++) if(Sl[i][j][0]!=(float)MISSVAL) SCEp[0]+=Sl[i][j][0]*Sl[i][j][0];
	for(l=1;l<=m;l++) for(i=1;i<n;i++)for(j=i+1;j<=n;j++) if(Sl[i][j][-l]!=(float)MISSVAL) SCEp[-l]+=Sl[i][j][-l]*Sl[i][j][-l]; //for jackknife
	SCEp[0]=SCEp[0]/m - S[0]*S[0]/(m*n*(n-1.)/2.);
	for(l=1;l<=m;l++) SCEp[-l]=SCEp[-l]/(m-1.) - S[-l]*S[-l]/((m-1.)*n*(n-1.)/2.);
	//SCEt
	SCEt[0]=T[0] - S[0]*S[0]/(m*n*(n-1.)/2.);
	for(l=1;l<=m;l++) SCEt[-l]=T[-l] - S[-l]*S[-l]/((m-1.)*n*(n-1.)/2.);
	//SCEe
	SCEe[0]=SCEt[0]-SCEl[0]-SCEp[0];
	CMl[0]=SCEl[0]/(m-1);
	CMp[0]=SCEp[0]/(n*(n-1.)/2.-1.);
	CMe[0]=SCEe[0]/((m-1) * (n*(n-1.)/2.-1.));
	Vp2[0]=(CMp[0]-CMe[0])/m;
	Vl[0]=(CMl[0]-CMe[0])/(n*(n-1.)/2.);
	Ve[0]=CMe[0];
	//intra-class correl (ICC(C,1) for ANOVA2 1 obs)
	Rl2[0]=Vp2[0]/(Vp2[0]+Ve[0]);
	for(l=1;l<=m;l++){
		SCEe[-l]=SCEt[-l]-SCEl[-l]-SCEp[-l];
		CMl[-l]=SCEl[-l]/((m-1)-1);
		CMp[-l]=SCEp[-l]/(n*(n-1.)/2.-1.);
		CMe[-l]=SCEe[-l]/(((m-1)-1) * (n*(n-1.)/2.-1.));
		Vp2[-l]=(CMp[-l]-CMe[-l])/(m-1);
		Vl[-l]=(CMl[-l]-CMe[-l])/(n*(n-1.)/2.);
		Ve[-l]=CMe[-l];
		//intra-class correl (ICC(C,1) for ANOVA2 1 obs)
		Rl2[-l]=Vp2[-l]/(Vp2[-l]+Ve[-l]);
	}


	//ANOVA1 à m obs
	SCEr[0]=SCEt[0]-SCEp[0];
	CMr[0]=SCEr[0]/((m-1)*n*(n-1.)/2.);
	Vp1[0]=(CMp[0]-CMr[0])/m;
	//intra-class correl (ICC(1) for ANOVA1 m obs)
	Rl1[0]=Vp1[0]/(Vp1[0]+CMr[0]);
	for(l=1;l<=m;l++){
		SCEr[-l]=SCEt[-l]-SCEp[-l];
		CMr[-l]=SCEr[-l]/(((m-1)-1)*n*(n-1.)/2.);
		Vp1[-l]=(CMp[-l]-CMr[-l])/(m-1);
		//intra-class correl (ICC(1) for ANOVA1 m obs)
		Rl1[-l]=Vp1[-l]/(Vp1[-l]+CMr[-l]);
	}

	//jackknife
	nJK=0;
//	for(l=1;l<=m;l++) if(Rll[0][l]!=(float)MISSVAL) nJK++;
	nJK=m;

	if(nJK>1){
		Vp1[1]=Vp2[1]=Rl1[1]=Rl2[1]=0.;
		for(l=1;l<=m;l++){
			if(Vp1[l]!=(float)MISSVAL){ //compute pseudovalues
				Vp1[-l]=nJK*Vp1[0]-(nJK-1)*Vp1[-l];
				Vp1[1]+=Vp1[-l];
				Vp2[-l]=nJK*Vp2[0]-(nJK-1)*Vp2[-l];
				Vp2[1]+=Vp2[-l];
				Rl1[-l]=nJK*Rl1[0]-(nJK-1)*Rl1[-l];
				Rl1[1]+=Rl1[-l];
				Rl2[-l]=nJK*Rl2[0]-(nJK-1)*Rl2[-l];
				Rl2[1]+=Rl2[-l];
			}
		}
		Vp1[1]/=nJK;  //mean JK estimate
		Vp2[1]/=nJK;
		Rl1[1]/=nJK;
		Rl2[1]/=nJK;

		Vp1[2]=Vp2[2]=Rl1[2]=Rl2[2]=0.;
		for(l=1;l<=m;l++)if(Vp1[l]!=(float)MISSVAL){ //compute variance JK estimate
			Vp1[2]+=(Vp1[-l]-Vp1[1])*(Vp1[-l]-Vp1[1]);
			Vp2[2]+=(Vp2[-l]-Vp2[1])*(Vp2[-l]-Vp2[1]);
			Rl1[2]+=(Rl1[-l]-Rl1[1])*(Rl1[-l]-Rl1[1]);
			Rl2[2]+=(Rl2[-l]-Rl2[1])*(Rl2[-l]-Rl2[1]);
		}
		Vp1[2]=sqrt((Vp1[2]/(nJK-1))/nJK); //standard error JK estimate
		Vp2[2]=sqrt((Vp2[2]/(nJK-1))/nJK);
		Rl1[2]=sqrt((Rl1[2]/(nJK-1))/nJK);
		Rl2[2]=sqrt((Rl2[2]/(nJK-1))/nJK);
	}
	else Vp1[1]=Vp1[2]=Vp2[1]=Vp2[2]=Rl1[1]=Rl1[2]=Rl2[1]=Rl2[2]=(float)MISSVAL;

	V[0][0]=(float)Vp1[0]; V[1][0]=(float)Rl1[0]; V[2][0]=(float)Vp2[0]; V[3][0]=(float)Vl[0]; V[4][0]=(float)Rl2[0]; // estimates
	V[0][1]=(float)Vp1[1]; V[1][1]=(float)Rl1[1]; V[2][1]=(float)Vp2[1]; V[3][1]=(float)Vl[1]; V[4][1]=(float)Rl2[1]; //mean jackknife estimates
	V[0][2]=(float)Vp1[2]; V[1][2]=(float)Rl1[2]; V[2][2]=(float)Vp2[2]; V[3][2]=(float)Vl[2]; V[4][2]=(float)Rl2[2]; //jackknife SE

	free_d3tensor(Sl,0,n,0,n,-m,m);
	free_dvector(Sp,-m,m);
	free_dvector(S,-m,m);
	free_dvector(T,-m,m);
	free_dvector(SCEl,-m,m);
	free_dvector(SCEp,-m,m);
	free_dvector(SCEt,-m,m);
	free_dvector(SCEe,-m,m);
	free_dvector(SCEr,-m,m);
	free_dvector(CMl,-m,m);
	free_dvector(CMp,-m,m);
	free_dvector(CMe,-m,m);
	free_dvector(CMr,-m,m);
	free_dvector(Ve,-m,m);
	free_dvector(Vl,-m,m);
	free_dvector(Vp1,-m,m);
	free_dvector(Vp2,-m,m);
	free_dvector(Rl1,-m,m);
	free_dvector(Rl2,-m,m);
/**/


}

/*************************************************************************************/

/*************************************************************************************/

void estimate_sigma_2D_kinship (int n,int m,double *xi,double *yi,double *zi,double **Mdij,
		int *sgi,float ***corrlij,float **corrlc,int ploidy,int Stat,int JKest,float density,float dwidth)
{
	int i,j,l,cneighb,niter,backiter,linit,lterm; 		
	double dij,logd;
	double SPvdlog,Sdlog,SSdlog,Svlog,SSvlog,val;
	int npairslog;
	double SSDvlog,SSDdlog,SPDvdlog,b,Nb[101];
	double dijmin, sigmaest[101],meansigma,meanNb,sigmasq,minNb,maxNb,minsigma,maxsigma;
	int converge=0,nJK;

	printf("\n\nEstimating dispersal distance for a 2D population at drift-dispersal equilibrium assuming an effective density = %G",density);
	printf("\nIter\t#pairs\tNb\tsigma");

	linit=lterm=0;
	if(JKest) linit=-m;
	if(m==1) linit=lterm=1;


	cneighb=1;


	for(l=linit;l<=lterm;l++){

		dijmin=sigmaest[0]=0.;
		niter=0;

		do{
			npairslog=0;
			Sdlog=Svlog=SSdlog=SSvlog=SPvdlog=0.;
			niter++;

			for(i=1;i<=n;i++){

				for(j=i+1;j<=n;j++)if(corrlij[l][i][j]!=(float)MISSVAL) {

					val=corrlij[l][i][j];

					if(Mdij[0][0]==1.){
						if(Mdij[i][j]!=(float)MISSVAL) dij=Mdij[i][j];
						else continue;
					}
					else if(Mdij[0][0]==-1.) dij=acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
					else dij=sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
				
 					if(Mdij[0][0]==0 && (sgi[i]==sgi[j])){ dij=-1.; cneighb=2;}	/*intra-group class*/

					/*compute correlation/regression btw pairwise autocorrelation coef & dist*/
					if(dij>0. && dij>=dijmin && (dij<=(dwidth*dijmin) || dijmin==0)){

						logd=log(dij);
						npairslog++;
						Svlog+=val;
						SSvlog+=val*val;
						Sdlog+=logd;
						SSdlog+=logd*logd;
						SPvdlog+=logd*val;
					}
				} //end loop j

			}/*end loop i*/

			/*reg with log d*/
			if(npairslog){

				SSDvlog=SSvlog-Svlog*Svlog/npairslog;
				if(SSDvlog<0.000000001) SSDvlog=0.;/*when val is constant, Vv might however take some value due to error of round*/
				SSDdlog=SSdlog-Sdlog*Sdlog/npairslog;
				if(SSDdlog<0.000000001) SSDdlog=0.;/*when val is constant, Vd might however take some value due to error of round*/

				if(SSDdlog && SSDvlog){
					SPDvdlog=SPvdlog-Svlog*Sdlog/npairslog;
					b=SPDvdlog/SSDdlog;
				}
				else b=MISSVAL;
			}
			else b=MISSVAL;

			if(Stat==1 || Stat==2 || Stat==11){ //for kinship estimators
				if(b<0. && b!=MISSVAL){
				//	Nb[niter]=1./(  (corrlc[l][cneighb]-1.)/b  );
				//	sigmasq=(1./Nb[niter])/(2*ploidy*3.1415*(double)density);
					Nb[niter]=(corrlc[l][cneighb]-1.)/b  ;
					sigmasq=Nb[niter]/(2*ploidy*3.1415*(double)density);
					
					sigmaest[niter]=sqrt(sigmasq);
					dijmin=sigmaest[niter];
				}
				else sigmaest[niter]=Nb[niter]=MISSVAL;
			}
			if(Stat==4){ //for Rousset estimators
				if(b>0. && b!=MISSVAL){
				//	Nb[niter]=1./(  1./b  );
				//	sigmasq=(1./Nb[niter])/(2*ploidy*3.1415*(double)density);
					Nb[niter]=1./b;
					sigmasq=Nb[niter]/(2*ploidy*3.1415*(double)density);
					sigmaest[niter]=sqrt(sigmasq);
					dijmin=sigmaest[niter];
				}
				else sigmaest[niter]=Nb[niter]=MISSVAL;
			}

			if(l>=0) printf("\n%i\t%i\t%G\t%G",niter,npairslog,Nb[niter],sigmaest[niter]);

			if(sigmaest[niter]!=MISSVAL){
				if(fabs((sigmaest[niter]-sigmaest[niter-1])/sigmaest[niter])<0.01) converge=1;
				else converge=0;
			}

		}while(sigmaest[niter]!=MISSVAL && converge==0 && niter<100);


		corrlc[l][-21]=corrlc[l][-22]=corrlc[l][-23]=(float)Nb[niter];
		corrlc[l][-25]=corrlc[l][-26]=corrlc[l][-27]=(float)sigmaest[niter];
		corrlc[l][-24]=(float)npairslog;

		if(niter==100){ //to get average over a cycle
			backiter=1;
			meansigma=minsigma=maxsigma=sigmaest[niter];
			meanNb=minNb=maxNb=Nb[niter];
			while((sigmaest[niter]!=sigmaest[niter-backiter]) && backiter<100){
				meansigma+=sigmaest[niter-backiter];
				meanNb+=Nb[niter-backiter];
				if(sigmaest[niter-backiter]<minsigma) minsigma=sigmaest[niter-backiter];
				if(sigmaest[niter-backiter]>maxsigma) maxsigma=sigmaest[niter-backiter];
				if(Nb[niter-backiter]<minNb) minNb=Nb[niter-backiter];
				if(Nb[niter-backiter]>maxNb) maxNb=Nb[niter-backiter];
				backiter++;
			}					
			corrlc[l][-21]=(float)meanNb/backiter;
			corrlc[l][-22]=(float)minNb;
			corrlc[l][-23]=(float)maxNb;
			corrlc[l][-25]=(float)meansigma/backiter;
			corrlc[l][-26]=(float)minsigma;
			corrlc[l][-27]=(float)maxsigma;
		}
	}

	//Jacknife estimates
	if(JKest){
		converge=1;
		for(l=1;l<=m;l++) if(corrlc[-l][-21]==(float)MISSVAL) { converge=0; break; }

		if(converge){
			nJK=m;
			corrlc[m+1][-21]=corrlc[m+2][-21]=corrlc[m+1][-25]=corrlc[m+2][-25]=0.;
			for(l=1;l<=m;l++)if(corrlc[-l][-21]!=(float)MISSVAL){ //compute pseudovalues
				corrlc[l][-21]=(float)(nJK*fabs((double)corrlc[0][-21])-(nJK-1)*fabs((double)corrlc[-l][-21]));
				corrlc[m+1][-21]+=corrlc[l][-21];
				corrlc[l][-25]=(float)(nJK*fabs((double)corrlc[0][-25])-(nJK-1)*fabs((double)corrlc[-l][-25]));
				corrlc[m+1][-25]+=corrlc[l][-25];
			}
			corrlc[m+1][-21]/=nJK;   //mean JK estimate
			corrlc[m+1][-25]/=nJK;   //mean JK estimate
			for(l=1;l<=m;l++)if(corrlc[-l][-21]!=(float)MISSVAL){ //compute variance of pseudoval
				corrlc[m+2][-21]+=(corrlc[l][-21]-corrlc[m+1][-21])*(corrlc[l][-21]-corrlc[m+1][-21]);
				corrlc[m+2][-25]+=(corrlc[l][-25]-corrlc[m+1][-25])*(corrlc[l][-25]-corrlc[m+1][-25]);
			}
			corrlc[m+2][-21]=(float)sqrt((corrlc[m+2][-21]/(nJK-1))/nJK); //standard error JK estimate
			corrlc[m+2][-25]=(float)sqrt((corrlc[m+2][-25]/(nJK-1))/nJK); //standard error JK estimate
		}
		else corrlc[m+1][-21]=corrlc[m+2][-21]=corrlc[m+1][-25]=corrlc[m+2][-25]=(float)MISSVAL;
	}
	else corrlc[m+1][-21]=corrlc[m+2][-21]=corrlc[m+1][-25]=corrlc[m+2][-25]=(float)MISSVAL;


}	/*end procedure estimate_sigma_2D_kinship*/


/*************************************************************************************/



void analyse_resampling(int m,int cinit,int nc,int Npermut,float ***corrlcp,
		struct resample_stat_type **r_statlc)
{
	int l,linit,c,p;
	int validp;
	float *val;

	val=vector(0,Npermut);
	
	if(m==1) linit=1;
	else linit=0;

	if(Npermut>0) for(l=linit;l<=m;l++){
		/*make stat analysis of the corrlcp and store it in r_statlc[l][c]*/
		for(c=cinit;c<=nc;c++){
			/*eliminate the bad values in frlcp*/
			validp=0;
			val[0]=corrlcp[l][c][0];
			for(p=1;p<=Npermut;p++){
				if(corrlcp[l][c][p]!=(float)MISSVAL){
  					validp++;
					val[validp]=corrlcp[l][c][p];
				}	
			}
			/*make the stat analysis*/
			r_statlc[l][c].n=validp;
			if(validp && val[0]!=(float)MISSVAL) resample_analyze(val,1,validp,&r_statlc[l][c],val[0]);
			else{
				r_statlc[l][c].nd=r_statlc[l][c].mode=0;
				r_statlc[l][c].obs=r_statlc[l][c].mean=r_statlc[l][c].sd=r_statlc[l][c].low95=r_statlc[l][c].high95=r_statlc[l][c].low99=r_statlc[l][c].high99=r_statlc[l][c].plow=r_statlc[l][c].phigh=r_statlc[l][c].pbil=(float)MISSVAL;	
			}
	/*		if(r_statlc[l][c].plow==r_statlc[l][c].phigh && r_statlc[l][c].phigh==1.)	r_statlc[l][c].obs=r_statlc[l][c].mean=r_statlc[l][c].low95=r_statlc[l][c].high95=r_statlc[l][c].low99=r_statlc[l][c].high99=r_statlc[l][c].plow=r_statlc[l][c].phigh=r_statlc[l][c].pbil=(float)MISSVAL;		*/
		}/*end of loop c*/

	} /*end of loop l,v*/
	else for(l=linit;l<=m;l++)for(c=cinit;c<=nc;c++){
		r_statlc[l][c].n=r_statlc[l][c].nd=r_statlc[l][c].mode=0;
		r_statlc[l][c].obs=r_statlc[l][c].mean=r_statlc[l][c].sd=r_statlc[l][c].low95=r_statlc[l][c].high95=r_statlc[l][c].low99=r_statlc[l][c].high99=r_statlc[l][c].plow=r_statlc[l][c].phigh=r_statlc[l][c].pbil=(float)MISSVAL;	
	}
 
	free_vector(val,0,Npermut);

}/*end of procedure analyse_resampling*/ 


/*************************************************************************************/

/*compute variance components for a 2-level nested ANOVA with unequal sample sizes, 
as described in Sokal and Rohlf, 1995, p. 294
Yijk : value of the variable for sample k (k=1 to Nij) of subgroup j (j=1 to bi) from group i (i=1 to a)
a : # groups 
bi : # subgroups in within group i
Nij : # samples in subgroup j from (nested in) group i

It gives the sum of squared deviates SS, mean squares MS and variance components s2
among groups (SS[1], MS[1], s2[1]),
among subgroups within groups (SS[2], MS[2], s2[2]),
among samples within subgroups (SS[3], MS[3], s2[3]), 
and the totals over all three levels of variation (SS[0], s2[0])*/

void NestedANOVA(int a,int *bi,int **Nij,double ***Yijk,double SS[4],double MS[4],double s2[4]){
	int i,j,k,r;				 /*counters*/
	double **mYij;	 /*averages: for i,j when i>0 and j>0; for i when j=0; grand average for i=j=0*/
	double 	S_Nij2, S_SNij2, S_SNij2_over_SNij, SNij2temp;	/*sum of Nij squared*/
	double n0,n0bis,nb0;	/*"average" samples*/
	int df[4];		/*degrees of freedom*/
	int bimax;

	bimax=0;
	for(i=1;i<=a;i++) if(bi[i]>bimax) bimax=bi[i];

	mYij=dmatrix(0,a,0,bimax);

	/*compute means and sample sizes*/
	Nij[0][0]=0;	   
	mYij[0][0]=0.;	   
	for(i=1;i<=a;i++){
		Nij[i][0]=0;   
		mYij[i][0]=0.; 
		for(j=1;j<=bi[i];j++){
			mYij[i][j]=0.;	
			Nij[i][0]+=Nij[i][j]; /*tot # samples for group i*/
			for(k=1;k<=Nij[i][j];k++) mYij[i][j]+=Yijk[i][j][k];
			mYij[i][0]+=mYij[i][j];
			if(Nij[i][j]) mYij[i][j]/=Nij[i][j];	/*mean of subgroup j from group i*/
		}
		Nij[0][0]+=Nij[i][0];	/*total # samples*/
		mYij[0][0]+=mYij[i][0];
		if(Nij[i][0]) mYij[i][0]/=Nij[i][0];	/*mean of group i*/
	}
	if(Nij[0][0]) mYij[0][0]/=Nij[0][0];		/*grand mean*/
	else {for(r=0;r<=3;r++) SS[r]=MS[r]=s2[r]=(float)MISSVAL; return;}


	/*compute sum of squared deviates among groups (SS[1]),
	among subgroups within groups (SS[2]), within subgroups - i.e. among samples 
	within subgroups - (SS[3]), and total (SS[0])*/
	SS[0]=SS[1]=SS[2]=SS[3]=0.;	   
	for(i=1;i<=a;i++) SS[1]+=Nij[i][0]*(mYij[i][0]-mYij[0][0])*(mYij[i][0]-mYij[0][0]);
	for(i=1;i<=a;i++)for(j=1;j<=bi[i];j++) SS[2]+=Nij[i][j]*(mYij[i][j]-mYij[i][0])*(mYij[i][j]-mYij[i][0]);
	for(i=1;i<=a;i++)for(j=1;j<=bi[i];j++)for(k=1;k<=Nij[i][j];k++) SS[3]+=(Yijk[i][j][k]-mYij[i][j])*(Yijk[i][j][k]-mYij[i][j]);
/**/for(i=1;i<=a;i++)for(j=1;j<=bi[i];j++)for(k=1;k<=Nij[i][j];k++) SS[0]+=(Yijk[i][j][k]-mYij[0][0])*(Yijk[i][j][k]-mYij[0][0]);
	SS[0]=SS[1]+SS[2]+SS[3];


	/*compute degrees of freedom*/
	df[1]=-1;
	for(i=1;i<=a;i++)if(Nij[i][0])df[1]++;
	df[2]=0;
	for(i=1;i<=a;i++)if(Nij[i][0]){
		for(j=1;j<=bi[i];j++)if(Nij[i][j]) df[2]++;
		df[2]--;
	}
	df[3]=0;
	for(i=1;i<=a;i++)for(j=1;j<=bi[i];j++)if(Nij[i][j]) df[3]+=Nij[i][j]-1;
/**/df[0]=df[1]+df[2]+df[3];
	df[0]=Nij[0][0]-1;

	/*compute Mean Squares*/
	if(df[1]) MS[1]=SS[1]/df[1];
	else MS[1]=(float)MISSVAL;
	if(df[2]) MS[2]=SS[2]/df[2];
	else MS[2]=(float)MISSVAL;
	if(df[3]) MS[3]=SS[3]/df[3];
	else MS[3]=(float)MISSVAL;

	/*compute "average" sample sizes*/
	S_Nij2=0.;
	for(i=1;i<=a;i++)for(j=1;j<=bi[i];j++) S_Nij2+=Nij[i][j]*Nij[i][j];

	S_SNij2=0.;
	for(i=1;i<=a;i++) S_SNij2+=Nij[i][0]*Nij[i][0];

	S_SNij2_over_SNij=0.;
	for(i=1;i<=a;i++)if(Nij[i][0]){
		SNij2temp=0.;
		for(j=1;j<=bi[i];j++) SNij2temp+=Nij[i][j]*Nij[i][j];
		S_SNij2_over_SNij+=SNij2temp/Nij[i][0];
	}

	n0bis=n0=nb0=0.;
	if(df[1]) n0bis=(S_SNij2_over_SNij-(S_Nij2/Nij[0][0]))/df[1];
	if(df[2]) n0=(Nij[0][0]-S_SNij2_over_SNij)/df[2];
	if(df[1]) nb0=(Nij[0][0]-(S_SNij2/Nij[0][0]))/df[1];

	/*compute variance components*/
	s2[3]=MS[3];
	if(MS[3]!=(float)MISSVAL){
		if(n0) s2[2]=(MS[2]-s2[3])/n0;
		else s2[2]=(float)MISSVAL;
	   	if(nb0) s2[1]=(MS[1]-s2[3]-n0bis*s2[2])/nb0;
		else s2[1]=(float)MISSVAL;
	}
	else{			/*one-way ANOVA*/
		s2[2]=MS[2];
		if(nb0) s2[1]=(MS[1]-s2[2])/nb0;
		else s2[1]=(float)MISSVAL;
	}

	s2[0]=0.;
	for(r=1;r<=3;r++) if(s2[r]!=(float)MISSVAL) s2[0]+=s2[r];

	free_dmatrix(mYij,0,a,0,bimax);
}

/************************************************************************************/
void permut_allelesizes
	(int m,int *Nallelel,int **allelesizela,int **allelesizelamix,long *seed)
{
	int a,l;

	for(l=1;l<=m;l++){
		for(a=1;a<=Nallelel[l];a++)allelesizelamix[l][a]=allelesizela[l][a];
		resample_shuffle(allelesizelamix[l],1,Nallelel[l],seed);
	}
}
/************************************************************************************/

void permut_allelesizes_among_2pop_alleles
	(int m,int *Nallelel,float ***Ppla,int pop1,int pop2,int **allelesizela,int **allelesizelamix,long *seed)
{
	int a,l,twopopa,size2popa[1000];

	for(l=1;l<=m;l++){

		for(a=1;a<=Nallelel[l];a++)allelesizelamix[l][a]=allelesizela[l][a];

		twopopa=0;
		for(a=1;a<=Nallelel[l];a++) if(Ppla[pop1][l][a] || Ppla[pop2][l][a]){
			twopopa++;
			size2popa[twopopa]=allelesizela[l][a];
		}
		resample_shuffle(size2popa,1,twopopa,seed);

		twopopa=0;
		for(a=1;a<=Nallelel[l];a++) if(Ppla[pop1][l][a] || Ppla[pop2][l][a]){
			twopopa++;
			allelesizelamix[l][a]=size2popa[twopopa];
		}
	}
}
/************************************************************************************/
void permut_genetic_distances
	(int m,int *Nallelel,float ***Mgdlaa,float ***Mgdlaamix,long *seed)
{
	int a,a1,a2,l,loca[MAXALLID];


	for(l=1;l<=m;l++){
		for(a=0;a<=Nallelel[l];a++) loca[a]=a;
		resample_shuffle(loca,1,Nallelel[l],seed);

		for(a1=1;a1<=Nallelel[l];a1++){
			for(a2=1;a2<=Nallelel[l];a2++) Mgdlaamix[l][a1][a2]=Mgdlaa[l][loca[a1]][loca[a2]];
		}
	}
}
/************************************************************************************/
void permut_genetic_distances_among_2pop_alleles
	(int m,int *Nallelel,float ***Ppla,int pop1,int pop2,float ***Mgdlaa,float ***Mgdlaamix,long *seed)
{
	int a,a1,a2,l,twopopa,twopopa1,twopopa2,loca[MAXALLID];


	for(l=1;l<=m;l++){

		twopopa=0;
		for(a=1;a<=Nallelel[l];a++) if(Ppla[pop1][l][a] || Ppla[pop2][l][a]){
			twopopa++;
			loca[twopopa]=a;
		}
		resample_shuffle(loca,1,twopopa,seed);

		twopopa1=0;
		for(a1=1;a1<=Nallelel[l];a1++) if(Ppla[pop1][l][a1] || Ppla[pop2][l][a1]){
			twopopa1++;
			twopopa2=0;
			for(a2=1;a2<=Nallelel[l];a2++) if(Ppla[pop1][l][a2] || Ppla[pop2][l][a2]){ 
				twopopa2++;
				Mgdlaamix[l][a1][a2]=Mgdlaa[l][loca[twopopa1]][loca[twopopa2]];
			}
		}
	}
}
/************************************************************************************/


void permut_locations
	(int n,double *x,double *y,double *z,double **Mdij,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed)
{
	int i,j,loci[NMAX];

	for(i=0;i<=n;i++) loci[i]=i;

	resample_shuffle(loci,1,n,seed);

	for(i=1;i<=n;i++){
		if(Mdij[0][0]==1.) for(j=1;j<=n;j++) Mdijmix[i][j]=Mdij[loci[i]][loci[j]];
		else{
			xmix[i]=x[loci[i]];
			ymix[i]=y[loci[i]];
			zmix[i]=z[loci[i]];
		}
	}
	
}
/************************************************************************************/
void permut_locations_of_groups
	(int n,double *x,double *y,double *z,double **Mdij,int *groupi,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed)
{
	int i,j,loci[NMAX];
	int locg[NMAX],group[NMAX],newgroup,ng,g;

	ng=0;
	for(i=1;i<=n;i++){
		newgroup=1;
		for(g=ng;g>=1;g--)if(groupi[i]==group[g]) {
			newgroup=0;
			break;
		}
		if(newgroup){
			ng++;
			group[ng]=groupi[i];
			locg[ng]=i;
		}
	}
	resample_shuffle(locg,1,ng,seed);

	for(i=1;i<=n;i++) loci[i]=locg[groupi[i]];

	for(i=1;i<=n;i++){
		if(Mdij[0][0]==1.) for(j=1;j<=n;j++) Mdijmix[i][j]=Mdij[loci[i]][loci[j]];
		else{
			xmix[i]=x[loci[i]];
			ymix[i]=y[loci[i]];
			zmix[i]=z[loci[i]];
		}
	}
	
}
/************************************************************************************/
void permut_locations_within_cat
	(int n,int *cati,int Ncat,double *x,double *y,double *z,double **Mdij,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed)
{
	int i,j,i2,j2,k,loci[NMAX];

	for(k=1;k<=Ncat;k++){
		j=0;
		for(i=1;i<=n;i++) if(cati[i]==k){
			j++;
			loci[j]=i;
		}
		resample_shuffle(loci,1,j,seed);
		j=0;
		for(i=1;i<=n;i++) if(cati[i]==k){
			j++;
			if(Mdij[0][0]!=1.){
				xmix[i]=x[loci[j]];
				ymix[i]=y[loci[j]];
				zmix[i]=z[loci[j]];
			}
			j2=0;
			if(Mdij[0][0]==1.) for(i2=1;i2<=n;i2++) if(cati[i2]==k){
				j2++;
				Mdijmix[i][i2]=Mdij[loci[j]][loci[j2]];
			}
		}
	}
}

/************************************************************************************/
void permut_locations_of_groups_within_cat
	(int n,int *cati,int Ncat,double *x,double *y,double *z,double **Mdij,int *groupi,
	double *xmix,double *ymix,double *zmix,double **Mdijmix,long *seed)
{
	int i,j,i2,j2,k,loci[NMAX];
	int locg[NMAX],group[NMAX],newgroup,ng,g;

	for(k=1;k<=Ncat;k++){
		ng=0;
		for(i=1;i<=n;i++) if(cati[i]==k){
			newgroup=1;
			for(g=ng;g>=1;g--)if(groupi[i]==group[g]) {
				newgroup=0;
				break;
			}
			if(newgroup){
				ng++;
				group[ng]=groupi[i];
				locg[ng]=i;
			}
		}
		resample_shuffle(locg,1,ng,seed);

		for(i=1;i<=n;i++) if(cati[i]==k) loci[i]=locg[groupi[i]];

		for(i=1;i<=n;i++) if(cati[i]==k){			
			if(Mdij[0][0]!=1.){
				xmix[i]=x[loci[i]];
				ymix[i]=y[loci[i]];
				zmix[i]=z[loci[i]];
			}
			if(Mdij[0][0]==1.) for(i2=1;i2<=n;i2++) if(cati[i2]==k){
				Mdijmix[i][i2]=Mdij[loci[i]][loci[i2]];
			}
		}
	}
}

/************************************************************************************/
void permut_genes_among_indiv
	(int n,int m,int ploidy,int *ploidyi,int ***gilc,int ***gilcmix,long *seed)
{
	int i,l,c,Nvalidgene;
	int *gene;

	gene=ivector(0,n*ploidy);
	for(l=1;l<=m;l++){
		Nvalidgene=0;
		for(i=1;i<=n;i++)for(c=0;c<ploidy;c++)if(gilc[i][l][c]){
			Nvalidgene++;
			gene[Nvalidgene]=gilc[i][l][c];
		}

		resample_shuffle(gene,1,Nvalidgene,seed);

		Nvalidgene=0;
		for(i=1;i<=n;i++)for(c=0;c<ploidy;c++){
			if(gilc[i][l][c]){
				Nvalidgene++;
				gilcmix[i][l][c]=gene[Nvalidgene];
			}
			else gilcmix[i][l][c]=0;
		}
	}

	free_ivector(gene,0,n*ploidy);
}
/************************************************************************************/
void permut_genes_among_indiv_within_pop
	(int n,int m,int ploidy,int *ploidyi,int *popi, int Npop,
	int ***gilc,int ***gilcmix,long *seed)
{
	int i,l,c,Nvalidgene,p;
	int *gene;

	gene=ivector(0,n*ploidy);
	for(l=1;l<=m;l++){
		for(p=1;p<=Npop;p++){
			Nvalidgene=0;
			for(i=1;i<=n;i++)if(popi[i]==p){
				for(c=0;c<ploidy;c++)if(gilc[i][l][c]){
					Nvalidgene++;
					gene[Nvalidgene]=gilc[i][l][c];
				}
			}
			resample_shuffle(gene,1,Nvalidgene,seed);

			Nvalidgene=0;
			for(i=1;i<=n;i++)if(popi[i]==p){
				for(c=0;c<ploidy;c++){
					if(gilc[i][l][c]){
						Nvalidgene++;
						gilcmix[i][l][c]=gene[Nvalidgene];
					}
					else gilcmix[i][l][c]=0;
				}
			}
		}
	}

	free_ivector(gene,0,n*ploidy);
}
/************************************************************************************/
void permut_indiv_among_pop
	(int n,int *popi, int Npop,int *popimix,long *seed)
{
	int i;
	int *ind;

	ind=ivector(0,n);

	for(i=1;i<=n;i++)ind[i]=i;
	resample_shuffle(ind,1,n,seed);
	for(i=1;i<=n;i++) popimix[i]=popi[ind[i]];

	free_ivector(ind,0,n);
}
/************************************************************************************/
void permut_indiv_among_pop_within_categ
	(int n,int *popi, int Npop,int *cati,int Ncat,int *popimix,long *seed)
{
	int i,k,Ni;
	int *ind;

	ind=ivector(0,n);

	for(k=1;k<=Ncat;k++){
		Ni=0;
		for(i=1;i<=n;i++)if(cati[i]==k){
			Ni++;
			ind[Ni]=i;
		}
		resample_shuffle(ind,1,Ni,seed);

		Ni=0;
		for(i=1;i<=n;i++)if(cati[i]==k){
			Ni++;
			popimix[i]=popi[ind[Ni]];
		}
	}

	free_ivector(ind,0,n);
}
/************************************************************************************/

