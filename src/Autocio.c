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

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#	include <unistd.h>
#endif
#ifdef HAVE_DIRECT_H
#	include <direct.h>
#	define chdir _chdir
#endif
#ifdef HAVE_LIBGEN_H
#	include <libgen.h>
#endif

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Autoccomp.h"
#include "Autocio.h" 
#include "Xatools.h"

#ifdef USE_GETCH
#	if defined(CURSES_HAVE_CURSES_H)
#		include <curses.h>
#	elif defined(CURSES_HAVE_NCURSES_H)
#		include <ncurses.h>
#	elif defined(CURSES_HAVE_NCURSES_NCURSES_H)
#		include <ncurses/ncurses.h>
#	elif defined(CURSES_HAVE_NCURSES_CURSES_H)
#		include <ncurses/curses.h>
#	endif
#	ifdef HAVE_CONIO_H
#		include <conio.h>
#	endif
#endif

#define write write_string

char errorfile[PATH_BUF_SIZE] = ERRORTXT;

void commom_errors()
{
	FILE *fp;

	while((fp=fopen(ERRORFILE,"wt"))==NULL){
		printf("\n\nWARNING: Cannot open error file \"%s\"\n  If it is being used by another application, close it first, then press RETURN\n  Close the window if you wish to stop the program now\n",ERRORFILE);
		wait_a_char();
	}

	fprintf(fp,"COMMON ERRORS IN DATA FILES");
	fprintf(fp,"\n- Be careful that the numbers of individuals and categories given (format numbers) are correct.");
	fprintf(fp,"\n- Use a point, not a coma, for decimal separator when defining distances or spatial coordinates.");
	fprintf(fp,"\n- Spatial coordinates must be given following an orthogonal system\n\t(degree of longitude / latitude will not be interpreted correctly, use an UTM system instead).");
	fprintf(fp,"\n- Alleles must be identified by numbers, not letters.");
	fprintf(fp,"\n- Empty columns are not allowed. \n\tMissing alleles or genotypes must be specified by 0 (or 00 or 000 according to the number of digits per allele).");
	fprintf(fp,"\n- If one allele of a genotype is undefined, put 0 (or 00, or 000) on the right side. \n\t(if put on the left side the genotype would be interpreted as of a lower ploidy level.)");
	fprintf(fp,"\n- The list of genotypes must be followed by a line with the word END.");
	fprintf(fp,"\n- When defining spatial or genetic distances, be careful that the matrix size or the number of lines (column format) is correct.");
	fprintf(fp,"\n- When defining pairwise spatial distances, be careful that individual or population names are correct (case matter!) and all are distinct.\n\tNote that population names vary according to the way populations are defined (check the manual).");
	fprintf(fp,"\n- When defining pairwise spatial distances, be careful that they are given for the correct level of analysis (individuals or populations).");
	fprintf(fp,"\n- When defining pairwise genetic distances between alleles, be careful that locus names are correct (case matter!).");
	fprintf(fp,"\n- If you get the message \"End of file unexpected on line ... from file ...\", check the first line of pairwise spatial distances (M... or C...) or genetic distances (locus names)");
	fprintf(fp,"\n- Statistics such as Nst or dm2 require a matrix with pairwise genetic distances between alleles.\n\tDo not select these statistics otherwise.");
	fprintf(fp,"\n- Some worksheet or text editors might produce a file not readable by SPAGeDi according to the saving format,\n\tin which case the error message is: \"WARNING: cannot read the first line containing the 6 format numbers...\"\n\t'Tab delimited text' or 'DOS text' saving formats should work.");
	
	fprintf(fp,"\n\nSPECIFIC ERRORS FOUND DURING THE LAST LAUNCH OF SPAGeDi\n");

	fclose(fp);
}

// remove newline character that fgets retains
char * chomp(char *str) {
	int i=0;
	if(str == NULL)
		return NULL;
	while(str[i] != '\0') { i++;}
	if(str[i-1] == '\n')
		str[i-1] = '\0';
	return str;
}

char *fgets_chomp(char * a, int b, FILE * c) {
	char * ret = fgets(a, b, c);
	chomp(a);
	return ret;
}

void copy_file_name(char *to, const char *from, const char *outdir) {
	const char sps[] = " \n\r\"\'";
	char filename[PATH_BUF_SIZE];
	size_t len,len2;
	
	len = strspn(from,sps); // some drag-n-drop strings can be quoted
	from += len;
	strlcpy(to,from,PATH_BUF_SIZE);
	len = strcspn(to,sps);  // skip non-sps chars
	while(to[len] != '\0') {
		len2 = len+strspn(to+len,sps); // skip sps chars
		if(to[len2] == '\0') {
			to[len] = '\0'; // shorten string
			break;
		}
		len = len2+strcspn(to+len2,sps);
	}
	if(outdir != NULL && to[0] != '/') {
		// Relative filename, so make it relative to outdir
		strlcpy(filename, outdir, PATH_BUF_SIZE);
		strlcat(filename, "/", PATH_BUF_SIZE);
		strlcat(filename, to, PATH_BUF_SIZE);
		strlcpy(to, filename, PATH_BUF_SIZE);
	}
}

/****************************************************************************************/

void get_input_output_file_names(int argc,char *argv[],char inputfile[],char outputfile[],char instrfile[])
{
	char smess[SMAX], filename[PATH_BUF_SIZE], outdir[PATH_BUF_SIZE], ch,*ptr, *ptrt;
	int check=1;
	FILE *fp;

	// printf("\n%s",argv[0]);
	// if(argc>1)printf("\n%s",argv[1]);
	// if(argc>2)printf("\n%s",argv[2]);
	
	printf("\n\n\n\nDATA / RESULTS FILE NAMES");
	
	// Determine Inputfile
	if( argc > 1 ) {
		strlcpy(inputfile, argv[1], PATH_BUF_SIZE);
		// Find base name
	} else { // Filename not specified on command line
		printf("\n\nEnter the name of the data file (with ext)\n"
		       "  or press RETURN for the default data file \"in.txt\"\n"
		       "  or enter a SPACE for importing a data file from FSTAT or GENEPOP format\n\n"
		       "  Data file: ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') {
			strlcpy(inputfile, "in.txt", PATH_BUF_SIZE);
			printf("%s\n",inputfile);
		} else if(smess[0]==' ')
			import_data_file(inputfile);
		else {
			copy_file_name(inputfile, smess, NULL);
		}
	}
	// Check inputfile
	while((fp=fopen(inputfile,"rt"))==NULL) {
		printf("\n\nWARNING: Cannot open data file \"%s\"\n"
		       "  If it is being used by another application, close it first, then press RETURN\n"
		       "  Otherwise enter a new name for the data file\n"
		       "  Press ctrl+c if you wish to stop the program now: ", inputfile);
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0') {
			copy_file_name(inputfile, smess, NULL);
		}
	}
	fclose(fp);

	// Find base name
	strlcpy(filename, basename(inputfile), PATH_BUF_SIZE);
	printf("\n\n  Data file: %s\n    full path: %s\n", filename, inputfile);
	
	// find directory of input file name, will be used as initial outdir
	// use filename as initial buffer in case dirname modifies it
	strlcpy(filename, inputfile, PATH_BUF_SIZE);
	strlcpy(outdir, dirname(filename), PATH_BUF_SIZE);
	
	// Determine Outputfile
	if( argc > 2) {
		strlcpy(outputfile, argv[2], PATH_BUF_SIZE);
		strlcpy(outdir, dirname(outputfile), PATH_BUF_SIZE);	
	} else { // Filename not specifed on command line
		printf("\n\nEnter the name of the results file (with ext)\n"
		       "or press RETURN for the default results file \"%s\"\n\n"
		       "  Results file: ", "out.txt");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') {
			// Create default output filename: outdir+"out.txt"
			strlcpy(outputfile, outdir, PATH_BUF_SIZE);
			strlcat(outputfile, "/out.txt", PATH_BUF_SIZE);
		} else {
			copy_file_name(outputfile, smess, outdir);
		}
	}
	while(check != 0) {
		check = 0;
		if((fp=fopen(outputfile,"r"))) {
			fclose(fp);
			printf("\n\nWARNING: The results file \"%s\" already exists.\n"
			       "  Enter 'a' or press RETURN to add the new results to the end of the file\n"
			       "  Enter 'e' to erase the present content of the file\n"
			       "  Otherwise, enter a new name for the results file: ", outputfile);
			fgets_chomp(smess, sizeof(smess), stdin);
			// Clear file
			if(!strcmp(smess,"e") || !strcmp(smess,"E")) {
				if((fp=fopen(outputfile,"wt"))!=NULL) 
					fclose(fp);
			}
			// New file specified?
			else if(strcmp(smess,"") && strcmp(smess,"a") && strcmp(smess,"A")) {
				check=1;
				copy_file_name(outputfile, smess, outdir);
			}
		}
 		while(check==0 && (fp=fopen(outputfile,"a"))==NULL ){
			printf("\nWARNING: Cannot open results file \"%s\".\n"
			       "If it is being used by another application, close it first. Then press RETURN.\n"
			       "  Otherwise, enter a new name for the results file: ", outputfile);
			fgets_chomp(smess, sizeof(smess), stdin);
			if(smess[0]!='\0'){
				copy_file_name(outputfile,smess,outdir);
				check=1;
			}
		}
		fclose(fp);
	}
	
	// Find base name
	strlcpy(filename, basename(outputfile), PATH_BUF_SIZE);
	printf("\n\n  Results file: %s\n    full path: %s\n", filename, outputfile);	
	
	// Use outdir of output
	strlcpy(outdir, dirname(outputfile), PATH_BUF_SIZE);
	
	// Create error.txt in outputdir
	strlcpy(errorfile, outdir, PATH_BUF_SIZE);
	strlcat(errorfile, "/", PATH_BUF_SIZE);
	strlcat(errorfile, ERRORTXT, PATH_BUF_SIZE);
}
/*end of get_input_output_file_names*/



/***************************************************************************************/




void import_data_file(char *inputfile)
{
	int line=0,line3=0,flag;
	int format,n,i,check;
	int m=0,Maxallele,ndigit,l;
	int Npop,npop,nump,pop,labels,numpop[501],popi;
	int genotypel[501];
	char s[SMAX], *s2,*s3,smess[SMAX],sbis[SMAX],*ptr;
	char ch,importfile[50],labelfile[50];
	char namelocus[MMAX][50],popname[501][50],indname[50];
	double X,Y;
	FILE *fp,*fp2,*fp3;

	/*define the data file to import*/
	format=0;
	do{
		printf("\n\nImport a data file in\n 1- FSTAT format\n 2- GENEPOP format with 2 digits per allele\n 3- GENEPOP format with 3 digits per allele\n");
		fgets_chomp(smess, sizeof(smess), stdin);
		sscanf(smess,"%i",&format);
	}while(format<1 && format>3);

	printf("\n\nEnter the name of the data file to import (with ext if necessary): ");
	fgets_chomp(importfile, sizeof(importfile), stdin);
	while((fp=fopen(importfile,"rt"))==NULL){
		printf("\nWARNING: Cannot open the file \"%s\".\n  If it is being used by another application, close it first.\n  If it does not exist in the right directory, bring it first and press RETURN.\n  Otherwise enter a new name for the data file to import\n  To stop the programm now, press Ctrl+c.\n",importfile);
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0') strncpy(importfile,smess,MAXNOM-1);
	}

	/*define the name of the data file in SPAGeDi format*/
	printf("\nEnter the name of the data file in format SPAGeDi: ");
	do{
		scanf("%s",inputfile);
	}while(inputfile[0]=='\0'); 
	do{
		check=0;
		if((fp2=fopen(inputfile,"r"))){
			fclose(fp2);
			printf("\nWARNING: File \"%s\" already exists.\n  Enter 'e' to erase the present content of the file\n  Otherwise, enter a new file name\n",inputfile);
			scanf("%s",smess);
			if(!strcmp(smess,"e") || !strcmp(smess,"E")) {
				if((fp2=fopen(inputfile,"wt"))!=NULL) fclose(fp2);
			}
			else {strncpy(inputfile,smess,MAXNOM-1);check=1;}
		}
	}while(check);

	while((fp2=fopen(inputfile,"a"))==NULL){
		printf("\nWARNING: Cannot open data file \"%s\".\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",inputfile);
		wait_a_char();
	}

	
	if(format==1){/*read Fstat data file*/

		/*read 4 format numbers*/
		readsfromfile(fp,s,importfile,&line);
		flag=sscanf(s,"%d%d%d%d",&Npop,&m,&Maxallele,&ndigit);
		if(flag!=4) {
			sprintf(smess,"\ncannot read # of samples, # of loci, highest # to label alleles, or # of digits to label alleles on line %d from file %s",line,importfile);
			write(ERRORFILE,smess);
			printf("\nPress any key to stop the program");
			wait_a_char();exit(1);
		}
		if(Npop > 500 || Npop < 1) {
			sprintf(smess,"\nNumber of samples =%d not allowed (it must be between 1 and 500)",Npop);
			write(ERRORFILE,smess);
			printf("\nPress any key to stop the program");
			wait_a_char();exit(1);
		}
		if(m > MMAX || m < 1) {
			sprintf(smess,"\nNumber of loci =%d not allowed (it must be between 1 and %d)",m,MMAX);
			write(ERRORFILE,smess);
			printf("\nPress any key to stop the program");
			wait_a_char();exit(1);
		}
		if(Maxallele > 999 || Maxallele < 2) {
			sprintf(smess,"\nMaximum # to label alleles =%d not allowed (it must be between 2 and 999)",Maxallele);
			write(ERRORFILE,smess);
			printf("\nPress any key to stop the program");
			wait_a_char();exit(1);
		}
		if(ndigit > 3 || ndigit < 1) {
			sprintf(smess,"\nNumber of digits to label alleles =%d not allowed (it must be between 1 and 3)",ndigit);
			write(ERRORFILE,smess);
			printf("\nPress any key to stop the program");
			wait_a_char();exit(1);
		}

		/*read locus names*/
		for(l=1;l<=m;l++){
			readsfromfile(fp,s,importfile,&line);
			readsfromstring(s,namelocus[l],MAXNOM-1,importfile,line);
 		}

		/*determine the number of individuals*/
		n=0;
		while(fgets(s,SMAX,fp) && sscanf(s,"%s",smess)) n++;
		fclose(fp);
		fp=fopen(importfile,"rt");
		for(l=0;l<=m;l++) fgets(s,SMAX,fp);


		/*write 4 first lines of data file*/
		sprintf(smess,"\" Data file in SPAGeDi format imported from file \"%s\" in FSTAT format",importfile);
		write_tofile_only(inputfile,smess);
		sprintf(smess,"\n%i\t%i\t0\t%i\t%i\t2\n0\nInd\tPop",n,Npop,m,ndigit);
		write_tofile_only(inputfile,smess);
		for(l=1;l<=m;l++){
			sprintf(smess,"\t%s",namelocus[l]);
			write_tofile_only(inputfile,smess);
		}


		/*Read labels*/
		printf("\nEnter the name the file of population labels in Fstat format\n (press RETURN if there is no such file): ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\r'){
			strncpy(labelfile,smess,MAXNOM-1);
			while((fp3=fopen(labelfile,"rt"))==NULL){
				printf("\nWARNING: Cannot open the labels file \"%s\".\nIf it is being used by another application, close it first.\nIf it does not exist in the specified directory, bring it first.\nThen press RETURN.\nPress Ctrl+c twice to stop the programm now.\n",labelfile);
				wait_a_char();
			}
			labels=1;
			for(pop=1;pop<=Npop;pop++){
				readsfromfile(fp3,s,labelfile,&line3);
				readsfromstring(s,popname[pop],MAXNOM-1,labelfile,line3);
			}
			fclose(fp3);
		}
		else labels=0;


		printf("\nConverting the file \"%s\" into \"%s\". Please, wait.",importfile,inputfile);

		/*read genotypes*/
		npop=1;
		for(i=1;i<=n;i++){
			fgets(s,SMAX,fp);
			/*get pop*/
			sscanf(s,"%d",&nump);
			if(i==1){
				numpop[1]=nump;
				if(!labels){sprintf(smess,"%i",nump); sscanf(smess,"%s",popname[1]);}
			}
			for(pop=npop;pop>=1;pop--) if(nump==numpop[pop]) break;
			if(pop<1){
				npop++;
				numpop[npop]=nump;
				if(!labels){sprintf(smess,"%i",nump); sscanf(smess,"%s",popname[npop]);}
				pop=npop;
			}
			/*get genotypes*/
			s2=&s[0];
			for(l=1;l<=m;l++){
				s2=strpbrk(s2,"0123456789");
				s2=strpbrk(s2," \t");
				sscanf(s2,"%d",&genotypel[l]);
			}

			/*write info in data file*/
			sprintf(smess,"\n%i\t%s",i,popname[pop]);
			write_tofile_only(inputfile,smess);
			for(l=1;l<=m;l++){
				sprintf(smess,"\t%i",genotypel[l]);
				write_tofile_only(inputfile,smess);
			}
		}
		sprintf(smess,"\nEND");
		write_tofile_only(inputfile,smess);

		printf("\n\nData from file \"%s\" succesfully imported into file \"%s\"",importfile,inputfile);
		printf("\n press RETURN to go on with SPAGeDi analyses (without defined distance classes)\n or press any other key to stop the program and edit the data file");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0') exit(0);											 
		printf("\n");
			
	}




	/*read GENEPOP data file*/
	if(format==2 || format==3){
		printf("\nConverting the file \"%s\" into \"%s\". Please, wait.",importfile,inputfile);
		/*read 1st line (comments)*/
		readsfromfile(fp,s,importfile,&line);
		/*write 2 comment lines*/
		sprintf(smess,"\" Data file in SPAGeDi format imported from file \"%s\" in GENEPOP format",importfile);
		write_tofile_only(inputfile,smess);
		sprintf(smess,"\n/* %s",s);
		write_tofile_only(inputfile,smess);

		/*read locus names*/
		readsfromfile(fp,s,importfile,&line);
		if(strchr(s,',')){	/*locus names written on a line and separated by comma*/
			s2=s3=&s[0];
			for(l=1;l<=MMAX;l++){
 				s2=strchr(s2,',');
				if(s2){
					*s2='\0';
					readsfromstring(s3,namelocus[l],MAXNOM-1,importfile,line);
					s3=s2+1;
				}
				else{
					readsfromstring(s3,namelocus[l],MAXNOM-1,importfile,line);
					m=l;
					readsfromfile(fp,s,importfile,&line);
					break;
				}
			}
		}
		else for(l=1;l<=MMAX;l++){ /*locus names written on successive lines*/
			readsfromstring(s,namelocus[l],MAXNOM-1,importfile,line);
			readsfromfile(fp,s,importfile,&line);
			if(!strncmp(s,"POP",3) || !strncmp(s,"Pop",3) || !strncmp(s,"pop",3)){ m=l; break;}
 		}

		/*determine the number of individuals and populations and read pop names*/
		n=0;
		Npop=1;
		for(i=1;i<=NMAX;i++){
			strcpy(smess,s);  
			if(!fgets(s,SMAX,fp)){
				s2=strrchr(smess,',');
				if(s2) *s2='\0';
				strncpy(popname[Npop],smess,MAXNOM-1);
				break;
			}
			if(!sscanf(s,"%s",sbis) || !strncmp(s,"POP",3) || !strncmp(s,"Pop",3) || !strncmp(s,"pop",3)){
				s2=strrchr(smess,',');
				if(s2) *s2='\0';
				strncpy(popname[Npop],smess,MAXNOM-1);
				if(!sscanf(s,"%s",sbis)) break;
				Npop++;
			}
			else n++;
		}
		fclose(fp);
		fp=fopen(importfile,"rt");
		do{ 
			fgets(s,SMAX,fp);
		}while(!(!strncmp(s,"POP",3) || !strncmp(s,"Pop",3) || !strncmp(s,"pop",3)));

		/*write format numbers and locus names of data file*/
		if(Npop==n) sprintf(smess,"%i\t0\t2\t%i\t%i\t2\n0\nInd\tX\tY",n,m,format);
		else sprintf(smess,"%i\t%i\t0\t%i\t%i\t2\n0\nInd\tPop",n,Npop,m,format);
		write_tofile_only(inputfile,smess);
		for(l=1;l<=m;l++){
			sprintf(smess,"\t%s",namelocus[l]);
			write_tofile_only(inputfile,smess);
		}


		/*read genotypes*/
		npop=1;
		for(i=1;i<=n;i++){
			fgets(s,SMAX,fp);
			if(!strncmp(s,"POP",3) || !strncmp(s,"Pop",3) || !strncmp(s,"pop",3)){ npop++; fgets(s,SMAX,fp);}
			s2=&s[0];
			/*get ind name*/
			s2=strrchr(s2,',');
			*s2='\0';
			s2++;
			strncpy(indname,s,MAXNOM-1);

			/*get genotypes*/
			s2=strpbrk(s2,"0123456789");
			sscanf(s2,"%d",&genotypel[1]);
			for(l=2;l<=m;l++){
				s2=strpbrk(s2," ");
				s2=strpbrk(s2,"0123456789");
				sscanf(s2,"%d",&genotypel[l]);
			}

			/*write info in data file*/
			if(Npop==n){
				sscanf(indname,"%le%le",&X,&Y);
				sprintf(smess,"\n%i\t%G\t%G",i,X,Y);
			}
			else sprintf(smess,"\n%s\t%s",indname,popname[npop]);
			write_tofile_only(inputfile,smess);
			for(l=1;l<=m;l++){
				sprintf(smess,"\t%i",genotypel[l]);
				write_tofile_only(inputfile,smess);
			}
		}
		sprintf(smess,"\nEND");
		write_tofile_only(inputfile,smess);

		printf("\n\nData from file \"%s\" succesfully imported into file \"%s\"",importfile,inputfile);
		printf("\n press RETURN to go on with SPAGeDi analyses (without defined distance classes)\n or press any other key to stop the program and edit the data file");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0') exit(0);											 
		printf("\n");
			
	}
	
	fclose(fp);
	fclose(fp2);



}/*end of import_data_file*/


/*************************************************************************************/
void export_data_file(int StatType,char inputfile[],int n,
	double *xi,double *yi,double *zi,double *xp,double *yp,double *zp,double **Mdij,
	int *popi,int Npop,int m,int ndigit,int *Nallelel,int ***gilc,int **allelesizela,
	struct name namei[],struct name namepop[],char namelocus[][MAXNOM])
{
	FILE *fp,*fp2,*fp3;
	int line=0,line3=0,flag;
	int format,i,check;
	int Maxallele,l;
	int pop,labels;
	int genotypel[201];
	char s[SMAX], *s2,*s3,smess[SMAX],sbis[SMAX];
	char ch,exportfile[50],labelfile[50];


	//export for PATRI
/*	while((fp=fopen("inPATRIb.txt","a"))==NULL){
		printf("\nWARNING: Cannot open data file inPATRI.txt.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n");
		wait_a_char();
	}

	for(i=1;i<=n;i++) if(strcmp(namepop[popi[i]].n,"M")==0){
		fprintf(fp,"\n  males in group %i = 1",i-1);
	}
	for(i=1;i<=n;i++) if(strcmp(namepop[popi[i]].n,"M")==0){
		fprintf(fp,"\ngroup %i = {\nm%i\n}\n",i-1,i);
	}
	for(i=1;i<=n;i++) if(strcmp(namepop[popi[i]].n,"M")==0){
		fprintf(fp,"\n\nmale: %s",namei[i].n);
		for(l=1;l<=m;l++) fprintf(fp,"\n%i\t%i",allelesizela[l][gilc[i][l][0]],allelesizela[l][gilc[i][l][1]]);
	}
	for(i=1;i<=n;i++) if(strcmp(namepop[popi[i]].n,"M")==0){
		fprintf(fp,"\n\nfemale: %s",namei[i].n);
		for(l=1;l<=m;l++) fprintf(fp,"\n%i\t%i",allelesizela[l][gilc[i][l][0]],allelesizela[l][gilc[i][l][1]]);
	}
	for(i=1;i<=n;i++) if(strcmp(namepop[popi[i]].n,"M")!=0){
		fprintf(fp,"\n\noffspring: %s",namei[i].n);
		fprintf(fp,"\nmother: %s",namepop[popi[i]].n);
		for(l=1;l<=m;l++) fprintf(fp,"\n%i\t%i",allelesizela[l][gilc[i][l][0]],allelesizela[l][gilc[i][l][1]]);
	}

	fclose(fp);
	exit(0);

*/



	/*define the export data file name*/
	format=0;
	do{
		printf("\n\nExport the data file \"%s\" in\n 1- GENEPOP format\n 2- FSTAT format\n",inputfile);
		fgets_chomp(smess, sizeof(smess), stdin);
		sscanf(smess,"%i",&format);
	}while(format<1 && format>2);

	if(format==1) printf("\nEnter the name of the new data file in format GENEPOP (without ext): ");
	if(format==2) printf("\nEnter the name of the new data file in format FSTAT (with ext .dat): ");
	fgets_chomp(exportfile, sizeof(exportfile), stdin);
	do{
		check=0;
		if((fp=fopen(exportfile,"r"))){
			fclose(fp);
			printf("\nWARNING: File \"%s\" already exists.\n  Enter 'e' to erase the present content of the file\n  Otherwise, enter a new file name\n",exportfile);
			fgets_chomp(smess, sizeof(smess), stdin);
			if(!strcmp(smess,"e") || !strcmp(smess,"E")) {
				if((fp=fopen(exportfile,"wt"))!=NULL) fclose(fp);
			}
			else {strcpy(exportfile,smess);check=1;}
		}
	}while(check);

	while((fp=fopen(exportfile,"a"))==NULL){
		printf("\nWARNING: Cannot open data file \"%s\".\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",exportfile);
		wait_a_char();
	}

	if(format==2){/*write Fstat data file*/
		Maxallele=0;
		for(l=1;l<=m;l++) if(Maxallele<allelesizela[l][Nallelel[l]]) Maxallele=allelesizela[l][Nallelel[l]];

		if(StatType>=2) fprintf(fp,"%i  %i  %i  %i",Npop,m,Maxallele,ndigit);
		if(StatType==1) fprintf(fp,"1  %i  %i  %i",m,Maxallele,ndigit);
		for(l=1;l<=m;l++) fprintf(fp,"\n%s",namelocus[l]);

		for(i=1;i<=n;i++){
			if(StatType>=2) fprintf(fp,"\n%i ",popi[i]);
			if(StatType==1) fprintf(fp,"\n1 ");
			if(ndigit==1) for(l=1;l<=m;l++) fprintf(fp," %.1d%.1d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
			if(ndigit==2) for(l=1;l<=m;l++) fprintf(fp," %.2d%.2d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
			if(ndigit==3) for(l=1;l<=m;l++) fprintf(fp," %.3d%.3d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
		}


		if(StatType>=2){
			printf("\nEnter the name of the file with population labels (with ext .lab)\n (or press RETURN if you don't want such a file): ");
			fgets_chomp(smess, sizeof(smess), stdin);
			if(smess[0]!='\0'){
				strcpy(labelfile,smess);
  				do{
					check=0;
					if((fp2=fopen(labelfile,"r"))){
						fclose(fp2);
						printf("\nWARNING: File \"%s\" already exists.\n  Enter 'e' to erase the present content of the file\n  Otherwise, enter a new file name\n",labelfile);
						fgets_chomp(smess, sizeof(smess), stdin);
						if(!strcmp(smess,"e") || !strcmp(smess,"E")) {
							if((fp2=fopen(labelfile,"wt"))!=NULL) fclose(fp2);
						}
						else {strcpy(labelfile,smess);check=1;}
					}
				}while(check);

				fp2=fopen(labelfile,"a");
				fprintf(fp2,"%s",namepop[1].n);
				for(pop=2;pop<=Npop;pop++) fprintf(fp2,"\n%s",namepop[pop].n);

				fclose(fp2);
			}
		}
	}



	if(format==1){/*write Genepop data file*/

		fprintf(fp,"Data file imported from \"%s\" (SPAGeDi format)",inputfile);
		for(l=1;l<=m;l++) fprintf(fp,"\n%s",namelocus[l]);

		if(StatType>=2) for(pop=1;pop<=Npop;pop++){
			fprintf(fp,"\nPop");
			for(i=1;i<=n;i++) if(popi[i]==pop){
				fprintf(fp,"\n%s , ",namei[i].n);
				if(ndigit<=2) for(l=1;l<=m;l++) fprintf(fp," %.2d%.2d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
		/*		if(ndigit==3) for(l=1;l<=m;l++) fprintf(fp," %.3d%.3d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);	   */
				if(ndigit==3) for(l=1;l<=m;l++) fprintf(fp," %.2d%.2d",gilc[i][l][1],gilc[i][l][0]);
			}
		}
		if(StatType==1){
			for(i=1;i<=n;i++){
				fprintf(fp,"\nPop\n%G %G , ",xi[i],yi[i]);
				if(ndigit<=2) for(l=1;l<=m;l++) fprintf(fp," %.2d%.2d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
				if(ndigit==3) for(l=1;l<=m;l++) fprintf(fp," %.3d%.3d",allelesizela[l][gilc[i][l][1]],allelesizela[l][gilc[i][l][0]]);
			}
		}
	}


	fclose(fp);
}
/*************************************************************************************/
void readbasicinfoF(char *inputfile,int *n,int *ncat,int *ncoord,int *mp,int *ndigit,int *ploidy,
					char catname[MAXNOM],char namecoord[] [MAXNOM],char namelocusp[] [MAXNOM],
					int *nc,double *maxc)
{
	FILE *fp;
	int i,flag;
	char s[SMAX], *s2,smess[SMAX],ch;
	int line=0;

	while((fp=fopen(inputfile,"rt"))==NULL){
		printf("\nWARNING: Cannot open data file %c%s%c.\nIf it is being used by another application, close it first.\nIf it does not exist in the specified directory, bring it first.\nThen press RETURN.\nPress Ctrl+c twice to stop the programm now.\n",'"',inputfile,'"');
		wait_a_char();
	}

	/*read first non comment line*/
	readsfromfile(fp,s,inputfile,&line);	
	while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  //skip over empty lines or lines beggining with   /* or // or "  

	flag=sscanf(s,"%d%d%d%d%d%d",n,ncat,ncoord,mp,ndigit,ploidy);
	if(flag!=6) {
		sprintf(smess,"\n\nWARNING: cannot read the first line containing the 6 format numbers (number of ind, number of categories, number of coordinates, number of loci, number of digits per allele, ploidy) on line %d from file %s",line,inputfile);
		write(ERRORFILE,smess);
		sprintf(smess,"\nVerify in the data file that comment lines begin with the characters   //");
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	if(*n > NMAX || *n < 2) {
		sprintf(smess,"\nNumber of individuals=%d higher than max allowed=%d, or less than 2",*n,NMAX);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	if(*ncoord > NCOORDMAX || (*ncoord < 0 && *ncoord != -2) ) {
		sprintf(smess,"\nNumber of spatial coordinates, %i, not in the range admitted (0 to 3 or -2)",*ncoord);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	if(*mp > MMAX || *mp < 1) {
		sprintf(smess,"\nNumber of loci=%d higher than max allowed=%d, or less than 1",*mp,MMAX);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	if(*ndigit > NDIGITMAX /*|| *ndigit < 1*/) {
		sprintf(smess,"\nNumber of digits per allele=%d higher than max allowed=%d"/*, or less than 1*/,*ndigit,NDIGITMAX);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}

	if(*ploidy > PLOIDYMAX || *ploidy < 1) {
		sprintf(smess,"\nPloidy=%d higher than max allowed=%d, or less than 1",*ploidy,PLOIDYMAX);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}

	if(*ndigit<=0 && *ploidy != 2) {
		sprintf(smess,"\nPloidy=%d, different from 2. When the number of digits per allele <= 0, data are assumed to hold for a dominant marker in a diploid. Hence, ploidy must be set to 2",*ploidy);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}

	 
	/*read second non comment line*/
	readsfromfile(fp,s,inputfile,&line);	
	while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  /*skip over lines beggining with comments   */

	flag=sscanf(s,"%d",nc);
	if(flag!=1) {
		sprintf(smess,"\nWARNING: cannot read number of distance intervals on line %d from file %s\nIf you do not need distance intervals, put 0 in the data file.",line,inputfile);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	if(abs(*nc) >= MAXINTERVALS) {
		sprintf(smess,"\nNumber of intervals=%d is higher than max allowed=%d",abs(*nc),MAXINTERVALS);
		write(ERRORFILE,smess);
		printf("\nPress any key to stop the program");
		wait_a_char();exit(1);
	}
	maxc[0]=0.0;
	if((*nc)>0){
		s2=nexttab(s,inputfile,line);
		readdoublefromstring(s2,&maxc[1],inputfile,line);
		for(i=2;i<=(*nc);i++) {
			s2=nexttab(s2,inputfile,line);
			readdoublefromstring(s2,&maxc[i],inputfile,line);
 		}  /*end of for, loop for max values of intervals*/
	}

	readsfromfile(fp,s,inputfile,&line);	/*read third non comment line*/
	while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  /*skip over lines beggining with comments   */

	s2=&s[0];
	if(*ncat){
		s2=nexttab(s2,inputfile,line);
		readsfromstring(s2,catname,MAXNOM-1,inputfile,line);
	}
	for(i=1;i<=abs(*ncoord);i++) {
		s2=nexttab(s2,inputfile,line);
		readsfromstring(s2,namecoord[i],MAXNOM-1,inputfile,line);
 	}  /*end of for, loop for loci names*/
	for(i=1;i<=*mp;i++) {
		s2=nexttab(s2,inputfile,line);
		readsfromstring(s2,namelocusp[i],MAXNOM-1,inputfile,line);
 	}  /*end of for, loop for loci names*/

	fclose(fp);

} /*end of readbasicinfoF*/


/***********************************************************************/

void readsecondinfoF(char *inputfile,int n,struct name namei[],int Ncat,struct name namecat[],
					 struct name namecati[],int *cati,int ncoord,double *xi,double *yi,double *zi,
					 int m,int ndigit,int ploidy,int ***gilc,int *Nallelel,int **allelesizela,
					 int *ploidyi,int Nip[])
{
	FILE *fp;
	int c,i,l,line=0,z,ncatobs,newcat,k,g,nsg;
	int gentemp,maxploidy,ploidyl;
	char s[SMAX], *s2,smess[SMAX],ch,*ptr;
	char newind[SMAX];
	int nallele,ok,coordproblem=0;
	float allelesizea[1000];
	double PI=3.14159265358979;
	
	while((fp=fopen(inputfile,"rt"))==NULL){
		printf("\nWARNING: Cannot open data file %c%s%c.\nIf it is being used by another application, close it first.\nIf it does not exist in the specified directory, bring it first.\nThen press RETURN.\nPress Ctrl+c twice to stop the programm now.\n",'"',inputfile,'"');
		wait_a_char();
	}

	/*skip over 3 first non comment lines*/
	for(i=0;i<3;i++){
		readsfromfile(fp,s,inputfile,&line);  
		while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  /*skip over lines beggining with comments   */
	}

	for(i=1;i<=n;i++) {  /*loop over total number of individuals in data file*/
		readsfromfile(fp,s,inputfile,&line);  /*read next line*/
		while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  /*skip over lines beggining with comments   */
		s2=&s[0];  

		/*read ind name*/
		readsfromstring(s2,newind,MAXNOM-1,inputfile,line);
		if(!strncmp(newind,"END",MAXNOM-1)) {	/*if individual is named END*/
			sprintf(smess,"\nEND message encountered after only %d individuals on line %d from file %s",i-1,line,inputfile);
			write(ERRORFILE,smess);
			sprintf(smess,"\nThe following individuals were successfully found");
			write_tofile_only(ERRORFILE,smess);
			for(z=1;z<=i;z++) {
				sprintf(smess,"\n%s",namei[z].n);
				write_tofile_only(ERRORFILE,smess);
			}
			printf("\nCheck in file %s the list of individuals successfully found.\nPress any key to stop the program.",ERRORFILE);
			wait_a_char();exit(1);					
		} 
		else strncpy(namei[i].n,newind,MAXNOM-1);
		
		/*read indiv category*/
		if(Ncat){  
			s2=nexttab(s2,inputfile,line);
			readsfromstring(s2,namecati[i].n,MAXNOM-1,inputfile,line);
		}
		else strncpy(namecati[i].n,"0",MAXNOM-1);

		/*indiv coordinates*/
		if(abs(ncoord)>0){
			s2=nexttab(s2,inputfile,line);
			if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading spatial coordinates on line %d of the data file: \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now",line);wait_a_char();exit(1);} 
			readdoublefromstring(s2,&xi[i],inputfile,line);
			if(ncoord==-2){
				if(xi[i]>90. || xi[i]<-90.) coordproblem=1;
				xi[i]=xi[i]/360.*2.*PI; //transformation into radian
			}
		}
		else xi[i]=0.;
		if(abs(ncoord)>1){
			s2=nexttab(s2,inputfile,line);
			if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading spatial coordinates on line %d of the data file: \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now",line);wait_a_char();exit(1);} 
			readdoublefromstring(s2,&yi[i],inputfile,line);
			if(ncoord==-2){
				if(yi[i]>180. || yi[i]<-180.) coordproblem=1;				
				yi[i]=yi[i]/360.*2.*PI; //transformation into radian
			}
		}
		else yi[i]=0.;
		if(abs(ncoord)>2){
			s2=nexttab(s2,inputfile,line);
			if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading spatial coordinates on line %d of the data file: \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now",line);wait_a_char();exit(1);} 
			readdoublefromstring(s2,&zi[i],inputfile,line);
		}
		else zi[i]=0.;
		
		/*read indiv genotype*/
		/*routine to read genotype in case of dominant markers*/
		if(ndigit<=0)for(l=1;l<=m;l++) {	/*loop loci*/
			s2=nexttab(s2,inputfile,line);   
			readintfromstring(s2,&gentemp,inputfile,line);
			if(ndigit==0){
				if(gentemp==0 || gentemp==1 || gentemp==2) gilc[i][l][0]=gentemp;
				else {sprintf(smess,"\nGenotype at locus %i from ind %i unexpected: it should be 0 for a missing value,\n1 for a recessive genotype, or 2 for a dominant genotype",l,i); write(ERRORFILE,smess); printf("\nPress any key to stop the program");wait_a_char();exit(551);}
			}
			else{
				if(gentemp==0) gilc[i][l][0]=1;			 /*recessive genotype*/
				else if(gentemp==1) gilc[i][l][0]=2;	 /*dominant genotype*/
				else if(gentemp==ndigit) gilc[i][l][0]=0;/*missing value*/
				else {sprintf(smess,"\nGenotype at locus %i from ind %i unexpected: it should be %i for a missing value,\n0 for a recessive genotype, or 1 for a dominant genotype",l,i,ndigit); write(ERRORFILE,smess); printf("\nPress any key to stop the program");wait_a_char();exit(552);}
			}
		}  /*end of for loop over loci*/

		/*this routine allows to read codominant genotypes where information per chromosome 
			is either ndigit long and without space between chromosomes, either
			no more than ndigit long and with any non numerical (except tab) 
			character to separate chromosomes*/
		else for(l=1;l<=m;l++) {	/*loop loci*/
			s2=nexttab(s2,inputfile,line);   
			ptr=strpbrk(s2,"\t\n\r");	 /*pointer to the first tab or return encountered in s2*/
			ptr--;
			if(isdigit(*ptr)==0){sprintf(smess,"\nGenotype at locus %i from ind %i unreadable (it ends with a non-numerical character)",l,i); write(ERRORFILE,smess); printf("\nPress any key to stop the program");wait_a_char();exit(55);}
			else for(c=0;c<ploidy;c++){
				k=1;
				while(isdigit(*(ptr-1)) && k<ndigit){ptr--;k++;};
				strncpy(smess,ptr,k);
				smess[k]='\0';
				readintfromstring(smess,&gentemp,inputfile,line);
				gilc[i][l][c]=gentemp;
				ptr--;
				while(isdigit(*ptr)==0 && (*ptr)!='\t') ptr--;
				if((*ptr)=='\t' && c<(ploidy-1)){
					for(k=c+1;k<ploidy;k++) gilc[i][l][k]=0;
					break;
				}
			} /*end of for loop over chromosomes*/
		}  /*end of for loop over loci*/
	}	/*end of for loop over individuals*/

	readsfromfile(fp,s,inputfile,&line);  /*read next line*/
	while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,inputfile,&line);  /*skip over lines beggining with comment   */

	readsfromstring(s,newind,MAXNOM-1,inputfile,line);
	if(strncmp(newind,"END",MAXNOM-1))  {	/*if END not encoutered after n individuals*/
		sprintf(smess,"\nEND message not encountered after %d individuals on line %d from file %s",n,line,inputfile);
		write(ERRORFILE,smess);
		write(ERRORFILE,"\nCheck the number of individuals in the data file");
		sprintf(smess,"\nThe following individuals were successfully found");
		write_tofile_only(ERRORFILE,smess);
		for(z=1;z<=n;z++) {
			sprintf(smess,"\n%s",namei[z].n);
			write_tofile_only(ERRORFILE,smess);
		}
		printf("\nCheck in the file %s the list of individuals successfully found.\nPress any key to stop the program.",ERRORFILE);
		wait_a_char();exit(1);					
	}
	if(coordproblem){
		write(ERRORFILE,"\nERROR: some coordinates are out of range for latitude (-90 to 90) or for\n\tlongitude (-180 to 180)");
		wait_a_char();exit(1);
	}



	
	/*check the ploidy of each individual, count the number of each ploidy (0 for individual without any valid genotype, ploidy+1 for individuals with different apparent ploidy levels according to locus) and count the # of ind with missing and incomplete genotypes for each locus*/
	for(g=0;g<=ploidy+1;g++) Nip[g]=0;
	for(i=1;i<=n;i++){
		maxploidy=ploidyi[i]=0;
		for(l=1;l<=m;l++){
			ploidyl=0;
			for(c=ploidy-1;c>=0;c--) if(gilc[i][l][c]) {ploidyl=c+1; break;}
			if(ploidyl){
				if(maxploidy && maxploidy!=ploidyl) ploidyi[i]=ploidy+1;
				else ploidyi[i]=maxploidy=ploidyl;
			}
		}
		Nip[ploidyi[i]]++;
	}

	/*replace allelesize by allele number (1,2,...) in gilc[][][] where numbers follow the order of allele lenght*/
	if(ndigit>0)for(l=1;l<=m;l++){
		nallele=0;
		for(i=1;i<=n;i++)for(c=0;c<ploidy;c++)if(gilc[i][l][c]){
			ok=0;
			for(k=1;k<=nallele;k++)if(gilc[i][l][c]==allelesizea[k]) ok=1;
			if(!ok){
				nallele++;
				allelesizea[nallele]=(float)gilc[i][l][c];
			}
		}
		shellsort(nallele,allelesizea);
		for(i=1;i<=n;i++)for(c=0;c<ploidy;c++)if(gilc[i][l][c]){
			for(k=1;k<=nallele;k++)if(gilc[i][l][c]==allelesizea[k]) gilc[i][l][c]=k;
		}
		Nallelel[l]=nallele;
		for(k=1;k<=nallele;k++) allelesizela[l][k]=(int)allelesizea[k];
		allelesizela[l][0]=0;
	}
	else for(l=1;l<=m;l++) Nallelel[l]=2;

	/*max # alleles*/
	Nallelel[0]=0;
	for(l=1;l<=m;l++) if(Nallelel[l]>Nallelel[0]) Nallelel[0]=Nallelel[l];

	fclose(fp);

}	/*end of readsecondinfoF*/


/****************
											displaybasicinfo()

		This function display the basic informations read from file data file by readbasicinfo.
	If errorcount is different from 0 then display the number of errors detected.

********/

void displaybasicinfoF(int argc,char *inputfilename,char *outputfilename,
					int n,int Ncat,int *Nik,int ncoord,int m,int ndigit,
					int ploidy,int *ploidyi,int Nip[],char namelocus[] [MAXNOM],
					char namecoord[] [MAXNOM],struct name *namecat,
					int nc,double *maxc,int Nsg,int *Nig,int Nskg,int *Niskg)
{
	int i,k,g,p;
	int Nikmin,Nikmax,Nigmin,Nigmax,Niskgmin,Niskgmax,Nk1i,Ng1i,Nskg1i;
	char ch;
	char smess[SMAX];
	FILE *fp;

	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file \"%s\".\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",outputfilename);
		wait_a_char();
	}
	fclose (fp);

	/*check min-max # ind per group + # of groups with only 1 ind*/
	Nigmin=Nigmax=Nig[1];
	Ng1i=0;
	for(g=1;g<=Nsg;g++){
		if(Nig[g]<Nigmin) Nigmin=Nig[g];
		if(Nig[g]>Nigmax) Nigmax=Nig[g];
		if(Nig[g]==1) Ng1i++;
	}
	Nikmin=Nikmax=Niskgmin=Niskgmax=Nk1i=Nskg1i=0;
	if(Ncat){
		Nikmin=Nikmax=Nik[1];
		for(k=1;k<=Ncat;k++) {
			if(Nik[k]<Nikmin) Nikmin=Nik[k]; 
			if(Nik[k]>Nikmax) Nikmax=Nik[k];
			if(Nik[k]==1) Nk1i++;
		}
	}
	if(Ncat>1 && Nsg>1){
		Niskgmin=Niskgmax=Niskg[1];
		for(g=1;g<=Nskg;g++){
			if(Niskg[g]<Niskgmin) Niskgmin=Niskg[g];
			if(Niskg[g]>Niskgmax) Niskgmax=Niskg[g];
			if(Niskg[g]==1) Nskg1i++;
		}
	}



	sprintf(smess,"\n%s - a program for Spatial Pattern Analysis of Genetic Diversity\n             Written by Olivier Hardy & Xavier Vekemans\n",VERSION);
	write_tofile_only(outputfilename,smess);
	printf("\n\n\n\n\n\nBASIC INFORMATION\n\nCheck the following informations from the data file \"%s\"\n",inputfilename);
	sprintf(smess,"\nInput file : \"%s\"",inputfilename);
	write_tofile_only(outputfilename,smess);
	sprintf(smess,"\nResults file : \"%s\"",outputfilename);
	write_tofile_only(outputfilename,smess);
	sprintf(smess,"\n%i individuals",n);
	write(outputfilename,smess);
	sprintf(smess,"\n%i categories: ",Ncat);
	write(outputfilename,smess);
	for(i=1;i<=Ncat;i++) {
		if(i==1)sprintf(smess,"%s",namecat[i].n); 
		else sprintf(smess,", %s",namecat[i].n);
		write(outputfilename,smess);
	}
	if(ncoord==-2) sprintf(smess,"\nLatitude-longitude coordinates (assumed to be in degrees; used to compute distances in km)"); 
	else sprintf(smess,"\n%i spatial coordinates",ncoord);
	write(outputfilename,smess);
	for(i=1;i<=abs(ncoord);i++) {
		if(i==1) sprintf(smess,": %s",namecoord[i]);
		sprintf(smess,", %s",namecoord[i]);
		write(outputfilename,smess);
	}
	sprintf(smess,"\n%i loci: %s",m,namelocus[1]);
	write(outputfilename,smess);
	for(i=2;i<=m;i++) {
		sprintf(smess,", %s",namelocus[i]);
		write(outputfilename,smess);
	}
	if(ndigit>0) sprintf(smess,"\n%i digits per allele",ndigit);
	else if(ndigit<0) sprintf(smess,"\nDominant marker: missing data=%i, recessive genotype=0, dominant genotype=1",ndigit);
	else sprintf(smess,"\nDominant marker: missing data=0, recessive genotype=1, dominant genotype=2");
	write(outputfilename,smess);
	sprintf(smess,"\nMax ploidy announced: %i",ploidy);
	write(outputfilename,smess);
	if(ndigit>0){ 
		sprintf(smess,"\n# indiv within each ploidy:  %ix=%i",ploidy,Nip[ploidy]);
		write(outputfilename,smess);
		for(p=ploidy-1;p>=1;p--)if(Nip[p]){sprintf(smess,",   %ix=%i",p,Nip[p]);write(outputfilename,smess);}
	}
	if(Nip[0]){	sprintf(smess,"\n  WARNING: There is(are) %i individual(s) without genotype at any locus",Nip[0]);write(outputfilename,smess); }
	if(ndigit>0)if(Nip[ploidy+1]){	
		sprintf(smess,"\n  WARNING: There is(are) %i individual(s) showing different ploidy levels according to the locus: indiv n° ",Nip[ploidy+1]);
		write(outputfilename,smess); 
		for(i=1;i<=n;i++) if(ploidyi[i]==(ploidy+1)){sprintf(smess,"%i  ",i); write(outputfilename,smess);}
		sprintf(smess,"\nYou must first resolve this problem. Note that 0's representing missing alleles of incomplete genotypes must be on the RIGHT, whereas 0's representing no information for individuals with a ploidy level inferior to that announced must be on the LEFT.\nPress any key to stop the program now.");
		write(outputfilename,smess); 
		wait_a_char();
		exit(56);
	}

	sprintf(smess,"\n\n# groups (min - max # indiv per group): ");
	write(outputfilename,smess);
	sprintf(smess,"\n  %i categorical groups",Ncat);
	write(outputfilename,smess);
	if(Ncat>1){sprintf(smess," (%i - %i)",Nikmin,Nikmax); write(outputfilename,smess);}
	if(Nk1i){sprintf(smess,"   WARNING: %i groups with 1 indiv",Nk1i);write(outputfilename,smess);}
	sprintf(smess,"\n  %i spatial groups",Nsg);
	write(outputfilename,smess);
	if(Nsg>1){sprintf(smess," (%i - %i)",Nigmin,Nigmax);write(outputfilename,smess);}
	if(Ng1i){sprintf(smess,"   WARNING: %i groups with 1 indiv",Ng1i);write(outputfilename,smess);}
	sprintf(smess,"\n  %i spatio-categorical groups",Nskg);
	write(outputfilename,smess);
	if(Nskg){sprintf(smess," (%i - %i)",Niskgmin,Niskgmax);write(outputfilename,smess);}
	if(Nskg1i){sprintf(smess,"   WARNING: %i groups with 1 indiv",Nskg1i);write(outputfilename,smess);}

	if(argc<4){
		printf("\n\nPress RETURN to continue, or enter 's' to stop now\n");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='s' || smess[0]=='S') exit(0);
	}

/*	if(nc<0) sprintf(smess,"\n\n%i distance intervals with equal # pairwise comparisons",abs(nc));
	else sprintf(smess,"\n\n%i distance intervals:",nc);
	write(outputfilename,smess);
	if(nc>0){
		sprintf(smess,"\nclass\tmax distance");
		write(outputfilename,smess);
		for(i=1;i<=nc;i++) {
			sprintf(smess,"\n  %d\t%G",i,maxc[i]);
			write(outputfilename,smess);
		}
	}		  */
}    /*end of displaybasicinfoF*/


/**************************************************************************************/

void define_analysisF(int argc,char *instrfile,int n,int ploidy,int ndigit,int m,int Ncat,int Nsg,int Nskg,int ncoord,
					  int *StatType,int *NS,int Stat[],int *TypeComp,int *cat1,int *cat2,
					  struct name *namecat,int *FreqRef,float *givenF,int *writeallelefreq,int *JKest,int *distmatrix,int Npermut[],
					  float *dijmin,float *dijmax,int *writeresampdistr,int *regdetails,int *varcoef,int *Rbtwloc,float *sigmaest,float *density,float *dwidth,
					  int *permutdetails,int *distm,char *inputfilename,char *distfilename,char *freqfilename,int *definealleledist,char *alleledistfilename,
					  int *export,long *seed)
{
	char ch,smess[SMAX];
	int stats,k,k2,Noptions,ok,otherstat=0,S;
	int permutloc=1,permutind=0,permutgenes=0,permutalleles=0,permutgil,reffreqpossible=0;
	float ram;

	FILE *stream;

	// Reassign "stdin" (=console) to the instruction file (instrfile) to read options from that file 
//	if(argc>=4) stream = freopen(instrfile, "r", stdin);


startagain:

	*cat1=*cat2=*FreqRef=0;
	*writeallelefreq=(*definealleledist)=0;
	if(ndigit>0 && (Ncat>0 || Nsg>1)){
		do{
			printf("\n\n\n\nLEVEL OF ANALYSES\n\nSpatial analyses carried out at the level of\n\n INDIVIDUALS (e.g. kinship coefficients,...) .....................1\n\n\n POPULATIONS (e.g. F-statistics,...) where one population is a");
			if(Ncat>0)			 printf("\n\n   categorical group         ->%4i pop ..........................2",Ncat);
			if(Nsg>1 && Ncat<1) printf("\n\n   spatial group             ->%4i pop ..........................2",Nsg);
			if(Nsg>1 && Ncat>0)  printf("\n\n   spatial group             ->%4i pop ..........................3\n\n   spatio-categorical group  ->%4i pop ..........................4",Nsg,Nskg);
			printf("\n");
			fgets_chomp(smess, sizeof(smess), stdin);
			ch = smess[0];
		}while( (ch!='1' && ch!='2' && ch!='3' && ch!='4') || ((ch=='3' || ch=='4')&&(Ncat<1 || Nsg<=1)) ); 
		if(ch=='1') *StatType=1;
		if(ch=='2') *StatType=2;
		if(ch=='3') *StatType=3;
		if(ch=='4') *StatType=4;
		if(ch=='2' && Ncat<1) *StatType=3;/*for the program, StatType always = 1 for individuals, = 2 for categorical groups, = 3 for spatial groups, = 4 for spatial groups of given categories*/
	}
	else *StatType=1; 


	if(*StatType==1){
		do{
			printf("\n\n\n\nSTATISTICS for individual level analyses\nSelect the statistics (you can select several, e.g. enter \"13A\"):");
			if(ndigit>0){ //codominant marker
				printf("\n\n Statistics based on allele identity / non-identity :");
				printf("\n   1- KINSHIP coefficient (Loiselle et al., 1995)");
				printf("\n   2- KINSHIP coefficient (Ritland, 1996)");
//				printf("\n");
				printf("\n   3- RELATIONSHIP coefficient (Moran's I on indiv allele freq)");
				if(ploidy==2){
			 		printf("\n   4- RELATIONSHIP coefficient (Queller and Goodnight, 1989)");		 
					printf("\n   5- RELATIONSHIP coefficient ('r' in Lynch & Ritland, 1999)");
					printf("\n   6- RELATIONSHIP coefficient ('r' in Wang, 2002)");
					printf("\n   7- RELATIONSHIP coefficient (Li et al., 1993)");
					printf("\n   8- FRATERNITY coefficient ('delta' in Lynch & Ritland, 1999)");
					printf("\n   9- FRATERNITY coefficient ('delta' in Wang, 2002)");
				}
				if(ploidy>1) printf("\n   A- Aij: Rousset's distance ('a' in Rousset, 2000)");
				printf("\n\n Statistic based on allele size (cf. microsatellites) :");
				printf("\n   R- Rij: Allele size correlation coefficient (I' in Streiff et al., 1998)");
				printf("\n\n Statistics based on distances between alleles (matrix to define) :");
				printf("\n   N- Nij: kinship coefficient for ordered alleles (OJ Hardy, unpublished)");
			}
			else{//dominant marker
				printf("\n\n 1- KINSHIP coefficient for dominant marker in diploids (Hardy, 2003)");
				printf("\n\n 2- RELATIONSHIP coefficient for dominant marker in diploids (Hardy, 2003)");
			}
			if(Ncat>0 || Nsg>1) printf("\n\nEnter 'R' to return to the first menu");
			printf("\n");

			do{	
				fgets_chomp(smess, sizeof(smess), stdin);
				if(smess[0]=='r' || smess[0]=='R')  goto startagain;
				(*NS)=strlen(smess);
			}while(smess[0]=='\0' || *NS>11 || *NS<1);

			for(k=1;k<=(*NS);k++){
				if(smess[k-1]=='1') Stat[k]=1;
				else if(smess[k-1]=='2') Stat[k]=2;
				else if(smess[k-1]=='3') Stat[k]=3;
				else if(smess[k-1]=='4') Stat[k]=8;
				else if(smess[k-1]=='5') Stat[k]=6;
				else if(smess[k-1]=='6') Stat[k]=9;
				else if(smess[k-1]=='7') Stat[k]=13;
				else if(smess[k-1]=='8') Stat[k]=7;
				else if(smess[k-1]=='9') Stat[k]=10;
				else if(smess[k-1]=='A' || smess[k]=='a') Stat[k]=4;
				else if(smess[k-1]=='R' || smess[k]=='r') Stat[k]=5;
				else if(smess[k-1]=='N' || smess[k]=='n') {Stat[k]=14; *definealleledist=1;}
				else Stat[k]=0;
				if(ndigit<=0){
					if(Stat[k]==1) Stat[k]=11;
					if(Stat[k]==2) Stat[k]=12;
				}
			}									  

			/*remove stat not allowed*/
			for(k=1;k<=(*NS);k++)if(Stat[k]<1 || Stat[k]>14 || (ploidy==1 && Stat[k]==4) || (ploidy!=2 && (Stat[k]==6 || Stat[k]==7 || Stat[k]==8 || Stat[k]==9 || Stat[k]==10 || Stat[k]==13)) || (ndigit<=0 && (Stat[k]>12 || Stat[k]<11))){
				(*NS)--;
				for(k2=k;k2<=(*NS);k2++) Stat[k2]=Stat[k2+1];
			}
		}while((*NS)<1);

		permutgenes=permutalleles=permutind=reffreqpossible=0;
		permutloc=1;
		if(ploidy>1 && ndigit>0)for(k=1;k<=(*NS);k++)if(Stat[k]==1 || Stat[k]==2 || Stat[k]==5) permutgenes=1;
		for(k=1;k<=(*NS);k++)if(Stat[k]==5 || Stat[k]==14) permutalleles=1;
		if(Nsg>1) permutind=1;
		for(k=1;k<=(*NS);k++)if(Stat[k]==1 || Stat[k]==2 || Stat[k]==5 || Stat[k]==6 || Stat[k]==7 || Stat[k]==8 || Stat[k]==9 || Stat[k]==10 || Stat[k]==11 || Stat[k]==12 || Stat[k]==13) reffreqpossible=1;

		if(ndigit<=0){
			do{
				printf("\n\nEnter the assumed inbreeding coefficient : ");
				fgets_chomp(smess, sizeof(smess), stdin);
				ok=sscanf(smess,"%f",givenF);
			}while(ok!=1 || *givenF<=-1. || *givenF>1.);
		}

	}
	
	if(*StatType>=2){
		do{
			otherstat=0;
			printf("\n\n\n\nSTATISTICS for population level analyses\nSelect the statistics (you can select several, e.g. enter \"15\"):");
			printf("\n\n Statistics based on allele identity / non-identity :");
			printf("\n\n   1- global F-statistics and pairwise Fst  (Weir and Cockerham 1984)");
			printf("\n   2- global F-statistics and pairwise Rho  (intra-class relationship coef)");
			printf("\n   3- global Gst and pairwise Gst  (Pons & Petit 1996)");
			printf("\n   4- global Gst and pairwise Gij  (inter-pop kinship coef)");
			printf("\n   5- global F-statistics and pairwise Ds   (Nei's standard distance, 1978)");
			printf("\n\n Statistics based on allele size (cf. microsatellites) :");
			printf("\n\n   6- global R-statistics and pairwise Rst  (ANOVA approach)");
			printf("\n   7- global R-statistics and pairwise dm2 distance (Goldstein 1995)");
			printf("\n\n Statistics based on distances between alleles (matrix to define) :");
			printf("\n\n   8- global Nst and pairwise Nst  (Pons & Petit 1996)");
			printf("\n   9- global Nst and pairwise Nij  (OJ Hardy, unpubl)");
			printf("\n\nEnter 'R' to return to the first menu\n\n");   

			do{	
				fgets_chomp(smess, sizeof(smess), stdin);
				if(smess[0]=='r' || smess[0]=='R')  goto startagain;
				ok=sscanf(smess,"%i",&stats);
				(*NS)=1+(int)log10(stats);
			}while(ok!=1 || *NS>9 || *NS<1);

			/*Chosen statistics, put in Stat[k] with k=1 to #different stat*/
			for(k=(*NS);k>0;k--){
				if(stats%10==1) {Stat[k]=1;	otherstat=1;}
				if(stats%10==2) {Stat[k]=3;	otherstat=1;}
				if(stats%10==3) Stat[k]=9;
				if(stats%10==4) Stat[k]=10;
				if(stats%10==5) Stat[k]=7;
				if(stats%10==6) {Stat[k]=5;	otherstat=1;}
				if(stats%10==7) Stat[k]=8;
				if(stats%10==8) {Stat[k]=11; *definealleledist=1;}
				if(stats%10==9) {Stat[k]=12; *definealleledist=1;}

				stats=(int)(stats/10);
			}
			for(k=1;k<=(*NS);k++)if(Stat[k]<1 || Stat[k]>12){
				(*NS)--;
				for(k2=k;k2<=(*NS);k2++) Stat[k2]=Stat[k2+1];
			}


		}while((*NS)<1);

		permutalleles=permutgenes=0;
		permutind=permutloc=1;
		if(ploidy>1) permutgenes=1;
		for(k=1;k<=(*NS);k++)if(Stat[k]==5 || Stat[k]==6 || Stat[k]==8 || Stat[k]==11 || Stat[k]==12) permutalleles=1;
	}

	if(*definealleledist){
		printf("\n\n\n\nFILE FOR PAIRWISE GENETIC DISTANCES BETWEEN ALLELES\n\nEnter the name of the file containing pairwise distances, \nor press RETURN if it is the data file \"%s\" :  ",inputfilename);
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') strncpy(alleledistfilename,inputfilename,MAXNOM);
		else strncpy(alleledistfilename,smess,MAXNOM);

	}
	
	/*Select computational options*/
	(*JKest)=Npermut[0]=(*distm)=(*TypeComp)=0;
	(*dijmin)=(*dijmax)=0.0f;
	printf("\n\n\n\n\n\n\n\nCOMPUTATIONAL OPTIONS\n\nSelect among the following options (you can select several) or press RETURN :");
	printf("\n\n 1- Use a matrix to define pairwise spatial distances");
	printf("\n\n 2- Make restricted regression analyses (i.e. over limited distance range)");
	printf("\n\n 3- Make permutation tests");
	if(m>1) printf("\n\n 4- Jackknife over loci");
	if(Ncat>1 && (*StatType==1 || *StatType==4)) printf("\n\n 5- Restrict pairwise comparisons within or among (selected) categories");
	if(*StatType>=2 && otherstat) printf("\n\n 6- Pairwise  Fst (or Rst or Rho) provided as Fst/(1-Fst) ratios");
	if(*StatType==1 && reffreqpossible) printf("\n\n 6- Define reference allele frequencies to compute relatedness coefficients");

	printf("\n\nEnter 'R' to return to the first menu\n\n");

	fgets_chomp(smess, sizeof(smess), stdin);
	if(smess[0]=='r' || smess[0]=='R') goto startagain; 
	if(smess[0]!='\0'){
		sscanf(smess,"%d",&stats);

		/*number of different options*/
		Noptions=1+(int)log10(stats);
		for(k=1;k<=Noptions;k++){
			if(stats%10==1) (*distm)=1;
			if(stats%10==2) (*dijmin)=1;
			if(stats%10==3) Npermut[0]=1;
			if(stats%10==4) (*JKest)=1;
			if(stats%10==5) (*TypeComp)=1;
			if(stats%10==6){
				if(*StatType>=2) for(S=1;S<=(*NS);S++)if(Stat[S]==1 || Stat[S]==3 || Stat[S]==5) Stat[S]++;
				if(*StatType==1 && reffreqpossible) *FreqRef=-1;
			}
			stats=(int)(stats/10);
		}
	}

	if(*TypeComp){
		do{
			printf("\n\n\n\n\n\n\n\nTYPE OF PAIRWISE COMPARISONS");
			if(*StatType>1) printf("\n\nSelect the type of pairs of spatio-categorical groups to compute the statistics");
			if(*StatType==1) printf("\n\nSelect the type of individual pairs to compute pairwise statistics:");
			printf("\n\n 1- all pairs (irrespective of categories = default option)");
			printf("\n\n 2- only pairs within category");
			printf("\n\n 3- only pairs among categories");
			printf("\n\n 4- only pairs within a specified category");
			printf("\n\n 5- only pairs between two specified categories\n\n");

			fgets_chomp(smess, sizeof(smess), stdin);
			sscanf(smess,"%c",&ch);
		}while(ch!='1' && ch!='2' &&  ch!='3' && ch!='4' &&  ch!='5');	
		if(ch=='1') *TypeComp=0;
		if(ch=='2') *TypeComp=1;
		if(ch=='3') *TypeComp=2;
		if(ch=='4') {*TypeComp=1; *cat1=1;}
		if(ch=='5') {*TypeComp=2; *cat1=*cat2=1;}

 		if(*cat1 && !*cat2){
			printf("\n");
			do{
				printf("\nGive the name of the category for within category analysis: ");
				fgets_chomp(smess, sizeof(smess), stdin);
				*cat1=0;
				for(k=1;k<=Ncat;k++) if(strcmp(smess,namecat[k].n)==0) *cat1=k;
				if(*cat1==0){
					printf("\n\"%s\" is not recognized as a category name; it must be one of the following name:\n",smess);
					for(k=1;k<=Ncat;k++) printf("   %s",namecat[k].n);
				}
			}while(*cat1==0);
		}
		if(*cat1 && *cat2){
			printf("\n");
			do{
				printf("\nGive the names of the two categories for between categories analysis: ");
				printf("\n  First category: ");
				fgets_chomp(smess, sizeof(smess), stdin);
				*cat1=0;
				for(k=1;k<=Ncat;k++) if(strcmp(smess,namecat[k].n)==0) *cat1=k;
				if(*cat1==0){
					printf("\n\"%s\" is not recognized as a category name; it must be one of the following name:\n",smess);
					for(k=1;k<=Ncat;k++) printf("   %s",namecat[k].n);
				}
				else{
					printf("  Second category: ");
					fgets_chomp(smess, sizeof(smess), stdin);
					*cat2=0;
					for(k=1;k<=Ncat;k++) if(strcmp(smess,namecat[k].n)==0) *cat2=k;
					if(*cat2==0){
						printf("\n\"%s\" is not recognized as a category name; it must be one of the following name:\n",smess);
						for(k=1;k<=Ncat;k++) printf("   %s",namecat[k].n);
					}
					if(*cat1 == *cat2) printf("\nError: the two categories must be different");
				}
			}while(*cat1==0 || *cat2==0 || *cat1==*cat2);
		}

		if(*TypeComp==1 && *StatType==1 && *FreqRef!=-1){
			printf("\n\nSelect the reference allele frequencies to compute pairwise statistics:");
			printf("\n\n 1- whole sample (enter '1' or press RETURN)");
			printf("\n\n 2- sample within category (enter '2')\n");
			fgets_chomp(smess, sizeof(smess), stdin);
			if(smess[0]=='2') *FreqRef=1;
			else *FreqRef=0;
		}

	}/*end of if(TypeComp)*/	
	
	
	/*if no pairwise dist computable, ask again*/
 	if(!ncoord && !(*distm)){
		printf("\n\n\nWARNING: There is no spatial coordinates in the data file. \nPress RETURN to go on without spatial information \nor enter the name of the file containing pairwise distances: ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0'){
			*distm=2;
			strncpy(distfilename,smess,MAXNOM);
		}
		else permutloc=0;
	}

	if(*distm==1){//pairwise spatial distances are defined in a matrix to specify
		printf("\n\n\n\nFILE FOR PAIRWISE SPATIAL DISTANCES\n\nEnter the name of the file containing pairwise distances, \nor press RETURN if it is the data file \"%s\" :  ",inputfilename);
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') strncpy(distfilename,inputfilename,MAXNOM);
		else strncpy(distfilename,smess,MAXNOM);
	}

	if(*FreqRef==-1){//reference allele frequencies given in a file to specify
		printf("\n\n\n\nFILE WITH REFERENCE ALLELE FREQUENCIES\n\nEnter the name of the file with reference allele frequencies, \nor press RETURN if it is the default file name \"freq.txt\" :   ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') strncpy(freqfilename,"freq.txt",MAXNOM);
		else strncpy(freqfilename,smess,MAXNOM);
	}

	/*define dist interval for partial regressions*/
	if(*dijmin){
		printf("\n\n\n\nDISTANCE RANGE FOR RESTRICTED REGRESSION\n\nEnter minimal and maximal distances for restricted spatial regressions\n(or press RETURN for an unbounded minimal or maximal distance)");
		printf("\n\nMinimal distance = ");
		k=0;
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') (*dijmin)=0.;
		else do { k=sscanf(smess,"%f",dijmin);
		}while(k!=1);
		if(k==0) printf("\n");
		printf("\nMaximal distance = ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='\0') (*dijmax)=0.;
		else do { k=sscanf(smess,"%f",dijmax);
		}while(k!=1);
	}

	/*define parameters for permutation tests*/
	(*writeresampdistr)=*permutdetails=*seed=0;
	if(Npermut[0]){
		printf("\n\n\n\n\n\n\n\nPERMUTATION OPTIONS\n\nSelect among the following options (you can select several) or press RETURN");
		
		if(permutalleles){
			printf("\n\n 1- Test of genetic structuring (permuting genes, individuals and/or locations)");
			printf("\n\n 2- Test of mutation effect on genetic structure (permuting alleles)");
			if(*StatType>=2) printf("\n\n 3- Test of mutation effect on genetic differentiation for each population pair");
			printf("\n\n");

			permutalleles=permutgil=0;
			fgets_chomp(smess, sizeof(smess), stdin);
			if(smess[0]!='\0'){	
				sscanf(smess,"%i",&stats);
				Noptions=1+(int)log10(stats);
				for(k=1;k<=Noptions;k++){
					if(stats%10==1) permutgil=1;
					if(stats%10==2 && permutalleles==0) permutalleles=1;
					if(stats%10==3) permutalleles=2;
					stats=(int)(stats/10);
				}
			}
  		}
		else permutgil=1;

		if(permutgil==0) permutloc=permutind=permutgenes=0;
		*permutdetails=2;

		if(permutalleles)printf("\n\nSelect among the following additional options (you can select several) or press RETURN");
		printf("\n\n 1- Report only P-values (otherwise details of permutation tests are reported)");
		if(permutgil) printf("\n\n 2- Define # of permutations for each randomised unit (otherwise same #)");
		printf("\n\n 3- Initialise random number generator (otherwise initialisation on clock)\n\n");

		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0'){	
			sscanf(smess,"%i",&stats);
			Noptions=1+(int)log10(stats);
			for(k=1;k<=Noptions;k++){
				if(stats%10==1) *permutdetails=0;
				if(stats%10==2 && permutgil) Npermut[0]=2;
				if(stats%10==3) *seed=1;
				stats=(int)(stats/10);
			}
		}


		/*define # permutations*/
		if(Npermut[0]==1){
			do{
				printf("\n\nEnter # permutations for tests (min=40, max=%i): ",NRESAMPLE-1);
				fgets_chomp(smess, sizeof(smess), stdin); 
				ok=sscanf(smess,"%i",&Npermut[0]);
			}while(ok!=1 || Npermut[0]>(NRESAMPLE-1) || (Npermut[0]>1 && Npermut[0]<40) || Npermut[0]<0);
			for(k=1;k<=4;k++) Npermut[k]=Npermut[0];
		}
		if(Npermut[0]==2){
			printf("\n\nEnter # permutations for each randomisation unit (>=40 & <=%i, 0 = no test):\n",NRESAMPLE-1);
			for(k=1;k<=4/*3*/;k++){
				if((!permutloc && k==1) || (!permutind && k==2) || (!permutgenes && k==3) || (!permutalleles && k==4)) continue;
				do{
					if(k==1)printf("\n LOCATION permutations (tests of spatial structure): ................ ");
					if(k==2 && (*StatType)==1) printf("\n INDIVIDUAL permutations (tests of intra-group coefficients):........ ");
					if(k==2 && (*StatType)>=2) printf("\n INDIVIDUAL permutations (tests of population differentiation):...... ");
					if(k==3)printf("\n GENE permutations (tests of inbreeding coefficients): .............. ");
					if(k==4)printf("\n ALLELE permutations (tests of mutation effect of genetic structure): ");
					fgets_chomp(smess, sizeof(smess), stdin); 
					ok=sscanf(smess,"%i",&Npermut[k]);
				}while(ok!=1 || abs(Npermut[k])>(NRESAMPLE-1) || (abs(Npermut[k])>1 && abs(Npermut[k])<40));
				if(abs(Npermut[k])>Npermut[0]) Npermut[0]=abs(Npermut[k]);
			}
		}
		if(!permutloc) Npermut[1]=0;
		if(!permutind) Npermut[2]=0;
		if(!permutgenes) Npermut[3]=0;
		if(!permutalleles) Npermut[4]=0;
		if(permutalleles==2) Npermut[4]*=-1; //negative # if tests for each pair of pop required 

		/*random number generators initialisation*/ 
		if(*seed){
			do{
				printf("\n\nEnter a number for the initial seed: ");
				fgets_chomp(smess, sizeof(smess), stdin);
				ok=sscanf(smess,"%li",seed);
				*seed=abs(*seed)*(-1);
			}while(ok!=1);
		}
		else{	/*initialisation on local time*/
			*seed=timeseed();
			*seed=(*seed)*(-1);
		}


/*		printf("Enter 'y' if you wish the values of the statistics for each permutation ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='y') *writeresampdistr=1;
		if(*writeresampdistr) printf("\nThese values will be written in a file called ""ResampDistr.txt"".\nWarning: if this file already exists, its content will be erased.");
*/	
	}	
   	

	/*Select output options*/
	(*distmatrix)=(*regdetails)=(*writeallelefreq)=(*varcoef)=(*export)=(*Rbtwloc)=0;
	(*sigmaest)=0.0f;
	printf("\n\n\n\n\n\n\n\nOUTPUT OPTIONS\n\nSelect among the following options (you can select several) or press RETURN :");
	if(*StatType==1 && Ncat) printf("\n\n 1- Report allele freq and diversity coef per categ (or only averages reported)");
	if(*StatType>=2) printf("\n\n 1- Report allele freq and diversity coef per pop (or only averages reported)");
	printf("\n\n 2- Report all stat of regression analyses (or only slopes reported)");
	printf("\n\n 3- Report matrices with pairwise spatial distances and genetic coefficients");
	if(ploidy==2) printf("\n\n 4- Convert data file into GENEPOP or FSTAT format");
	if(*StatType==1 && *TypeComp==0){
		for(k=1;k<=(*NS);k++) if(Stat[k]==1 || Stat[k]==2 || Stat[k]==4 || Stat[k]==11) (*sigmaest)=1.0f;
		if((*sigmaest)){ printf("\n\n 5- Estimate gene dispersal sigma"); (*sigmaest)=0.0f;}
	}
	if(m>1){
		printf("\n\n 6- Report actual variance of pairwise genetic coefficients (Ritland 2000)");
//		printf("\n\n 7- Report inter-locus correlation for pairwise genetic coefficients");
	}

	printf("\n\nEnter 'R' to return to the first menu\n\n");

	fgets_chomp(smess, sizeof(smess), stdin);
	if(smess[0]=='r' || smess[0]=='R')  goto startagain;
	if(smess[0]!='\0'){	
		sscanf(smess,"%i",&stats);

		/*number of different options*/
		Noptions=1+(int)log10(stats);
		for(k=1;k<=Noptions;k++){
			if(stats%10==1) (*writeallelefreq)=1;
			if(stats%10==2) (*regdetails)=1;
			if(stats%10==3) (*distmatrix)=1;
			if(stats%10==4) (*export)=1;
			if(stats%10==5) (*sigmaest)=1.;
			if(stats%10==6) (*varcoef)=1;
			if(stats%10==7) (*Rbtwloc)=1;
			stats=(int)(stats/10);
		}
	}


	if(*sigmaest!=0.){
		if(ncoord==-2 && (*distm)==0) printf("\n\nEnter the effective population density for gene dispersal sigma estimation \n(density in # individuals per squared km): ");
		else printf("\n\nEnter the effective population density for gene dispersal sigma estimation \n(density in # individuals per squared distance unit, using same distance unit \nas for the spatial coordinates or the distance matrix): ");
		do{
			fgets_chomp(smess, sizeof(smess), stdin);
			ok=sscanf(smess,"%f",density);
		}while (ok!=1 || (*density)<=0.);
		printf("\nDefine X of the distance range sigma to X.sigma to compute regression slope\n or press RETURN for the default value X=20 : ");
		do{
			fgets_chomp(smess, sizeof(smess), stdin);
			if(smess[0]!='\0') 	ok=sscanf(smess,"%f",dwidth);
			else {(*dwidth)=20.0f; ok=1;}
		}while (ok!=1 || (*dwidth)<=1.);
	}


	if(*StatType>1)if(!ncoord && !(*distm) && !(*distmatrix)){
		printf("\n\n\nWARNING: In the absence of defined pairwise spatial distances, pairwise \n         statistics are not computed unless you ask to report matrices with \n         pairwise spatial and genetic distances\n         Do you want such matrices ? (y/n): ");
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]=='y' || smess[0]=='Y') *distmatrix=1;
	}

	
	if(*distmatrix){
		printf("\n\n\n\nFORMAT FOR PAIRWISE SPATIAL AND GENETIC DISTANCES\n\nPairwise spatial distances and genetic coefficients given for");
		printf("\n\n 1  : multilocus estimates only (columnar form)");
		printf("\n\n 2  : each locus and multilocus estimates (columnar form)");
		printf("\n\n 3  : multilocus estimates only (matrix and columnar forms)");
		printf("\n\n 4  : each locus and multilocus estimates (matrix and columnar forms)");
		printf("\n\n 5  : multilocus estimates, written in Phylip format\n\n");
		do{
			fgets_chomp(smess, sizeof(smess), stdin);
			*distmatrix=0;
			if(smess[0]=='1') *distmatrix=1;
			if(smess[0]=='2') *distmatrix=2;
			if(smess[0]=='3') *distmatrix=-1;
			if(smess[0]=='4') *distmatrix=-2;
			if(smess[0]=='5') *distmatrix=3;
		}while(*distmatrix<-2 || *distmatrix>3);
	}


}



/**************************************************************************************/



void read_pairwise_spatial_distances(char *inputfile,char *distfile,int n,
		struct name *namei,double **Mdij)
{
	int i,j,msize,line=0;
	int c,l,*iCol,*iLine,i1,i2;
	char s[SMAX],*s2,format,name[MAXNOM],smess[SMAX];
	double d;
	int npdefined;

	FILE *fp; 

	printf("\n\nReading pairwise spatial distances"); 

	/*check firstly that all names are different*/
	for(i=1;i<n;i++)for(j=i+1;j<=n;j++)if(!strncmp(namei[i].n,namei[j].n,MAXNOM-1)){
		printf("\n\nWARNING: Not all individual/population names are different.\nTo define each pairwise distance, all names must be different.\nCheck in the file \"%s\" the list of names found.\nPress any key to stop the program.",ERRORFILE);
		for(i=1;i<=n;i++) {
			sprintf(smess,"\n%s",namei[i].n);
			write(ERRORFILE,smess);
		}
		wait_a_char(); exit(0);
	}


	while((fp=fopen(distfile,"rt"))==NULL){
		printf("\nWARNING: Cannot open data file %c%s%c.\nIf it is being used by another application, close it first.\nIf it does not exist in the specified directory, bring it first.\nThen press RETURN.\nPress Ctrl+c twice to stop the programm now.\n",'"',distfile,'"');
		wait_a_char();
	}

	if(!strncmp(inputfile,distfile,MAXNOM-1)) do{	//if the distfile is the datafile, skip over lines until "END" is encountered
		readsfromfile(fp,s,distfile,&line);	
	}while(strncmp(s,"END",3)!=0);  

	fgets(s,SMAX,fp); line++;
	while(!feof(fp) && (strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1)) {fgets(s,SMAX,fp); line++;}  //skip over comment or empty lines

 	if(sscanf(s,"%c%i",&format,&msize)!=2){
		printf("\nCannot read matrix format type (M or C) or the number of elements of data \nfor pairwise distances in file \"%s\"",distfile);
		wait_a_char(); exit(0);
		
	}
	if(format!='C' && format!='c' && format!='M' && format!='m'){
		printf("\nIn data for pairwise distances in \"%s\", the format type \"%c\" is incorrect, it must be \"M\" (Matrix format) or \"C\" (Column format)",distfile,format); 
		wait_a_char(); exit(0);
	}
/*	if(msize!=n && (format=='M' || format=='m')){
		printf("\nIn data for pairwise distances in \"%s\", the number of elements announced (%i) does not match that expected (%i)",msize,n); 
		wait_a_char(); exit(0);
	}	   */

	/*init*/
	for(i=0;i<=n;i++)for(j=0;j<=n;j++) Mdij[i][j]=MISSVAL;

	if(format=='m' || format=='M'){	/*read data in matrix format*/
		iCol=ivector(0,n);
		iLine=ivector(0,n);
		s2=&s[0];
		for(c=1;c<=msize;c++){ /*get names heading columns*/
			s2=nexttab(s2,distfile,line);
			readsfromstring(s2,name,MAXNOM-1,distfile,line);
			for(i=1;i<=n;i++) if(!strncmp(name,namei[i].n,MAXNOM-1)) break;
			if(i>n){
				printf("\nIn data for pairwise distances in \"%s\", the name \"%s\" heading column %i is not recognized\nCheck in the error file the list of names expected",distfile,name,c);
				for(i=1;i<=n;i++) {
					sprintf(smess,"\n%s",namei[i].n);
					write_tofile_only(ERRORFILE,smess);
				}
				wait_a_char(); exit(0);
			}
			iCol[c]=i;
		}
		for(l=1;l<=msize;l++){
			readsfromfile(fp,s,distfile,&line);
			while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,distfile,&line);  /*skip over lines beggining with comments   */
	
			readsfromstring(s,name,MAXNOM-1,distfile,line);
			if(!strncmp(name,"END",MAXNOM-1)) {	/*if END is encountered*/
				printf("\nWARNING: END message encountered after only %d elements of the distance matrix on line %d from file %s",l-1,line,distfile);
				printf("\nPress 'g' if you want to go on anyway, otherwise press any other key to stop the program.");
				if(wait_a_char()!='g') exit(1);					
			} 
			for(i=1;i<=n;i++) if(!strncmp(name,namei[i].n,MAXNOM-1)) break;
			if(i>n){
				printf("\nIn data for pairwise distances in \"%s\", the name \"%s\" beggining line %i is not recognized\nCheck in the error file the list of names expected",distfile,name,l);
				for(i=1;i<=n;i++) {
					sprintf(smess,"\n%s",namei[i].n);
					write_tofile_only(ERRORFILE,smess);
				}
				wait_a_char(); exit(0);
			}
			iLine[l]=i;
			s2=&s[0];
			for(c=1;c<=msize;c++){
				if(!strpbrk(s2,"\t")) break;
				s2=nexttab_blank(s2,distfile,line);
				if(*s2=='\n') break;
				if(*s2=='\t') continue;
				if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading spatial distances on line %d of file %s \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now",line,distfile);wait_a_char();exit(1);} 
				readdoublefromstring(s2,&d,distfile,line);
				if(iLine[l]>=iCol[c]) {i=iCol[c];j=iLine[l];}
				else {j=iCol[c];i=iLine[l];}
				if(Mdij[i][j]==MISSVAL) Mdij[i][j]=Mdij[j][i]=d;
				else if(Mdij[i][j]!=d){
					printf("\nIn data for pairwise distances in \"%s\", the distance between \"%s\" and \"%s\" is redefine: first value= %f, second value= %f",distfile,namei[i].n,namei[j].n,Mdij[i][j],d);
					wait_a_char(); exit(0);
				}
			}
		}
		
		readsfromfile(fp,s,distfile,&line);
		while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,distfile,&line);  /*skip over lines beggining with comments */
		readsfromstring(s,name,MAXNOM-1,distfile,line);
		if(strncmp(name,"END",MAXNOM-1))  {	/*if END not encoutered after n individuals*/
			sprintf(smess,"\nEND message not encountered after %d elements of the matrix of pairwise distances on line %d from file %s",msize,line,distfile);
			write(ERRORFILE,smess);
			printf("\nPress 'g' if you want to go on anyway, otherwise press any other key to stop the program.");
			if(wait_a_char()!='g') exit(1);					
		}

		free_ivector(iCol,0,n);
		free_ivector(iLine,0,n);
	}	 /*end read matrix in M format*/

	if(format=='c' || format=='C'){	/*read data in column format*/
		for(l=1;l<=msize;l++){
			readsfromfile(fp,s,distfile,&line);
			while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,distfile,&line);  /*skip over lines beggining with comments */
			readsfromstring(s,name,MAXNOM-1,distfile,line);
			if(!strncmp(name,"END",MAXNOM-1)) {	/*if END is encountered*/
				printf("\nWARNING: END message encountered after only %d elements of the distance matrix on line %d from file %s",l-1,line,distfile);
				printf("\nPress 'g' if you want to go on anyway, otherwise press any other key to stop the program.");
				if(wait_a_char()!='g') exit(1);					
			} 

			for(i=1;i<=n;i++) if(!strncmp(name,namei[i].n,MAXNOM-1)) break;
			if(i>n){
				printf("\nIn data for pairwise distances in \"%s\", the first name \"%s\" on line %i is not recognized\nCheck in the error file the list of names expected",distfile,name,l);
				for(i=1;i<=n;i++) {
					sprintf(smess,"\n%s",namei[i].n);
					write_tofile_only(ERRORFILE,smess);
				}
				wait_a_char(); exit(0);
			}
			i1=i;

			s2=nexttab(s,distfile,line);
			readsfromstring(s2,name,MAXNOM-1,distfile,line);
			for(i=1;i<=n;i++) if(!strncmp(name,namei[i].n,MAXNOM-1)) break;
			if(i>n){
				printf("\nIn data for pairwise distances in \"%s\", the second name \"%s\" on line %i is not recognized\nCheck in the error file the list of names expected",distfile,name,l);
				for(i=1;i<=n;i++) {
					sprintf(smess,"\n%s",namei[i].n);
					write_tofile_only(ERRORFILE,smess);
				}
				wait_a_char(); exit(0);
			}
			i2=i;

			s2=nexttab(s2,distfile,line);
			if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading spatial distances on line %d of file %s \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now.",line,distfile);wait_a_char();exit(1);} 
			readdoublefromstring(s2,&d,distfile,line);
			if(i1>=i2) {i=i2;j=i1;}
			else {j=i2;i=i1;}
			if(Mdij[i][j]==MISSVAL) Mdij[i][j]=Mdij[j][i]=d;
			else if(Mdij[i][j]!=d){
				printf("\nIn data for pairwise distances in \"%s\", the distance between \"%s\" and \"%s\" is redefine: first value= %f, second value= %f",distfile,namei[i].n,namei[j].n,Mdij[i][j],d);
				wait_a_char(); exit(0);
			}
		}
		
		readsfromfile(fp,s,distfile,&line);
		while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,distfile,&line);  /*skip over lines beggining with comments  */
		readsfromstring(s,name,MAXNOM-1,distfile,line);
		if(strncmp(name,"END",MAXNOM-1))  {	/*if END not encoutered after n individuals*/
			sprintf(smess,"\nEND message not encountered after %d elements of the matrix of pairwise distances on line %d from file %s",msize,line,distfile);
			write(ERRORFILE,smess);
			printf("\nPress 'g' if you want to go on anyway, otherwise press any other key to stop the program.");
			if(wait_a_char()!='g') exit(1);					
		}

	}  /*end read matrix in C format*/

	npdefined=0;
	for(i=1;i<n;i++)for(j=i+1;j<=n;j++) if(Mdij[i][j]!=MISSVAL) npdefined++;
	if(npdefined!=((n*(n-1))/2)){
		printf("\nThere are only %i pairwise distances defined for %i potential pairwise comparisons (excluding self comparisons)\nPress any key to go on",npdefined,(n*(n-1))/2 );
		wait_a_char();
	}

/*	nselfdefined=0;		
	for(i=1;i<=n;i++) if(Mdij[i][i]!=MISSVAL) nselfdefined++;
	if(nselfdefined!=n && nselfdefined!=0){
		printf("\nThere are %i self pairwise distances defined for %i potential self comarisons\nPress any key to go on",nselfdefined,n);
		wait_a_char();
	}					 */

	fclose(fp);

	printf("\nPairwise spatial distances read successfully"); 

}



/**************************************************************************************/




void read_pairwise_genetic_distances(char *inputfile,char *distfile,int m,
		char namelocusp[] [MAXNOM],int *Nallelel,int **allelenamela,float ***Mgdlaa)
{
	int a1=0,a2,loc,loc2,line=0;
	int nameallele,nlocusread,nallelerecogn,nallelegiven;
	int l,c,*aCol,*aLine;
	char s[SMAX],*s2,nameloc[MAXNOM],smess[SMAX];
	float d;

	FILE *fp; 


	printf("\n\nReading pairwise genetic distances between alleles"); 


	//check that all locus names are different
	for(loc=1;loc<m;loc++)for(loc2=loc+1;loc2<=m;loc2++)if(!strncmp(namelocusp[loc],namelocusp[loc2],MAXNOM-1)){
		printf("\n\nWARNING: Not all locus names are different.\nTo define pairwise genetic distances, all locus names must be different.\nCheck in the file \"%s\" the list of names found.\nPress any key to stop the program.",ERRORFILE);
		for(loc=1;loc<=m;loc++) {
			sprintf(smess,"\n%s",namelocusp[loc]);
			write(ERRORFILE,smess);
		}
		wait_a_char(); exit(0);
	}




	while((fp=fopen(distfile,"rt"))==NULL){
		printf("\nWARNING: Cannot open data file %c%s%c.\nIf it is being used by another application, close it first.\nIf it does not exist in the specified directory, bring it first.\nThen press RETURN.\nClose the window to stop the programm now.\n",'"',distfile,'"');
		wait_a_char();
	}

	if(!strncmp(inputfile,distfile,MAXNOM-1)) do{	//if the distfile is the datafile, skip over lines until "END" is encountered
		readsfromfile(fp,s,distfile,&line);	
	}while(strncmp(s,"END",3)!=0);  




	//init
	for(loc=1;loc<=m;loc++)for(a1=0;a1<=Nallelel[loc];a1++)for(a2=0;a2<=Nallelel[loc];a2++) Mgdlaa[loc][a1][a2]=(float)MISSVAL;
	aCol=ivector(0,1000);
	aLine=ivector(0,1000);
 
	nlocusread=0;
 	do{

		//get locus name and find it (skip over data for loci undefined in the data file)
		do{
			fgets(s,SMAX,fp); line++;
			if(sscanf(s,"%s",nameloc)>0){
				for(loc=1;loc<=m;loc++) if(!strncmp(nameloc,namelocusp[loc],MAXNOM-1)) break;
			}
		}while(loc>m && !feof(fp));

		//check if end of file encountered before all loci read
		if(feof(fp) && nlocusread!=m){	  
			printf("\n\nGenetic distances between alleles in file \"%s\" found for %i loci whereas \nexpected for %i loci from the data file\nPress any key to stop",distfile,nlocusread,m);
			wait_a_char(); exit(0);
		}

		nlocusread++;
		s2=&s[0];
		//get allele names heading columns
		nallelerecogn=0;
		c=0;
		for(;;){
			if(!(s2=strpbrk(s2,"\t"))) break;
			if(s2[1]=='\t')	break; 
			else s2++;
			c++;
			readintfromstring(s2,&nameallele,distfile,line);
			for(a1=1;a1<=Nallelel[loc];a1++) if(nameallele==allelenamela[loc][a1]) break;
			if(a1>Nallelel[loc]) aCol[c]=aLine[c]=0;	 //allele not recognized
			else{ 
				aCol[c]=aLine[c]=a1;
				nallelerecogn++;
			}
		}
		nallelegiven=c;

		//check that all alleles from data file were encountered
		if(nallelerecogn!=Nallelel[loc]){
			printf("\n\nGenetic distances between alleles at locus \"%s\" in file \"%s\" not found for all alleles given in the data file\nPress any key to stop",namelocusp[loc],distfile);
			wait_a_char(); exit(0);
		}
  

		//read data in matrix format
		for(l=1;l<=nallelegiven;l++){
			readsfromfile(fp,s,distfile,&line);
			while(strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1) readsfromfile(fp,s,distfile,&line);  /*skip over lines beggining with comments   */

			if(strncmp(s,"END",3)==0 || feof(fp)) {	//if END is encountered or the end of the file
				printf("\nERROR: END message or end of file encountered after only %d alleles of the genetic distance matrix between alleles at locus \"%s\" from file %s for %d expected\nCheck that you defined a squared matrix with the same orders of alleles in lines and columns\nPress any key to stop the program",a1-1,namelocusp[loc],distfile,nallelegiven);
				wait_a_char(); exit(0);					
			} 

			//skip over alleles not given in the data file
			if(aLine[l]==0) continue;

			//read allele name and check it is the expected one
			readintfromstring(s,&nameallele,distfile,line);
			if(nameallele!=allelenamela[loc][aLine[l]]){
				printf("\nFATAL ERROR: The line for allele %d of the genetic distance matrix between alleles at locus \"%s\" in file %s was expected to be for allele %d\nCheck that you defined a squared matrix with the same orders of alleles in lines and columns\nPress any key to stop the program",nameallele,namelocusp[loc],distfile,allelenamela[loc][aLine[l]]);
				wait_a_char(); exit(0);					
			} 

			//read pairwise genetic distance value btw alleles
  			s2=&s[0];
			for(c=1;c<=nallelegiven;c++){
				if(!(s2=strpbrk(s2,"\t"))) break;
				else s2++;
				if(*s2=='\n') break;
				if(*s2=='\t') continue;
				if(aCol[c]==0) continue;   // skip over alleles not given in the data file

				if(strchr(s2,','))if(strchr(s2,',')<strpbrk(s2,"\t\n\r")){printf("\nWARNING: problem reading genetic distance between alleles on line %d of file %s \n\tdecimals must be indicated by a point, not a coma.\nPress any key to stop the program now.",line,distfile);wait_a_char();exit(1);} 
				readfloatfromstring(s2,&d,distfile,line);
				if(aLine[l]>=aCol[c]) {a1=aCol[c];a2=aLine[l];}
				else {a2=aCol[c];a1=aLine[l];}
				if(Mgdlaa[loc][a1][a2]==MISSVAL) Mgdlaa[loc][a1][a2]=Mgdlaa[loc][a2][a1]=d;
				else if(Mgdlaa[loc][a1][a2]!=d){
					printf("\nERROR: In data for pairwise genetic distances between alleles for locus \"%s\" in file \"%s\", the distance between alleles %d and %d is redefined: first value= %g, second value= %g",namelocusp[loc],distfile,allelenamela[loc][a1],allelenamela[loc][a2],Mgdlaa[loc][a1][a2],d);
					wait_a_char(); exit(0);
				}
			}
		}
		for(a1=1;a1<Nallelel[loc];a1++)for(a2=a1+1;a2<=Nallelel[loc];a2++) if(Mgdlaa[loc][a1][a2]==MISSVAL){
			printf("\nERROR: The genetic distance between alleles %d and %d at locus \"%s\" in file \"%s\" was undefined\nPress any key to stop the program",allelenamela[loc][a1],allelenamela[loc][a2],namelocusp[loc],distfile);
			wait_a_char(); exit(0);					
		} 
		for(a1=1;a1<=Nallelel[loc];a1++){
			if(Mgdlaa[loc][a1][a1]==MISSVAL) Mgdlaa[loc][a1][a1]=0.;
			if(Mgdlaa[loc][a1][a1]!=0.){
				printf("\nERROR: The genetic distance between allele %d and itself at locus \"%s\" in file \"%s\" was set as %g whereas it must be set as 0 (or undefined)\nPress any key to stop the program",allelenamela[loc][a1],namelocusp[loc],distfile,Mgdlaa[loc][a1][a1]);
				wait_a_char(); exit(0);	
			}
		} 

	}while(nlocusread<m); //end loop on nlocusread

	free_ivector(aCol,0,1000);
	free_ivector(aLine,0,1000);


	fclose(fp);

	printf("\nPairwise genetic distances between alleles read successfully"); 


}	//end of read_pairwise_genetic_distances

/***************************************************************************/
void read_allele_frequencies(char *freqfile,int m,char namelocus[][MAXNOM],int *Nallelel,
							 int **allelesizela, float **givenPla,int *Ngivenallelel)
{
	int line=0,a,allele,allelename,nallelel[MMAX],Nallelemax;
	int locus[MMAX],loc,loc2,l,Nloc;
	int error=0,error1,error2;
	float freq;
	char s[SMAX],smess[SMAX],locusname[MAXNOM],*s2;

	FILE *fp;


	printf("\n\nReading reference allele frequencies."); 

	while((fp=fopen(freqfile,"rt"))==NULL){
		printf("\n\nWARNING: Cannot open reference allele frequency file \"%s\"\n  If it is being used by another application, close it first, then press RETURN\n  Otherwise enter the correct file name\n  Close the window if you wish to stop the program now\n",freqfile);
		fgets_chomp(smess, sizeof(smess), stdin);
		if(smess[0]!='\0') strncpy(freqfile,smess,MAXNOM-1);
	}

	//read first line with locus names and number of alleles
	fgets(s,SMAX,fp); line++;
	while(!feof(fp) && (strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1)) {fgets(s,SMAX,fp); line++;}  //skip over comment or empty lines

	for(l=1;l<=m;l++) Ngivenallelel[l]=0;
	loc=0;
	s2=&s[0];
	while(s2){
		sscanf(s2,"%s",locusname);
		loc++;
		for(l=1;l<=m;l++) if(strcmp(namelocus[l],locusname)==0) break;
		if((s2=strpbrk(s2,"\t"))){
			s2++;
			if(l<=m){sscanf(s2,"%d",&Ngivenallelel[l]); locus[loc]=l;}
			else locus[loc]=0;
		}
		if((s2=strpbrk(s2,"\t"))) s2++;
	}
	Nloc=loc;//number of loci in the allele freq file

	//check that all loci of the data file are represented
	error=0;
	for(l=1;l<=m;l++) if(Ngivenallelel[l]==0){
		printf("\nERROR: Reference allele frequencies not given for locus \"%s\".",namelocus[l]);
	}
	if(error){
		printf("\nPress any key to stop the program now.\n");
		wait_a_char(); exit(1);
	}


	Nallelemax=0;
	for(l=1;l<=m;l++) if(Ngivenallelel[l]>Nallelemax)	Nallelemax=Ngivenallelel[l];

	for(l=1;l<=m;l++) for(a=0;a<=Nallelemax+1;a++) givenPla[l][a]=0.;

	for(l=1;l<=m;l++) nallelel[l]=Nallelel[l];

	//read next lines with allele names and frequencies
	for(allele=1;allele<=Nallelemax;allele++){
		fgets(s,SMAX,fp); line++;
		while(!feof(fp) && (strncmp(s,"/*",2)==0 || strncmp(s,"//",2)==0 || s[0]=='"' || sscanf(s,"%s",smess)<1)) {fgets(s,SMAX,fp); line++;}  //skip over comment or empty lines
		s2=&s[0];
		for(loc=1;loc<=Nloc;loc++){
			if(locus[loc])if(allele<=Ngivenallelel[locus[loc]]){
				sscanf(s2,"%d%g",&allelename,&freq);
				for(a=1;a<=nallelel[locus[loc]];a++) if(allelename==allelesizela[locus[loc]][a]) break;
				if(a>nallelel[locus[loc]]){ //add alleles non existant in the data file
					allelesizela[locus[loc]][a]=allelename;
					nallelel[locus[loc]]++;
				}
				givenPla[locus[loc]][a]=freq;
			}
			for(loc2=loc+1;loc2<=Nloc;loc2++) if(allele<=Ngivenallelel[locus[loc2]]){
				s2=strpbrk(s2,"\t");s2++;s2=strpbrk(s2,"\t");s2++;
				break;
			}
		}
	}
	fclose(fp);

	for(l=1;l<=m;l++) Ngivenallelel[l]=nallelel[l];


	//test if sum of allele frequencies equal one
	error1=0;
	for(l=1;l<=m;l++){
		freq=0.;
		for(a=1;a<=Ngivenallelel[l];a++) freq+=givenPla[l][a];
		if(freq<0.999 || freq>1.001){
			printf("\nERROR: Sum of reference allele frequencies for locus \"%s\" different\nfrom one (=%g).",namelocus[l],freq);
			error1=1;
		}
	}

	//check that all alleles from the data file have non zero frequency in the reference population
	error2=0;
	for(l=1;l<=m;l++){
		for(a=1;a<=Nallelel[l];a++) if(givenPla[l][a]==0.){
			if(error2==0) printf("\nWARNING: the following alleles found in the data file have a zero frequency in the reference allele frequency file. \nThis currently not allowed.");
			printf("\n\tlocus \"%s\", allele \"%d\"",namelocus[l],allelesizela[l][a]);
			error2=1;
		}
	}

	if(error1 || error2){
		printf("\nPress any key to stop the program now.\n");
		wait_a_char();
	}

}


/**************************************************************************************/
void write_allele_freq(char *outputfilename,int n,int m,char namelocus[][MAXNOM],int ndigit,
		int Npop,int *Nip,struct name namepop[],int StatType,int NS,int Stat[],int **allelesizela,
		int **Nallelepl,float ***Ppla,float **Hepl,float **hTpl,float **vTpl,float **Dmpl,float **Dwmpl,float **Masizepl,float **Vasizepl,
		int **Nmissinggenotpl,int **Nincompletegenotpl,int **Nvalgenpl,int printallelefreq,int FreqRef,int *Ngivenallelel,float **givenPla)
{
	int l,S,a,p,allelesize,orderedalleles,Ntab,t;
	FILE *fp;

	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}

	allelesize=0;
	for(S=1;S<=NS;S++)if((StatType==1 && Stat[S]==5)||(StatType>1 && (Stat[S]==5 || Stat[S]==6 || Stat[S]==8))) allelesize=1;
	orderedalleles=0;
	for(S=1;S<=NS;S++)if((StatType>1 && (Stat[S]==11 || Stat[S]==12))) orderedalleles=1;

	fprintf(fp,"\n\nALLELE FREQUENCIES");
	if(ndigit>0){
		Ntab=6;
		fprintf(fp,"\nLocus\t# missing genotypes (%%)\t# incomplete genotypes (%%)\t# of defined genes\t# alleles\tHe (gene diversity corrected for sample size, Nei 1978)");
		if(StatType>1){
			Ntab++;
			fprintf(fp,"\th (gene diversity with UNORDERED alleles, Pons & Petit 1996)");
			if(orderedalleles){	fprintf(fp,"\tv (non-standardised gene diversity with ORDERED alleles, Pons & Petit 1996)\tDm (unweighted mean genetic distance between existing alleles)\tDwm (weighted mean genetic distance between alleles, weights = allele frequencies)");Ntab+=3;}
		}
		if(allelesize){	fprintf(fp,"\tMean allele size\tVariance of allele size"); Ntab+=2;}
	//	fprintf(fp,"\t\tAllele names\tAllele frequencies");
		Ntab++;
	}
	else {fprintf(fp,"\nLocus\t# missing genotypes\tFrequency of dominant genotype"); Ntab=3;}
	for(p=0;p<=Npop;p++){
		if(Npop && p==0){
			if(StatType==1) fprintf(fp,"\nAll categories confounded (%i individuals)",Nip[p]);
			else fprintf(fp,"\nAll populations confounded (%i individuals)",Nip[p]);
		}
		if(printallelefreq==0 && p>0) break;
		if(Npop && p){
			if(StatType==1) fprintf(fp,"\nCategory %s (%i individuals)",namepop[p].n,Nip[p]);
			else fprintf(fp,"\nPopulation %s (%i individuals)",namepop[p].n,Nip[p]);
		}

		for(l=1;l<=m;l++){
			if(ndigit>0) fprintf(fp,"\n%s\t%i (%.1f%%)\t%i (%.1f%%)\t%i",namelocus[l],Nmissinggenotpl[p][l],(float)(Nmissinggenotpl[p][l]*100./Nip[p]),Nincompletegenotpl[p][l],(float)(Nincompletegenotpl[p][l]*100./Nip[p]),Nvalgenpl[p][l]);
 			else fprintf(fp,"\n%s\t%i (%.1f%%)",namelocus[l],Nmissinggenotpl[p][l],(float)(Nmissinggenotpl[p][l]*100./Nip[p]));

			if(ndigit>0){
				if(Nvalgenpl[p][l]){
					fprintf(fp,"\t%i\t%.4f",Nallelepl[p][l],Hepl[p][l]);
					if(StatType>1){
						fprintf(fp,"\t%.4f",hTpl[p][l]);
						if(orderedalleles){
							fprintf(fp,"\t%.4f",vTpl[p][l]);
							if(Dmpl[p][l]!=MISSVAL) fprintf(fp,"\t%.4f",Dmpl[p][l]);
							else fprintf(fp,"\t");
							if(Dwmpl[p][l]!=MISSVAL) fprintf(fp,"\t%.4f",Dwmpl[p][l]);
							else fprintf(fp,"\t");
						}
					}
					if(allelesize) fprintf(fp,"\t%.1f\t%.1f",Masizepl[p][l],Vasizepl[p][l]);
					fprintf(fp,"\t");
				}
			}

			//write allele names and freq
		/*	if(ndigit>0){
				if(Nvalgenpl[p][l]){
					fprintf(fp,"\t%i",allelesizela[l][1]);
					for(a=2;a<=Nallelepl[0][l];a++) fprintf(fp,", %i",allelesizela[l][a]);
 					for(a=1;a<=Nallelepl[0][l];a++) fprintf(fp,"\t%g",Ppla[p][l][a]);
				}
			}		*/
			if(ndigit>0){
				if(Nvalgenpl[p][l]){
					fprintf(fp,"\tallele");
					for(a=1;a<=Nallelepl[0][l];a++) fprintf(fp,"\t%i",allelesizela[l][a]);
					fprintf(fp,"\n");
					for(t=1;t<Ntab;t++) fprintf(fp,"\t");
					fprintf(fp,"\tfrequency");
 					for(a=1;a<=Nallelepl[0][l];a++) fprintf(fp,"\t%g",Ppla[p][l][a]);
				}
			}		
			else if(Nvalgenpl[p][l]) fprintf(fp,"\t%g",Ppla[p][l][2]);
		}
	}

	if(FreqRef==-1){
		fprintf(fp,"\n\nREFERENCE ALLELE FREQUENCIES");

		for(l=1;l<=m;l++){
			fprintf(fp,"\n%s",namelocus[l]);
			/*write allele names and freq*/
			if(ndigit>0){
				for(a=1;a<=Ngivenallelel[l];a++) fprintf(fp,"\t%i",allelesizela[l][a]);
				fprintf(fp,"\n");
 				for(a=1;a<=Ngivenallelel[l];a++) fprintf(fp,"\t%g",givenPla[l][a]);
			}
			else fprintf(fp,"\t%g",givenPla[l][2]);
		}
	}


	fclose (fp);
}	/*end of write_allele_freq*/

/**************************************************************************************/
void displaydist(int argc,char *outputfilename,int nc,double *maxc,double *mdc,double *mlndc,
				 int *npc,float **indexpartic)
{
	int c;
	char ch,smess[SMAX];
	FILE *fp;

	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}
	fclose (fp);

	printf("\n\nDISTANCE INTERVALS (%i) :",nc);
 /*	fprintf(fp,"\n\n%i distance intervals:",nc);
	write(outputfilename,smess);   */
	printf("\ninterval   max d     mean d   mean ln(d)   # pairs   %% partic  CV #partic");
/*	fprintf(fp,"\ninterval\tmax dist\tmean dist\tmean ln(dist)\t# pairs\t%% participation\tCV #participation");
	write_tofile_only(outputfilename,smess);	 */
	for(c=1;c<=nc;c++) {
		if(npc[c]) printf("\n%4d %11.2f %10.2f %9.3f %12d %9.1f %9.2f",c,maxc[c],mdc[c],mlndc[c],npc[c],100*indexpartic[c][0],indexpartic[c][1]);
		else  printf("\n%4d   %9.2f                         %9d",c,maxc[c],npc[c]);
	/*	if(npc[c]) fprintf(fp,"\n%d\t%G\t%G\t%G\t%ld\t%.1f\t%.2f",c,maxc[c],mdc[c],mlndc[c],npc[c],100*indexpartic[c][0],indexpartic[c][1]);
		else  fprintf(fp,"\n%d\t%G\t\t\t%ld",c,maxc[c],npc[c]);
		write_tofile_only(outputfilename,smess);   */
		if(nc>20 && (c%20)==0 && argc<4){printf("\nPress any key to see the next intervals");wait_a_char();}
	}
 	printf("\n\n");
/*	printf("\nCheck the values and press RETURN to continue, or enter 's' to stop now\n");
	if((ch=wait_a_char())=='s' || ch=='S') exit(0);	  */
	
}    /*end of displaydist*/


/**************************************************************************************/

void writeIndStatresults(char *outputfilename,int n,int Nsg,int m,char namelocus[][MAXNOM],int nc,double *maxc,
		int *npc,float **indexpartic,double *mdc,double *mlndc,float dijmin,float dijmax,float givenF,
		int TypeComp,int cat1,int cat2,struct name *namecat,int FreqRef,int NS,int Stat[12],float **corrSlc[],float density,float dwidth,
		int JKest,int regdetails,int varcoef,int Rbtwloc,float ***RSll[12],float ***V[12],float **R2pl[12])
{
	int c,l,S,a;
	char smess[SMAX];
	FILE *fp;

	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}

	fprintf(fp,"\n\n\nVALUES OF THE COMPUTED STATISTICS");	
	fprintf(fp,"\nGenetic analyses at INDIVIDUAL level (%i individuals)",n);
	
	if(TypeComp){
		if(TypeComp==0) fprintf(fp," for all pairs of individuals");
		if(TypeComp==1 && !cat1) fprintf(fp," for pairs of individuals WITHIN category");
		if(TypeComp==2 && !cat1 && !cat2) fprintf(fp," for pairs of individuals AMONG categories");
		if(TypeComp==1 && cat1 && !cat2) fprintf(fp," for pairs of individuals WITHIN category \"%s\"",namecat[cat1].n);
		if(TypeComp==2 && cat1 && cat2) fprintf(fp," for pairs of individuals BETWEEN categories \"%s\" and \"%s\"",namecat[cat1].n,namecat[cat2].n);		
	}
	if(TypeComp){
		if(FreqRef==0) fprintf(fp," (reference sample (allele frequencies) = WHOLE SAMPLE)");
		if(FreqRef==1) fprintf(fp," (reference sample (allele frequencies) = SAMPLE WITHIN CATEGORY)");		
	}
 	if(FreqRef==-1) fprintf(fp," (relative to GIVEN REFERENCE ALLELE FREQUENCIES (except for Moran's I relationship coef and Rousset's 'a' coef))");

	fprintf(fp,"\nDist classes\t");	
	for(c=1;c<=nc;c++) fprintf(fp,"\t%i",c);
	fprintf(fp,"\nMax distance\t");	
	for(c=1;c<=nc;c++) fprintf(fp,"\t%.4f",maxc[c]);
	fprintf(fp,"\nNumber of pairs\t");	
	for(c=1;c<=nc;c++) fprintf(fp,"\t%i",npc[c]);
	fprintf(fp,"\n%% partic");	
	for(c=0;c<=nc;c++){if(npc[c])fprintf(fp,"\t%.1f",100*indexpartic[c][0]);else fprintf(fp,"\t");}
	fprintf(fp,"\nCV partic");	
	for(c=0;c<=nc;c++){if(npc[c])fprintf(fp,"\t%.2f",indexpartic[c][1]);else fprintf(fp,"\t");}
	fprintf(fp,"\nMean distance\t");	
	for(c=1;c<=nc;c++){if(npc[c])fprintf(fp,"\t%.4f",mdc[c]);else fprintf(fp,"\t");}
	fprintf(fp,"\nMean ln(distance)\t");	
	for(c=1;c<=nc;c++){if(npc[c] && mdc[c]>0.)fprintf(fp,"\t%.4f",mlndc[c]);else fprintf(fp,"\t");}
		
	for(S=1;S<=NS;S++){
		if(Stat[S]==1) fprintf(fp,"\n\nPairwise KINSHIP coefficients (Loiselle et al., 1995)");
		if(Stat[S]==2) fprintf(fp,"\n\nPairwise KINSHIP coefficients (Ritland, 1996)");
	/*	if(Stat[S]==3) fprintf(fp,"\n\nKinship according to Ritland (modified)");*/
		if(Stat[S]==3) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
		if(Stat[S]==4) fprintf(fp,"\n\nPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
		if(Stat[S]==5) fprintf(fp,"\n\nPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
		if(Stat[S]==6) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
		if(Stat[S]==7) fprintf(fp,"\n\nPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
		if(Stat[S]==8) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
		if(Stat[S]==9) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
		if(Stat[S]==10) fprintf(fp,"\n\nPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
		if(Stat[S]==11) fprintf(fp,"\n\nPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
		if(Stat[S]==12) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
		if(Stat[S]==13) fprintf(fp,"\n\nPairwise RELATIONSHIP coefficients (Li et al., 1993)");
		if(Stat[S]==14) fprintf(fp,"\n\nPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");
		
		for(c=0;c<=nc+2;c++) fprintf(fp,"\t");
		fprintf(fp,"distance range for regression analyses");
		
		if(regdetails) fprintf(fp,"\tRegression with linear distance\t\t\t\t\t\t\t\t\tRegression with ln(distance)");			

		fprintf(fp,"\nLocus\tintra-individual (inbreeding coef)");		
		if(Nsg>1) fprintf(fp,"\tintra-group");
		else fprintf(fp,"\t1");
		
		for(c=2;c<=nc;c++) fprintf(fp,"\t%i",c);
		if(nc>1) fprintf(fp,"\taverage");
		else fprintf(fp,"\t");		

		if(regdetails) fprintf(fp,"\t%G-%G\tb (slope)\ta (intercept)\tr2 (coef of determination)\tN (# pairs)\tMd (mean distance)\tVd (variance of distance)\tMv (mean statistic)\tVv (variance of statistic)\t\tb (slope)\ta (intercept)\tr2 (coef of determination)\tN (# pairs)\tMd (mean ln(distance))\tVd (variance of ln(distance))\tMv (mean statistic)\tVv (variance of statistic)\t",dijmin,dijmax); 
		else fprintf(fp,"\t%G-%G\tb-lin(slope linear dist)\tb-log(slope log dist)",dijmin,dijmax);
		
		for(l=0;l<=m+4;l++){
			if(m==1 && l==0) l=1;
			if(l==0) fprintf(fp,"\nALL LOCI");
			if(l>0 && l<=m) fprintf(fp,"\n%s",namelocus[l]);
			if(m==1 && l>m) break;
			if(JKest==0 && l>m && l!=m+3) continue;
			if(varcoef==0 && l>m+2) break;
			if(l==m+1) fprintf(fp,"\nJackknifed estimators (over loci)\nMean");
			if(l==m+2) fprintf(fp,"\nSE");
			if(l==m+3) fprintf(fp,"\nActual variance of pairwise coefficients (method of Ritland 2000)\nEstimate");
			if(l==m+4) fprintf(fp,"\nSE (jackknife over loci)");			

			for(c=0;c<=nc+1;c++){
				if(((Stat[S]!=1 && Stat[S]!=2 && Stat[S]!=5 && Stat[S]!=14) && c==0) || corrSlc[S][l][c]==MISSVAL) fprintf(fp,"\t");
				else if(l<=m+2) fprintf(fp,"\t%.4f",corrSlc[S][l][c]);
				else fprintf(fp,"\t%G",corrSlc[S][l][c]);
				if(nc==1 && c==nc+1) fprintf(fp,"\t");				
			}
			fprintf(fp,"\t");
			
			if(l<=m+2) for(c=-20;c<=-1;c++){
				if(regdetails && (c==-12 || c>=-2)) continue;
				if(!regdetails && c<-2) c=-2;
				if(corrSlc[S][l][c]==MISSVAL || (l>m && (c==-8 ||c==-7 ||c==-6 ||c==-5 ||c==-3 ||c==-18 ||c==-17 ||c==-16 ||c==-15 ||c==-13)) ) fprintf(fp,"\t");
				else fprintf(fp,"\t%G",corrSlc[S][l][c]);				
			}
		}/*end of loop l*/

		if(density!=0. && (Stat[S]==1 || Stat[S]==2 || Stat[S]==11 || Stat[S]==4) ){ //estimates of gene dispersal distances
			if(m>1) l=0;
			else l=1;
			fprintf(fp,"\nEstimated gene dispersal parameters for an assumed effective pop density = %G  (#individuals per square distance unit) in a 2-dimensional population at drift-dispersal equilibrium under isotropic dispersal. Estimation based on the regression slope between sigma and %.1fsigma following an iterative procedure.",density,dwidth);
			fprintf(fp,"\n\tNeighbourhood size Nb:");
			if(corrSlc[S][l][-21]==(float)MISSVAL) fprintf(fp,"\t\tNo convergence");
			else  fprintf(fp,"\t\t%G",corrSlc[S][l][-21]);
			if(corrSlc[S][l][-22]!=corrSlc[S][l][-23]) fprintf(fp,"\t(mean over iterations cycling in the range %G to %G)",corrSlc[S][l][-22],corrSlc[S][l][-23]);
			else fprintf(fp,"\t");
			if(JKest){
				fprintf(fp,"\tSE (standard error by jackknifying over loci):");
				if(corrSlc[S][m+2][-21]!=(float)MISSVAL) fprintf(fp,"\t%G",corrSlc[S][m+2][-21]);
				else fprintf(fp,"\tConvergence not achieved when removing some loci");
			}
			fprintf(fp,"\n\tSigma (half the mean square parent-offspring distance):");
			if(corrSlc[S][l][-25]==(float)MISSVAL) fprintf(fp,"\t\tNo convergence");
			else  fprintf(fp,"\t\t%G",corrSlc[S][l][-25]);
			if(corrSlc[S][l][-26]!=corrSlc[S][l][-27]) fprintf(fp,"\t(mean over iterations cycling in the range %G to %G)",corrSlc[S][l][-26],corrSlc[S][l][-27]);
			else fprintf(fp,"\t");
			if(JKest){
				fprintf(fp,"\tSE (standard error by jackknifying over loci):");
				if(corrSlc[S][m+2][-25]!=(float)MISSVAL) fprintf(fp,"\t%G",corrSlc[S][m+2][-25]);
				else fprintf(fp,"\tConvergence not achieved when removing some loci");
			}
		}

		if(Rbtwloc){
			fprintf(fp,"\n\n\t\tMean inter-locus correlation coefficient for pairwise coefficients (Pairon & Hardy, unpublished)\n\t\t%G\tSE (Jackknife over loci)=\t%G\t\tVp1=\t%G\tSE(Vp1)=\t%G\tRl1=\t%G\tSE(Rl1)=\t%G\t\tVp2=\t%G\tSE(Vp2)=\t%G\tRl2=\t%G\tSE(Rl2)=\t%G\tVl=\t%G\tSE(Vl)=\t%G",RSll[S][0][0][0],RSll[S][0][2][0],V[S][0][0][0],V[S][0][0][2],V[S][0][1][0],V[S][0][1][2],V[S][0][2][0],V[S][0][2][2],V[S][0][3][0],V[S][0][3][2],V[S][0][4][0],V[S][0][4][2]);			
			fprintf(fp,"\n\t\tMean inter-locus correlation coefficient for residuals of pairwise coefficients after accounting for means per class (Pairon & Hardy, unpublished)\n\t\t%G\tSE (Jackknife over loci)=\t%G\t\tVp1=\t%G\tSE(Vp1)=\t%G\tRl1=\t%G\tSE(Rl1)=\t%G\t\tVp2=\t%G\tSE(Vp2)=\t%G\tRl2=\t%G\tSE(Rl2)=\t%G\tVl=\t%G\tSE(Vl)=\t%G",RSll[S][1][0][0],RSll[S][1][2][0],V[S][1][0][0],V[S][1][0][2],V[S][1][1][0],V[S][1][1][2],V[S][1][2][0],V[S][1][2][2],V[S][1][3][0],V[S][1][3][2],V[S][1][4][0],V[S][1][4][2]);			

			fprintf(fp,"\n\n\t\tMean correlation coefficient between 2 random pool of loci for pairwise coefficients (Pairon & Hardy, unpublished)\n\t\t%G\tSE (Jackknife over loci)=\t%G",R2pl[S][0][0],R2pl[S][0][2]);			
			fprintf(fp,"\n\t\tMean correlation coefficient between 2 random pool of loci for residuals of pairwise coefficients after accounting for means per class (Pairon & Hardy, unpublished)\n\t\t%G\tSE (Jackknife over loci)=\t%G",R2pl[S][1][0],R2pl[S][1][2]);						
		}
		
	}/*end of loop S*/

	fclose (fp);

}	/*end of procedure writeFrresults*/


/**************************************************************************************/
void writePopStatresults(char *outputfilename,int Npop,int m,char namelocus[][MAXNOM],int nc,double *maxc,
		int *npc,float **indexpartic,double *mdc,double *mlndc,float dijmin,float dijmax,
		int StatType,int TypeComp,int cat1,int cat2,struct name *namecat,int NS,int Stat[12],float **corrSlc[],float **FstatSlr[],
		int JKest,int regdetails,int varcoef,int Rbtwloc,float ***RSll[12],float ***V[12],float **R2pl[12],int PWstat)
{
	int c,l,S,a,r;
	char smess[SMAX];
	FILE *fp;

	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}

	fprintf(fp,"\n\n\nVALUES OF THE COMPUTED STATISTICS");
	
	if(StatType==2) fprintf(fp,"\nGenetic analyses at POPULATION level (1 pop = 1 CATEGORICAL GROUP -> %i populations)",Npop);
	if(StatType==3) fprintf(fp,"\nGenetic analyses at POPULATION level (1 pop = 1 SPATIAL GROUP -> %i populations)",Npop);
	if(StatType==4) fprintf(fp,"\nGenetic analyses at POPULATION level (1 pop = 1 SPATIAL GROUP of given CATEGORY -> %i populations)",Npop);
	
	if(StatType==4){
		if(TypeComp==0) fprintf(fp," for all pairs of populations");
		if(TypeComp==1 && !cat1) fprintf(fp," for pairs of populations WITHIN category");
		if(TypeComp==2 && !cat1 && !cat2) fprintf(fp," for pairs of populations AMONG categories");
		if(TypeComp==1 && cat1 && !cat2) fprintf(fp," for pairs of populations WITHIN category \"%s\"",namecat[cat1].n);
		if(TypeComp==2 && cat1 && cat2) fprintf(fp," for pairs of populations BETWEEN categories \"%s\" and \"%s\"",namecat[cat1].n,namecat[cat2].n);
		
	}

	if(PWstat>1){
		fprintf(fp,"\nDist classes\t\t\t\t\t");
		
		for(c=1;c<=nc;c++)fprintf(fp,"\t%i",c);
		fprintf(fp,"\nMax distance\t\t\t\t\t");
		
		for(c=1;c<=nc;c++)fprintf(fp,"\t%G",maxc[c]);
		fprintf(fp,"\nNumber of pairs\t\t\t\t\t");
		
		for(c=1;c<=nc;c++)fprintf(fp,"\t%i",npc[c]);
		fprintf(fp,"\n%% partic\t\t\t\t\t");
		
		for(c=1;c<=nc;c++){if(npc[c])fprintf(fp,"\t%.1f",100*indexpartic[c][0]);else fprintf(fp,"\t");}
		fprintf(fp,"\nCV partic\t\t\t\t\t");
		
		for(c=1;c<=nc;c++){if(npc[c])fprintf(fp,"\t%.2f",indexpartic[c][1]);else fprintf(fp,"\t");}
		fprintf(fp,"\nMean distance\t\t\t\t\t");
		
		for(c=1;c<=nc;c++){if(npc[c])fprintf(fp,"\t%G",mdc[c]);else fprintf(fp,"\t");}
		fprintf(fp,"\nMean ln(distance)\t\t\t\t\t");
		
		for(c=1;c<=nc;c++){if(npc[c] && mdc[c])fprintf(fp,"\t%G",mlndc[c]);else fprintf(fp,"\t");}
	}
	
	for(S=1;S<=NS;S++){
		if(PWstat){
			if(Stat[S]==1) fprintf(fp,"\n\n\tGlobal F-statistics \t\t\t\t\tPairwise Fst (ANOVA approach)");
			if(Stat[S]==2) fprintf(fp,"\n\n\tGlobal F-statistics \t\t\t\t\tPairwise Fst/(1-Fst) (ANOVA approach)");
			if(Stat[S]==3) fprintf(fp,"\n\n\tGlobal F-statistics \t\t\t\t\tPairwise Rho (ANOVA approach)");
			if(Stat[S]==4) fprintf(fp,"\n\n\tGlobal F-statistics \t\t\t\t\tPairwise Rho/(1-Rho) (ANOVA approach)");
			if(Stat[S]==5) fprintf(fp,"\n\n\tGlobal R-statistics \t\t\t\t\tPairwise Rst (ANOVA approach)");
			if(Stat[S]==6) fprintf(fp,"\n\n\tGlobal R-statistics \t\t\t\t\tPairwise Rst/(1-Rst) (ANOVA approach)");
			if(Stat[S]==7) fprintf(fp,"\n\n\tGlobal F-statistics \t\t\t\t\tPairwise Ds (Nei's 1978 standard distance)");
			if(Stat[S]==8) fprintf(fp,"\n\n\tGlobal R-statistics \t\t\t\t\tPairwise dm2 (Goldstein's 1995 distance)");
			if(Stat[S]==9) fprintf(fp,"\n\n\tGlobal Gst \t\t\t\t\tPairwise Gst (Pons & Petit 1995)");
			if(Stat[S]==10) fprintf(fp,"\n\n\tGlobal Gst \t\t\t\t\tPairwise Gij (Hardy)");
			if(Stat[S]==11) fprintf(fp,"\n\n\tGlobal Nst \t\t\t\t\tPairwise Nst (Pons & Petit 1996)");
			if(Stat[S]==12) fprintf(fp,"\n\n\tGlobal Nst \t\t\t\t\tPairwise Nij (Hardy)");
		}
		else{
			if(Stat[S]==1) fprintf(fp,"\n\n\tGlobal F-statistics");
			if(Stat[S]==2) fprintf(fp,"\n\n\tGlobal F-statistics");
			if(Stat[S]==3) fprintf(fp,"\n\n\tGlobal F-statistics");
			if(Stat[S]==4) fprintf(fp,"\n\n\tGlobal F-statistics");
			if(Stat[S]==5) fprintf(fp,"\n\n\tGlobal R-statistics");
			if(Stat[S]==6) fprintf(fp,"\n\n\tGlobal R-statistics");
			if(Stat[S]==7) fprintf(fp,"\n\n\tGlobal F-statistics");
			if(Stat[S]==8) fprintf(fp,"\n\n\tGlobal R-statistics");
			if(Stat[S]==9) fprintf(fp,"\n\n\tGlobal Gst");
			if(Stat[S]==10) fprintf(fp,"\n\n\tGlobal Gst");
			if(Stat[S]==11) fprintf(fp,"\n\n\tGlobal Nst");
			if(Stat[S]==12) fprintf(fp,"\n\n\tGlobal Nst");
		}
		

		if(PWstat>1){
			for(c=1;c<=nc+1;c++){fprintf(fp,"\t");}
			fprintf(fp,"distance range for regression analyses");			
			if(regdetails) fprintf(fp,"\tRegression with LINEAR DISTANCE\t\t\t\t\t\t\t\t\tRegression with LN(DISTANCE)");			
		}
		
		if(Stat[S]<=4 || Stat[S]==7) fprintf(fp,"\nLocus\tFit\tFis\tFst\t");
		if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8) fprintf(fp,"\nLocus\tRit\tRis\tRst\t");
		if(Stat[S]==9 || Stat[S]==10) fprintf(fp,"\nLocus\t\t\tGst\t");
		if(Stat[S]==11 || Stat[S]==12) fprintf(fp,"\nLocus\t\t\tNst\t");
		
		if(Stat[S]==3 || Stat[S]==4) fprintf(fp,"Rho\t");
		else fprintf(fp,"\t");
		

		if(PWstat>1){
			for(c=1;c<=nc;c++) fprintf(fp,"\t%i",c);
			if(nc>1) fprintf(fp,"\taverage");
			else fprintf(fp,"\t");
			
			if(regdetails) fprintf(fp,"\t%G-%G\tb (slope)\ta (intercept)\tr2 (coef of determination)\tN (# pairs)\tMd (mean distance)\tVd (variance of distance)\tMv (mean statistic)\tVv (variance of statistic)\t\tb (slope)\ta (intercept)\tr2 (coef of determination)\tN (# pairs)\tMd (mean ln(distance))\tVd (variance of ln(distance))\tMv (mean statistic)\tVv (variance of statistic)\t",dijmin,dijmax);
			else fprintf(fp,"\t%G-%G\tb-lin(slope linear dist)\tb-log(slope log dist)",dijmin,dijmax);
			
		}
		else if(PWstat==1) fprintf(fp,"\taverage");
			
		for(l=0;l<=m+4;l++){
			if(m==1 && l==0) l=1;
			if(l==0) fprintf(fp,"\nALL LOCI");
			if(l>0 && l<=m) fprintf(fp,"\n%s",namelocus[l]);
			if(m==1 && l>m) break;
			if(JKest==0 && l>m && l!=m+3) continue;
			if(l>m+2 && (varcoef==0 || Stat[S]==7 || Stat[S]==8)) break;
			if(l==m+1) fprintf(fp,"\nJackknifed estimators (over loci)\nMean");
			if(l==m+2) fprintf(fp,"\nSE");
			if(l==m+3) fprintf(fp,"\nActual variance of pairwise differentiation (method of Ritland 2000)\nEstimate");
			if(l==m+4) fprintf(fp,"\nSE (jackknife over loci)");
			

			if(l<=m+2){
				for(r=1;r<=4;r++){
					if(FstatSlr[S][l][r]==MISSVAL) fprintf(fp,"\t");
					else fprintf(fp,"\t%.4f",FstatSlr[S][l][r]);					
				}
				fprintf(fp,"\t");				
			}
			else{
				fprintf(fp,"\t\t\t\t\t");
				
			}

			if(PWstat==2){
				for(c=1;c<=nc+1;c++){
					if(nc==1 && c==nc+1) fprintf(fp,"\t");
					else{
						if(corrSlc[S][l][c]==MISSVAL) fprintf(fp,"\t");
						else if(l<=m+2) fprintf(fp,"\t%.4f",corrSlc[S][l][c]);
						else fprintf(fp,"\t%G",corrSlc[S][l][c]);
					}				
				}
				fprintf(fp,"\t");
				
				if(l<=m+2) for(c=-20;c<=-1;c++){
					if(regdetails && (c==-12 || c>=-2)) continue;
					if(!regdetails && c<-2) c=-2;
					if(corrSlc[S][l][c]==MISSVAL || (l>m && (c==-8 ||c==-7 ||c==-6 ||c==-5 ||c==-3 ||c==-18 ||c==-17 ||c==-16 ||c==-15 ||c==-13)) ) fprintf(fp,"\t");
					else fprintf(fp,"\t%G",corrSlc[S][l][c]);
					
				}
			}
			else if(PWstat==1){
				if(corrSlc[S][l][nc+1]==MISSVAL) fprintf(fp,"\t");
				else fprintf(fp,"\t%.4f",corrSlc[S][l][nc+1]);
				
			}
		}/*end of loop l*/

		if(Rbtwloc){
			fprintf(fp,"\n\n\t\t\t\t\tMean inter-locus correlation coefficient for pairwise coefficients (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G\t\tVp1=\t%G\tSE(Vp1)=\t%G\tRl1=\t%G\tSE(Rl1)=\t%G\t\tVp2=\t%G\tSE(Vp2)=\t%G\tRl2=\t%G\tSE(Rl2)=\t%G\tVl=\t%G\tSE(Vl)=\t%G",RSll[S][0][0][0],RSll[S][0][2][0],V[S][0][0][0],V[S][0][0][2],V[S][0][1][0],V[S][0][1][2],V[S][0][2][0],V[S][0][2][2],V[S][0][4][0],V[S][0][4][2],V[S][0][3][0],V[S][0][3][2]);			
			fprintf(fp,"\n\t\t\t\t\tMean inter-locus correlation coefficient for residuals of pairwise coefficients after accounting for means per class (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G\t\tVp1=\t%G\tSE(Vp1)=\t%G\tRl1=\t%G\tSE(Rl1)=\t%G\t\tVp2=\t%G\tSE(Vp2)=\t%G\tRl2=\t%G\tSE(Rl2)=\t%G\tVl=\t%G\tSE(Vl)=\t%G",RSll[S][1][0][0],RSll[S][1][2][0],V[S][1][0][0],V[S][1][0][2],V[S][1][1][0],V[S][1][1][2],V[S][1][2][0],V[S][1][2][2],V[S][1][4][0],V[S][1][4][2],V[S][1][3][0],V[S][1][3][2]);			

			fprintf(fp,"\n\n\t\t\t\t\tMean correlation coefficient between 2 random pool of loci for pairwise coefficients (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G",R2pl[S][0][0],R2pl[S][0][2]);			
			fprintf(fp,"\n\t\t\t\t\tMean correlation coefficient between 2 random pool of loci for residuals of pairwise coefficients after accounting for means per class (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G",R2pl[S][1][0],R2pl[S][1][2]);			
		
		}

/*		if(Rbtwloc){
			fprintf(fp,"\n\n\t\t\t\t\tMean inter-locus correlation coefficient for pairwise coefficients (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G",RSll[S][0][0][0],RSll[S][0][2][0]);			
			fprintf(fp,"\n\t\t\t\t\tMean inter-locus correlation coefficient for residuals of pairwise coefficients after accounting for means per class (Pairon & Hardy, unpublished)\n\t\t\t\t\t%G\tSE (Jackknife over loci)=\t%G",RSll[S][1][0][0],RSll[S][1][2][0]);			
		}
*/
		
	}/*end of loop S*/

	fclose (fp);

}	/*end of writePopStatresults*/


/**************************************************************************************/



void writedistmatrices (char *outputfilename,int n,int m,float givenF,int TypeComp, int *cati,
		int printdistmatrix,double *xi,double *yi,double *zi,double **Mdij,int *sgi,
		int StatType,int NS,int Stat[12],float ***corrSlij[],struct name namei[],
		char namelocus[][MAXNOM])
{	
	int i,j,S,l,indivF,intraPop;
	char smess[SMAX];
	char ch;
	float dij;
	FILE *fp;
	
	
	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}
	
	if(printdistmatrix<0){ /*write pairwise results in matrix form*/
		printf("\nWriting pairwise distances in matrix form. Please wait.");
		fprintf(fp,"\n\nPAIRWISE SPATIAL AND GENETIC DISTANCES written as matrices");
		fprintf(fp,"\n\nMatrix of pairwise spatial distances");
		if(Mdij[0][0]==-1.) fprintf(fp," (km)");
		if(StatType==1) fprintf(fp," between INDIVIDUALS (a -1 value is given for pairs of individuals belonging to a same spatial group)");
		if(StatType==2) fprintf(fp," between POPULATIONS (1 pop = 1 CATEGORICAL GROUP)");
		if(StatType==3) fprintf(fp," between POPULATIONS (1 pop = 1 SPATIAL GROUP)");
		if(StatType==4) fprintf(fp," between POPULATIONS (1 pop = 1 SPATIAL GROUP of given CATEGORY)");
		fprintf(fp,"\n");

		for(i=1;i<=n;i++) fprintf(fp,"\t%s",namei[i].n);
		for(j=1;j<=n;j++){
			fprintf(fp,"\n%s",namei[j].n);
			for(i=1;i<=n;i++){
				if(Mdij[0][0]==1.) dij=(float)Mdij[i][j];
				else{
					if(Mdij[0][0]==-1.) dij=(float)acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
					else dij=(float)sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
					if(StatType==1 && sgi[i]==sgi[j]) dij=-1.0f;
				}
				fprintf(fp,"\t%G",dij);
			}
		}
		
		for(S=1;S<=NS;S++){
			if(StatType==1){
				if(Stat[S]==1) sprintf(smess,"\n\nPairwise KINSHIP coefficients (Loiselle et al., 1995)");
				if(Stat[S]==2) sprintf(smess,"\n\nPairwise KINSHIP coefficients (Ritland, 1996)");
				if(Stat[S]==3) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
				if(Stat[S]==4) sprintf(smess,"\n\nPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
				if(Stat[S]==5) sprintf(smess,"\n\nPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
				if(Stat[S]==6) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
				if(Stat[S]==7) sprintf(smess,"\n\nPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
				if(Stat[S]==8) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
				if(Stat[S]==9) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
				if(Stat[S]==10) sprintf(smess,"\n\nPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
				if(Stat[S]==11) sprintf(smess,"\n\nPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==12) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==13) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Li et al., 1993)");
				if(Stat[S]==14) sprintf(smess,"\n\nPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");
			}
			else{
				if(Stat[S]==1) sprintf(smess,"\n\nPairwise Fst (ANOVA approach)");
				if(Stat[S]==2) sprintf(smess,"\n\nPairwise Fst/(1-Fst) (ANOVA approach)");
				if(Stat[S]==3) sprintf(smess,"\n\nPairwise Rho (ANOVA approach)");
				if(Stat[S]==4) sprintf(smess,"\n\nPairwise Rho/(1-Rho) (ANOVA approach)");
				if(Stat[S]==5) sprintf(smess,"\n\nPairwise Rst (ANOVA approach)");
				if(Stat[S]==6) sprintf(smess,"\n\nPairwise Rst/(1-Rst) (ANOVA approach)");
				if(Stat[S]==7) sprintf(smess,"\n\nPairwise Ds (Nei's 1978 standard distance)");
				if(Stat[S]==8) sprintf(smess,"\n\nPairwise dm2 (Goldstein's 1995 distance)");
				if(Stat[S]==9) sprintf(smess,"\n\nPairwise Gst (Pons & Petit 1995)");
				if(Stat[S]==10) sprintf(smess,"\n\nPairwise Gij (Hardy)");
				if(Stat[S]==11) sprintf(smess,"\n\nPairwise Nst (Pons & Petit 1996)");
				if(Stat[S]==12) sprintf(smess,"\n\nPairwise Nij (OJ Hardy, unpublished)");

			}
			fprintf(fp,"%s",smess);

			for(l=0;l<=m;l++){
				if(m==1) l=1;
				if(abs(printdistmatrix)==1 && m>1 && l>0) break;
				if(l==0) sprintf(smess,"\nALL LOCI");
				if(l>0) sprintf(smess,"\n%s",namelocus[l]);
				fprintf(fp,"\n%s",smess);

				for(i=1;i<=n;i++) fprintf(fp,"\t%s",namei[i].n);

				for(j=1;j<=n;j++){
					fprintf(fp,"\n%s",namei[j].n);

					for(i=1;i<=n;i++){
						if(i==j) sprintf(smess,"\t");
						if(i<j){
							if(corrSlij[S][l][i][j]!=MISSVAL) sprintf(smess,"\t%.4f",corrSlij[S][l][i][j]);
							else sprintf(smess,"\t");
						}
						if(i>j){ 
							if(corrSlij[S][l][j][i]!=MISSVAL) sprintf(smess,"\t%.4f",corrSlij[S][l][j][i]);
							else sprintf(smess,"\t");
						}
						fprintf(fp,"%s",smess);
					}
				}
			}/*end loop l*/

		}/*end of loop S*/
	}/*end of if(printdistmatrix<0)*/

	

	/*write pairwise results in column form*/
	printf("\nWriting pairwise distances in column form. Please wait.");
	fprintf(fp,"\n\n\nPAIRWISE SPATIAL AND GENETIC DISTANCES written in column form");
	fprintf(fp,"\nName i\tName j\tN°i\tN°j\tSpatial dist");
	if(Mdij[0][0]==-1.) fprintf(fp," (km)");
	for(S=1;S<=NS;S++){
		if(StatType==1){
			if(Stat[S]==1) sprintf(smess,"\tPairwise KINSHIP coefficients (Loiselle et al., 1995)");
			if(Stat[S]==2) sprintf(smess,"\tPairwise KINSHIP coefficients (Ritland, 1996)");
			if(Stat[S]==3) sprintf(smess,"\tPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
			if(Stat[S]==4) sprintf(smess,"\tPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
			if(Stat[S]==5) sprintf(smess,"\tPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
			if(Stat[S]==6) sprintf(smess,"\tPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
			if(Stat[S]==7) sprintf(smess,"\tPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
			if(Stat[S]==8) sprintf(smess,"\tPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
			if(Stat[S]==9) sprintf(smess,"\tPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
			if(Stat[S]==10) sprintf(smess,"\tPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
			if(Stat[S]==11) sprintf(smess,"\tPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
			if(Stat[S]==12) sprintf(smess,"\tPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
			if(Stat[S]==13) sprintf(smess,"\tPairwise RELATIONSHIP coefficients (Li et al., 1993)");
			if(Stat[S]==14) sprintf(smess,"\tPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");
		}
		if(StatType>=2){
			if(Stat[S]==1) sprintf(smess,"\tFst");
			if(Stat[S]==2) sprintf(smess,"\tFst/(1-Fst)");
			if(Stat[S]==3) sprintf(smess,"\tRho");
			if(Stat[S]==4) sprintf(smess,"\tRho/(1-Rho)");
			if(Stat[S]==5) sprintf(smess,"\tRst");
			if(Stat[S]==6) sprintf(smess,"\tRst/(1-Rst)");
			if(Stat[S]==7) sprintf(smess,"\tDs (Nei's 1978 standard distance)");
			if(Stat[S]==8) sprintf(smess,"\tdm2 (Goldstein's 1995 distance)");
			if(Stat[S]==9) sprintf(smess,"\tGst");
			if(Stat[S]==10) sprintf(smess,"\tGij (OJ Hardy, unpublished)");
			if(Stat[S]==11) sprintf(smess,"\tNst");
			if(Stat[S]==12) sprintf(smess,"\tNij (OJ Hardy, unpublished)");

		}
		fprintf(fp,"%s",smess);
		if(m>1) fprintf(fp,"\tALL LOCI");
		if(m==1 || abs(printdistmatrix)==2) for(l=1;l<=m;l++) fprintf(fp,"\t%s",namelocus[l]);
	}/*end of loop S*/
/*	for(i=1;i<n;i++)for(j=i+1;j<=n;j++){	*/
	indivF=intraPop=0;
	if(StatType==1) for(S=1;S<=NS;S++)if(Stat[S]==1 || Stat[S]==2 || Stat[S]==5) indivF=1;
	if(StatType>1) for(S=1;S<=NS;S++)if(Stat[S]==10 || Stat[S]==12) intraPop=1;
	for(i=1;i<=n;i++)for(j=i;j<=n;j++){

		if(i==j)if(indivF==0 && intraPop==0) continue;
		if(TypeComp==1)if(cati[i]!=cati[j]) continue;
		if(TypeComp==2)if(cati[i]==cati[j]) continue;

		if(Mdij[0][0]==1.) dij=(float)Mdij[i][j];
		else{
			if(Mdij[0][0]==-1.) dij=(float)acos(sin(xi[i])*sin(xi[j])+cos(xi[i])*cos(xi[j])*cos(yi[j]-yi[i]))*6371;
			else dij=(float)sqrt( (xi[i]-xi[j])*(xi[i]-xi[j])+(yi[i]-yi[j])*(yi[i]-yi[j])+(zi[i]-zi[j])*(zi[i]-zi[j]) );
			if(StatType==1 && sgi[i]==sgi[j]) dij=-1.0f;
		}

		if(i==j && StatType==1) fprintf(fp,"\n%s\t%s\t%i\t%i\tinbreeding coef",namei[i].n,namei[j].n,i,j);
		else fprintf(fp,"\n%s\t%s\t%i\t%i\t%G",namei[i].n,namei[j].n,i,j,dij);

		for(S=1;S<=NS;S++){
			fprintf(fp,"\t");
			for(l=0;l<=m;l++){
				if(m==1) l=1;
				if(abs(printdistmatrix)==1 && m>1 && l>0) break;
				if(i==j && ((StatType==1 && Stat[S]!=1 && Stat[S]!=2 && Stat[S]!=5)||(StatType>1 && Stat[S]!=10 && Stat[S]!=12)) ) fprintf(fp,"\t");
				else {
					if(corrSlij[S][l][i][j]!=MISSVAL) fprintf(fp,"\t%.6f",corrSlij[S][l][i][j]);
					else fprintf(fp,"\t");
				}
			}
		}
	}



	if(printdistmatrix==3){ /*write pairwise results in matrix in Phylip format*/
		printf("\nWriting pairwise distances in Phylip format. Please wait.");
		fprintf(fp,"\n\nPAIRWISE GENETIC DISTANCES written as matrices in PHYLIP format\n   note 1: the distance between an element and itself is always set to zero\n   note 2: negative values are replaced by -0.0000");

		for(S=1;S<=NS;S++){
			if(StatType==1){
				if(Stat[S]==1) sprintf(smess,"\n\nPairwise KINSHIP coefficients (Loiselle et al., 1995)");
				if(Stat[S]==2) sprintf(smess,"\n\nPairwise KINSHIP coefficients (Ritland, 1996)");
				if(Stat[S]==3) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
				if(Stat[S]==4) sprintf(smess,"\n\nPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
				if(Stat[S]==5) sprintf(smess,"\n\nPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
				if(Stat[S]==6) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
				if(Stat[S]==7) sprintf(smess,"\n\nPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
				if(Stat[S]==8) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
				if(Stat[S]==9) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
				if(Stat[S]==10) sprintf(smess,"\n\nPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
				if(Stat[S]==11) sprintf(smess,"\n\nPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==12) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==13) sprintf(smess,"\n\nPairwise RELATIONSHIP coefficients (Li et al., 1993)");
				if(Stat[S]==14) sprintf(smess,"\n\nPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");
			}
			else{
				if(Stat[S]==1) sprintf(smess,"\n\nPairwise Fst (ANOVA approach)");
				if(Stat[S]==2) sprintf(smess,"\n\nPairwise Fst/(1-Fst) (ANOVA approach)");
				if(Stat[S]==3) sprintf(smess,"\n\nPairwise Rho (ANOVA approach)");
				if(Stat[S]==4) sprintf(smess,"\n\nPairwise Rho/(1-Rho) (ANOVA approach)");
				if(Stat[S]==5) sprintf(smess,"\n\nPairwise Rst (ANOVA approach)");
				if(Stat[S]==6) sprintf(smess,"\n\nPairwise Rst/(1-Rst) (ANOVA approach)");
				if(Stat[S]==7) sprintf(smess,"\n\nPairwise Ds (Nei's 1978 standard distance)");
				if(Stat[S]==8) sprintf(smess,"\n\nPairwise dm2 (Goldstein's 1995 distance)");
				if(Stat[S]==9) sprintf(smess,"\n\nPairwise Gst (Pons & Petit 1995)");
				if(Stat[S]==10) sprintf(smess,"\n\nPairwise Gij (OJ Hardy, unpublished)");
				if(Stat[S]==11) sprintf(smess,"\n\nPairwise Nst (Pons & Petit 1996)");
				if(Stat[S]==12) sprintf(smess,"\n\nPairwise Nij (OJ Hardy, unpublished)");

			}
			fprintf(fp,"%s",smess);
			fprintf(fp,"\n\n%d",n);
			l=0;
			if(m==1) l=1;
			if(l==0) sprintf(smess,"\nALL LOCI");

			for(j=1;j<=n;j++){
				fprintf(fp,"\n%-10.10s",namei[j].n);

				for(i=1;i<=n;i++){
					if(i==j) sprintf(smess,"%10.4f",0.0);
					if(i<j){
						if(corrSlij[S][l][i][j]==MISSVAL) sprintf(smess,"%10s","undefined");
						else if(corrSlij[S][l][i][j]>=0.0) sprintf(smess,"%10.4f",corrSlij[S][l][i][j]);
						else sprintf(smess,"   -0.0000");
					}
					if(i>j){ 
						if(corrSlij[S][l][j][i]==MISSVAL) sprintf(smess,"%10s","undefined");
						else if(corrSlij[S][l][j][i]>=0.0) sprintf(smess,"%10.4f",corrSlij[S][l][j][i]);
						else sprintf(smess,"   -0.0000");
					}
					fprintf(fp,"%s",smess);
				}
			}

		}/*end of loop S*/
	}/*end of if(printdistmatrix==3)*/

	fclose(fp);

}
/***********************************************************************/
/*************************************************************************************/



void WriteIndPermutRes(char *outputfilename,int n,int ploidy,int Ncat,
		double *maxc,int m,char namelocus[][MAXNOM],int TypeComp,float givenF,
		int nc,int Nsg,int NS,int Stat[],int Npermut[],int permutalleles,
		int permutdetails,struct resample_stat_type **r_statSlc[],long seedinit)

{
	int c,l=0,cinit,S;
	float Pval;
	char smess[SMAX];
	char *permutmode[70];

	FILE *fp;

	permutmode[0]="";
	permutmode[10]="ASaAwL";	/*allele sizes among alleles within locus*/
	permutmode[11]="RCoDMbA";	/*permutation of rows and columns of distance matrices between alleles*/
	permutmode[20]="GaI";		/*genes among all individuals*/
	permutmode[21]="GaIwC";		/*genes among individuals within category*/
	permutmode[22]="GaIwP";		/*genes among individuals within population*/
	permutmode[30]="IaSG";		/*individuals among spatial groups*/
	permutmode[31]="IaSGwC";	/*individuals among spatial groups within categories*/
	permutmode[32]="IaP";		/*individuals among all populations*/
	permutmode[33]="IaPwC";		/*individuals among populations within category*/
	permutmode[40]="ILaI";		/*individual locations among all individuals*/
	permutmode[41]="ILaIwC";	/*individual locations among individuals within category*/
	permutmode[50]="SGLaSG";	/*spatial group locations among all spatial groups*/
	permutmode[51]="SGLaSGwC";	/*spatial group locations among spatial groups within category*/
	permutmode[60]="PLaP";		/*population locations among all populations*/
	permutmode[61]="PLaPwC";	/*population locations among population within category*/
	
	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}

	if(permutalleles) fprintf(fp,"\n\nALLELES PERMUTATION TESTS -> Ho: obs=exp (exp = mean value after permutation)  Note: all alleles found in the data file are considered for permutation");
	else fprintf(fp,"\n\nLOCATIONS, INDIVIDUALS and/or GENES PERMUTATION TESTS -> Ho: obs=exp (exp = mean value after permutation)");
 	
	if(TypeComp==2 && !permutalleles) fprintf(fp,"\nWARNING: Significant tests for spatial structuring (location permutations) in the case of pairwise comparisons among categories cannot be interpreted as a demonstration that gene flow occurs among categories (Hardy and Vekemans 2001)");		

	fprintf(fp,"\nInitial seed= %i",abs(seedinit));
 	
	if(permutdetails==0){
		if(permutalleles) fprintf(fp,"\nP values of 2-sided tests (obs<>exp) after %i random permutations of allele sizes",Npermut[4]);
		else fprintf(fp,"\nP values of 2-sided tests (obs<>exp) after %i random permutations of locations, %i random permutations of individuals, and %i random permutations of genes",Npermut[1],Npermut[2],Npermut[3]);	 	
	}

	for(S=1;S<=NS;S++){
		if(permutalleles)if(Stat[S]!=5 && Stat[S]!=14) continue;
		if(permutdetails==0 || (permutdetails==1 && l!=0)){
			fprintf(fp,"\nSTATISTIC: ");			
			if(Stat[S]==1) fprintf(fp,"\tPairwise KINSHIP coefficients (Loiselle et al., 1995)");
			if(Stat[S]==2) fprintf(fp,"\tPairwise KINSHIP coefficients (Ritland, 1996)");
			if(Stat[S]==3) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
			if(Stat[S]==4) fprintf(fp,"\tPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
			if(Stat[S]==5) fprintf(fp,"\tPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
			if(Stat[S]==6) fprintf(fp,"\tPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
			if(Stat[S]==7) fprintf(fp,"\tPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
			if(Stat[S]==8) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
			if(Stat[S]==9) fprintf(fp,"\tPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
			if(Stat[S]==10) fprintf(fp,"\tPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
			if(Stat[S]==11) fprintf(fp,"\tPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
			if(Stat[S]==12) fprintf(fp,"\tPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
			if(Stat[S]==13) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Li et al., 1993)");
			if(Stat[S]==14) fprintf(fp,"\tPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");

			for(c=1;c<=nc+2;c++){fprintf(fp,"\t");}
			fprintf(fp,"\tSlopes of regression analyses");
			if(Nsg>1) fprintf(fp,"\nLocus\tintra-individual\tintra-group");
			else fprintf(fp,"\nLocus\tintra-individual\t1");	
			for(c=2;c<=nc;c++) fprintf(fp,"\t%i",c);
			if(Nsg>1) cinit=2;
			else cinit=1;
 			fprintf(fp,"\t\t\tb-lin (slope linear dist)\tb-log (slope ln(dist))");
			
		}

		for(l=0;l<=m;l++){
			if(m==1) l=1;

			if(permutdetails==0 || (permutdetails==1 && l!=0)){
				if(l==0) fprintf(fp,"\nALL LOCI");
				if(l>0) fprintf(fp,"\n%s",namelocus[l]);
				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].pbil==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%.4f",r_statSlc[S][l][c].pbil);					
				}
				fprintf(fp,"\t\t");
				
				for(c=-2;c<=-1;c++){
					if(Stat[S]==4) Pval=r_statSlc[S][l][c].phigh;
					else Pval=r_statSlc[S][l][c].plow;
					if(Pval==MISSVAL) fprintf(fp,"\t");
					else fprintf(fp,"\t%.4f",Pval);					
				}
			}

			/*write complete analysis*/
			if(permutdetails==2 || (permutdetails==1 && l==0)){

				if(l==0) fprintf(fp,"\n\nALL LOCI");
				if(l>0) fprintf(fp,"\n%s",namelocus[l]);
				
				if(Stat[S]==1) fprintf(fp,"\tPairwise KINSHIP coefficients (Loiselle et al., 1995)");
				if(Stat[S]==2) fprintf(fp,"\tPairwise KINSHIP coefficients (Ritland, 1996)");
				if(Stat[S]==3) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)");
				if(Stat[S]==4) fprintf(fp,"\tPairwise Rousset's distance between individuals ('a' in Rousset, 2000)");
				if(Stat[S]==5) fprintf(fp,"\tPairwise correlation coefficients of allele size (I' in Streiff et al., 1998)");
				if(Stat[S]==6) fprintf(fp,"\tPairwise RELATIONSHIP coefficients ('r' in Lynch & Ritland, 1999)");
				if(Stat[S]==7) fprintf(fp,"\tPairwise FRATERNITY coefficients ('delta' in Lynch & Ritland, 1999)");
				if(Stat[S]==8) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Queller & Goodnight, 1989)");
				if(Stat[S]==9) fprintf(fp,"\tPairwise RELATIONSHIP coefficients ('r' in Wang, 2002)");
				if(Stat[S]==10) fprintf(fp,"\tPairwise FRATERNITY coefficients ('delta' in Wang, 2002)");
				if(Stat[S]==11) fprintf(fp,"\tPairwise KINSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==12) fprintf(fp,"\tPairwise RELATIONSHIP coefficients for a DOMINANT marker, assuming inbreeding coef=%.4f (Hardy, 2003)",givenF);
				if(Stat[S]==13) fprintf(fp,"\tPairwise RELATIONSHIP coefficients (Li et al., 1993)");
				if(Stat[S]==14) fprintf(fp,"\tPairwise KINSHIP coefficients for ORDERED alleles (OJ Hardy, unpublished)");
				
				for(c=1;c<=nc+2;c++) fprintf(fp,"\t");
				fprintf(fp,"\tSlopes of regression analyses");

				if(Nsg>1) fprintf(fp,"\n\tintra-individual\tintra-group");
				else fprintf(fp,"\n\tintra-individual\t1");
				
				for(c=2;c<=nc;c++) fprintf(fp,"\t%i",c);
					
				if(Nsg>1) cinit=2;
				else cinit=1;
 				fprintf(fp,"\t\t\tb-lin (slope linear dist)\tb-log (slope ln(dist))");
							
				fprintf(fp,"\nObject permuted");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][0][c].mode==0) fprintf(fp,"\t");
					else fprintf(fp,"\t%s",permutmode[r_statSlc[S][0][c].mode]);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nN valid permut");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].n==0) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%i",r_statSlc[S][l][c].n);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nN different permut val");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].nd==0) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%i",r_statSlc[S][l][c].nd);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nObs val");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].obs==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%G",r_statSlc[S][l][c].obs);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nMean permut val");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].mean==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%G",r_statSlc[S][l][c].mean);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nSD permut val");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].sd==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%G",r_statSlc[S][l][c].sd);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\n95%%CI-inf");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].low95==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%G",r_statSlc[S][l][c].low95);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\n95%%CI-sup");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].high95==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%G",r_statSlc[S][l][c].high95);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nP(1-sided test, H1: obs<exp)");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].plow==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%.4f",r_statSlc[S][l][c].plow);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nP(1-sided test, H1: obs>exp)");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].phigh==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%.4f",r_statSlc[S][l][c].phigh);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}
				fprintf(fp,"\nP(2-sided test, H1: obs<>exp)");				
				for(c=0;c<=nc;c++){
					if(r_statSlc[S][l][c].pbil==MISSVAL) fprintf(fp,"\t"); 
					else fprintf(fp,"\t%.4f",r_statSlc[S][l][c].pbil);					
					if(c==nc){c=-3; fprintf(fp,"\t\t"); }
					if(c==-1) break;
				}			
				fprintf(fp,"\n");				
			}
		}/*end of loop l for variables*/
		fprintf(fp,"\n");
		
	}/*end of loop S*/

	fclose (fp);

}/*end of procedure WriteIndPermutRes*/

/*************************************************************************************/

void WritePopPermutRes(char *outputfilename,int n,int ploidy,
		double *maxc,int m,char namelocus[][MAXNOM],int TypeComp,
		int nc,int NS,int Stat[],int Npermut[],int permutalleles,
		int permutdetails,struct resample_stat_type **r_statFSlr[],
		struct resample_stat_type **r_statSlc[],long seedinit,int PWstat)

{
	int c,l=0,cinit,S,r,ri,rf;
	float Pval;
	char smess[SMAX];
	char *permutmode[70];

	FILE *fp;

	permutmode[0]="";
	permutmode[10]="ASaAwL";	/*allelic sizes within locus*/
	permutmode[11]="RCoDMbA";	/*rows and columns of distance matrices between alleles*/
	permutmode[20]="GaI";		/*genes among all individuals*/
	permutmode[21]="GaIwC";		/*genes among individuals within category*/
	permutmode[22]="GaIwP";		/*genes among individuals within population*/
	permutmode[30]="IaSG";		/*individuals among spatial groups*/
	permutmode[31]="IaSGwC";	/*individuals among spatial groups within categories*/
	permutmode[32]="IaP";		/*individuals among all populations*/
	permutmode[33]="IaPwC";		/*individuals among populations within category*/
	permutmode[40]="ILaI";		/*individual locations among all individuals*/
	permutmode[41]="ILaIwC";	/*individual locations among individuals within category*/
	permutmode[50]="SGLaSG";	/*spatial group locations among all spatial groups*/
	permutmode[51]="SGLaSGwC";	/*spatial group locations among spatial groups within category*/
	permutmode[60]="PLaP";		/*population locations among all populations*/
	permutmode[61]="PLaPwC";	/*population locations among population within category*/
	
	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file %c%s%c.\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",'"',outputfilename,'"');
		wait_a_char();
	}
	fclose (fp);

	printf("\nWriting permutation results. Please, wait.");

	if(permutalleles) sprintf(smess,"\n\nALLELES PERMUTATION TESTS -> Ho: obs=exp (exp = mean value after permutation)   Note: all alleles found in the data file are considered for permutation (for Nst, the alleles from the distance matrix that are not found in the data file are not considered)");
	else sprintf(smess,"\n\nLOCATIONS, INDIVIDUALS and/or GENES PERMUTATION TESTS -> Ho: obs=exp (exp = mean value after permutation)");
 	write_tofile_only(outputfilename,smess);
	if(TypeComp==2 && !permutalleles){
		sprintf(smess,"\nWARNING: Significant tests for spatial structuring (location permutations) in the case of pairwise comparisons among categories cannot be interpreted as a demonstration that gene flow occurs among categories (Hardy and Vekemans 2001)");
 		write_tofile_only(outputfilename,smess);
	}
	sprintf(smess,"\nInitial seed= %i",abs(seedinit));
 	write_tofile_only(outputfilename,smess);
	if(permutdetails==0){
		if(permutalleles) sprintf(smess,"\nP values of 2-sided tests (obs<>exp) after %i random permutations of allele sizes",Npermut[4]);
		else sprintf(smess,"\nP values of 2-sided tests (obs<>exp) after %i random permutations of locations, %i random permutations of individuals, and %i random permutations of genes",Npermut[1],Npermut[2],Npermut[3]);
	 	write_tofile_only(outputfilename,smess);
	}

	for(S=1;S<=NS;S++){
		if(permutalleles)if(Stat[S]!=5 && Stat[S]!=6 && Stat[S]!=8 && Stat[S]!=11 && Stat[S]!=12) continue;
		
		if(permutdetails==0 || (permutdetails==1 && l!=0)){
			sprintf(smess,"\nSTATISTICS:");
			write_tofile_only(outputfilename,smess);
			if(PWstat){
				if(Stat[S]==1) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Fst");
				if(Stat[S]==2) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Fst/(1-Fst)");
				if(Stat[S]==3) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Rho");
				if(Stat[S]==4) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Rho/(1-Rho)");
				if(Stat[S]==5) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Rst");
				if(Stat[S]==6) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Rst/(1-Rst)");
				if(Stat[S]==7) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Nei's (1978) Ds standard distance");
				if(Stat[S]==8) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Goldstein's (1995) dm2 distance");
				if(Stat[S]==9) sprintf(smess,"\tGlobal Gst \t\t\t\t\tPairwise Gst");
				if(Stat[S]==10) sprintf(smess,"\tGlobal Gst \t\t\t\t\tPairwise Gij (OJ Hardy, unpublished)");
				if(Stat[S]==11) sprintf(smess,"\tGlobal Nst \t\t\t\t\tPairwise Nst");
				if(Stat[S]==12) sprintf(smess,"\tGlobal Nst \t\t\t\t\tPairwise Nij (OJ Hardy, unpublished)");
			}
			else{
				if(Stat[S]==1) sprintf(smess,"\tGlobal F-statistics");
				if(Stat[S]==2) sprintf(smess,"\tGlobal F-statistics");
				if(Stat[S]==3) sprintf(smess,"\tGlobal F-statistics");
				if(Stat[S]==4) sprintf(smess,"\tGlobal F-statistics");
				if(Stat[S]==5) sprintf(smess,"\tGlobal R-statistics");
				if(Stat[S]==6) sprintf(smess,"\tGlobal R-statistics");
				if(Stat[S]==7) sprintf(smess,"\tGlobal F-statistics");
				if(Stat[S]==8) sprintf(smess,"\tGlobal R-statistics");
				if(Stat[S]==9) sprintf(smess,"\tGlobal Gst");
				if(Stat[S]==10) sprintf(smess,"\tGlobal Gst");
				if(Stat[S]==11) sprintf(smess,"\tGlobal Nst");
				if(Stat[S]==12) sprintf(smess,"\tGlobal Nst");
			}
			write_tofile_only(outputfilename,smess);
			if(Stat[S]<=2 || Stat[S]==7) sprintf(smess,"\nLocus\tFit\tFis\tFst\t\t");
			if(Stat[S]==3 || Stat[S]==4) sprintf(smess,"\nLocus\tFit\tFis\tFst\tRho\t");
			if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8) sprintf(smess,"\nLocus\tRit\tRis\tRst\t\t");
			if(Stat[S]==9 || Stat[S]==10) sprintf(smess,"\nLocus\t\t\tGst\t\t");
			if(Stat[S]==11 || Stat[S]==12) sprintf(smess,"\nLocus\t\t\tNst\t\t");

			write_tofile_only(outputfilename,smess);
			if(PWstat){
				for(c=1;c<=nc;c++){
					sprintf(smess,"\t%i",c);
					write_tofile_only(outputfilename,smess);
				}
 				sprintf(smess,"\taverage\t\tb-lin (slope linear dist)\tb-log (slope ln(dist))");
				write_tofile_only(outputfilename,smess);
			}
		}

		
		for(l=0;l<=m;l++){
			if(m==1) l=1;
	
			if(permutdetails==0 || (permutdetails==1 && l!=0)){
				if(l==0) sprintf(smess,"\nALL LOCI");
				if(l>0) sprintf(smess,"\n%s",namelocus[l]);
				write_tofile_only(outputfilename,smess);

				for(r=1;r<=3;r++){
					if(r_statFSlr[S][l][r].pbil!=MISSVAL) sprintf(smess,"\t%.4f",r_statFSlr[S][l][r].pbil);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				if(Stat[S]==3 || Stat[S]==4) sprintf(smess,"\t");
				else sprintf(smess,"\t\t"); 
				write_tofile_only(outputfilename,smess);
				if(PWstat)for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].pbil==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][c].pbil);
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat)for(c=-2;c<=-1;c++){
					if(r_statSlc[S][l][c].phigh==MISSVAL) sprintf(smess,"\t");
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][c].phigh);
					write_tofile_only(outputfilename,smess);
				}
			}

			/*write complete analysis*/
			if(permutdetails==2 || (permutdetails==1 && l==0)){
				if(l==0) sprintf(smess,"\n\nALL LOCI");
				if(l>0) sprintf(smess,"\n%s",namelocus[l]);
				write_tofile_only(outputfilename,smess);

				if(PWstat){
					if(Stat[S]==1) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Fst");
					if(Stat[S]==2) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Fst/(1-Fst)");
					if(Stat[S]==3) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Rho");
					if(Stat[S]==4) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Rho/(1-Rho)");
					if(Stat[S]==5) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Rst");
					if(Stat[S]==6) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Rst/(1-Rst)");
					if(Stat[S]==7) sprintf(smess,"\tGlobal F-statistics \t\t\t\t\tPairwise Nei's (1978) Ds standard distance");
					if(Stat[S]==8) sprintf(smess,"\tGlobal R-statistics \t\t\t\t\tPairwise Goldstein's (1995) dm2 distance");
					if(Stat[S]==9) sprintf(smess,"\tGlobal Gst \t\t\t\t\tPairwise Gst");
					if(Stat[S]==10) sprintf(smess,"\tGlobal Gst \t\t\t\t\tPairwise Gij (OJ Hardy, unpublished)");
					if(Stat[S]==11) sprintf(smess,"\tGlobal Nst \t\t\t\t\tPairwise Nst");
					if(Stat[S]==12) sprintf(smess,"\tGlobal Nst \t\t\t\t\tPairwise Nij (OJ Hardy, unpublished)");
				}
				else{
					if(Stat[S]==1) sprintf(smess,"\tGlobal F-statistics");
					if(Stat[S]==2) sprintf(smess,"\tGlobal F-statistics");
					if(Stat[S]==3) sprintf(smess,"\tGlobal F-statistics");
					if(Stat[S]==4) sprintf(smess,"\tGlobal F-statistics");
					if(Stat[S]==5) sprintf(smess,"\tGlobal R-statistics");
					if(Stat[S]==6) sprintf(smess,"\tGlobal R-statistics");
					if(Stat[S]==7) sprintf(smess,"\tGlobal F-statistics");
					if(Stat[S]==8) sprintf(smess,"\tGlobal R-statistics");
					if(Stat[S]==9) sprintf(smess,"\tGlobal Gst");
					if(Stat[S]==10) sprintf(smess,"\tGlobal Gst");
					if(Stat[S]==11) sprintf(smess,"\tGlobal Nst");
					if(Stat[S]==12) sprintf(smess,"\tGlobal Nst");
				}
				write_tofile_only(outputfilename,smess);


				if(Stat[S]<=2 || Stat[S]==7) sprintf(smess,"\n\tFit\tFis\tFst\t\t");
				if(Stat[S]==3 || Stat[S]==4) sprintf(smess,"\n\tFit\tFis\tFst\tRho\t");
				if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8) sprintf(smess,"\n\tRit\tRis\tRst\t\t");
				if(Stat[S]==9 || Stat[S]==10) sprintf(smess,"\n\t\t\tGst\t\t");
				if(Stat[S]==11 || Stat[S]==12) sprintf(smess,"\n\t\t\tNst\t\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat){
					for(c=1;c<=nc;c++){
						sprintf(smess,"\t%i",c);
						write_tofile_only(outputfilename,smess);
					}
 					sprintf(smess,"\taverage\t\tb-lin (slope linear dist)\tb-log (slope ln(dist))");
					write_tofile_only(outputfilename,smess);
				}

				
				sprintf(smess,"\nObject permuted");
				write_tofile_only(outputfilename,smess);
				sprintf(smess,"\t%s\t%s\t%s\t",permutmode[r_statFSlr[S][0][1].mode],permutmode[r_statFSlr[S][0][2].mode],permutmode[r_statFSlr[S][0][3].mode]);
				write_tofile_only(outputfilename,smess);
				if(Stat[S]==3 || Stat[S]==4) sprintf(smess,"%s\t",permutmode[r_statFSlr[S][0][4].mode]);
				else sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][0][c].mode==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%s",permutmode[r_statSlc[S][0][c].mode]);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				sprintf(smess,"\nN valid permut");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].n) sprintf(smess,"\t%i",r_statFSlr[S][l][r].n);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].n==0) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%i",r_statSlc[S][l][c].n);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				sprintf(smess,"\nN different permut val");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].nd) sprintf(smess,"\t%i",r_statFSlr[S][l][r].nd);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].nd==0) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%i",r_statSlc[S][l][c].nd);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				sprintf(smess,"\nObs val");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].obs!=MISSVAL) sprintf(smess,"\t%G",r_statFSlr[S][l][r].obs);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].obs==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%G",r_statSlc[S][l][c].obs);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}
			
				sprintf(smess,"\nMean permut val");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].mean!=MISSVAL) sprintf(smess,"\t%G",r_statFSlr[S][l][r].mean);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].mean==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%G",r_statSlc[S][l][c].mean);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				sprintf(smess,"\nSD permut val");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].sd!=MISSVAL) sprintf(smess,"\t%G",r_statFSlr[S][l][r].sd);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].sd==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%G",r_statSlc[S][l][c].sd);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				sprintf(smess,"\n95%%CI-inf");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].low95!=MISSVAL) sprintf(smess,"\t%G",r_statFSlr[S][l][r].low95);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].low95==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%G",r_statSlc[S][l][c].low95);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}
			
				sprintf(smess,"\n95%%CI-sup");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].high95!=MISSVAL) sprintf(smess,"\t%G",r_statFSlr[S][l][r].high95);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].high95==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%G",r_statSlc[S][l][c].high95);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}
			
				sprintf(smess,"\nP(1-sided test, H1: obs<exp)");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].plow!=MISSVAL) sprintf(smess,"\t%.4f",r_statFSlr[S][l][r].plow);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].plow==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][c].plow);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}
			
				sprintf(smess,"\nP(1-sided test, H1: obs>exp)");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].phigh!=MISSVAL) sprintf(smess,"\t%.4f",r_statFSlr[S][l][r].phigh);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].phigh==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][c].phigh);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}
			
				sprintf(smess,"\nP(2-sided test, H1: obs<>exp)");
				write_tofile_only(outputfilename,smess);
				for(r=1;r<=4;r++){
					if(r_statFSlr[S][l][r].pbil!=MISSVAL) sprintf(smess,"\t%.4f",r_statFSlr[S][l][r].pbil);
					else sprintf(smess,"\t");
					write_tofile_only(outputfilename,smess);
				}
				sprintf(smess,"\t");
				write_tofile_only(outputfilename,smess);
				if(PWstat) for(c=1;c<=nc+1;c++){
					if(r_statSlc[S][l][c].pbil==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][c].pbil);
					write_tofile_only(outputfilename,smess);
					if(c==nc+1){c=-3; sprintf(smess,"\t"); write_tofile_only(outputfilename,smess);}
					if(c==-1) break;
				}

				
				sprintf(smess,"\n");
				write_tofile_only(outputfilename,smess);
			}
		}/*end of loop l for variables*/
		sprintf(smess,"\n");
		write_tofile_only(outputfilename,smess);
	}/*end of loop S*/

}/*end of procedure WritePopPermutRes*/

/*******************************************************************************/
void WriteAllelesPermutResForPopPair(char *outputfilename,int m,char namelocus[][MAXNOM],
		int p1,int p2,struct name *namepop,int NS,int Stat[],
		struct resample_stat_type **r_statSlc[],int permutdetails,int comp)

{
	int l,linit,S,k;
	char smess[SMAX];
	char *permutmode[70];

	FILE *fp;

	permutmode[0]="";
	permutmode[10]="ASwL";	/*allelelic sizes within locus*/
	permutmode[11]="RCoDMbA";	/*rows and columns of distance matrices between alleles*/
	permutmode[20]="GaI";		/*genes among all individuals*/
	permutmode[21]="GaIwC";		/*genes among individuals within category*/
	permutmode[22]="GaIwP";		/*genes among individuals within population*/
	permutmode[30]="IaSG";		/*individuals among spatial groups*/
	permutmode[31]="IaSGwC";	/*individuals among spatial groups within categories*/
	permutmode[32]="IaP";		/*individuals among all populations*/
	permutmode[33]="IaPwC";		/*individuals among populations within category*/
	permutmode[40]="ILaI";		/*individual locations among all individuals*/
	permutmode[41]="ILaIwC";		/*individual locations among individuals within category*/
	permutmode[50]="SGLaSG";	/*spatial group locations among all spatial groups*/
	permutmode[51]="SGLaSGwC";	/*spatial group locations among spatial groups within category*/
	permutmode[60]="PLaP";		/*population locations among all populations*/
	permutmode[61]="PLaPwC";	/*population locations among population within category*/
	
	while((fp=fopen(outputfilename,"a"))==NULL){
		printf("\nWARNING: Cannot open results file \"%s\".\nIf it is being used by another application, close it first. Then press RETURN.\nPress ctrl+c if you wish to stop the program now.\n",outputfilename);
		wait_a_char();
	}
	fclose (fp);

	linit=0;
	if(m==1) linit=1;

	if(comp==1){  /*write 3 first lines for the first call*/
		sprintf(smess,"\nALLELES PERMUTATION TESTS FOR EACH PAIR OF POPULATIONS -> Ho: obs=exp (exp = mean value after permutation)   Note: for Rst and Nst, only the alleles found within the compared populations are permuted; for Nij, all the alleles found in the data file are permuted (the alleles from the distance matrix that are not found in the data file are not considered)");
	 	write_tofile_only(outputfilename,smess);
		sprintf(smess,"\n\tStatistic ->\t\t");
		write_tofile_only(outputfilename,smess);
		for(S=1;S<=NS;S++) if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
			if(Stat[S]==5) sprintf(smess,"\tRst");
			if(Stat[S]==6) sprintf(smess,"\tRst/(1-Rst)");
			if(Stat[S]==8) sprintf(smess,"\tdm2");
			if(Stat[S]==11) sprintf(smess,"\tNst");
			if(Stat[S]==12) sprintf(smess,"\tNij");
			for(l=linit;l<=m;l++){
				if(permutdetails==2 || (permutdetails==1 && l==0)) for(k=1;k<=11;k++) write_tofile_only(outputfilename,smess);
				if(permutdetails==0 || (permutdetails==1 && l!=0)) write_tofile_only(outputfilename,smess);
			}
		}
		sprintf(smess,"\n\tLocus ->\t\t");
		write_tofile_only(outputfilename,smess);
		for(S=1;S<=NS;S++) if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
			for(l=linit;l<=m;l++){
				if(l==0) sprintf(smess,"\tALL LOCI");
				if(l>0) sprintf(smess,"\t%s",namelocus[l]);
				if(permutdetails==2 || (permutdetails==1 && l==0)) for(k=1;k<=11;k++) write_tofile_only(outputfilename,smess);
				if(permutdetails==0 || (permutdetails==1 && l!=0)) write_tofile_only(outputfilename,smess);
			}
		}
		sprintf(smess,"\nName i\tName j\tN°i\tN°j");
		write_tofile_only(outputfilename,smess);
		for(S=1;S<=NS;S++) if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
			for(l=linit;l<=m;l++){
				if(permutdetails==2 || (permutdetails==1 && l==0)) sprintf(smess,"\tObject permuted\tN valid permut\tN different permut val\tObs val\tMean permut val\tSD permut val\t95%%CI-inf\t95%%CI-sup\tP(1-sided test, H1: obs<exp)\tP(1-sided test, H1: obs>exp)\tP(2-sided test, H1: obs<>exp)");
				if(permutdetails==0 || (permutdetails==1 && l!=0)){
					if(Stat[S]==12) sprintf(smess,"\tP(2-sided test, H1: obs<>exp)");
					else sprintf(smess,"\tP(1-sided test, H1: obs>exp)");
				}
				write_tofile_only(outputfilename,smess);
			}
		}
	} /*end of if(comp==1)*/



	sprintf(smess,"\n%s\t%s\t%i\t%i",namepop[p1].n,namepop[p2].n,p1,p2);
	write_tofile_only(outputfilename,smess);
	for(S=1;S<=NS;S++)if(Stat[S]==5 || Stat[S]==6 || Stat[S]==8 || Stat[S]==11 || Stat[S]==12){
		for(l=linit;l<=m;l++){
			/*write limited analysis*/
			if(permutdetails==0 || (permutdetails==1 && l!=0)){
				if(Stat[S]==12){ 
					if(r_statSlc[S][l][0].pbil==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][0].pbil);
				}
				else{
					if(r_statSlc[S][l][0].phigh==MISSVAL) sprintf(smess,"\t"); 
					else sprintf(smess,"\t%.4f",r_statSlc[S][l][0].phigh);
				}
				write_tofile_only(outputfilename,smess);
			}
			/*write complete analysis*/
			if(permutdetails==2 || (permutdetails==1 && l==0)){
				
				if(r_statSlc[S][0][0].mode<=0 || r_statSlc[S][0][0].mode>=70) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%s",permutmode[r_statSlc[S][0][0].mode]);
				write_tofile_only(outputfilename,smess);

				if(r_statSlc[S][l][0].n==0) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%i",r_statSlc[S][l][0].n);
				write_tofile_only(outputfilename,smess);

				if(r_statSlc[S][l][0].nd==0) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%i",r_statSlc[S][l][0].nd);
				write_tofile_only(outputfilename,smess);

				if(r_statSlc[S][l][0].obs==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%G",r_statSlc[S][l][0].obs);
				write_tofile_only(outputfilename,smess);
			
				if(r_statSlc[S][l][0].mean==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%G",r_statSlc[S][l][0].mean);
				write_tofile_only(outputfilename,smess);
			
				if(r_statSlc[S][l][0].sd==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%G",r_statSlc[S][l][0].sd);
				write_tofile_only(outputfilename,smess);
			
				if(r_statSlc[S][l][0].low95==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%G",r_statSlc[S][l][0].low95);
				write_tofile_only(outputfilename,smess);
			
				if(r_statSlc[S][l][0].high95==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%G",r_statSlc[S][l][0].high95);
				write_tofile_only(outputfilename,smess);
			
				if(r_statSlc[S][l][0].plow==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%.4f",r_statSlc[S][l][0].plow);
				write_tofile_only(outputfilename,smess);

				if(r_statSlc[S][l][0].phigh==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%.4f",r_statSlc[S][l][0].phigh);
				write_tofile_only(outputfilename,smess);

				if(r_statSlc[S][l][0].pbil==MISSVAL) sprintf(smess,"\t"); 
				else sprintf(smess,"\t%.4f",r_statSlc[S][l][0].pbil);
				write_tofile_only(outputfilename,smess);
			}
		}/*end of loop l for variables*/
	}/*end of loop S*/

}/*end of procedure WriteAllelesPermutResForPopPair*/


/*******************************************************************************/

char *nexttab(char *s, char *inputfile, int line)
/*find next tab in string s and increment pointer s to just after the tab*/
{
    char *sn1;
    char smess[SMAX];
    
    sn1=strpbrk(s,"\t");
    if(!sn1) {
		sprintf(smess,"\nERROR: one column missing on line %d from file %s (tab not found)",line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
    	exit(1);
    }
    if(sn1[1]=='\t') {
		sprintf(smess,"\nERROR: two consecutive tabs on line %d from file %s",line,inputfile);
		write(ERRORFILE,smess);
 		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
    }
    return ++sn1;
}
/***********************************************************************/
	
	
char *nexttab_blank(char *s, char *inputfile, int line)
/*find next tab in string s and increment pointer s to just after the tab
same as nexttab but do not warn if two consecutive tabs are found*/
{
    char *sn1;
    char smess[SMAX];
    
    sn1=strpbrk(s,"\t");
    if(!sn1) {
		sprintf(smess,"\ntab not found one column missing on line %d from file %s",line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
    	exit(1);
    }
/*    if(sn1[1]=='\t') {
		sprintf(smess,"\ntwo consecutive tabs on line %d from file %s",line,inputfile);
		write(ERRORFILE,smess);
    	exit(1);
    }
 */
    return ++sn1;
}

/***********************************************************************/	

void readsfromfile(FILE *fp, char *s, char *inputfile, int *line)
{
	char *flag,smess[SMAX]; 

	(*line)++;
	flag = fgets(s,SMAX,fp);
	if(! flag) {
		if (feof(fp)) {
			sprintf(smess,"\nEnd of file unexpected on line %d from file %s",*line,inputfile);
			write(ERRORFILE,smess);
			printf("\npress any key to stop the program now");
			wait_a_char();
			exit(1);
		} else  {
			sprintf(smess,"\ncannot read line %d from file %s",*line,inputfile);
			write(ERRORFILE,smess);
			printf("\npress any key to stop the program now");
			wait_a_char();
			exit(1);
		}
	}
} /*end of readsfromfile*/

/***********************************************************************/	

void readsfromfile_no_end_of_file_check(FILE *fp, char *s, char *inputfile, int *line)
{
	char *flag,smess[SMAX]; 

	(*line)++;
	flag = fgets(s,SMAX,fp);
/*	if(! flag) {
		if (feof(fp)) {
			sprintf(smess,"\nEnd of file unexpected on line %d from file %s",*line,inputfile);
			write(ERRORFILE,smess);
			printf("\npress any key to stop the program now");
			wait_a_char();
			exit(1);
		} else  {
			sprintf(smess,"\ncannot read line %d from file %s",*line,inputfile);
			write(ERRORFILE,smess);
			printf("\npress any key to stop the program now");
			wait_a_char();
			exit(1);
		}
	}	   */
} /*end of readsfromfile*/

/***********************************************************************/	


void readsfromstring(char *ins, char *outs,int nchar, char *inputfile, int line) 
{
	int flag;
	char smess[SMAX],stemp[SMAX];
	
	flag = sscanf(ins,"%s",stemp);
	if (flag == EOF) {
		sprintf(smess,"\ncannot read string from \"%s\" on line %d from file \"%s\"",ins,line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
	}
	strncpy(outs,stemp,nchar);
	outs[nchar]='\0';	/*attention this asks for a string with length = nchar+1 at least*/
}  /*end readsfromstring*/

/***********************************************************************/	

void readintfromstring(char *ins, int *outi, char *inputfile, int line)
{
	int flag,temp;
	char smess[SMAX];
	
	flag = sscanf(ins,"%d",&temp);
	if (flag == EOF) {
		sprintf(smess,"\ncannot read int from %s on line %d from file %s",ins,line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
	}
	/*if(((temp/(MAXALLID+1))>MAXALLID)||((temp % (MAXALLID+1))>MAXALLID)) {
		sprintf(smess,"\nallele ID from genotype %d is higher than MAXALLID on line %d from file %s",temp,line,inputfile);
		write(ERRORFILE,smess);
		exit(1);
	}*/
	*outi = temp;
}	/*end of readintfromstring*/	

/***********************************************************************/	

void readintfromstring_bis(char *ins, int *outi, char *inputfile, int line)
/*same but does not test value of temp*/
{
	int flag,temp;
	char smess[SMAX];
	
	flag = sscanf(ins,"%d",&temp);
	if (flag == EOF) {
		sprintf(smess,"\ncannot read int from %s on line %d from file %s",ins,line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
	}
	*outi = temp;
}	/*end of readintfromstring_bis*/	

/***********************************************************************/

void readfloatfromstring(char *ins, float *outf, char *inputfile, int line)
{
	int flag;
	float temp;
	char smess[SMAX];
	
	flag = sscanf(ins,"%f",&temp);
	if (flag == EOF) {
		sprintf(smess,"\ncannot read float from %s on line %d from file %s",ins,line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
	}
	*outf = temp;
}	/*end of readfloatfromstring*/	

/***********************************************************************/

void readdoublefromstring(char *ins, double *outf, char *inputfile, int line)
{
	int flag;
	double temp;
	char smess[SMAX];
	
	flag = sscanf(ins,"%le",&temp);
	if (flag == EOF) {
		sprintf(smess,"\ncannot read float from %s on line %d from file %s",ins,line,inputfile);
		write(ERRORFILE,smess);
		printf("\npress any key to stop the program now");
		wait_a_char();
		exit(1);
	}
	*outf = temp;
}	/*end of readfloatfromstring*/	




/*****************************************************************************/
//  NON ANSI CODE

// Function waiting for a keystroke from the console and returning the character
// For Windows operating systems

char wait_a_char()
{
#ifdef NO_MESSAGE_PAUSE
	return 'x';
#endif
#ifdef USE_GETCH
	return (char)getch();
#else
	char ch = 'x';
	scanf("%c", &ch);
	return ch;
#endif
}

