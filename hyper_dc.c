/*
This SW is provided under a NASA Open Source Software Usage Agreement for case GSC-15992-1.

Copyright © 2009 United States Government as represented by the Administrator of the
National Aeronautics and Space Administration. No copyright is claimed in the United States
under Title 17, U.S. Code. All Other Rights Reserved.


Hyper_dc.c, is based on a hyper-spectral preprocessor and uses the CCSDS121B as the entropy coder using szip code 
		in '-ec' (entropy coding) mode 
Note that the mapping of the prediction result needs work, current is using simple mapping, by:
		(ek >= 0? 2*ek : -(2*ek)-1);
		==> Note that this does NOT guarantee values to be within 16-bit.
This version uses BSQ format as input.

The hyperspectral prediction is based on:  
	paper 1: "Low-Complexity Lossless Compression of Hyperspectral Imagery Via Adaptive Filtering", M. Klimesh, Sept, 2005
	paper 2: "CCSDS Concept Paper:  Fast Lossless Compression of Multispectral Imagery", M. Klimesh, Sept, 2007
	paper 3: "CCSDS Concept Paper:  Fast Lossless Compressiom of Multispectral Imagery", M. Klimesh and A. Kiely, June, 2009

hyper_dc.c (main.c here) version 1.0 uses 3 bands and all the local neigbhorhood processing
   Usage: hyper_dc -w [ ] -h [ ] -nbit [ ] -nband [ ] [-s [ ] ]-fout Resultfilename filename
          -w: width of the 2d spatial
          -h: height of the 2d spatial
          -nbit: input dynamic range of 3d cube, limited to 16. Default to 8.

          filename: band-sequential (BSQ) input file name
          mapped predictor error saved as:  filename.perr
		  -fout: stats of results saved in Resultfilename
		  [-s [ ] ]: option to enter the "-s" paramter for szip to avoid excessive padding and short lines.

    Output unsigned i2 array of predictor error values as input file for szip.     
;    The Wk dimension is set at 6, but the valid number of elements are:
;    Band 0: [1/3, 1/3, 1/3, 0,0,0]
;    Band 1: [1/4, 1/4, 1/4, 1/4, 0, 0]
;    Band 2: [1/5, 1/5, 1/5, 1/5, 1/5, 0]
;    Others: [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
; 
; Boundary conditions such as first line of each band is taken care of:
;    First line in each band: by setting appropriate Si values to be 0, thus 
;    Band 0: first pixel is saved as the ek value; other pixels use previous pixel value to fill
;     the needed Si values
;    Other bands: appropriately fill in with previous values
;    Last pixel on each line: the S4 value is filled with S3 value.
; 
; in same size as input file
; Assumption: width: >3, height: >2, bands: >= 4
  Output max/min prediction error values/band to filename.log

	3/12/09, psy, use only float data type
	3/24/09, psy, change the way boundary conditions are used, output first pix as ref. then when
		x=0, or x=width-1, or when y=0, setting appropriate si to be zero,and s0p, s5p, s10p and s15p
		with averages on available values.

	10/7/09, use mu_start and mu_end (as in the June 2009 paper) and write out the compressed file size in a 
	given filename to to make it easier later for doing batch job by using -fout option
	3/8/10, switching to the publicly available szip (for CCSDS 121B) codes; also do the simple prediction error
	mapping according to what is in the Klimesh paper.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#define ROUND(a)	( (a >= 0) ? ( int) (a+0.5) : ( int) (a-0.5) )
#define MAX(a,b)	( (a >= b) ? a : b)
#define MIN(a,b)	( (a < b) ? a : b)
int main(int argc, char *argv[])
{   char filename[128],f_result[128],filename_perr[128],widthchar[64],nbitchar[64],filename_cmp[128];
	FILE *fp_in,*fp_out,*fp_log,*fp_result,*fp_cmp;
	int *buf0,*buf1,*buf2,*buf3,*pt0,*pt1,*pt2,*pt3,s[20],stmp0,stmp1,stmp2,stmp3,arraymax;
	int *error,*pterror;
	unsigned short int *errori2,*pterrori2;
	int width,height,nbit,nband,emax,emin,spredmax,spredmin;
	int emaxfile,eminfile,spredmaxfile,spredminfile,s0predictor;
	int ek=0;
	int f_out_flag=1;   /* normal case, do entropy coding and collect the compressed file size into Result file */

	int *s_pred,*pt_s_pred,sp;
	register i,iline,ipix,iband;
	int mu_schedule=2;
	int s_uses=0;			/* 10/11/09 added this to handle the case when width is too small */

	double wk[6]={0.0,0.0,0.0,0.0,0.0,0.0};
	double mu = 0.0;
/*	double c = 0.0; */
	double dk,dkp;
	char command1[128],command2[128];			/* 10/3/09 */
	struct stat stbuf; 
	long int int_f_size,filesize_total;
	double f_size=0.0;
	double mu_start,mu_end;		/* 10/07/09 */

	void get_ek(int *s, double *wk, double mu, int arraymax, int *ek, double *dk, double *dkp, 
				int *s0predictor);
  
	void rpic(int nbit,int *inbuf,int npix, int nline, FILE *fp_in, int byteswap);
	void wpic(int nbit, int *inbuf, int npix, int nline,FILE *fp_out, int byteswap);

  	nbit=8;	   /* default to 8bit to avoid random values */
	width=512;	/* default value  */
	height=512;	/* default value  */
    nband = 256; 
	

    stmp0=stmp1=stmp2=stmp3=0;
	emaxfile=eminfile=0; 	/* find out the max and min value of e */
	spredmaxfile=spredminfile=0; 
	s0predictor=0;
	
    
   	if(argc < 12)
	{   printf("hyper_dc version 1.0 \n");
		printf("hyper_dc  -w [width] -h [height] -nbit [ ] -nband [ ]  -fout Resultfile filename\n");
		  exit(1);
	}
  /* getting parameters */
	for(i=1;i<argc;i+=2) 
	{	if(strcmp(argv[i],"-w")==0) 
			width=atoi(argv[i+1]);
		
		if(strcmp(argv[i],"-h")==0) 
			height=atoi(argv[i+1]);	

		if(strcmp(argv[i],"-nbit")==0) 
			nbit=atoi(argv[i+1]);
			
		if(strcmp(argv[i],"-nband")==0) 
			nband=atoi(argv[i+1]);	

		if(strcmp(argv[i],"-s")==0)
			s_uses = atoi(argv[i+1]);

		if(strcmp(argv[i],"-fout")==0)
		{	f_out_flag=1;
			strcpy(f_result,argv[i+1]);
		}
	}
	arraymax = (int) pow(2,nbit)-1;
	mu_start = pow(2,(-nbit-2)); 

	mu_end = pow(2, (-nbit-6));
	mu_schedule = ROUND(1024./width);		/* change mu every mu_schedule line,slightly different from paper 3, 2009 */
	if (mu_schedule == 0) mu_schedule=1;
	printf("array dynamic range max for %d bit data is: %d \n",nbit,arraymax);
	strcpy(filename,argv[argc-1]);
	if((fp_in = fopen(filename,"rb")) == NULL) 
    {
		printf("Failed to open input file \n");
		exit(1);
	}
	

	strcpy(filename_perr,argv[argc-1]);			/* 3/8/10 */
	strcat(filename_perr,".perr");		/* 3/8/10 mapped prediction error */

	if( (fp_out = fopen(filename_perr,"wb")) == NULL) 
	{	
	 	printf("Failed to open output file \n");
		exit(1);
	}
	strcpy(filename,argv[argc-1]);
 	strcat(filename,".log"); 
	if( (fp_log = fopen(filename,"w")) == NULL) {	 
	 	printf("Failed to open log file \n");
		exit(1);
	} 
	if(f_out_flag)
	{  if( (fp_result = fopen(f_result,"a")) == NULL)
		{
		printf("Failed to open result file to append compressed file size: \n");
		exit(1);
		}
	}

	fprintf(fp_log,"%s width: %d , height: %d , bands: %d \n ",filename,width,height,nband);
/* memory allocation */
	buf0 = (int *) calloc(width*height,sizeof(int));
	buf1 = (int *) calloc(width*height,sizeof(int));
	buf2 = (int *) calloc(width*height,sizeof(int));
	buf3 = (int *) calloc(width*height,sizeof(int));
	error =(int *) calloc(width*height,sizeof(int));
	errori2=(unsigned short int *) calloc(width*height,sizeof(short int));		/* 3/8/10 */
	s_pred = (int *) calloc(width*height,sizeof(int));
			
   for(iband=0;iband<nband;++iband)
   {  emax=0;
	  emin=0;
	  spredmax=0;
	  spredmin=0;
      for(i=0;i<20;++i) s[i]=0;                             
      pterror = error;   
	  pterrori2 = errori2;  /* 3/8/10 */
	  pt_s_pred = s_pred;
      if (iband == 0) 
        for(i=0;i<3;++i) wk[i]= (float) 0.3333;
      else if (iband == 1) 
        for(i=0;i<4;++i) wk[i]=0.25;
      else if (iband == 2) 
        for(i=0;i<5;++i) wk[i]= (float) 0.2;
      else 
/*      for(i=0;i<6;++i) wk[i]= (float) 0.16667;  */
	  {     wk[0]=0.1;
			wk[1]=0.1;
			wk[2]=0.1;
			wk[3]=0.5;
			wk[4]=0.1;
			wk[5]=0.1;
	  }

      /* ring buffer assignment. always read into pt0 */
      switch (iband%4)
      { 
             case 0: 
                  pt0=buf0;
                  pt1=buf3;
                  pt2=buf2;
                  pt3=buf1;
                  break;
             case 1:
                  pt0=buf1;
                  pt1=buf0;
                  pt2=buf3;
                  pt3=buf2;
                  break;
             case 2:
                  pt0=buf2;
                  pt1=buf1;
                  pt2=buf0;
                  pt3=buf3;
                  break;
             case 3:
                  pt0=buf3;
                  pt1=buf2;
                  pt2=buf1;
                  pt3=buf0;
                  break;    
      }
      rpic(nbit,pt0,width,height,fp_in,0);
      
	  mu=mu_start;				/* 10/07 09 */
      for (iline=0;iline<height;iline++)
      {
/*         if (iline <= 9 ) mu =  c;
          else mu = c*alpha;   
*/
		 if ((mu > mu_end) && (iline%mu_schedule == 0) && (iline > 0))  mu /= 2;  /* 10/11/09, reduce mu every 2 lines */
          
          if(iline==0)
          {             
             for(ipix=0;ipix<width;ipix++)   
             { 
               /* s0 to s4 */                                
               s[0] = *pt0;
               if (ipix > 0)
               { 
                 s[1] = *(pt0-1);
                 s[2] = s[1];
                 s[3] = s[1];
                 s[4] = s[3];
               }
               pt0++;    /* current band */
              /* s5 to s9 */
              s[5] = *pt1;    /* same pixel from previous band */
              if (ipix > 0)
              {  s[6] = *(pt1-1);
                 s[7] = s[6];
                 s[8] = s[6];
                 s[9] = s[8]; 
              }
              pt1++;
              /* s10 to s14 */
              s[10] = *pt2;
              if (ipix > 0)
              {   s[11] = *(pt2-1);
                  s[12] = s[11];
                  s[13] = s[11];
                  s[14] = s[13];
              }
              pt2++;
              /* s15 to s19 */
              s[15] = *pt3;
              if (ipix > 0)
              {  s[16] = *(pt3-1);
                 s[17] = s[16];
                 s[18] = s[16];
                 s[19] = s[18];
              }
              pt3++;
              /* now process to get ek */
              if (ipix == 0)
              {  ek=s[0]-s[5];         /* only uses previous band for pix 0 */
                 dk = (float) s[0];
                 dkp = (float) s[5];
				 /* sp = s[5]; */
				 s0predictor = s[5];  /* 10/5/09 */
              }              
              else   /* ipix >0 for the iline=0 case */
              { 
		/*		  if(iband==33 && iline==0 && ipix <= 30)
                { 
                  printf("iline:%d, ipix: %d s values:",iline,ipix);
                  for(i=0;i<20;++i) printf("%d ",s[i]);
				  printf("before computing ek: \n");
                  printf("wk, mu: %f %f %f %f %f %f %f \n",wk[0],wk[1],wk[2],wk[3],wk[4],wk[5],mu);
                  
                } 
         */
				  get_ek(s,wk,mu,arraymax,&ek,&dk,&dkp,&s0predictor);
				  emax = MAX(emax,ek);
				  emin = MIN(emin,ek);
				  spredmax = MAX(spredmax,s0predictor);
				  spredmin = MIN(spredmin,s0predictor);
				  
         /*       if(iband==33 && iline==0 && ipix <= 30)
                { 
				  printf("ek: %d,  dk: %f,  dkp: %f  \n",ek,dk,dkp);
                  printf("after ek, wk: %f %f %f %f %f %f \n \n",wk[0],wk[1],wk[2],wk[3],wk[4],wk[5]);
                } 
		*/		
			  }
              *pterror++ = ek;    
			  *pt_s_pred++ = s0predictor;
			  *pterrori2++ = (unsigned short int) (ek >= 0? 2*ek : -(2*ek)-1); /* 3/8/10 */
            } /* end of ipix loop */  
          }  /* now pointer at first pixel of 2nd line */
          else      /* all other lines */
          {  
  
             for(ipix=0;ipix<width;ipix++)
             {  /* s0 to s4 */   
                stmp0 = *(pt0-width);  /* pixel above as s1*/
                s[0] = *pt0;
                if(ipix==0)
                { s[1]=stmp0;
                  s[2]=stmp0;
                  }
                else  /* other pixels on line */
                {
                   s[1] = *(pt0-1);
                   s[2] = *(pt0-1-width);
                   }
                s[3]= stmp0;  /* pixel above */
                if( ipix == width -1) s[4]=s[3];   /*last pixel on line */
                else s[4]= *(pt0 - width +1);              
                pt0++;     /* increment pointer */
                /* s5 to s9 */
                stmp1 = *(pt1-width);
                s[5] = *pt1;
                if(ipix == 0)
                {  s[6] = stmp1;
                   s[7] = s[6];
                }
                else  /*other pixels on line */
                {  s[6] = *(pt1-1);
                   s[7] = *(pt1-1-width);
                }
                s[8] = stmp1;
                if( ipix == width-1) s[9] = s[8];
                else s[9] = *(pt1-width+1);
                pt1++;
                /* s10 to s14 */
                stmp2 = *(pt2-width);
                s[10] = *pt2;
                if(ipix == 0)
                {  s[11]=stmp2;
                   s[12] = s[11];
                }
                else
                {
                    s[11] = *(pt2-1);
                    s[12] = *(pt2-1-width);
                }
                s[13] = stmp2;
                if (ipix == width -1) s[14] = s[13];
                else s[14] = *(pt2-width+1);
                pt2++;
                /* s15 to s19 */
                stmp3 = *(pt3-width);
                s[15] = *pt3;
                if(ipix ==0)
                {  s[16] = stmp3;
                   s[17] = s[16];
                }
                else
                {  s[16] = *(pt3-1);
                   s[17] = *(pt3-1-width);
                }
                s[18] = stmp3;
                if (ipix == width-1) s[19] = s[18];
                else s[19] = *(pt3-width+1);
                pt3++;  
                /* now process to get ek and ek_map;  */ 
       
                get_ek(s,wk,mu,arraymax,&ek,&dk,&dkp,&s0predictor);
				emax = MAX(emax,ek);
				emin = MIN(emin,ek);
				spredmax = MAX(spredmax,s0predictor);
				spredmin = MIN(spredmin,s0predictor);

    
                *pterror++ = ek;     
				*pterrori2++ = (short int) (ek >= 0? 2*ek : -(2*ek)-1);  /* 3/8/10 */
				*pt_s_pred++ =  s0predictor;
              }   /* end of ipix loop */  
           } /* end of else for lines */   
      } /* end of iline loop in current band */          
   /* output  for each plane to file, currently as 2byte data */
   /*  fwrite(s_pred,sizeof(unsigned short), width*height,fp_out); */
/*	wpic(nbit,s_pred,width,height,fp_out,0); */ /* comment out 3/8/10. No byteswap for szip */
	 
	 
   fwrite(errori2,sizeof(unsigned short int),width*height,fp_out);  /* 3/8/10 */
   fprintf(fp_log,"iband: %d , emax: %d ,  emin: %d  spredmax: %d  spredmin: %d \n",
		iband,emax,emin,spredmax,spredmin);
		emaxfile =MAX(emaxfile,emax);
		eminfile = MIN(eminfile,emin);
		spredmaxfile=MAX(spredmaxfile,spredmax);
		spredminfile=MIN(spredminfile,spredmin);
   }  /* end of iband loop */ 
   fclose(fp_in);
   fclose(fp_out);
   fprintf(fp_log,"file emax: %d ,  emin: %d \n",emaxfile,eminfile);
   fprintf(fp_log,"file s_pred_max: %d ,  s_pred_max: %d \n",spredmaxfile,spredminfile);
   if(f_out_flag)		/* do compression and write result file size to file f_result */
	{	
/* 10/11/09 added to avoid excessive padding for uses */
		if(s_uses >0 ) width =s_uses;			/* so that the entered s parameter is used for szip */
		sprintf(widthchar,"%d ",width);
		sprintf(nbitchar,"%d ",nbit);
		sprintf(command2,"C:\\psyeh\\Cprograms\\szip\\szip -ec -ki -s ");   /* keep file to be compressed */
		strcat(command2,widthchar);
		strcat(command2," -n ");
		strcat(command2,nbitchar);
		strcat(command2,"  ");
		strcat(command2,filename_perr);
		system(command2);
		strcpy(filename_cmp,filename_perr);
		strcat(filename_cmp,".sz");
		if(!stat(filename_cmp,&stbuf))
        {     f_size = (double)stbuf.st_size; 
                                                /* get the file size, now is float array file */
              int_f_size = (int) f_size;
        }
		fprintf(fp_result,"%s     %d  \n",filename_cmp,int_f_size);
   }  /* end of f_out_flag = 1 case */


   fclose(fp_log);  
   if(f_out_flag) 
	   fclose(fp_result);
	   
   return 0;
	
}
  


void get_ek(int *s, double *wk, double mu, int arraymax, int *ek, double *dk, double *dkp, int *s_pred)
{
     double s0p,s5p,s10p,s15p,uk[6];
	 int sgn_ek,sptmp;
     register i;
     
     s0p = (double) (s[1]+s[2]+s[3]+s[4])/4.0;
     s5p = (double) (s[6]+s[7]+s[8]+s[9])/4.0;
     s10p= (double) (s[11]+s[12]+s[13]+s[14])/4.0;
     s15p= (double) (s[16]+s[17]+s[18]+s[19])/4.0;
     uk[0]=s[1] -s0p;
     uk[1]=s[2] -s0p;
     uk[2]=s[3] -s0p;
     uk[3]=s[5]- s5p;
     uk[4]=s[10]-s10p;
     uk[5]=s[15]-s15p;
     
 /*    *dk = s[0] - s0p; 10/5/09 */
	 *dk = s0p - s[0];   /* 10/05/09 */
     *dkp = 0.0;
     for(i=0;i<6;++i)
     {  *dkp += uk[i]*wk[i];
     }
  /*   *ek = ROUND(*dk) - ROUND(*dkp);  */
	 if(s0p < 0) printf("local mean < 0: %d \n",ROUND(s0p));
/*	 sptmp = ROUND(s0p)+ROUND(*dkp); 10/5/09 */
	 sptmp = ROUND(s0p + *dkp);		/*10/5/09 */
	 if(sptmp <0) sptmp = 0;
	 if(sptmp > arraymax) sptmp =  arraymax;
	 *s_pred = sptmp;
/*	 *ek = s[0] - sptmp;   10/05/09 */
	 *ek = sptmp - s[0];   /* 10/05/09 */
     sgn_ek = ( (*ek) >= 0? 1 : -1);
     for(i=0;i<6;++i)
     { /*  wk[i] += (mu*uk[i]*sgn_ek); 10/05/09 */
		wk[i] -= (mu*uk[i]*sgn_ek);
     }
     
}    
