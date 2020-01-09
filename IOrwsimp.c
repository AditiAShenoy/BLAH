/************************************************************************
This SW is provided under a NASA Open Source Software Usage Agreement for case GSC-15992-1.

Copyright © 2009 United States Government as represented by the Administrator of the
National Aeronautics and Space Administration. No copyright is claimed in the United States
under Title 17, U.S. Code. All Other Rights Reserved.

July 2009
***************************************************************************/
/*	This IOrws.c is a simplied version of IOroutine.c. It reads/writes in either byte or unsigned int2 data only.
    It does not pad data or lines, no bias adjustment. Input data will be put into int4 array.
    But it does offer to do byteswap

    7/17/09
*/


/* rpic.c reads nbit data into int buffer, does byteswap if set to do so
	7/17/09, psy
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void rpic(int nbit,int *inbuf,int npix, int nline, FILE *fp_in, int byteswap)
{

	unsigned char *charbuf;
	register i,iline;
	int *pt_inbuf;
	unsigned short int *inbuf_line;

    if(nbit <=8) charbuf=(unsigned char *) malloc(npix);
	 else inbuf_line = (unsigned short int *) malloc(npix*sizeof(short int));


	pt_inbuf=inbuf;
	for(iline=0;iline<nline;++iline)
	{
		if(nbit <=8)
		{ fread(charbuf,1,npix,fp_in);
			for(i=0;i<npix;++i) *(pt_inbuf+i) = (int) *(charbuf+i);
			pt_inbuf += npix;
		}
    	else  /* larger than 8 bits, unsigned  */
		{	fread(inbuf_line,sizeof(short int),npix,fp_in);
			if(nbit < 16)
				for(i=0;i<npix;++i) *(pt_inbuf+i) = (int) *(inbuf_line+i);
			else  /* 16-bit data unsigned */
			{
				if(byteswap)
				{ for(i=0;i<npix;++i) *(pt_inbuf+i) = (int) ( (( *(inbuf_line+i)&0xff00)>>8)
													| (( *(inbuf_line+i) &0x00ff)<<8) );

				}
				else /* no swapping */
				/*  for(i=0;i<npix;++i) *(pt_inbuf+i) = (int) *(inbuf_line+i);*/
				for(i=0;i<npix;++i)
				{
					if( *(inbuf_line+i) <0) printf("input data < 0 \n");
					*(pt_inbuf+i) = (int) *(inbuf_line+i);
					if ( *(pt_inbuf+i) <0 ) printf("i4 input data <0  \n");
				}
			}
			pt_inbuf += npix;
		}
	}  /* end of nline */

}

/*  wpic.c writes int buffer data into file, either as char or as unsigned int2
	depending on nbit, also clips the data at max and min value, and byteswapping if selected
	7/17/09, psy
*/

void wpic(int nbit, int *inbuf, int npix, int nline,FILE *fp_out, int byteswap)
{
	unsigned char *charbuf,*pt_char;
	register i,iline;
	int *pt_inbuf,max_val,min_val;
	unsigned short int *pt_i2,*int2buf,temp;


	int2buf = (unsigned short int *) malloc(npix*sizeof(unsigned short int));

	max_val=(int) pow(2.0,nbit)-1;
	min_val=0;

	if(nbit <=8) charbuf=(unsigned char *) malloc(npix);

	pt_inbuf=inbuf;
	if(nbit <= 8)
	{
		pt_char=charbuf;
		pt_inbuf=inbuf;
		for(iline=0;iline<nline;++iline)
		{
		    for(i=0;i<npix;++i)
		    {	if( *pt_inbuf > max_val) *(pt_char+i) = (unsigned char)	max_val;
				else if ( *pt_inbuf < min_val)
					*(pt_char+i) = (unsigned char) min_val;
				else	*(pt_char+i) = (unsigned char) *pt_inbuf;
				pt_inbuf++ ;
		    }
		    fwrite(charbuf,sizeof(unsigned char),npix,fp_out);
		}
	}
	if(nbit > 8)
	{
		pt_inbuf=inbuf;
		for(iline=0;iline<nline;++iline)
		{
		   for(i=0;i<npix;++i)
		    {	if( *(pt_inbuf+i) >max_val) *(pt_inbuf+i)=max_val;
			if( *(pt_inbuf+i) <min_val) *(pt_inbuf+i)=min_val;

		    }

		    for(i=0;i<npix;i++) *(int2buf+i) = (unsigned short int) (*(pt_inbuf+i)
		    					 & 0x00ffff);
		    if(byteswap)
				{ for(i=0;i<npix;++i)
					{	temp = *(int2buf+i);
						*(int2buf+i) =  ((temp>>8) & 0x00ff) | ((temp<<8)&0xff00);

					}
				}


		    pt_i2 = int2buf;
		    fwrite(pt_i2,sizeof(short int),npix,fp_out);
		    pt_inbuf += npix;
		}
	}
}


