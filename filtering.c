#include <string.h>
#include <math.h>
#include <omp.h>
#include "main.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "txtdump.h"
#include "gmx_fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "assert.h"
#include "smalloc.h"
#include "names.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"
#include "vec.h"
#include "statutil.h"
//#include "gp_memory.h"
#include "index.h"
#include "confio.h"
#include "pbc.h"
//#include <lap.h>
#include "pdbio.h"
#include "gstat.h"
#include <fftw3.h>


void filter_curve_abs(real *grid, real *out, int binx, int biny, real bin_sizex, real bin_sizey,
		real filt_low, real filt_high, int fr_num, gmx_bool filter_verbose)
{
        //////// estimate and report number of modes in both directions /////////////
        int mq_high = floor(filt_high*(binx*bin_sizex)/(2.0*M_PI));
        int mq_low = floor(filt_low*(binx*bin_sizex)/(2.0*M_PI));
        int nq_high = floor(filt_high*(biny*bin_sizey)/(2.0*M_PI));
        int nq_low = floor(filt_low*(biny*bin_sizey)/(2.0*M_PI));

	if(filt_high < 0)
	{
		filt_high = 9999.99;
		mq_high = 9999.99;
		nq_high = 99999.99;
	}

	if(filter_verbose)
	{
		printf("\n******* Number of modes selected for filtering ******\n");
		printf("mq_low: %d \n",mq_low);
		printf("mq_high: %d \n",mq_high);
		printf("nq_low: %d \n",nq_low);
		printf("nq_high: %d \n",nq_high);
		printf("***************************\n\n");
	}
        //////// estimate and report number of modes in both directions /////////////

	filt_low = pow(filt_low,2); //because will be calculated squared radius
	filt_high = pow(filt_high,2);

	int size = binx*biny;
	real rsize = (real)size;
	int shiftX, shiftY;
	int sx=ceil(binx/2.0), sy=ceil(biny/2.0);
	real ssx=(real)2.0*M_PI/(binx*bin_sizex), ssy=(real)2.0*M_PI/(biny*bin_sizey); //these are for scaling
//printf("interesting: %d %f %f\n",sx,bin_sizex,ssx);
	real rad, frame_num;
	frame_num = (real) fr_num;

	int i,j,aux_ind,shift_aux_ind;
	fftw_complex *complex_grid, *complexfoo, *outfoo;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	//it is easier to convert real grid to a complex grid
	complexfoo = fftw_malloc(sizeof(fftw_complex) * size);
	outfoo = fftw_malloc(sizeof(fftw_complex) * size);
	for(i=0; i<size; i++)
	{
		complexfoo[i][0] = grid[i]/frame_num; //real
	    complexfoo[i][1] = 0.0f; //imaginary
	}

	//run FFTW
	complex_grid = fftw_malloc(sizeof(fftw_complex) * size);
	plan_forward = fftw_plan_dft_2d(biny, binx, complexfoo, complex_grid, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

/*    FILE *mat,*fft,*shift,*ifft;
    mat = fopen("mat.dat","w");
    fft = fopen("fft.dat","w");
    shift = fopen("shift.dat","w");
    ifft = fopen("ifft.dat","w");*/

	//filtering
	for(j=biny-1; j>=0; j--) //go over rows
	{
		for(i=0; i<binx; i++) //go over columns
		{
			  aux_ind = get_ind(i,j,binx);

			  //shift
			  if(i<sx && j<sy)
			  { shiftX=i+sx; shiftY=j+sy; }
			  else if(i>=sx && j<sy)
			  { shiftX=i-sx; shiftY=j+sy; }
			  else if(i<sx && j>=sy)
			  { shiftX=i+sx; shiftY=j-sy; }
			  else if(i>=sx && j>=sy)
			  { shiftX=i-sx; shiftY=j-sy; }
			  shift_aux_ind = get_ind(shiftX,shiftY,binx);

			  //calculate R from the mid-point
			  rad = ( pow((shiftX-sx)*ssx,2) + pow((shiftY-sy)*ssy,2) );

			  //filter

			  //*** brick-wall ****//
			  if(rad<filt_low || rad>filt_high)
			  {
				  complex_grid[aux_ind][0] = 0.0;
				  complex_grid[aux_ind][1] = 0.0;
			  }

			  //*** Butterworth ***//
			  /*if(rad<filt_low || rad>filt_high)
			  {
				  complex_grid[aux_ind][0] = 0.0;
				  complex_grid[aux_ind][1] = 0.0;
			  }
			  else
			  {
				  complex_grid[aux_ind][0] = 1/(pow(1+rad/filt_high,2));
				  complex_grid[aux_ind][1] = 1/(pow(1+rad/filt_high,2));
			  }*/

//			  fprintf(mat,"%f ",grid[aux_ind]/frame_num);
//			  fprintf(fft,"%f ",complex_grid[aux_ind][0]);
		}
//		fprintf(mat,"\n");
//		fprintf(fft,"\n");
	}
//	fclose(mat);
//	fclose(fft);

	//inverse FFTW
	plan_backward = fftw_plan_dft_2d(biny, binx, complex_grid, outfoo, FFTW_BACKWARD, FFTW_ESTIMATE);
	//plan_backward = fftw_plan_dft_c2r_2d ( binx, biny, complex_grid, dout, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );

	//save real part
	for(i=0; i<size; i++)
	{
		out[i] = outfoo[i][0]/rsize; //real
	}

	//output
	for(j=biny-1; j>=0; j--) //go over rows
	{
		for(i=0; i<binx; i++) //go over columns
		{
			  //shift
			  if(i<sx && j<sy)
			  { shiftX=i+sx; shiftY=j+sy; }
			  else if(i>=sx && j<sy)
			  { shiftX=i-sx; shiftY=j+sy; }
			  else if(i<sx && j>=sy)
			  { shiftX=i+sx; shiftY=j-sy; }
			  else if(i>=sx && j>=sy)
			  { shiftX=i-sx; shiftY=j-sy; }
			  shift_aux_ind = get_ind(shiftX,shiftY,binx);

			  aux_ind = get_ind(i,j,binx);
//                          fprintf(shift,"%f ",complex_grid[shift_aux_ind][0]);
//                          fprintf(ifft,"%f ",out[aux_ind]);
		}
//		fprintf(shift,"\n");
//		fprintf(ifft,"\n");
	}
//	fclose(shift);
//	fclose(ifft);

	fftw_destroy_plan( plan_forward );
	fftw_destroy_plan( plan_backward );
	sfree(complex_grid);
	sfree(complexfoo);
	sfree(outfoo);
}


void filter_curve_rel(real *grid, real *out, int binx, int biny, real bin_sizex, real bin_sizey, real filt_low, real filt_high, int fr_num, gmx_bool filter_verbose)
{
        //////// estimate and report filter radius in reciprocal space /////////////
        real q_x_high = M_PI*filt_high/bin_sizex;
        real q_x_low = M_PI*filt_low/bin_sizex;
        real q_y_high = M_PI*filt_high/bin_sizey;
        real q_y_low = M_PI*filt_low/bin_sizey;

	if(filt_high < 0)
	{
		filt_high = 9999.99;
		q_x_high = 9999.99;
		q_y_high = 99999.99;
	}

	if(filter_verbose)
	{
		printf("\n******* Relative number of modes was selected for filtering ******\n");
		printf("********* the relative numbers correspond to the following radii in the reciprocal space *****\n");
		printf("q_x_low: %f nm^-1\n",q_x_low);
		printf("q_x_high: %f nm^-1\n",q_x_high);
		printf("q_y_low: %f nm^-1\n",q_y_low);
		printf("q_y_high: %f nm^-1\n",q_y_high);
		printf("***************************\n\n");
	}
        //////// estimate and report filter radius in reciprocal space /////////////

	filt_low = pow(filt_low/2.0,2); //because will be calculated squared radius
	filt_high = pow(filt_high/2.0,2);

	int size = binx*biny;
	real rsize = (real)size;
	int shiftX, shiftY;
	int sx=ceil(binx/2.0), sy=ceil(biny/2.0);
	real ssx=(real)1.0/binx, ssy=(real)1.0/biny;
	real rad, frame_num;
	frame_num = (real) fr_num;

	int i,j,aux_ind,shift_aux_ind;
	fftw_complex *complex_grid, *complexfoo, *outfoo;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	//it is easier to convert real grid to a complex grid
	complexfoo = fftw_malloc(sizeof(fftw_complex) * size);
	outfoo = fftw_malloc(sizeof(fftw_complex) * size);
	for(i=0; i<size; i++)
	{
		complexfoo[i][0] = grid[i]/frame_num; //real
	    complexfoo[i][1] = 0.0f; //imaginary
	}

	//run FFTW
	complex_grid = fftw_malloc(sizeof(fftw_complex) * size);
	plan_forward = fftw_plan_dft_2d(biny, binx, complexfoo, complex_grid, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

/*    FILE *mat,*fft,*shift,*ifft;
    mat = fopen("mat.dat","w");
    fft = fopen("fft.dat","w");
    shift = fopen("shift.dat","w");
    ifft = fopen("ifft.dat","w");*/

	//filtering
	for(j=biny-1; j>=0; j--) //go over rows
	{
		for(i=0; i<binx; i++) //go over columns
		{
			  aux_ind = get_ind(i,j,binx);

			  //shift
			  if(i<sx && j<sy)
			  { shiftX=i+sx; shiftY=j+sy; }
			  else if(i>=sx && j<sy)
			  { shiftX=i-sx; shiftY=j+sy; }
			  else if(i<sx && j>=sy)
			  { shiftX=i+sx; shiftY=j-sy; }
			  else if(i>=sx && j>=sy)
			  { shiftX=i-sx; shiftY=j-sy; }
			  shift_aux_ind = get_ind(shiftX,shiftY,binx);

			  //calculate R from the mid-point
			  rad = ( pow((shiftX-sx)*ssx,2) + pow((shiftY-sy)*ssy,2) );

			  //filter

			  //*** brick-wall ****//
			  if(rad<filt_low || rad>filt_high)
			  {
				  complex_grid[aux_ind][0] = 0.0;
				  complex_grid[aux_ind][1] = 0.0;
			  }

			  //*** Butterworth ***//
			  /*if(rad<filt_low || rad>filt_high)
			  {
				  complex_grid[aux_ind][0] = 0.0;
				  complex_grid[aux_ind][1] = 0.0;
			  }
			  else
			  {
				  complex_grid[aux_ind][0] = 1/(pow(1+rad/filt_high,2));
				  complex_grid[aux_ind][1] = 1/(pow(1+rad/filt_high,2));
			  }*/

//			  fprintf(mat,"%f ",grid[aux_ind]/frame_num);
//			  fprintf(fft,"%f ",complex_grid[aux_ind][0]);
		}
//		fprintf(mat,"\n");
//		fprintf(fft,"\n");
	}
//	fclose(mat);
//	fclose(fft);

	//inverse FFTW
	plan_backward = fftw_plan_dft_2d(biny, binx, complex_grid, outfoo, FFTW_BACKWARD, FFTW_ESTIMATE);
	//plan_backward = fftw_plan_dft_c2r_2d ( binx, biny, complex_grid, dout, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );

	//save real part
	for(i=0; i<size; i++)
	{
		out[i] = outfoo[i][0]/rsize; //real
	}

	//output
	for(j=biny-1; j>=0; j--) //go over rows
	{
		for(i=0; i<binx; i++) //go over columns
		{
			  //shift
			  if(i<sx && j<sy)
			  { shiftX=i+sx; shiftY=j+sy; }
			  else if(i>=sx && j<sy)
			  { shiftX=i-sx; shiftY=j+sy; }
			  else if(i<sx && j>=sy)
			  { shiftX=i+sx; shiftY=j-sy; }
			  else if(i>=sx && j>=sy)
			  { shiftX=i-sx; shiftY=j-sy; }
			  shift_aux_ind = get_ind(shiftX,shiftY,binx);

			  aux_ind = get_ind(i,j,binx);
//                          fprintf(shift,"%f ",complex_grid[shift_aux_ind][0]);
//                          fprintf(ifft,"%f ",out[aux_ind]);
		}
//		fprintf(shift,"\n");
//		fprintf(ifft,"\n");
	}
//	fclose(shift);
//	fclose(ifft);

	fftw_destroy_plan( plan_forward );
	fftw_destroy_plan( plan_backward );
	sfree(complex_grid);
	sfree(complexfoo);
	sfree(outfoo);
}


