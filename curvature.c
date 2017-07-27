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



void curvature(int dirx, int diry, int dirz, int stepx, int stepy, real bin_size_x, real bin_size_y, int binx, int biny, int fr_num,
		real *grid, real *gausCurve, real *meanCurve, int up_down)
{
	int i=0, j=0;
	int grid_size = binx*biny;
	int aux_ind_current = 0;
	int aux_ind_prev = 0;
	int aux_ind_next = 0;
	real denom_gaus = 0.0;
	real denom_mean = 0.0;
	real frame_num;
	frame_num = (real) fr_num;

	//derivatives
	real *gradX = NULL;
	rvec *gradXvec = NULL;
	real *gradY = NULL;
	rvec *gradYvec = NULL;
	snew(gradX,grid_size);
	snew(gradY,grid_size);
	snew(gradXvec,grid_size);
	snew(gradYvec,grid_size);
	real gradXX=0.0, gradYY=0.0, gradXY=0.0;

	//first parameters
	real paramE=0.0, paramF=0.0, paramG=0.0, paramp=0.0;
	rvec paramm = {0.0, 0.0, 0.0};
	rvec paramn = {0.0, 0.0, 0.0};

	//second parameters
	real paramL=0.0, paramM=0.0, paramN=0.0;

	//first order derivatives
	for(j=biny-1; j>=0; j--)
	{
		for(i=0; i<binx; i++)
		{
			  aux_ind_current = get_ind(i,j,binx);

			  //gradient in Y
			  if(j==0) //first row
			  {
				  aux_ind_prev = get_ind(i,biny-1,binx);
				  aux_ind_next = get_ind(i,j+1,binx);
				  gradY[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
				  gradYvec[aux_ind_current][dirx] = 0.0;
				  gradYvec[aux_ind_current][diry] = bin_size_y;
				  gradYvec[aux_ind_current][dirz] = gradY[aux_ind_current];
			  }
			  else if(j==biny-1) //last row
					{
						aux_ind_prev = get_ind(i,j-1,binx);
						aux_ind_next = get_ind(i,0,binx);
						gradY[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
						gradYvec[aux_ind_current][dirx] = 0.0;
						gradYvec[aux_ind_current][diry] = bin_size_y;
						gradYvec[aux_ind_current][dirz] = gradY[aux_ind_current];
					}
			  else //middle rows
			  {
					if(j+1>biny-1)
					{
						aux_ind_next = get_ind(i,0,binx);
					}
					else
					{
						aux_ind_next = get_ind(i,j+1,binx);
					}

					if(j-1<0)
					{
						aux_ind_prev = get_ind(i,biny-1,binx);
					}
					else
					{
						aux_ind_prev = get_ind(i,j-1,binx);
					}
					gradY[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
					gradYvec[aux_ind_current][dirx] = 0.0;
					gradYvec[aux_ind_current][diry] = bin_size_y;
					gradYvec[aux_ind_current][dirz] = gradY[aux_ind_current];
			  }

			  //gradient in X
			  if(i==0) //first column
			  {
				  aux_ind_prev = get_ind(binx-1,j,binx);
				  aux_ind_next = get_ind(i+1,j,binx);
				  gradX[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
				  gradXvec[aux_ind_current][dirx] = bin_size_x;
				  gradXvec[aux_ind_current][diry] = 0.0;
				  gradXvec[aux_ind_current][dirz] = gradX[aux_ind_current];
			  }
			  else if(i==binx-1) //last column
					{
						aux_ind_prev = get_ind(i-1,j,binx);
						aux_ind_next = get_ind(0,j,binx);
						gradX[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
						gradXvec[aux_ind_current][dirx] = bin_size_x;
						gradXvec[aux_ind_current][diry] = 0.0;
						gradXvec[aux_ind_current][dirz] = gradX[aux_ind_current];
					}
			  else //middle columns
			  {
					if(i+1>binx-1)
					{
						aux_ind_next = get_ind(0,j,binx);
					}
					else
					{
						aux_ind_next = get_ind(i+1,j,binx);
					}

					if(i-1<0)
					{
						aux_ind_prev = get_ind(binx-1,j,binx);
					}
					else
					{
						aux_ind_prev = get_ind(i-1,j,binx);
					}
					gradX[aux_ind_current] = up_down*(grid[aux_ind_next]/frame_num-grid[aux_ind_prev]/frame_num)/(2.0);
					gradXvec[aux_ind_current][dirx] = bin_size_x;
					gradXvec[aux_ind_current][diry] = 0.0;
					gradXvec[aux_ind_current][dirz] = gradX[aux_ind_current];
			  }
		}
	}

	//second order derivatives
	for(j=biny-1; j>=0; j--)
	{
		for(i=0; i<binx; i++)
		{
			  aux_ind_current = get_ind(i,j,binx);

			  //gradient in YY and XY
			  if(j==0) //first row
			  {
				  aux_ind_prev = get_ind(i,biny-1,binx);
				  aux_ind_next = get_ind(i,j+1,binx);
				  gradYY = (gradY[aux_ind_next]-gradY[aux_ind_prev])/(2.0);
				  gradXY = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
			  }
			  else if(j==biny-1) //last row
					{
						aux_ind_prev = get_ind(i,j-1,binx);
						aux_ind_next = get_ind(i,0,binx);
						gradYY = (gradY[aux_ind_next]-gradY[aux_ind_prev])/(2.0);
						gradXY = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
					}
			  else //middle rows
			  {
					if(j+1>biny-1)
					{
						aux_ind_next = get_ind(i,j+1-(biny-1)-1,binx);
					}
					else
					{
						aux_ind_next = get_ind(i,j+1,binx);
					}

					if(j-1<0)
					{
						aux_ind_prev = get_ind(i,biny-1+(j-1+1),binx);
					}
					else
					{
						aux_ind_prev = get_ind(i,j-1,binx);
					}
					gradYY = (gradY[aux_ind_next]-gradY[aux_ind_prev])/(2.0);
					gradXY = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
			  }


			  //gradients in XX
			  if(i==0) //first column
			  {
				  aux_ind_prev = get_ind(binx-1,j,binx);
				  aux_ind_next = get_ind(i+1,j,binx);
				  gradXX = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
			  }
			  else if(i==binx-1)
					{
						aux_ind_prev = get_ind(i-1,j,binx);
						aux_ind_next = get_ind(0,j,binx);
						gradXX = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
					}
			  else
			  {
					if(i+1>binx-1)
					{
						aux_ind_next = get_ind(i+1-(binx-1)-1,j,binx);
					}
					else
					{
						aux_ind_next = get_ind(i+1,j,binx);
					}

					if(i-1<0)
					{
						aux_ind_prev = get_ind(binx-1+(i-1+1),j,binx);
					}
					else
					{
						aux_ind_prev = get_ind(i-1,j,binx);
					}
					gradXX = (gradX[aux_ind_next]-gradX[aux_ind_prev])/(2.0);
			  }


			  //First fundamental Coefficients of the surface E, F, G, n, p
			  //paramE = sq_bin_size_x + gradX[aux_ind_current]*gradX[aux_ind_current];
			  paramE = gradXvec[aux_ind_current][dirx]*gradXvec[aux_ind_current][dirx] 
			          +gradXvec[aux_ind_current][diry]*gradXvec[aux_ind_current][diry]
			          +gradXvec[aux_ind_current][dirz]*gradXvec[aux_ind_current][dirz];

//			  paramF = gradX[aux_ind_current]*gradY[aux_ind_current];
			  paramF = gradXvec[aux_ind_current][dirx]*gradYvec[aux_ind_current][dirx]
			          +gradXvec[aux_ind_current][diry]*gradYvec[aux_ind_current][diry]
			          +gradXvec[aux_ind_current][dirz]*gradYvec[aux_ind_current][dirz];

//			  paramG = sq_bin_size_y + gradY[aux_ind_current]*gradY[aux_ind_current];
			  paramG =  gradYvec[aux_ind_current][dirx]*gradYvec[aux_ind_current][dirx]	//x
			          +gradYvec[aux_ind_current][diry]*gradYvec[aux_ind_current][diry]
			          +gradYvec[aux_ind_current][dirz]*gradYvec[aux_ind_current][dirz];


//			  cprod(gradXvec[aux_ind_current],gradYvec[aux_ind_current],paramm);
			  paramm[dirx] = gradXvec[aux_ind_current][diry]*gradYvec[aux_ind_current][dirz]-gradXvec[aux_ind_current][dirz]*gradYvec[aux_ind_current][diry];
			  paramm[diry] = gradXvec[aux_ind_current][dirz]*gradYvec[aux_ind_current][dirx]-gradXvec[aux_ind_current][dirx]*gradYvec[aux_ind_current][dirz];
			  paramm[dirz] = gradXvec[aux_ind_current][dirx]*gradYvec[aux_ind_current][diry]-gradXvec[aux_ind_current][diry]*gradYvec[aux_ind_current][dirx];

			  paramp = sqrt(iprod(paramm,paramm));
			  paramn[dirx] = paramm[dirx]/paramp;
			  paramn[diry] = paramm[diry]/paramp;
			  paramn[dirz] = paramm[dirz]/paramp;


			  //Second fundamental Coeffecients of the surface (L,M,N)
			  paramL = gradXX*paramn[dirz];
			  paramM = gradXY*paramn[dirz]; //x
			  paramN = gradYY*paramn[dirz];


			  //Gaussian Curvature
			  denom_gaus = paramE*paramG-paramF*paramF;
			  if(denom_gaus!=0.0)
			  {
				  gausCurve[aux_ind_current] = (paramL*paramN-paramM*paramM)/denom_gaus;
			  }
			  else
			  {
				  gausCurve[aux_ind_current] = 0.0;
			  }

			  //Mean Curvature
			  denom_mean = 2.0*(paramE*paramG - paramF*paramF);
			  if(denom_mean!=0.0)
			  {
				  meanCurve[aux_ind_current]	= (paramE*paramN + paramG*paramL - 2.0*paramF*paramM)/denom_mean;
			  }
			  else
			  {
				  meanCurve[aux_ind_current] = 0.0;
			  }


				//Principal Curvatures
		}
	}

	sfree(gradX);
	sfree(gradY);
	sfree(gradXvec);
	sfree(gradYvec);
}

