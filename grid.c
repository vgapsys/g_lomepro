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
#include "headers.h"


void fill_grid(int *is_cell_prot1, int *is_cell_prot2,
		int dirx, int diry, int dirz, rvec a1, rvec a2, int binx, int nliptop, int nlipbot,
		t_pbc pbc, int i, int j, real bin_sizex, real bin_sizey, real min_dist, int *aux_ind,
		int lip_ind, int k, int l, int *top_ind, rvec *lipidCOM,
		rvec dx, real dist, real *grid_up, real *grid_down, int *top_index,
		gmx_bool is_prot, int prot_ind, int *bot_ind, rvec *framex, real *height1, real *height2,
		atom_id *idlip, int nprot_top, int *ptop_ind, int *pbot_ind, int *bottom_index,
		int nprot_bot, real left_x, real left_y)
{
	  *is_cell_prot1=0;
	  *is_cell_prot2=0;

	  a1[dirx] = left_x+i*bin_sizex;
	  a1[diry] = left_y+j*bin_sizey;

	  //top leaflet
	  min_dist=9999.9;
	  *aux_ind = get_ind(i,j,binx);
	  lip_ind = 0;
	  for(k=0; k<nliptop; k++)
	  {
		  lip_ind = top_ind[k];

		  a2[dirx] = lipidCOM[lip_ind][dirx];
		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  grid_up[*aux_ind] = lipidCOM[lip_ind][dirz];
			  *height1 = grid_up[*aux_ind];
			  min_dist = dist;
			  *top_index=idlip[lip_ind];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(k=0; k<nprot_top; k++)
		  {
			  prot_ind = ptop_ind[k];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  real dist = norm2(dx);

			  if(dist < min_dist)
			  {
				  grid_up[*aux_ind] = framex[prot_ind][dirz];
				  min_dist = dist;
				  *is_cell_prot1=1;
				  *top_index=-1;
			  }
		  }
	  }

	  //bottom leaflet
	  min_dist=9999.9;
	  lip_ind = 0;
	  for(l=0; l<nlipbot; l++)
	  {
		  lip_ind = bot_ind[l];

		  a2[dirx] = lipidCOM[lip_ind][dirx];
		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  grid_down[*aux_ind] = lipidCOM[lip_ind][dirz];
			  *height2 = grid_down[*aux_ind];
			  min_dist = dist;
			  *bottom_index=idlip[lip_ind];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(l=0; l<nprot_bot; l++)
		  {
			  prot_ind = pbot_ind[l];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  dist = norm2(dx);

			  if(dist < min_dist)
			  {
				  grid_down[*aux_ind] = framex[prot_ind][dirz];
				  min_dist = dist;
				  *is_cell_prot2=1;
				  *bottom_index=-1;
			  }
		  }
	  }
}

//special grid for order parameters
void fill_grid_order(int dirx, int diry, int dirz,
		int binx, int nliptop, int nlipbot,
		t_pbc pbc, int i, int j, real bin_sizex, real bin_sizey, int aux_ind,
		int lip_ind, int k, int l, int *top_ind,
		real ***grid_up, real ***grid_down,
		gmx_bool is_prot, int *bot_ind, rvec *framex,
		atom_id *idlip, int *nprot_top_order, int **ptop_ind_order, int **pbot_ind_order,
		int *nprot_bot_order, atom_id *idorder,
		int order_atom_num, int acyl_atom,
		real **order_lip,rvec *lipidCOM,real z_mid,real order_val,
		real *grid_up_order, real *grid_down_order, real scale, real left_x, real left_y,
		int **order_count_up, int **order_count_down)
{
	  int acyl_atom_ind=0;
	  rvec a1 = {0.0, 0.0, 0.0};
	  rvec a2 = {0.0, 0.0, 0.0};
	  a1[dirz] = z_mid;
	  a2[dirz] = z_mid;

	  rvec dx = {0.0, 0.0, 0.0};
	  real dist = 0.0;
	  real best_avg = 0.0;
	  real best_sd = 0.0;

	  a1[dirx] = left_x+i*bin_sizex;
	  a1[diry] = left_y+j*bin_sizey;

	  //variables for z values
	  int is_cell_prot1=0, is_cell_prot2=0;
	  real foo_up=0.0, foo_down=0.0;

	  //top leaflet
	  real min_dist=9999.9;
	  int prot_ind = 0;
	  for(k=0; k<nliptop; k++)
	  {
		  acyl_atom_ind = idorder[acyl_atom+top_ind[k]*order_atom_num];

		  a2[dirx] = framex[acyl_atom_ind][dirx];
		  a2[diry] = framex[acyl_atom_ind][diry];
//		  a2[dirx] = lipidCOM[lip_ind][dirx];
//		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  best_avg = order_lip[top_ind[k]][acyl_atom-1];
			  best_sd = pow(order_lip[top_ind[k]][acyl_atom-1],2);
			  min_dist = dist;
			  foo_up = framex[acyl_atom_ind][dirz];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(k=0; k<nprot_top_order[acyl_atom-1]; k++)
		  {
			  prot_ind = ptop_ind_order[acyl_atom-1][k];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  dist = norm2(dx);

			  if(dist < min_dist)
			  {
//				  best_avg = order_val;
//				  best_sd = order_val*order_val;
				  min_dist = dist;
				  is_cell_prot1=1;
			  }
		  }
	  }

	  if(is_cell_prot1==0)
	  {
		  grid_up[acyl_atom-1][aux_ind][0] += best_avg;
		  grid_up[acyl_atom-1][aux_ind][1] += best_sd;
		  order_count_up[acyl_atom-1][aux_ind]++;
	  }


	  //bottom leaflet
	  min_dist=9999.9;
	  for(l=0; l<nlipbot; l++)
	  {
		  acyl_atom_ind = idorder[acyl_atom+bot_ind[l]*order_atom_num];

		  a2[dirx] = framex[acyl_atom_ind][dirx];
		  a2[diry] = framex[acyl_atom_ind][diry];
//		  a2[dirx] = lipidCOM[lip_ind][dirx];
//		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  best_avg = order_lip[bot_ind[l]][acyl_atom-1];
			  best_sd = pow(order_lip[bot_ind[l]][acyl_atom-1],2);
			  min_dist = dist;
			  foo_down = framex[acyl_atom_ind][dirz];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(l=0; l<nprot_bot_order[acyl_atom-1]; l++)
		  {
			  prot_ind = pbot_ind_order[acyl_atom-1][l];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  dist = norm2(dx);

			  if(dist < min_dist)
			  {
//				  best_avg = order_val;
//				  best_sd = order_val*order_val;
				  min_dist = dist;
				  is_cell_prot2=1;
			  }
		  }
	  }
	  if(is_cell_prot2==0)
	  {
		  grid_down[acyl_atom-1][aux_ind][0] += best_avg;
		  grid_down[acyl_atom-1][aux_ind][1] += best_sd;
		  order_count_down[acyl_atom-1][aux_ind]++;
	  }

	  //z values for the grid
	  if(is_cell_prot1==0 && is_cell_prot2==0)		//lipids on both layers
	  {
		  grid_up_order[aux_ind] += foo_up;
		  grid_down_order[aux_ind] += foo_down;
	  }
	  else
		  if(is_cell_prot1+is_cell_prot2==1)	//lipid in one layer, protein in another
		  {
			  grid_up_order[aux_ind] += z_mid+scale*(foo_up-z_mid);
			  grid_down_order[aux_ind] += z_mid+scale*(foo_down-z_mid);
		  }
	  else
		  if(is_cell_prot1==1 && is_cell_prot2==1)	//proteins on both layers
		  {
			  grid_up_order[aux_ind] += z_mid;
			  grid_down_order[aux_ind] += z_mid;
		  }
}


void fill_grid_diffus(int *is_cell_prot1, int *is_cell_prot2,
		int dirx, int diry, int dirz, rvec a1, rvec a2, int binx, int nliptop, int nlipbot,
		t_pbc pbc, int i, int j, real bin_sizex, real bin_sizey, real min_dist, int *aux_ind,
		int lip_ind, int k, int l, int *top_ind, rvec *lipidCOM,
		rvec dx, real dist, real *grid_up, real *grid_down, int *top_index,
		gmx_bool is_prot, int prot_ind, int *bot_ind, rvec *framex, real *height1, real *height2,
		atom_id *idlip, int nprot_top, int *ptop_ind, int *pbot_ind, int *bottom_index,
		int nprot_bot, real left_x, real left_y,
		int diffus_steps, rvec **diffus_grid_pos_up, rvec **diffus_grid_pos_down)
{
	  *is_cell_prot1=0;
	  *is_cell_prot2=0;

	  a1[dirx] = left_x+i*bin_sizex;
	  a1[diry] = left_y+j*bin_sizey;

/////////////////////////////////////////
	  //top leaflet
	  min_dist=9999.9;
	  *aux_ind = get_ind(i,j,binx);
	  lip_ind = 0;
	  for(k=0; k<nliptop; k++)
	  {
		  lip_ind = top_ind[k];

		  a2[dirx] = lipidCOM[lip_ind][dirx];
		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  grid_up[*aux_ind] = lipidCOM[lip_ind][dirz];
			  *height1 = grid_up[*aux_ind];
			  min_dist = dist;
			  *top_index=idlip[lip_ind];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(k=0; k<nprot_top; k++)
		  {
			  prot_ind = ptop_ind[k];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  real dist = norm2(dx);

			  if(dist < min_dist)
			  {
				  grid_up[*aux_ind] = framex[prot_ind][dirz];
				  min_dist = dist;
				  *is_cell_prot1=1;
				  *top_index=-1;
			  }
		  }
	  }

///// diffusion stuff ////////////
	  real foo=0.0;
	  for(k=diffus_steps-1; k<=0; k--)
	  {
		  if(k==0)
		  {
			  if(*top_index != -1)
			  {
				  //diffus_grid_pos_up[aux_ind][k][XX] = lipidCOM[lip_ind][dirx];
				  //diffus_grid_pos_up[aux_ind][k][YY] = lipidCOM[lip_ind][dirx];
				  //diffus_grid_pos_down[aux_ind][k][XX] = lipidCOM[lip_ind][dirx];
				  //diffus_grid_pos_down[aux_ind][k][YY] = lipidCOM[lip_ind][dirx];
			  }
			  else
			  {

			  }
		  }
		  else
		  {
			  //diffus_grid_pos_up[aux_ind][k] = diffus_grid_pos_up[aux_ind][k-1];
			  //diffus_grid_pos_down[aux_ind][k] = diffus_grid_pos_down[aux_ind][k-1];
		  }
	  }



////////////////////////////////////////////
	  //bottom leaflet
	  min_dist=9999.9;
	  lip_ind = 0;
	  for(l=0; l<nlipbot; l++)
	  {
		  lip_ind = bot_ind[l];

		  a2[dirx] = lipidCOM[lip_ind][dirx];
		  a2[diry] = lipidCOM[lip_ind][diry];

		  pbc_dx(&pbc,a1,a2,dx);
		  dist = norm2(dx);

		  if(dist < min_dist)
		  {
			  grid_down[*aux_ind] = lipidCOM[lip_ind][dirz];
			  *height2 = grid_down[*aux_ind];
			  min_dist = dist;
			  *bottom_index=idlip[lip_ind];
		  }
	  }
	  if(is_prot)
	  {
		  prot_ind = 0;
		  for(l=0; l<nprot_bot; l++)
		  {
			  prot_ind = pbot_ind[l];

			  a2[dirx] = framex[prot_ind][dirx];
			  a2[diry] = framex[prot_ind][diry];

			  pbc_dx(&pbc,a1,a2,dx);
			  dist = norm2(dx);

			  if(dist < min_dist)
			  {
				  grid_down[*aux_ind] = framex[prot_ind][dirz];
				  min_dist = dist;
				  *is_cell_prot2=1;
				  *bottom_index=-1;
			  }
		  }
	  }
}


