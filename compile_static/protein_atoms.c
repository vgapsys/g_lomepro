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

void protein_atoms(int nprot, real z_mid, rvec *framex, rvec *lipidCOM,
		int dirx, int diry, int dirz, atom_id *idprot, int nliptop, int nlipbot,
		t_pbc pbc, real pr2, atom_id *top_ind, atom_id *bot_ind,
		atom_id *ptop_ind, atom_id *pbot_ind, int *nprot_top, int *nprot_bot)
{
	int i=0, j=0, counter=0, switch1=0, switch2=0, save_time=0;
	rvec a1 = {0.0, 0.0, 0.0};
	rvec a2 = {0.0, 0.0, 0.0};
	a1[dirz] = z_mid;
	a2[dirz] = z_mid;

	rvec dx = {0.0, 0.0, 0.0};
	real dist = 0.0;
	for(i=0; i<nprot; i++)
	{
		//top leaflet
		j=0;
		counter=0;
		switch1=0;
		switch2=0;
		save_time=0;
		a1[dirx] = framex[idprot[i]][dirx];
		a1[diry] = framex[idprot[i]][diry];

		for(j=0; j<nliptop; j++)
		{
			a2[dirx] = lipidCOM[top_ind[j]][dirx];
			a2[diry] = lipidCOM[top_ind[j]][diry];
			pbc_dx(&pbc,a1,a2,dx);
			dist = norm2(dx);

			if(dist<=pr2)
			{
				if(counter==2 || framex[idprot[i]][dirz] == lipidCOM[top_ind[j]][dirz])
				{
					counter=2;
					break;
				}
				if(switch1==0 && framex[idprot[i]][dirz] < lipidCOM[top_ind[j]][dirz])
				{
					counter++;
					switch1++;
				}
				else
				 if(switch2==0 && framex[idprot[i]][dirz] > lipidCOM[top_ind[j]][dirz])
				 {
					 counter++;
					 switch2++;
				 }
			 }
		}
		if(counter==2)
		{
			ptop_ind[*nprot_top] = idprot[i];
			*nprot_top = *nprot_top + 1;
			save_time=1;
		}

		//bottom leaflet
		counter=0;
		switch1=0;
		switch2=0;
		if(save_time==0)
		{
			for(j=0; j<nlipbot; j++)
			{
				a2[dirx] = lipidCOM[bot_ind[j]][dirx];
				a2[diry] = lipidCOM[bot_ind[j]][diry];
				pbc_dx(&pbc,a1,a2,dx);
				dist = norm2(dx);

				if(dist<=pr2)
				{
					if(counter==2 || framex[idprot[i]][dirz] == lipidCOM[bot_ind[j]][dirz])
					{
						counter=2;
						break;
					}
					if(switch1==0 && framex[idprot[i]][dirz] < lipidCOM[bot_ind[j]][dirz])
					{
						counter++;
						switch1++;
					}
					else
					  if(switch2==0 && framex[idprot[i]][dirz] > lipidCOM[bot_ind[j]][dirz])
					  {
						  counter++;
						  switch2++;
					  }
					}
			}
			if(counter==2)
			{
				pbot_ind[*nprot_bot] = idprot[i];
				*nprot_bot = *nprot_bot + 1;
			}
		}
	}
}


//special function to detect protein atoms for order parameters
void protein_atoms_order(int nprot, real z_mid, rvec *framex,
                int dirx, int diry, int dirz, atom_id *idprot, int nliptop, int nlipbot,
                t_pbc pbc, real pr2, atom_id *top_ind, atom_id *bot_ind,
                int order_atom_num, int acyl_atom,
                atom_id **ptop_ind, atom_id **pbot_ind, int *nprot_top, int *nprot_bot,
                atom_id *idorder)
{
	int i=0, j=0, counter=0, switch1=0, switch2=0, save_time=0;
	int acyl_atom_ind;
	rvec a1 = {0.0, 0.0, 0.0};
	rvec a2 = {0.0, 0.0, 0.0};
	a1[dirz] = z_mid;
	a2[dirz] = z_mid;

	rvec dx = {0.0, 0.0, 0.0};
	real dist = 0.0;
	for(i=0; i<nprot; i++)
	{
		//top leaflet
		j=0;
		counter=0;
		switch1=0;
		switch2=0;
		save_time=0;
		a1[dirx] = framex[idprot[i]][dirx];
		a1[diry] = framex[idprot[i]][diry];

		for(j=0; j<nliptop; j++)
		{
			acyl_atom_ind = idorder[acyl_atom+top_ind[j]*order_atom_num];
			a2[dirx] = framex[acyl_atom_ind][dirx];
			a2[diry] = framex[acyl_atom_ind][diry];
			pbc_dx(&pbc,a1,a2,dx);
			dist = norm2(dx);

			if(dist<=pr2)
			{
				if(counter==2 || framex[idprot[i]][dirz] == framex[acyl_atom_ind][dirz])
				{
					counter=2;
					break;
				}
				if(switch1==0 && framex[idprot[i]][dirz] < framex[acyl_atom_ind][dirz])
				{
					counter++;
					switch1++;
				}
				else
				 if(switch2==0 && framex[idprot[i]][dirz] > framex[acyl_atom_ind][dirz])
				 {
					 counter++;
					 switch2++;
				 }
			 }
		}
		if(counter==2)
		{
			ptop_ind[acyl_atom-1][nprot_top[acyl_atom-1]] = idprot[i];
			nprot_top[acyl_atom-1]++;
			save_time=1;
		}

		//bottom leaflet
		counter=0;
		switch1=0;
		switch2=0;
		if(save_time==0)
		{
			for(j=0; j<nlipbot; j++)
			{
				acyl_atom_ind = idorder[acyl_atom+bot_ind[j]*order_atom_num];
				a2[dirx] = framex[acyl_atom_ind][dirx];
				a2[diry] = framex[acyl_atom_ind][diry];
				pbc_dx(&pbc,a1,a2,dx);
				dist = norm2(dx);

				if(dist<=pr2)
				{
					if(counter==2 || framex[idprot[i]][dirz] == framex[acyl_atom_ind][dirz])
					{
						counter=2;
						break;
					}
					if(switch1==0 && framex[idprot[i]][dirz] < framex[acyl_atom_ind][dirz])
					{
						counter++;
						switch1++;
					}
					else
					  if(switch2==0 && framex[idprot[i]][dirz] > framex[acyl_atom_ind][dirz])
					  {
						  counter++;
						  switch2++;
					  }
					}
			}
			if(counter==2)
			{
				pbot_ind[acyl_atom-1][nprot_bot[acyl_atom-1]] = idprot[i];
				nprot_bot[acyl_atom-1]++;
			}
		}
	}
}


