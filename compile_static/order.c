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

real calc_unsat(int id_i, int id_min_i, int id_plus_i, rvec *x, int axis, int found_first_unsat)
{
    rvec dist, Sx,Sy,Sz, foo,bar;
    real length=0.0, result=0.0;
    real gamma=0.0, phi=0.0, alpha=0.0;

    if(found_first_unsat==0) //first atom
    {
        //first get Sz, the vector from Cn to Cn+1
        rvec_sub(x[id_i], x[id_plus_i], Sz);
		svmul(1.0/norm(Sz),Sz,Sz);
		//calculate Sx
		rvec_sub(x[id_plus_i],x[id_i],foo);
		rvec_sub(x[id_min_i],x[id_i],bar);
		cprod(foo,bar,Sx);
		svmul(1.0/norm(Sx),Sx,Sx);
		//calculate Sy
		cprod(Sz,Sx,Sy);
		svmul(1.0/norm(Sy),Sy,Sy);

		//calculation of the order parameters
		svmul(1.0/norm(bar), bar, bar);
        gamma = acos(iprod(bar,Sz));
        svmul(1.0/norm(foo),foo,foo);
        phi = acos(iprod(foo,bar));
        alpha = -M_PI_2 + phi/2.0 + gamma;
        result = (cos(alpha)*cos(alpha)*0.5*(3.0*sqr(Sy[axis])-1.0)
                 + sin(alpha)*sin(alpha)*0.5*(3.0*sqr(Sz[axis])-1.0)
                 - 2.0*cos(alpha)*sin(alpha)*0.5*(3.0*Sz[axis]*Sy[axis]));
    }
    else if(found_first_unsat==1) //second atom
    {
        //first get Sz, the vector from Cn to Cn+1
        rvec_sub(x[id_min_i], x[id_i], Sz);  //for testing try reversing the direction
        svmul(1.0/norm(Sz), Sz, Sz);
		//calculate Sx
		rvec_sub(x[id_min_i],x[id_i],foo);
		rvec_sub(x[id_plus_i],x[id_i],bar);
//		cprod(foo,bar,Sx); //try swapping for testing
		cprod(bar,foo,Sx);
		svmul(1.0/norm(Sx),Sx,Sx);
		//calculate Sy
		cprod(Sz,Sx,Sy);
		svmul(1.0/norm(Sy),Sy,Sy);

		//calculation of the order parameters
		svmul(1.0/norm(bar),bar,bar);
		svmul(-1.0/norm(foo),foo,foo);
        gamma = acos(iprod(bar,foo));
        svmul(-1.0,foo,foo);
        phi = acos(iprod(foo,bar));
        alpha = -M_PI_2 + phi/2.0 + gamma;
        result = (cos(alpha)*cos(alpha)*0.5*(3.0*sqr(Sy[axis])-1.0)
                 + sin(alpha)*sin(alpha)*0.5*(3.0*sqr(Sz[axis])-1.0)
                 + 2.0*cos(alpha)*sin(alpha)*0.5*(3.0*Sz[axis]*Sy[axis]));
    }

	result *= (-1.0); //-Scd
    return(result);
}

void order_param(int order_atom_num1, int order_atom_num2, int norder1, int norder2,
		atom_id *idorder1, atom_id *idorder2, real **order_lip1, real **order_lip2,
		int lipid_num, rvec *x, int axis, real *order_sum1, real *order_sum2,
		real *order_sum1_sd, real *order_sum2_sd,
		int unsat, int nunsat1, atom_id *idunsat1, int nunsat2, atom_id *idunsat2)
{
	//order_lip[lipid_ID][atom_in_acyl_chain]
	int i=0, j=0, m=0;
	rvec Sx,Sy,Sz,foo,bar,cossum;
	int id_i, id_min_i, id_plus_i;
	int unsat_in_chain1 = 0;
	int *found_first_unsat;
	snew(found_first_unsat,lipid_num);

	//acyl chain sn-1
	for(i=1; i<order_atom_num1-1; i++) //go over all the atoms in one chain
	{
		for(j=0; j<lipid_num; j++) //go over all the chains
		{
			id_i = idorder1[i+j*order_atom_num1]; //atom i ID
			id_min_i = idorder1[i+j*order_atom_num1-1]; //atom i-1 ID
			id_plus_i = idorder1[i+j*order_atom_num1+1]; //atom i+1 ID
			if( (unsat>0) && (id_i==idunsat1[found_first_unsat[j]+j*2]) )
			{
				order_lip1[j][i-1] = calc_unsat(id_i,id_min_i,id_plus_i,x,axis,found_first_unsat[j]);
				found_first_unsat[j]=1;
				unsat_in_chain1=1;
			}
			else
			{
				//calculate Sz
				rvec_sub(x[id_plus_i],x[id_min_i],Sz); //vector Sz(Cn-1,Cn+1)
				svmul(1/norm(Sz),Sz,Sz);
				//calculate Sx
				rvec_sub(x[id_plus_i],x[id_i],foo);
				rvec_sub(x[id_min_i],x[id_i],bar);
				cprod(foo,bar,Sx);
				svmul(1/norm(Sx),Sx,Sx);
				//calculate Sy
				cprod(Sz,Sx,Sy);
				svmul(1/norm(Sy),Sy,Sy);
				//cosine and order param
				cossum[XX] = sqr(Sx[axis]);
				cossum[YY] = sqr(Sy[axis]);
				cossum[ZZ] = sqr(Sz[axis]);
				for (m = 0; m < DIM; m++)
				{
					foo[m] = 0.5 * (3 * cossum[m] - 1);
				}
				order_lip1[j][i-1]=(-1)*(0.6667*foo[XX] + 0.333*foo[YY]);
			}
			order_sum1[i-1] += order_lip1[j][i-1];
			order_sum1_sd[i-1] += pow(order_lip1[j][i-1],2);
//			printf("ID: %d %d    Scd: %f\n",i,j,order_lip1[j][i-1]);
		}
//		printf("ID: %d    Scd: %f\n",i,order_sum1[i-1]/lipid_num);
	}


	sfree(found_first_unsat);
	snew(found_first_unsat,lipid_num);

	//acyl chain sn-2
	for(i=1; i<order_atom_num2-1; i++) //go over all the atoms in one chain
	{
		for(j=0; j<lipid_num; j++) //go over all the chains
		{
			id_i = idorder2[i+j*order_atom_num2]; //atom i ID
			id_min_i = idorder2[i+j*order_atom_num2-1]; //atom i-1 ID
			id_plus_i = idorder2[i+j*order_atom_num2+1]; //atom i+1 ID

			if( (unsat==1) && (unsat_in_chain1==0) )
			{
				if(id_i==idunsat1[found_first_unsat[j]+j*2])
				{
					order_lip2[j][i-1] = calc_unsat(id_i,id_min_i,id_plus_i,x,axis,found_first_unsat[j]);
					found_first_unsat[j]=1;
				}
			}
			if( (unsat==2) && (id_i==idunsat2[found_first_unsat[j]+j*2]) )
			{
				order_lip2[j][i-1] = calc_unsat(id_i,id_min_i,id_plus_i,x,axis,found_first_unsat[j]);
				found_first_unsat[j]=1;
			}
			else
			{
				//calculate Sz
				rvec_sub(x[id_plus_i],x[id_min_i],Sz); //vector Sz(Cn-1,Cn+1)
				svmul(1/norm(Sz),Sz,Sz);
				//calculate Sx
				rvec_sub(x[id_plus_i],x[id_i],foo);
				rvec_sub(x[id_min_i],x[id_i],bar);
				cprod(foo,bar,Sx);
				svmul(1/norm(Sx),Sx,Sx);
				//calculate Sy
				cprod(Sz,Sx,Sy);
				svmul(1/norm(Sy),Sy,Sy);
				//cosine and order param
				cossum[XX] = sqr(Sx[axis]);
				cossum[YY] = sqr(Sy[axis]);
				cossum[ZZ] = sqr(Sz[axis]);
				for (m = 0; m < DIM; m++)
				{
					foo[m] = 0.5 * (3 * cossum[m] - 1);
				}
				order_lip2[j][i-1]=(-1)*(0.6667*foo[XX] + 0.333*foo[YY]);
			}
			order_sum2[i-1] += order_lip2[j][i-1];
			order_sum2_sd[i-1] += pow(order_lip2[j][i-1],2);
		}
//		printf("ID: %d    Scd: %f\n",i,order_sum2[i-1]/lipid_num);
	}

	sfree(found_first_unsat);
}

