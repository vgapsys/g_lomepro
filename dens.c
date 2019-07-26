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



// assign density to a grid cell
void assign_density(real ***dens_grid, real **apl_grid, int lipidGroup_num, int *nlipGroup, int **dictLipNum, int **dictLipNumInv, int aux_ind)
{
    int i,j,k;
    real densVal = 0.0;
    int lipidID = (int)apl_grid[aux_ind][0]; // ID of a lipid atom assigned to this grid cell

    for(i=0; i<lipidGroup_num; i++) // investigate each lipid group
    {
        densVal = 0.0;

        if(dictLipNum[lipidID][0]==i)
        {
            densVal = 1.0/apl_grid[aux_ind][1];
        }

        dens_grid[i][aux_ind][0] += densVal;
        dens_grid[i][aux_ind][1] += densVal*densVal;
    }
}

