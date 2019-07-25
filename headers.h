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
//#include <fftw3.h>

int get_ind(int i, int j, int Nx);

void protein_atoms(int nprot, real z_mid, rvec *framex, rvec *lipidCOM,
                int dirx, int diry, int dirz, atom_id *idprot, int nliptop, int nlipbot,
                t_pbc pbc, real pr2, atom_id *top_ind, atom_id *bot_ind,
                atom_id *ptop_ind, atom_id *pbot_ind, int *nprot_top, int *nprot_bot);

void protein_atoms_order(int nprot, real z_mid, rvec *framex,
                int dirx, int diry, int dirz, atom_id *idprot, int nliptop, int nlipbot,
                t_pbc pbc, real pr2, atom_id *top_ind, atom_id *bot_ind,
                int order_atom_num, int acyl_atom,
                atom_id **ptop_ind, atom_id **pbot_ind, int *nprot_top, int *nprot_bot,
                atom_id *idorder);

void filter_curve_abs(real *grid, real *out, int binx, int biny, real bin_sizex, real bin_sizey, real filt_low, real filt_high, int fr_num, gmx_bool filter_verbose);

void filter_curve_rel(real *grid, real *out, int binx, int biny, real bin_sizex, real bin_sizey, real filt_low, real filt_high, int fr_num, gmx_bool filter_verbose);

void curvature(int dirx, int diry, int dirz, int stepx, int stepy, real bin_size_x, real bin_size_y, int binx, int biny, int frame_num,
                real *grid, real *gausCurve, real *meanCurve, int up_down);

void fill_grid(int *is_cell_prot1, int *is_cell_prot2,
                int dirx, int diry, int dirz, rvec a1, rvec a2, int binx, int nliptop, int nlipbot,
                t_pbc pbc, int i, int j, real bin_sizex, real bin_sizey, real min_dist, int *aux_ind,
                int lip_ind, int k, int l, int *top_ind, rvec *lipidCOM,
                rvec dx, real dist, real *grid_up, real *grid_down, int *top_index,
                gmx_bool is_prot, int prot_ind, int *bot_ind, rvec *framex, real *height1, real *height2,
                atom_id *idlip, int nprot_top, int *ptop_ind, int *pbot_ind, int *bottom_index,
                int nprot_bot, real left_x, real left_y);

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
                int **order_count_up, int **order_count_down);

void fill_grid_diffus(int *is_cell_prot1, int *is_cell_prot2,
		int dirx, int diry, int dirz, rvec a1, rvec a2, int binx, int nliptop, int nlipbot,
		t_pbc pbc, int i, int j, real bin_sizex, real bin_sizey, real min_dist, int *aux_ind,
		int lip_ind, int k, int l, int *top_ind, rvec *lipidCOM,
		rvec dx, real dist, real *grid_up, real *grid_down, int *top_index,
		gmx_bool is_prot, int prot_ind, int *bot_ind, rvec *framex, real *height1, real *height2,
		atom_id *idlip, int nprot_top, int *ptop_ind, int *pbot_ind, int *bottom_index,
		int nprot_bot, real left_x, real left_y,
		int diffus_steps, rvec **diffus_grid_pos_up, rvec **diffus_grid_pos_down);

real calc_unsat(int id_i, int id_min_i, int id_plus_i, rvec *x, int axis, int second_unsat);

void order_param(int order_atom_num1, int order_atom_num2, int norder1, int norder2,
                atom_id *idorder1, atom_id *idorder2, real **order_lip1, real **order_lip2,
                int lipid_num, rvec *x, int axis, real *order_sum1, real *order_sum2,
                real *order_sum1_sd, real *order_sum2_sd,
                int unsat, int nunsat1, atom_id *idunsat1, int nunsat2, atom_id *idunsat2);

void assign_density(real ***dens_grid, real **apl_grid, int lipidGroup_num, int *nlipGroup, atom_id **idlipGroup, int aux_ind);

void get_xyz_for_pdb( real *x, real *y, real *z, real left_x, real left_y, real bin_sizex, real bin_sizey, int normal, gmx_bool swapxy, real grid_aux_ind, int i, int j);

