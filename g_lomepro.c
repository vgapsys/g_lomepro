/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

/*  #define _VERBOSE  */
#define true 1
#define false 0

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
//#include <fftw3.h>



int main(int argc,char *argv[])
{
  const char *desc[] = {
		  "[PAR]",
		    "________________[PAR]",
		    "Basic Options.[PAR]",
		    "----------------[PAR]",

    "g_lomepro is a grid based tool to calculate Local Membrane Properties ",
    "of molecular dynamics trajectories [REF].",
    "Thickness (-thick), area per lipid (APL) (-apl), Curvature (-curve) and ",
    "SCD Order parameters (-order) can be calculated.[PAR]",

    "The user should specify a reference structure file (option -s) and a trajectoy file (option -f).",
    "An index group (option -n) together with the option -num_lip are used to specify the representative lipid atom(s)",
    "(Phosphorous, headgroups, whole lipid, etc.) and the total number of lipids.[PAR]",

    "The number of grid elements in x and y axes can be specified (default 10) with -binx and -biny.[PAR]",

    "With the option -prot a group of atoms (i.e. Protein) embedded in the membrane can be",
    "explicitly considered for membrane property calculation.",
    "If relevant, the center of mass displacement removal of the Protein group in the",
    "trajectory file should be done *before* running g_lomepro (i.e. 'trjconv -fit').[PAR]",

    "-precision option (default 1.5) defines the radius (nm) used to search for lipids when considering protein atoms.[PAR]",

    "_________________________[PAR]",
    "Property Specific Options[PAR]",
    "-------------------------[PAR]",

    //"...........................................................................[PAR]",
    "Thickness (-thick):[PAR]",

    "output: time averages (avg) and standard deviation (sd)",
    "as matrix (dat) files and as pdb files with the thickness values in the B-factor field.[PAR]",

    "-prot_val (default 0) specifies the thickness value for the grid elements occupied by the protein group on both leaflets.[PAR]",

    "-scale factor (default 0.75) scales the thickness value in a grid element",
    "in the situations where a grid cell is assigned a protein atom on one leaflet and the corresponding",
    "grid cell on another leaflet has a lipid assigned[PAR]",


    "...........................................................................[PAR]",
    "APL option outputs (-apl):[PAR]",

    "output: time averages (avg) and standard deviation (sd)",
    "as matrix (dat) files and as pdb files with the APL values in the B-factor field.",
    "In addition APL values for every lipid/protein averaged over trajectory are printed out,",
    "as well as the APL value for every lipid/protein at every time step. [PAR]",

    "...........................................................................[PAR]",
    "Curvature (-curve):[PAR]",

    "output: time averaged mean and gaussian curvature for the top and bottom leaflets separately",
    "(in pdb and matrix.dat formats)[PAR]",

    "-r_filter_low (default 0) relatively set lower bound for filtering in",
    "frequency space (after FFT) controls high pass or band pass filters[PAR]",

    "-r_filter_high (default -1) relatively set higher bound for filtering in",
    "frequency space (after FFT) controls low pass or band pass filters[PAR]",

    "-q_filter_low  (default 0) absolute value of lower bound (nm^-1) for filtering in",
    "frequency space (after FFT) controls high pass or band pass filters[PAR]",

    "-q_filter_high (default 100000) absolute value of higher bound (nm^-1) for ",
    "filtering in frequency space (after FFT) controls low pass or band pass filters[PAR]",

    "-scale_mcurve (default 1) mean curvature values may be small, apply scaling",
    " when printing pdb (no scaling for matrix output)[PAR]",

    "-scale_gcurve (default 1) gaussian curvature values may be small, ",
    "apply scaling when printing pdb (no scaling for matrix output)[PAR]",

    "-inv_mean_curve (default FALSE) allows inverting the signs of the mean curvature.",
    "In the original g_lomepro paper, the sign of the mean curvature was determined from a point of view of an observer looking down onto a bilayer leaflet.",
    "This way the positive mean curvature was assigned to the surface which bent towards the bilayer center.",
    "Setting this flag to TRUE inverts the sign: the mean curvature becomes positive when the surface bends away from the bilayer midpoint.",


    "...........................................................................[PAR]",
    "Scd Order parameter (-order):[PAR]",

    "output: time averaged pdb and data matrices containing -Scd order parameters are written for",
    "every acyl chain atom for the top and bottom leaflets separately. ",
    "Also the average -Scd values for every lipid are printed out.[PAR]",

    "NOTE: many output files are generated (acyl_chain_atoms x top_down_leaflets)[PAR]",

    "Index file is required to contain separate groups for sn1 and sn2 chains with the acyl chain carbon atoms.[PAR]",

    "-unsat option (default 0) defines the number of unsaturated acyl chains [0|1|2].",
    "Currently one double bond is allowed per acyl chain.",
    "The carbons around a double bond need to be defined in a separate group in the index file. [PAR]",

    "In combination with -prot, -order_val controls the Scd values (default -1) asigned to the grid points occupied by protein atoms.[PAR]",

    "-nt (number of threads) option can be used to speed up the Scd calculation by running on a number of threads in parallel.[PAR]",


    "________________[PAR]",
    "Other Options.[PAR]",
    "----------------[PAR]",

    "-normal flag allows to selecting axis normal to the bilayer (default 2, z-axis).[PAR]",

    "If -breath option is set the grid is modified according to the box at every step.",
    "This may be useful for highly fluctuating boxes.[PAR]",

    "-swapxy may be useful in combination with the -normal flag, when the normal axis is not matching Oz.",
    "If the flag is set to TRUE, Ox and Oy axes are swapped. [PAR]",

    "-nonflat option can be used to treat highly curved bilayers.",
    "If the option is set, a group of lipid tail atoms will be required in an index file.",
    "The tail atoms will be used to split the leaflets into top/bottom layers.[PAR]",

    "-mov_pdb and -mov_dat can be used to export multiple-model pdb files or movie matrix files.",
    "Since pdb files may quickly increase in size, only one pdb file is generated, irrespective of the number of selected properties.",
    "The pdb movie file contains, the coordinates of a grid. In case -thick option is selected, the pdb movie file",
    "also contains membrane thickness in the B-factor field (for the other options the field is left blank).",
    "The supplementary Perl scripts distributed together with the g_lomepro can be used to conveniently",
    "merge matrix movie output with the pdb movie.",
    "-mov_dat prints data matrix for -thick, -apl and -order options.",
    "For -order a separate file for every acyl chain carbon is generated.",
    "A movie output option for the -curve is currently disabled. [PAR]",

    "The -smooth option defines the number of frames over which a running average will be calculated",
    "to smoothen the movie output.[PAR]",

    "---------\n"
    " Citation:\n",
    "---------\n",
    "Vytautas Gapsys, Bert L. de Groot, Rodolfo Briones\n",
    "Computational analysis of local membrane properties\n",
    "Journal of Computer-Aided Molecular Design, 2013\n",
    "doi:10.1007/s10822-013-9684-0\n",
    "software version: 1.0.2\n",
    "---------\n",
    "[PAR]"
  };

  t_filenm fnm[] = {
      { efTRX, "-f", NULL, ffREAD },
      { efNDX, "-n", "index.ndx",    ffREAD},
      { efTPS, "-s", "tpr.tpr",    ffREAD  },
      { efOUT, "-thick", "thickness",    ffOPTWR  },
      { efOUT, "-apl", "apl",  ffOPTWR},
      { efOUT, "-curve", "curvature",  ffOPTWR},
      { efOUT, "-order", "order",  ffOPTWR},
      //{ efOUT, "-diffus", "diffusion",  ffOPTWR},
      { efTRX, "-mov_pdb", "mov_pdb.pdb",  ffOPTWR},
      { efOUT, "-mov_mat", "mov_mat",  ffOPTWR}
  };

#define NFILE asize(fnm)

  /* Command line options */
  const char *in_file;
  const char *index_file;
//  const char *tpr_file;
  static int lipid_num=0;
  static gmx_bool is_prot = FALSE;
  static int binx = 10;
  static int biny = 10;
  static real prot_val = 0.0;
  static real order_val = -1.0;
  static real precision = 1.5;
  static real scale = 0.75;
  static int curve_step_x = 0;
  static int curve_step_y = 0;
  static real r_filter_low = 0.0;
  static real r_filter_high = -1.0;
  static real q_filter_low = 0.0;
  static real q_filter_high = 99999.99;
  static gmx_bool nonflat=FALSE;
  static gmx_bool inv_mean_curve=FALSE;
//  static gmx_bool rm_pbc=FALSE;
  static int normal=2;
  static gmx_bool swapxy=FALSE;
  static gmx_bool breath=FALSE;
  static int unsat=0;
  static int nt=1, smooth=1, i=0, j=0;
  static real mcurve_scale = 1.0, gcurve_scale = 1.0;
  int foo=0; //use this for a temporary integer
  output_env_t oenv;

  t_pargs pa[] = {
      { "-lip_num", FALSE, etINT, {&lipid_num},
    		  "Number of lipids" },
      { "-breath", FALSE,  etBOOL, {&breath},
    	      "If set, the grid is modified according to the box at every step" },
     { "-prot", FALSE,  etBOOL, {&is_prot},
             "Put this flag if there is a protein in your membrane" },
     { "-unsat", FALSE,  etINT, {&unsat},
             "Number of unsaturated acyl chains [0|1|2] (important for order parameters)" },
     { "-prot_val", FALSE,  etREAL, {&prot_val},
             "Value of the thickness in a grid element in case protein atoms are on both layers" },
     { "-order_val", FALSE,  etREAL, {&order_val},
            "Value of the order parameter Scd in a grid element in case it is assigned to a protein" },
     { "-scale", FALSE,  etREAL, {&scale},
             "Value of the thickness in a grid element in case protein atom is on one layer and lipid atom is on another layer" },
      { "-binx", FALSE, etINT, {&binx},
             "Number of bins in x direction"
          },
      { "-biny", FALSE, etINT, {&biny},
                  "Number of bins in y direction"
      },
      { "-precision", FALSE, etREAL, {&precision},
                     "Precision value: radius (nm) to search for the lipids when considering protein atoms"
         },
//     { "-curvature_step_x", FALSE, etINT, {&curve_step_x},
//                        "step size in x direction for curvature calculation. By default step is 0, thus making "
//                        "the step size to be equal 1 bin size in x direction"
//     },
//      { "-curvature_step_y", FALSE, etINT, {&curve_step_y},
//					  "step size in y direction for curvature calculation. By default step is 0, thus making "
//					  "the step size to be equal 1 bin size in y direction"
//       },
      { "-r_filter_low", FALSE, etREAL, {&r_filter_low},
                          "Relatively set lower bound for filtering in frequency space (after FFT)"
						  " controls high pass or band pass filters"
      },
      { "-r_filter_high", FALSE, etREAL, {&r_filter_high},
                          "Relatively set higher bound for filtering in frequency space (after FFT)"
						  " controls low pass or band pass filters"
      },
      { "-q_filter_low", FALSE, etREAL, {&q_filter_low},
                          "Absolute value of lower bound (nm^-1) for filtering in frequency space (after FFT)"
                                                  " controls high pass or band pass filters"
      },
      { "-q_filter_high", FALSE, etREAL, {&q_filter_high},
                          "Absolute value of higher bound (nm^-1) for filtering in frequency space (after FFT)"
                                                  " controls low pass or band pass filters"
      },
      { "-scale_mcurve", FALSE, etREAL, {&mcurve_scale},
                          "Mean curvature values may be small,"
						  "apply scaling when printing pdb (no scaling for matrix output)"
      },
      { "-scale_gcurve", FALSE, etREAL, {&gcurve_scale},
						  "Gaussian curvature values may be small,"
						  "apply scaling when printing pdb (no scaling for matrix output)"
      },
      { "-inv_mean_curve", FALSE,  etBOOL, {&inv_mean_curve},
               "Invert the sign of the mean curvature" },
      { "-nonflat", FALSE,  etBOOL, {&nonflat},
               "Put this flag if the membrane is highly curved" },
      { "-normal", FALSE,  etINT, {&normal},
                  "Choose axis normal to the surface of the membrane: 0 - x, 1 - y or 2 - z" },
      { "-swapxy", FALSE,  etBOOL, {&swapxy},
                  "If the normal to the bilayer is not parallel to Oz axis and the -normal flag was set to 0 or 1, it may be necessary to swap Ox and Oy axes" },
      { "-smooth", FALSE,  etINT, {&smooth},
                   "Number of frames over which to run averaging for the -mov_mat flag" },
      { "-nt", FALSE,  etINT, {&nt},
              "Number of threads to use (for -order option only)" },
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_BE_NICE ,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL,&oenv);


  in_file=         opt2fn("-f",NFILE,fnm);
  index_file=      opt2fn("-n",NFILE,fnm);
  //tpr_file=      opt2fn("-s",NFILE,fnm);

  //deal with the topology file
  t_topology top;
  int        ePBC;
  char       title[STRLEN];
  rvec       *xtop;
  matrix     box;
  t_pbc pbc;

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  set_pbc(&pbc,ePBC,box);


  int nlip, nprot, ntail, norder1, norder2, nunsat1, nunsat2; //number of atoms
  atom_id *idlip, *idprot, *idtail, *idorder1, *idorder2, *idunsat1, *idunsat2;   //IDs of the atoms
  char *namelip, *nameprot, *nametail, *nameorder; //names of the groups


  //check if lipid_num was specified
  if(lipid_num == 0)
  {
	  printf("\nMust specify the number of lipids: -lip_num \n");
	  return(0);
  }

  //deal with index files
  printf("\nChoose the lipid group\n");
  rd_index(index_file,1,&nlip,&idlip,&namelip);

  if(is_prot)
  {
	  printf("\nChoose the protein group:\n");
	  rd_index(index_file,1,&nprot,&idprot,&nameprot);
  }

  if(nonflat)
  {
	  printf("\nChoose the group of a lipid tail end atom:\n");
	  rd_index(index_file,1,&ntail,&idtail,&nametail);
  }

  int mean_curve_sign_up = 1;
  int mean_curve_sign_down = -1;
  if(inv_mean_curve)
  {
      mean_curve_sign_up = -1;
      mean_curve_sign_down = 1;
  }

  gmx_bool order = FALSE;
  order = opt2bSet("-order",NFILE,fnm);
  if(order)
  {
	  printf("\nChoose the group of atoms in sn-1 acyl chain for order parameter calculation:\n");
	  rd_index(index_file,1,&norder1,&idorder1,&nameorder);

	  printf("\nChoose the group of atoms in sn-2 acyl chain for order parameter calculation:\n");
	  rd_index(index_file,1,&norder2,&idorder2,&nameorder);

	  if(unsat==1)
	  {
		  printf("\nChoose the group of atoms for unsaturated lipids (two atoms around double bond for every lipid):\n");
		  rd_index(index_file,1,&nunsat1,&idunsat1,&nameorder);
	  }
	  else if(unsat==2)
	  {
		  printf("\nChoose the group of atoms for unsaturated lipids from the chain sn-1 (two atoms around double bond for every lipid):\n");
		  rd_index(index_file,1,&nunsat1,&idunsat1,&nameorder);
		  printf("\nChoose the group of atoms for unsaturated lipids from the chain sn-2 (two atoms around double bond for every lipid):\n");
		  rd_index(index_file,1,&nunsat2,&idunsat2,&nameorder);
	  }
  }


  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/
  gmx_bool pdb = opt2bSet("-mov_pdb",NFILE,fnm);
  gmx_bool mat = opt2bSet("-mov_mat",NFILE,fnm);


  /*------------ THICKNESS -------------------*/
  /*------------ THICKNESS -------------------*/
  char *thick_name_avg_dat, *thick_name_sd_dat, *thick_name_avg_pdb, *thick_name_sd_pdb;
  FILE *thick_fp_avg_dat, *thick_fp_sd_dat, *thick_fp_avg_pdb, *thick_fp_sd_pdb;
  char *thick_name_mov_pdb, *thick_name_mov_mat;
  FILE *fp_mov_pdb_thick, *fp_mov_mat_thick;

  //for lipid index and thickness of each
  char *thick_name_lipids_up, *thick_name_lipids_down, *thick_name_over_time;
  FILE *thick_fp_lipids_up, *thick_fp_lipids_down, *thick_fp_over_time;

  gmx_bool thick = opt2bSet("-thick",NFILE,fnm);
  if(thick)
  {
	  const char *foo_name = opt2fn("-thick",NFILE,fnm);

	  //first deal with lipid index files
	  snew(thick_name_lipids_up,strlen(foo_name)+20);
	  strcpy(thick_name_lipids_up,foo_name);
	  strcat(thick_name_lipids_up,"_up_lipids.dat");
	  snew(thick_name_lipids_down,strlen(foo_name)+20);
	  strcpy(thick_name_lipids_down,foo_name);
	  strcat(thick_name_lipids_down,"_down_lipids.dat");
	  snew(thick_name_over_time,strlen(foo_name)+20);
	  strcpy(thick_name_over_time,foo_name);
	  strcat(thick_name_over_time,"_over_time.dat");
	  //lipid index
	  thick_fp_lipids_up = fopen(thick_name_lipids_up,"w");
	  thick_fp_lipids_down = fopen(thick_name_lipids_down,"w");
	  thick_fp_over_time = fopen(thick_name_over_time,"w");

	  // then the others
	  snew(thick_name_avg_dat,strlen(foo_name)+20);
	  strcpy(thick_name_avg_dat,foo_name);
	  strcat(thick_name_avg_dat,"_avg.dat");

	  snew(thick_name_sd_dat,strlen(foo_name)+20);
	  strcpy(thick_name_sd_dat,foo_name);
	  strcat(thick_name_sd_dat,"_sd.dat");

	  snew(thick_name_avg_pdb,strlen(foo_name)+20);
	  strcpy(thick_name_avg_pdb,foo_name);
	  strcat(thick_name_avg_pdb,"_avg.pdb");

	  snew(thick_name_sd_pdb,strlen(foo_name)+20);
	  strcpy(thick_name_sd_pdb,foo_name);
	  strcat(thick_name_sd_pdb,"_sd.pdb");

	  thick_fp_sd_dat = fopen(thick_name_sd_dat,"w");
	  thick_fp_avg_pdb = fopen(thick_name_avg_pdb,"w");
	  thick_fp_sd_pdb = fopen(thick_name_sd_pdb,"w");
	  thick_fp_avg_dat = fopen(thick_name_avg_dat,"w");

	  if(pdb)
	  {
		  foo_name = opt2fn("-mov_pdb",NFILE,fnm);
		  snew(thick_name_mov_pdb,strlen(foo_name)+20);
		  strcpy(thick_name_mov_pdb,foo_name);
		  strcat(thick_name_mov_pdb,"_thickness.pdb");
		  fp_mov_pdb_thick=ffopen(thick_name_mov_pdb,"w");
	  }
	  if(mat)
	  {
		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(thick_name_mov_mat,strlen(foo_name)+20);
		  strcpy(thick_name_mov_mat,foo_name);
		  strcat(thick_name_mov_mat,"_thickness.dat");
		  fp_mov_mat_thick=ffopen(thick_name_mov_mat,"w");
	  }
  }
  else //pdb movie may be outputted without thickness being declared
	  if(pdb)
	  {
		  const char *foo_name = opt2fn("-mov_pdb",NFILE,fnm);
		  fp_mov_pdb_thick=ffopen(foo_name,"w");
	  }
  /*------------ THICKNESS -------------------*/
  /*------------ THICKNESS -------------------*/


  /*------------- APL ---------------*/
  /*------------- APL ---------------*/
  char *apl_up_name_avg_dat,*apl_up_name_sd_dat,*apl_down_name_avg_dat,*apl_down_name_sd_dat;
  char *apl_name_avg_pdb,*apl_name_sd_pdb;
  FILE *apl_up_fp_avg_dat,*apl_up_fp_sd_dat,*apl_down_fp_avg_dat,*apl_down_fp_sd_dat;
  FILE *apl_fp_avg_pdb,*apl_fp_sd_pdb;
  char *apl_up_name_mov_mat, *apl_down_name_mov_mat;
  FILE *fp_mov_mat_apl_up, *fp_mov_mat_apl_down;

  //for lipid index and apl of each
  char *apl_name_lipids_up, *apl_name_lipids_down, *apl_name_over_time;
  FILE *apl_fp_lipids_up, *apl_fp_lipids_down, *apl_fp_over_time;


  gmx_bool apl = opt2bSet("-apl",NFILE,fnm);
  if(apl)
  {
	  const char *foo_name =	opt2fn("-apl",NFILE,fnm);

	  //first deal with lipid index files
	  snew(apl_name_lipids_up,strlen(foo_name)+20);
	  strcpy(apl_name_lipids_up,foo_name);
	  strcat(apl_name_lipids_up,"_up_lipids.dat");
	  snew(apl_name_lipids_down,strlen(foo_name)+20);
	  strcpy(apl_name_lipids_down,foo_name);
	  strcat(apl_name_lipids_down,"_down_lipids.dat");
	  snew(apl_name_over_time,strlen(foo_name)+20);
	  strcpy(apl_name_over_time,foo_name);
	  strcat(apl_name_over_time,"_over_time.dat");

	  //then the others
	  snew(apl_up_name_avg_dat,strlen(foo_name)+20);
	  strcpy(apl_up_name_avg_dat,foo_name);
	  strcat(apl_up_name_avg_dat,"_up_avg.dat");
	  snew(apl_down_name_avg_dat,strlen(foo_name)+20);
	  strcpy(apl_down_name_avg_dat,foo_name);
	  strcat(apl_down_name_avg_dat,"_down_avg.dat");

	  snew(apl_up_name_sd_dat,strlen(foo_name)+20);
	  strcpy(apl_up_name_sd_dat,foo_name);
	  strcat(apl_up_name_sd_dat,"_up_sd.dat");
	  snew(apl_down_name_sd_dat,strlen(foo_name)+20);
	  strcpy(apl_down_name_sd_dat,foo_name);
	  strcat(apl_down_name_sd_dat,"_down_sd.dat");

	  snew(apl_name_avg_pdb,strlen(foo_name)+20);
	  strcpy(apl_name_avg_pdb,foo_name);
	  strcat(apl_name_avg_pdb,"_avg.pdb");

	  snew(apl_name_sd_pdb,strlen(foo_name)+20);
	  strcpy(apl_name_sd_pdb,foo_name);
	  strcat(apl_name_sd_pdb,"_sd.pdb");

	  //lipid index
	  apl_fp_lipids_up = fopen(apl_name_lipids_up,"w");
	  apl_fp_lipids_down = fopen(apl_name_lipids_down,"w");
	  apl_fp_over_time = fopen(apl_name_over_time,"w");
	  //others
	  apl_up_fp_avg_dat = fopen(apl_up_name_avg_dat,"w");
	  apl_down_fp_avg_dat = fopen(apl_down_name_avg_dat,"w");
	  apl_up_fp_sd_dat = fopen(apl_up_name_sd_dat,"w");
	  apl_down_fp_sd_dat = fopen(apl_down_name_sd_dat,"w");
	  apl_fp_avg_pdb = fopen(apl_name_avg_pdb,"w");
	  apl_fp_sd_pdb = fopen(apl_name_sd_pdb,"w");

	  if(mat)
	  {
		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(apl_up_name_mov_mat,strlen(foo_name)+20);
		  strcpy(apl_up_name_mov_mat,foo_name);
		  strcat(apl_up_name_mov_mat,"_up_apl.dat");
		  fp_mov_mat_apl_up=ffopen(apl_up_name_mov_mat,"w");

		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(apl_down_name_mov_mat,strlen(foo_name)+20);
		  strcpy(apl_down_name_mov_mat,foo_name);
		  strcat(apl_down_name_mov_mat,"_down_apl.dat");
		  fp_mov_mat_apl_down=ffopen(apl_down_name_mov_mat,"w");
	  }
  }
  /*------------- APL ---------------*/
  /*------------- APL ---------------*/



  /*------------- CURVATURE ---------------*/
  /*------------- CURVATURE ---------------*/
  char *mcurve_up_name_avg_dat,*gcurve_up_name_avg_dat,*mcurve_down_name_avg_dat,*gcurve_down_name_avg_dat;
  char *mcurve_name_avg_pdb,*gcurve_name_avg_pdb;
  FILE *mcurve_up_fp_avg_dat,*gcurve_up_fp_avg_dat,*mcurve_down_fp_avg_dat,*gcurve_down_fp_avg_dat;
  FILE *mcurve_fp_avg_pdb,*gcurve_fp_avg_pdb;
  char *mcurve_up_name_mov_mat, *mcurve_down_name_mov_mat;
  char *gcurve_up_name_mov_mat, *gcurve_down_name_mov_mat;
  FILE *fp_mov_mat_mcurve_up, *fp_mov_mat_mcurve_down;
  FILE *fp_mov_mat_gcurve_up, *fp_mov_mat_gcurve_down;

  gmx_bool curve= opt2bSet("-curve",NFILE,fnm);
  if(curve)
  {
	  const char *foo_name = opt2fn("-curve",NFILE,fnm);

	  //file names
	  snew(mcurve_up_name_avg_dat,strlen(foo_name)+20);
	  strcpy(mcurve_up_name_avg_dat,foo_name);
	  strcat(mcurve_up_name_avg_dat,"_mean_curve_up_avg.dat");
	  snew(mcurve_down_name_avg_dat,strlen(foo_name)+20);
	  strcpy(mcurve_down_name_avg_dat,foo_name);
	  strcat(mcurve_down_name_avg_dat,"_mean_curve_down_avg.dat");

	  snew(gcurve_up_name_avg_dat,strlen(foo_name)+20);
	  strcpy(gcurve_up_name_avg_dat,foo_name);
	  strcat(gcurve_up_name_avg_dat,"_gauss_curve_up_avg.dat");
	  snew(gcurve_down_name_avg_dat,strlen(foo_name)+20);
	  strcpy(gcurve_down_name_avg_dat,foo_name);
	  strcat(gcurve_down_name_avg_dat,"_gauss_curve_down_avg.dat");

	  snew(mcurve_name_avg_pdb,strlen(foo_name)+20);
	  strcpy(mcurve_name_avg_pdb,foo_name);
	  strcat(mcurve_name_avg_pdb,"_mean_curve_avg.pdb");

	  snew(gcurve_name_avg_pdb,strlen(foo_name)+20);
	  strcpy(gcurve_name_avg_pdb,foo_name);
	  strcat(gcurve_name_avg_pdb,"_gauss_curve_avg.pdb");

	  //file pointers
	  mcurve_up_fp_avg_dat = fopen(mcurve_up_name_avg_dat,"w");
	  mcurve_down_fp_avg_dat = fopen(mcurve_down_name_avg_dat,"w");
	  gcurve_up_fp_avg_dat = fopen(gcurve_up_name_avg_dat,"w");
	  gcurve_down_fp_avg_dat = fopen(gcurve_down_name_avg_dat,"w");
	  mcurve_fp_avg_pdb = fopen(mcurve_name_avg_pdb,"w");
	  gcurve_fp_avg_pdb = fopen(gcurve_name_avg_pdb,"w");

	  if(mat)
	  {
		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(mcurve_up_name_mov_mat,strlen(foo_name)+20);
		  strcpy(mcurve_up_name_mov_mat,foo_name);
		  strcat(mcurve_up_name_mov_mat,"_up_mean_curve.dat");
		  fp_mov_mat_mcurve_up=ffopen(mcurve_up_name_mov_mat,"w");

		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(mcurve_down_name_mov_mat,strlen(foo_name)+20);
		  strcpy(mcurve_down_name_mov_mat,foo_name);
		  strcat(mcurve_down_name_mov_mat,"_down_mean_curve.dat");
		  fp_mov_mat_mcurve_down=ffopen(mcurve_down_name_mov_mat,"w");

		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(gcurve_up_name_mov_mat,strlen(foo_name)+20);
		  strcpy(gcurve_up_name_mov_mat,foo_name);
		  strcat(gcurve_up_name_mov_mat,"_up_gauss_curve.dat");
		  fp_mov_mat_gcurve_up=ffopen(gcurve_up_name_mov_mat,"w");

		  foo_name = opt2fn("-mov_mat",NFILE,fnm);
		  snew(gcurve_down_name_mov_mat,strlen(foo_name)+20);
		  strcpy(gcurve_down_name_mov_mat,foo_name);
		  strcat(gcurve_down_name_mov_mat,"_down_gauss_curve.dat");
		  fp_mov_mat_gcurve_down=ffopen(gcurve_down_name_mov_mat,"w");
	  }
  }
  /*------------- CURVATURE ---------------*/
  /*------------- CURVATURE ---------------*/



  /*------------- ORDER PARAMETERS ---------------*/
  /*------------- ORDER PARAMETERS ---------------*/
  //number of atoms in a single chain
  int order_atom_num1=0, order_atom_num2=0;
  //up, sn-1 and sn-2 acyl chains
  char *order_up_name_avg_dat_sn1;
  char *order_up_name_avg_dat_sn2;
  //down, sn-1 and sn-2
  char *order_down_name_avg_dat_sn1;
  char *order_down_name_avg_dat_sn2;
  //pdb file names
  char *order_name_avg_dat_sn1_pdb;
  char *order_name_avg_dat_sn2_pdb;
  //average over all lipids (as g_order), sn-1 and sn-2 acyl chains
  char *order_AVG_name_sn1;
  char *order_AVG_name_sn2;
  //up, sn-1 and sn-2
  FILE **order_up_fp_avg_dat_sn1;
  FILE **order_up_fp_avg_dat_sn2;
  //down, sn-1 and sn2
  FILE **order_down_fp_avg_dat_sn1;
  FILE **order_down_fp_avg_dat_sn2;
  //pdb files, sn-1 and sn2
  FILE **order_fp_avg_pdb_sn1;
  FILE **order_fp_avg_pdb_sn2;
  //average over all lipids (as g_order), sn-1 and sn-2 acyl chains
  FILE *order_fp_AVG_sn1;
  FILE *order_fp_AVG_sn2;
  /**** movie files *******/
  char *order_up_name_mov_mat1, *order_up_name_mov_mat2;
  char *order_down_name_mov_mat1, *order_down_name_mov_mat2;
  FILE **fp_mov_mat_order_up1, **fp_mov_mat_order_up2;
  FILE **fp_mov_mat_order_down1, **fp_mov_mat_order_down2;

  if(order)
  {
  	  order_atom_num1 = norder1/lipid_num;
  	  order_atom_num2 = norder2/lipid_num;

  	  snew(order_up_fp_avg_dat_sn1,order_atom_num1-2);
  	  snew(order_down_fp_avg_dat_sn1,order_atom_num1-2);
  	  snew(order_up_fp_avg_dat_sn2,order_atom_num2-2);
  	  snew(order_down_fp_avg_dat_sn2,order_atom_num2-2);
  	  /*pdb file pointer*/
  	  snew(order_fp_avg_pdb_sn1,order_atom_num1-2);
  	  snew(order_fp_avg_pdb_sn2,order_atom_num2-2);
  	  /**** movie files *****/
  	  snew(fp_mov_mat_order_up1,order_atom_num1-2);
  	  snew(fp_mov_mat_order_down1,order_atom_num1-2);
  	  snew(fp_mov_mat_order_up2,order_atom_num2-2);
  	  snew(fp_mov_mat_order_down2,order_atom_num2-2);


	  const char *foo_name =	opt2fn("-order",NFILE,fnm);
	  const char *foo_mov_name=NULL;
	  if(mat)
	  {
		  foo_mov_name = opt2fn("-mov_mat",NFILE,fnm);
	  }

  	  /*** average over all lipids (as g_order) ***/
	  //sn1
	  snew(order_AVG_name_sn1,strlen(foo_name)+20);
	  strcpy(order_AVG_name_sn1,foo_name);
	  strcat(order_AVG_name_sn1,"_AVG_sn1.dat");
	  order_fp_AVG_sn1=ffopen(order_AVG_name_sn1,"w");
	  //sn2
	  snew(order_AVG_name_sn2,strlen(foo_name)+20);
	  strcpy(order_AVG_name_sn2,foo_name);
	  strcat(order_AVG_name_sn2,"_AVG_sn2.dat");
	  order_fp_AVG_sn2=ffopen(order_AVG_name_sn2,"w");

	  //sn1
	  for(i=2; i<order_atom_num1; i++)
	  {
		  if(mat)
		  {
			  //up mov_mat sn1
			  snew(order_up_name_mov_mat1,strlen(foo_mov_name)+20);
			  strcpy(order_up_name_mov_mat1,foo_mov_name);
			  strcat(order_up_name_mov_mat1,"_up_sn1_atom");
			  sprintf(order_up_name_mov_mat1,"%s%d",order_up_name_mov_mat1,i);
			  strcat(order_up_name_mov_mat1,".dat");
			  //down mov_mat sn1
			  snew(order_down_name_mov_mat1,strlen(foo_mov_name)+20);
			  strcpy(order_down_name_mov_mat1,foo_mov_name);
			  strcat(order_down_name_mov_mat1,"_down_sn1_atom");
			  sprintf(order_down_name_mov_mat1,"%s%d",order_down_name_mov_mat1,i);
			  strcat(order_down_name_mov_mat1,".dat");
			  //file pointers
			  fp_mov_mat_order_up1[i-2] = fopen(order_up_name_mov_mat1,"w"); sfree(order_up_name_mov_mat1);
			  fp_mov_mat_order_down1[i-2] = fopen(order_down_name_mov_mat1,"w"); sfree(order_down_name_mov_mat1);
		  }

		  //up avg sn1
		  snew(order_up_name_avg_dat_sn1,strlen(foo_name)+20);
		  strcpy(order_up_name_avg_dat_sn1,foo_name);
		  strcat(order_up_name_avg_dat_sn1,"_up_avg_sn1_atom");
		  sprintf(order_up_name_avg_dat_sn1,"%s%d",order_up_name_avg_dat_sn1,i);
		  strcat(order_up_name_avg_dat_sn1,".dat");
			  /*pdb avg sn1*/
			  snew(order_name_avg_dat_sn1_pdb,strlen(foo_name)+20);
			  strcpy(order_name_avg_dat_sn1_pdb,foo_name);
			  strcat(order_name_avg_dat_sn1_pdb,"_avg_sn1_atom");
			  sprintf(order_name_avg_dat_sn1_pdb,"%s%d",order_name_avg_dat_sn1_pdb,i);
			  strcat(order_name_avg_dat_sn1_pdb,".pdb");
		  //down avg sn1
		  snew(order_down_name_avg_dat_sn1,strlen(foo_name)+20);
		  strcpy(order_down_name_avg_dat_sn1,foo_name);
		  strcat(order_down_name_avg_dat_sn1,"_down_avg_sn1_atom");
		  sprintf(order_down_name_avg_dat_sn1,"%s%d",order_down_name_avg_dat_sn1,i);
		  strcat(order_down_name_avg_dat_sn1,".dat");
		  //file pointers
		  order_up_fp_avg_dat_sn1[i-2] = fopen(order_up_name_avg_dat_sn1,"w"); sfree(order_up_name_avg_dat_sn1);
//		  order_up_fp_sd_dat_sn1[i-2] = fopen(order_up_name_sd_dat_sn1,"w"); sfree(order_up_name_sd_dat_sn1);
		  order_down_fp_avg_dat_sn1[i-2] = fopen(order_down_name_avg_dat_sn1,"w"); sfree(order_down_name_avg_dat_sn1);
//		  order_down_fp_sd_dat_sn1[i-2] = fopen(order_down_name_sd_dat_sn1,"w"); sfree(order_down_name_sd_dat_sn1);
		  	  /*pdb file pointers*/
		  	  order_fp_avg_pdb_sn1[i-2] = fopen(order_name_avg_dat_sn1_pdb,"w"); sfree(order_name_avg_dat_sn1_pdb);
//		  	  order_fp_sd_pdb_sn1[i-2] = fopen(order_name_sd_dat_sn1_pdb,"w"); sfree(order_name_sd_dat_sn1_pdb);
	  }
	  //sn2
	  for(i=2; i<order_atom_num2; i++)
	  {
		  if(mat)
		  {
			  //up mov_mat sn2
			  snew(order_up_name_mov_mat2,strlen(foo_mov_name)+20);
			  strcpy(order_up_name_mov_mat2,foo_mov_name);
			  strcat(order_up_name_mov_mat2,"_up_sn2_atom");
			  sprintf(order_up_name_mov_mat2,"%s%d",order_up_name_mov_mat2,i);
			  strcat(order_up_name_mov_mat2,".dat");
			  //down mov_mat sn2
			  snew(order_down_name_mov_mat2,strlen(foo_mov_name)+20);
			  strcpy(order_down_name_mov_mat2,foo_mov_name);
			  strcat(order_down_name_mov_mat2,"_down_sn2_atom");
			  sprintf(order_down_name_mov_mat2,"%s%d",order_down_name_mov_mat2,i);
			  strcat(order_down_name_mov_mat2,".dat");
			  //file pointers
			  fp_mov_mat_order_up2[i-2] = fopen(order_up_name_mov_mat2,"w"); sfree(order_up_name_mov_mat2);
			  fp_mov_mat_order_down2[i-2] = fopen(order_down_name_mov_mat2,"w"); sfree(order_down_name_mov_mat2);
		  }

		  //up avg sn2
		  snew(order_up_name_avg_dat_sn2,strlen(foo_name)+20);
		  strcpy(order_up_name_avg_dat_sn2,foo_name);
		  strcat(order_up_name_avg_dat_sn2,"_up_avg_sn2_atom");
		  sprintf(order_up_name_avg_dat_sn2,"%s%d",order_up_name_avg_dat_sn2,i);
		  strcat(order_up_name_avg_dat_sn2,".dat");
			  /*pdb avg sn1*/
			  snew(order_name_avg_dat_sn2_pdb,strlen(foo_name)+20);
			  strcpy(order_name_avg_dat_sn2_pdb,foo_name);
			  strcat(order_name_avg_dat_sn2_pdb,"_avg_sn2_atom");
			  sprintf(order_name_avg_dat_sn2_pdb,"%s%d",order_name_avg_dat_sn2_pdb,i);
			  strcat(order_name_avg_dat_sn2_pdb,".pdb");
		  //down avg sn2
		  snew(order_down_name_avg_dat_sn2,strlen(foo_name)+20);
		  strcpy(order_down_name_avg_dat_sn2,foo_name);
		  strcat(order_down_name_avg_dat_sn2,"_down_avg_sn2_atom");
		  sprintf(order_down_name_avg_dat_sn2,"%s%d",order_down_name_avg_dat_sn2,i);
		  strcat(order_down_name_avg_dat_sn2,".dat");
		  //file pointers
		  order_up_fp_avg_dat_sn2[i-2] = fopen(order_up_name_avg_dat_sn2,"w"); sfree(order_up_name_avg_dat_sn2);
//		  order_up_fp_sd_dat_sn2[i-2] = fopen(order_up_name_sd_dat_sn2,"w"); sfree(order_up_name_sd_dat_sn2);
		  order_down_fp_avg_dat_sn2[i-2] = fopen(order_down_name_avg_dat_sn2,"w"); sfree(order_down_name_avg_dat_sn2);
//		  order_down_fp_sd_dat_sn2[i-2] = fopen(order_down_name_sd_dat_sn2,"w"); sfree(order_down_name_sd_dat_sn2);
	  	  	  /*pdb file pointers*/
	  	  	  order_fp_avg_pdb_sn2[i-2] = fopen(order_name_avg_dat_sn2_pdb,"w"); sfree(order_name_avg_dat_sn2_pdb);
//	  	  	  order_fp_sd_pdb_sn2[i-2] = fopen(order_name_sd_dat_sn2_pdb,"w"); sfree(order_name_sd_dat_sn2_pdb);
	  }
  }
  /*------------- ORDER PARAMETERS ---------------*/
  /*------------- ORDER PARAMETERS ---------------*/



  /*------------ DIFFUSION -------------------*/
  /*------------ DIFFUSION -------------------*/
  char *diffus_name_up_dat, *diffus_name_down_dat, *diffus_name_pdb_avg;
  FILE *diffus_fp_up_dat, *diffus_fp_down_dat, *diffus_fp_pdb_avg;

  gmx_bool diffus = FALSE; //opt2bSet("-diffus",NFILE,fnm);
  if(diffus)
  {
	  const char *foo_name = opt2fn("-diffus",NFILE,fnm);

	  snew(diffus_name_up_dat,strlen(foo_name)+20);
	  strcpy(diffus_name_up_dat,foo_name);
	  strcat(diffus_name_up_dat,"_up.dat");

	  snew(diffus_name_down_dat,strlen(foo_name)+20);
	  strcpy(diffus_name_down_dat,foo_name);
	  strcat(diffus_name_down_dat,"_down.dat");

	  snew(diffus_name_pdb_avg,strlen(foo_name)+20);
	  strcpy(diffus_name_pdb_avg,foo_name);
	  strcat(diffus_name_pdb_avg,"_avg.pdb");

	  diffus_fp_up_dat = fopen(diffus_name_up_dat,"w");
	  diffus_fp_down_dat = fopen(diffus_name_down_dat,"w");
	  diffus_fp_pdb_avg = fopen(diffus_name_pdb_avg,"w");

	  if(mat)
	  {

	  }
  }
  /*------------ DIFFUSION -------------------*/
  /*------------ DIFFUSION -------------------*/


  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/
  /************************ File names and booleans ***************************/


  real left_x = 0.0, left_y = 0.0;

  int dirz=0,dirx=1,diry=2;
  if(normal==0)
  {
	  dirz=0;dirx=1;diry=2;
	  if(swapxy)
	  {
		  dirz=0;dirx=2;diry=1;
	  }
  }
  if(normal==1)
  {
	  dirz=1;dirx=0;diry=2;
	  if(swapxy)
	  {
		  dirz=1;dirx=2;diry=0;
	  }
  }
  if(normal==2)
  {
	  dirz=2;dirx=0;diry=1;
	  if(swapxy)
	  {
		  dirz=2;dirx=1;diry=0;
	  }
  }

  //get the box dimensions: x, y and z
  real grid_y = box[diry][diry];
  real grid_x = box[dirx][dirx];
  real grid_z = box[dirz][dirz];

  //start reading the trajectory
  t_trxstatus *trxhandle;
  t_trxframe frame;



  if (!read_first_frame(oenv, &trxhandle,in_file,&frame,TRX_READ_X || TRX_READ_V || TRX_READ_F))
	gmx_fatal(FARGS,"Could not read first frame from trajectory %s",in_file);

  int nlip_group = nlip/lipid_num; //number of atoms in one lipid group


  int grid_size = binx*biny;
  real area_of_cell = grid_x*grid_y/grid_size;
  real pr2 = precision*precision;
  real *grid_up_avg=NULL;	//z-values up
  real *grid_down_avg=NULL; //z-values down
  snew(grid_up_avg,grid_size);
  snew(grid_down_avg,grid_size);
  int frame_num=0;
  real bin_sizex = 0.0;
  real bin_sizey = 0.0;

  /***** if -mat_mov declared *****/
  int counter_smooth=0; //tells when to print


  /*****************THICKNESS***********/
  real *grid=NULL;	//stores thickness
  real *grid_sd=NULL;	//stores thickness standard deviation
  real **grid_smooth_frames; //saves the last grids of the last #smooth frames
  real *grid_smooth_avg; //grid averaged over the last #smooth frames
  real **z_smooth_frames_up; //saves the last z-coordinates of the last #smooth frames
  real **z_smooth_frames_down;
  real *z_smooth_avg_up; //z-coord averaged over the last #smooth frames
  real *z_smooth_avg_down;
  real **thick_lip_up=NULL; //thickness by lipid index: thick_lip_up[lipid_num]=[idlip,thick,thick_sum,thick_sum^2,tmp_counter_of_cells_per_lipid]
  real **thick_lip_down=NULL; //thickness by lipid index
  if(thick)
  {
	  snew(grid,grid_size);
	  snew(grid_sd,grid_size);

	  if(mat || pdb)
	  {
		  snew(grid_smooth_avg,grid_size);
		  snew(grid_smooth_frames,smooth);
		  for(foo=0; foo<smooth; foo++)
		  {
			  snew(grid_smooth_frames[foo],grid_size);
		  }
	  }
       
          // thickness by lipid index
	  if(is_prot)
	  {
		  snew(thick_lip_up,lipid_num+1);
		  snew(thick_lip_down,lipid_num+1);
		  for(i=0;i<lipid_num+1;i++)
		  {
			  snew(thick_lip_up[i],5);
			  snew(thick_lip_down[i],5);
			  if(i<lipid_num)
			  {
				  thick_lip_up[i][0]=idlip[i];
				  thick_lip_down[i][0]=idlip[i];
			  }
			  else
			  {
				  thick_lip_up[i][0]=-1;
				  thick_lip_down[i][0]=-1;
			  }
		  }
	  }
	  else
	  {
		  snew(thick_lip_up,lipid_num);
		  snew(thick_lip_down,lipid_num);
		  int i=0;
		  for(i=0;i<lipid_num;i++)
		  {
			  snew(thick_lip_up[i],5);
			  thick_lip_up[i][0]=idlip[i];
			  snew(thick_lip_down[i],5);
			  thick_lip_down[i][0]=idlip[i];
		  }
	  }
  }

  if( pdb || (curve && mat) )
  {
	  snew(z_smooth_avg_up,grid_size);
	  snew(z_smooth_avg_down,grid_size);
	  snew(z_smooth_frames_up,smooth);
  	  snew(z_smooth_frames_down,smooth);
	  for(foo=0; foo<smooth; foo++)
	  {
		  snew(z_smooth_frames_up[foo],grid_size);
		  snew(z_smooth_frames_down[foo],grid_size);
	  }
  }
  /*****************THICKNESS***********/


  /*****************APL***********/
  real **apl_grid_up=NULL;
  real **apl_grid_down=NULL;
  real **apl_lip_up=NULL;
  real **apl_lip_down=NULL;
  real **apl_smooth_up_frames; //saves the last grids of the last #smooth frames
  real **apl_smooth_down_frames; //saves the last grids of the last #smooth frames
  real *apl_smooth_up_avg; //grid averaged over the last #smooth frames
  real *apl_smooth_down_avg; //grid averaged over the last #smooth frames
  real *apl_smooth_down_avg_Xinv;

  if(apl)
  {
	  snew(apl_grid_up,grid_size);
	  snew(apl_grid_down,grid_size);
	  int i=0;
	  for(i=0;i<grid_size;i++)
	  {
		  snew(apl_grid_up[i],3);
		  snew(apl_grid_down[i],3);
	  }
	  if(is_prot)
	  {
		  snew(apl_lip_up,lipid_num+1);
		  snew(apl_lip_down,lipid_num+1);
		  for(i=0;i<lipid_num+1;i++)
		  {
			  snew(apl_lip_up[i],4);
			  snew(apl_lip_down[i],4);
			  if(i<lipid_num)
			  {
				  apl_lip_up[i][0]=idlip[i];
				  apl_lip_down[i][0]=idlip[i];
			  }
			  else
			  {
				  apl_lip_up[i][0]=-1;
				  apl_lip_down[i][0]=-1;
			  }
		  }
	  }
	  else
	  {
		  snew(apl_lip_up,lipid_num);
		  snew(apl_lip_down,lipid_num);
		  int i=0;
		  for(i=0;i<lipid_num;i++)
		  {
			  snew(apl_lip_up[i],4);
			  apl_lip_up[i][0]=idlip[i];
			  snew(apl_lip_down[i],4);
			  apl_lip_down[i][0]=idlip[i];
		  }
	  }

	  if(mat)
	  {
		  snew(apl_smooth_up_avg,grid_size);
		  snew(apl_smooth_down_avg,grid_size);
		  snew(apl_smooth_up_frames,smooth);
		  snew(apl_smooth_down_frames,smooth);
		  snew(apl_smooth_down_avg_Xinv,binx);
		  for(foo=0; foo<smooth; foo++)
		  {
			  snew(apl_smooth_up_frames[foo],grid_size);
			  snew(apl_smooth_down_frames[foo],grid_size);
		  }
	  }
  }
  /*****************APL***********/


  /*****************Order Parameters***********/
  real **order_lip1=NULL;	//order_lip[lipid_ID][atom_in_acyl_chain]
  real **order_lip2=NULL;
  real ***order_grid_up_sn1=NULL;	//order_grid[atom_in_acyl_chain][grid_element]
  //0 means avg, 1 standard deviation
  real ***order_grid_down_sn1=NULL;
  real ***order_grid_up_sn2=NULL;
  real ***order_grid_down_sn2=NULL;
  real *order_sum1=NULL; //order_sum1[atom_ID] sum of the order parameters over all lipids for one atom (as g_order output)
  real *order_sum2=NULL;
  real *order_sum1_sd=NULL; //order_sum1_sd[atom_ID]
  real *order_sum2_sd=NULL;
  int **ptop_ind_order1=NULL;   //protein indeces for each acyl atom in top leaflet
  int **ptop_ind_order2=NULL;
  int **pbot_ind_order1=NULL;   //protein indeces for each acyl atom in bottom leaflet
  int **pbot_ind_order2=NULL;
  int *nprot_top_order1=NULL; //number of protein atoms in top leaflet
  int *nprot_top_order2=NULL;
  int *nprot_bot_order1=NULL; //number of protein atoms in bottom leaflet
  int *nprot_bot_order2=NULL;
  real order_smooth_up1, order_smooth_up2; //in case -mov_mat is selected
  real order_smooth_down1, order_smooth_down2;
  real ***order_smooth_up_frames1, ***order_smooth_up_frames2; //[atom][smooth][grid]
  real ***order_smooth_down_frames1, ***order_smooth_down_frames2;
//  rvec *x_rm_pbc=NULL; //order parameters need x with rm_pbc
//  int order_atom_num1=0, order_atom_num2=0; //number of atoms in a single chain
  real **grid_up_order1=NULL;	//z-values up: grid_up_order[atom][grid_element]
  real **grid_down_order1=NULL; //z-values down
  real **grid_up_order2=NULL;
  real **grid_down_order2=NULL;
  int **order_count_sn1_up; //counting number of frames when lipids were assigned to the grids, needed when protein is considered
  int **order_count_sn1_down;
  int **order_count_sn2_up; //order_count_sn[atom_id][grid_cell]
  int **order_count_sn2_down;
  real **order_smooth_down1_Xinv;
  real **order_smooth_down2_Xinv;

  if(order)
  {
//	  snew(x_rm_pbc,top.atoms.nr);

	  //grid elements
	  snew(order_grid_up_sn1,order_atom_num1-2);
	  snew(order_grid_down_sn1,order_atom_num1-2);
	  snew(order_grid_up_sn2,order_atom_num2-2);
	  snew(order_grid_down_sn2,order_atom_num2-2);
	  snew(grid_up_order1,order_atom_num1-2);
	  snew(grid_down_order1,order_atom_num1-2);
	  snew(grid_up_order2,order_atom_num2-2);
	  snew(grid_down_order2,order_atom_num2-2);
	  snew(order_count_sn1_up,order_atom_num1-2);
	  snew(order_count_sn1_down,order_atom_num1-2);
	  snew(order_count_sn2_up,order_atom_num2-2);
	  snew(order_count_sn2_down,order_atom_num2-2);
	  for(i=0;i<order_atom_num1-2;i++)
	  {
		  snew(order_grid_up_sn1[i],grid_size);
		  snew(order_grid_down_sn1[i],grid_size);
		  for(j=0;j<grid_size;j++)
		  {
			  snew(order_grid_up_sn1[i][j],2);
			  snew(order_grid_down_sn1[i][j],2);
		  }
		  snew(grid_up_order1[i],grid_size);
		  snew(grid_down_order1[i],grid_size);
		  snew(order_count_sn1_up[i],grid_size);
		  snew(order_count_sn1_down[i],grid_size);
	  }
	  for(i=0;i<order_atom_num2-2;i++)
	  {
		  snew(order_grid_up_sn2[i],grid_size);
		  snew(order_grid_down_sn2[i],grid_size);
		  for(j=0;j<grid_size;j++)
		  {
			  snew(order_grid_up_sn2[i][j],2);
			  snew(order_grid_down_sn2[i][j],2);
		  }
		  snew(grid_up_order2[i],grid_size);
		  snew(grid_down_order2[i],grid_size);
		  snew(order_count_sn2_up[i],grid_size);
		  snew(order_count_sn2_down[i],grid_size);
	  }

	  //lipid indeces
	  snew(order_lip1,lipid_num);
	  snew(order_lip2,lipid_num);
	  for(i=0;i<lipid_num;i++)
	  {
		  snew(order_lip1[i],order_atom_num1-2);
		  snew(order_lip2[i],order_atom_num2-2);
	  }

	  //sum over lipids
	  snew(order_sum1,order_atom_num1-2);
	  snew(order_sum2,order_atom_num2-2);
	  snew(order_sum1_sd,order_atom_num1-2);
	  snew(order_sum2_sd,order_atom_num2-2);

	  //protein indeces
	  if(is_prot)
	  {
		  snew(nprot_top_order1,order_atom_num1-2);
		  snew(nprot_top_order2,order_atom_num2-2);
		  snew(nprot_bot_order1,order_atom_num1-2);
		  snew(nprot_bot_order2,order_atom_num2-2);
		  snew(ptop_ind_order1,order_atom_num1-2);
		  snew(ptop_ind_order2,order_atom_num2-2);
		  snew(pbot_ind_order1,order_atom_num1-2);
		  snew(pbot_ind_order2,order_atom_num2-2);
		  for(i=0;i<order_atom_num1-2;i++)
		  {
			  snew(ptop_ind_order1[i],nprot);
			  snew(pbot_ind_order1[i],nprot);
		  }
		  for(i=0;i<order_atom_num2-2;i++)
		  {
			  snew(ptop_ind_order2[i],nprot);
			  snew(pbot_ind_order2[i],nprot);
		  }
	  }
	  //flag -mov_mat
	  if(mat)
	  {
		  snew(order_smooth_up_frames1,order_atom_num1-2);
		  snew(order_smooth_down_frames1,order_atom_num1-2);
		  snew(order_smooth_up_frames2,order_atom_num2-2);
		  snew(order_smooth_down_frames2,order_atom_num2-2);
		  snew(order_smooth_down1_Xinv,order_atom_num1-2);
		  snew(order_smooth_down2_Xinv,order_atom_num2-2);
		  for(i=0;i<order_atom_num1-2;i++)
		  {
			  snew(order_smooth_up_frames1[i],smooth);
			  snew(order_smooth_down_frames1[i],smooth);
			  for(j=0;j<smooth;j++)
			  {
				  snew(order_smooth_up_frames1[i][j],grid_size);
				  snew(order_smooth_down_frames1[i][j],grid_size);
			  }
			  snew(order_smooth_down1_Xinv[i],binx);
		  }
		  for(i=0;i<order_atom_num2-2;i++)
		  {
			  snew(order_smooth_up_frames2[i],smooth);
			  snew(order_smooth_down_frames2[i],smooth);
			  for(j=0;j<smooth;j++)
			  {
				  snew(order_smooth_up_frames2[i][j],grid_size);
				  snew(order_smooth_down_frames2[i][j],grid_size);
			  }
			  snew(order_smooth_down2_Xinv[i],binx);
		  }
	  }
  }
  /*****************Order Parameters***********/


  /***************** Curvature ***********/
  real *mcurve_grid_up=NULL;
  real *mcurve_grid_down=NULL;
  real *gcurve_grid_up=NULL;
  real *gcurve_grid_down=NULL;
  real *mat_low_mcurve, *mat_low_gcurve;
  real *filtered_up = NULL;
  real *filtered_down = NULL;
  gmx_bool filter_verbose = TRUE;
  int curve_mat_frame_num = 10; //a constant

  if(curve)
  {
	  snew(mcurve_grid_up,grid_size);
	  snew(mcurve_grid_down,grid_size);
	  snew(gcurve_grid_up,grid_size);
	  snew(gcurve_grid_down,grid_size);
	  snew(mat_low_mcurve,binx);
	  snew(mat_low_gcurve,binx);

      if( (q_filter_low>0.0 || q_filter_high<99999.99) || (r_filter_low>0.0 || r_filter_high<1.0) )
      {
          snew(filtered_up,grid_size);
          snew(filtered_down,grid_size);
      }
  }
  /***************** Curvature ***********/


  /***************** DIFFUSION ***********/
  int diffus_steps = 10; //user provides this
  real *diffus_grid_offset_up, *diffus_grid_offset_down; //stores offset for every grid cell
  // the following variables have the form of: diffus[grid_cell][step]
  rvec **diffus_grid_pos_up, **diffus_grid_pos_down; //stores lipid positions for every grid cell
  real **diffus_grid_dist_up, **diffus_grid_dist_down; //stores lipid travelled distance for every grid cell
  //int **

  if(diffus)
  {
	  snew(diffus_grid_offset_up,grid_size);
	  snew(diffus_grid_offset_down,grid_size);

	  snew(diffus_grid_pos_up,grid_size);
	  snew(diffus_grid_pos_down,grid_size);
	  snew(diffus_grid_dist_up,grid_size);
	  snew(diffus_grid_dist_down,grid_size);

	  for(i=0; i<grid_size; i++)
	  {
		  snew(diffus_grid_pos_up[i],diffus_steps);
		  snew(diffus_grid_pos_down[i],diffus_steps);
		  snew(diffus_grid_dist_up[i],diffus_steps);
		  snew(diffus_grid_dist_down[i],diffus_steps);
	  }
  }
  /***************** DIFFUSION ***********/






  char pdbform[128];
  strcpy(pdbform,"%-6s%5u %-4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f");
  strcat(pdbform,"%6.2f%6.2f\n");
  char chA = 'A';
  char chB = 'B';
  int low_i=0;

  if(pdb)
  {
	  fprintf(fp_mov_pdb_thick,"TITLE     MEMBRANE\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
  }

  // these 6 parameters are for -mov_pdb flag (for both leaflets separately)
  real mov_pdb_x1=0.0, mov_pdb_y1=0.0, mov_pdb_z1=0.0;
  real mov_pdb_x2=0.0, mov_pdb_y2=0.0, mov_pdb_z2=0.0;

///////////////////////////////////////////////////////////////////
///////////////////// HERE STARTS the MAIN LOOP ///////////////////
///////////////////////////////////////////////////////////////////
  do
  {
	  set_pbc(&pbc,frame.ePBC,frame.box);

	  real z_min=9999.9, z_max=-9999.9, z_mid=0.0, right_x = -9999.9, right_y = -9999.9;
	  if(breath)
	  {
		  left_x = 9999.9; left_y = 9999.9;
	  }
	  rvec *lipidCOM=NULL;
	  snew(lipidCOM,lipid_num);

	  if(pdb || mat)
	  {
		  counter_smooth++;

		  if(pdb)
		  {
			  if(frame_num >= smooth-1)
			  {
				  fprintf(fp_mov_pdb_thick,"MODEL        %d\n",frame_num);
			  }
		  }
		  if(mat)
		  {

			  if(frame_num >= smooth-1)
			  {
				  if(thick)
				  {
					  fprintf(fp_mov_mat_thick,"FRAME %d\n",frame_num);
				  }
				  if(apl)
				  {
					  fprintf(fp_mov_mat_apl_up,"FRAME %d\n",frame_num);
					  fprintf(fp_mov_mat_apl_down,"FRAME %d\n",frame_num);
				  }
				  if(order)
				  {
					  for(i=2; i<order_atom_num1; i++)
					  {
						  fprintf(fp_mov_mat_order_up1[i-2],"FRAME %d\n",frame_num);
						  fprintf(fp_mov_mat_order_down1[i-2],"FRAME %d\n",frame_num);
					  }
					  for(i=2; i<order_atom_num2; i++)
					  {
						  fprintf(fp_mov_mat_order_up2[i-2],"FRAME %d\n",frame_num);
						  fprintf(fp_mov_mat_order_down2[i-2],"FRAME %d\n",frame_num);
					  }
				  }
				  if(curve)
				  {
					  fprintf(fp_mov_mat_mcurve_up,"FRAME %d\n",frame_num);
					  fprintf(fp_mov_mat_mcurve_down,"FRAME %d\n",frame_num);
					  fprintf(fp_mov_mat_gcurve_up,"FRAME %d\n",frame_num);
					  fprintf(fp_mov_mat_gcurve_down,"FRAME %d\n",frame_num);
				  }
			  }
		  }
	  }


	  /**************====0000000000====******************************/
	  //calculate order parameters
	  if(order)
	  {
		  // be aware that broken molecules will not be made whole
		  // process them prior starting the analysis
		  order_param(order_atom_num1,order_atom_num2,norder1,norder2,
				  idorder1,idorder2,order_lip1,order_lip2,lipid_num,frame.x,normal,
				  order_sum1,order_sum2,order_sum1_sd,order_sum2_sd,
				  unsat,nunsat1,idunsat1,nunsat2,idunsat2);
	  }


	  /**************====1111111111====******************************/
	  //go over the lipids for the first time: get COMs and z_mid
	  int i=0, mod=0, lip_count=-1;
	  real total_mass=0.0;
	  for(i=0; i<nlip; i++)
	  {
		  mod = i % nlip_group;
		  if(mod == 0)
		  {
			  if(i>0)
			  {
				  lipidCOM[lip_count][dirx] /= total_mass;
				  lipidCOM[lip_count][diry] /= total_mass;
				  lipidCOM[lip_count][dirz] /= total_mass;

				  if(lipidCOM[lip_count][dirz]>z_max)
				  {
					  z_max = lipidCOM[lip_count][dirz];
				  }
				  if(lipidCOM[lip_count][dirz]<z_min)
				  {
					  z_min = lipidCOM[lip_count][dirz];
				  }
				  if(breath)
				  {
                      if(lipidCOM[lip_count][dirx]<left_x)
                      {
                              left_x = lipidCOM[lip_count][dirx];
                      }
                      if(lipidCOM[lip_count][diry]<left_y)
                      {
                              left_y = lipidCOM[lip_count][diry];
                      }
                      if(lipidCOM[lip_count][dirx]>right_x)
                      {
                              right_x = lipidCOM[lip_count][dirx];
                      }
                      if(lipidCOM[lip_count][diry]>right_y)
                      {
                              right_y = lipidCOM[lip_count][diry];
                      }
				  }
			  }
			  lip_count++;
			  lipidCOM[lip_count][dirx] = 0.0;
			  lipidCOM[lip_count][diry] = 0.0;
			  lipidCOM[lip_count][dirz] = 0.0;
			  total_mass = 0.0;
		  }
		  lipidCOM[lip_count][dirx] += top.atoms.atom[idlip[i]].m * frame.x[idlip[i]][dirx];
		  lipidCOM[lip_count][diry] += top.atoms.atom[idlip[i]].m * frame.x[idlip[i]][diry];
		  lipidCOM[lip_count][dirz] += top.atoms.atom[idlip[i]].m * frame.x[idlip[i]][dirz];
		  total_mass += top.atoms.atom[idlip[i]].m;
	  }
	  lipidCOM[lip_count][dirx] /= total_mass;
	  lipidCOM[lip_count][diry] /= total_mass;
	  lipidCOM[lip_count][dirz] /= total_mass;
	  if(lipidCOM[lip_count][dirz]>z_max)
	  {
		  z_max = lipidCOM[lip_count][dirz];
	  }
	  if(lipidCOM[lip_count][dirz]<z_min)
	  {
		  z_min = lipidCOM[lip_count][dirz];
	  }
	  z_mid = (z_max+z_min)/2;
	  if(breath)
	  {
          if(lipidCOM[lip_count][dirx]<left_x)
          {
                  left_x = lipidCOM[lip_count][dirx];
          }
          if(lipidCOM[lip_count][diry]<left_y)
          {
                  left_y = lipidCOM[lip_count][diry];
          }
          if(lipidCOM[lip_count][dirx]>right_x)
          {
                  right_x = lipidCOM[lip_count][dirx];
          }
          if(lipidCOM[lip_count][diry]>right_y)
          {
                  right_y = lipidCOM[lip_count][diry];
          }
          area_of_cell = fabs(right_x-left_x)*fabs(right_y-left_y)/grid_size;
	  }
	  else
	  {
		  right_x = frame.box[dirx][dirx];
		  right_y = frame.box[diry][diry];
//                  printf("\ndirx %f %f %f\n",frame.box[dirx][0],frame.box[dirx][1],frame.box[dirx][2]);
//                 printf("diry %f %f %f\n",frame.box[diry][0],frame.box[diry][1],frame.box[diry][2]);
//                  printf("dirz %f %f %f\n",frame.box[dirz][0],frame.box[dirz][1],frame.box[dirz][2]);
	  }

	  /**************====22222222====******************************/
	  //create the grid
	  bin_sizex = (right_x-left_x)/binx;
	  bin_sizey = (right_y-left_y)/biny;

	  area_of_cell = bin_sizex*bin_sizey;
//          printf("CELL SIZE: %f\n",area_of_cell);
//          printf("right_x: %f\n",right_x);
//          printf("right_y: %f\n",right_y);
//          printf("left_x: %f\n",left_x);
//          printf("left_y: %f\n",left_y);
//          printf("bin_x: %d\n",binx);
//          printf("bin_y: %d\n",biny);
//          exit(0);

	  real *grid_up=NULL;	//stores z-values,
	  snew(grid_up,grid_size);

	  real *grid_down=NULL;	//stores z-values,
	  snew(grid_down,grid_size);

	  real *grid_thick=NULL;	//stores thickness
	  snew(grid_thick,grid_size);

	  int *top_ind=NULL;	//stores top leaflet lipid indeces
	  snew(top_ind,lipid_num);

	  int *bot_ind=NULL;	//stores bottom leaflet lipid indeces
	  snew(bot_ind,lipid_num);

	  int *ptop_ind=NULL;	//stores top leaflet protein indeces
	  if(is_prot)
	  {
		  snew(ptop_ind,nprot);
	  }

	  int *pbot_ind=NULL;	//stores bottom leaflet protein indeces
	  if(is_prot)
	  {
		  snew(pbot_ind,nprot);
	  }



	  /**************====3333333333====******************************/
	  //go over the lipids for the second time: divide them into top and bottom leaflets
//	  real z_max_top=-9999.9, z_min_top=9999.9, z_max_bot=-9999.9, z_min_bot=9999.9;
	  int nliptop = 0;
	  int nlipbot = 0;
	  for(i=0; i<lipid_num; i++)
	  {

		  if(nonflat)
		  {
			  if(lipidCOM[i][dirz]>=frame.x[idtail[i]][dirz])	//top leaflet
			  {
				  top_ind[nliptop] = i;
				  nliptop++;
			  }
			  else	//bottom leaflet
			  {
				  bot_ind[nlipbot] = i;
				  nlipbot++;
			  }
		  }
		  else
		  {
			  if(lipidCOM[i][dirz]>=z_mid)	//top leaflet
			  {
				  top_ind[nliptop] = i;
				  nliptop++;
			  }
			  else	//bottom leaflet
			  {
				  bot_ind[nlipbot] = i;
				  nlipbot++;
			  }
		  }
	  }



	  /**************====4444444444====******************************/
	  //go over the protein atoms
	  int nprot_top=0;
	  int nprot_bot=0;
	  rvec a1 = {0.0, 0.0, 0.0};
	  rvec a2 = {0.0, 0.0, 0.0};
	  a1[dirz] = z_mid;
	  a2[dirz] = z_mid;

	  rvec dx = {0.0, 0.0, 0.0};

	  if(is_prot)
	  {
		  protein_atoms(nprot,z_mid,frame.x,lipidCOM,dirx,diry,dirz,idprot,nliptop,nlipbot,
				  pbc,pr2,top_ind,bot_ind,ptop_ind,pbot_ind,&nprot_top,&nprot_bot);

		  if(order)
		  {
#pragma omp parallel num_threads(nt)
{ //starting parallel
#pragma omp for
			  for(foo=1; foo<order_atom_num1-1; foo++)
			  {
				  nprot_top_order1[foo-1] = 0;
				  nprot_bot_order1[foo-1] = 0;
				  protein_atoms_order(nprot,z_mid,frame.x,dirx,diry,dirz,idprot,nliptop,nlipbot,
					  pbc,pr2,top_ind,bot_ind, order_atom_num1, foo,
					  ptop_ind_order1, pbot_ind_order1, nprot_top_order1, nprot_bot_order1,
					  idorder1);
			  }
#pragma omp for
			  for(foo=1; foo<order_atom_num2-1; foo++)
			  {
				  nprot_top_order2[foo-1] = 0;
				  nprot_bot_order2[foo-1] = 0;
				  protein_atoms_order(nprot,z_mid,frame.x,dirx,diry,dirz,idprot,nliptop,nlipbot,
					  pbc,pr2,top_ind,bot_ind, order_atom_num2, foo,
					  ptop_ind_order2, pbot_ind_order2, nprot_top_order2, nprot_bot_order2,
					  idorder2);
			  }
}//ending parallel
		  }

	  }

	  /**************====5555555555====******************************/
	  //go over the grid elements and fill them
	  int k=0, l=0;
	  real height1=0.0, height2=0.0;
	  int top_index=0;
	  int bottom_index=0;

	  int is_cell_prot1=0, is_cell_prot2=0, aux_ind=0, lip_ind=0, prot_ind=0;
	  int time_saver_up=0, time_saver_down=0;
	  real min_dist=9999.9, dist=0.0;

//#pragma omp parallel num_threads(nt)
//{ //starting parallel
//#pragma omp for private(is_cell_prot1,is_cell_prot2,j,a1,a2,min_dist,aux_ind,lip_ind,k,dist,height1,height2,top_index,bottom_index,prot_ind,dx,l,time_saver_up,time_saver_down)
//#pragma omp for private(k,a1,a2,is_cell_prot1,is_cell_prot2,j,min_dist,aux_ind,lip_ind,dist,height1,height2,top_index,bottom_index,prot_ind,dx,l,time_saver_up,time_saver_down)
	  for(j=biny-1; j>=0; j--)
	  {
		  for(i=0; i<binx; i++)
		  {
			  if(diffus)
			  {
				  fill_grid_diffus(&is_cell_prot1,&is_cell_prot2,
						dirx,diry,dirz,a1,a2,binx,nliptop,nlipbot,
						pbc,i,j,bin_sizex,bin_sizey,min_dist,&aux_ind,
						lip_ind,k,l,top_ind,lipidCOM,
						dx,dist,grid_up,grid_down,&top_index,
						is_prot,prot_ind,bot_ind,frame.x,&height1,&height2,
						idlip,nprot_top,ptop_ind,pbot_ind,&bottom_index,
						nprot_bot,left_x,left_y,
						diffus_steps,diffus_grid_pos_up,diffus_grid_pos_down);
			  }
			  else
			  {
				  fill_grid(&is_cell_prot1,&is_cell_prot2,
						dirx,diry,dirz,a1,a2,binx,nliptop,nlipbot,
						pbc,i,j,bin_sizex,bin_sizey,min_dist,&aux_ind,
						lip_ind,k,l,top_ind,lipidCOM,
						dx,dist,grid_up,grid_down,&top_index,
						is_prot,prot_ind,bot_ind,frame.x,&height1,&height2,
						idlip,nprot_top,ptop_ind,pbot_ind,&bottom_index,
						nprot_bot,left_x,left_y);
			  }


/************************************ THICKNESS *******************************************/
			  if(is_cell_prot1==0 && is_cell_prot2==0)		//lipids on both layers
			  {
				  grid_thick[aux_ind] = grid_up[aux_ind] - grid_down[aux_ind];
			  }
			  else
				  if(is_cell_prot1+is_cell_prot2==1)	//lipid in one layer, protein in another
				  {
					  grid_thick[aux_ind] = scale*(grid_up[aux_ind] - grid_down[aux_ind]);
					  height1 = z_mid+scale*(height1-z_mid);
					  height2 = z_mid+scale*(height2-z_mid);
				  }
			  else
				  if(is_cell_prot1==1 && is_cell_prot2==1)	//proteins on both layers
				  {
					  grid_thick[aux_ind] = prot_val;
					  height1 = z_mid;
					  height2 = z_mid;
				  }

			  if(thick) //thickness
			  {
				  // thickness for every lipid by id
				  time_saver_up=0;
				  time_saver_down=0;
				  if(is_prot) //system with protein
				  {
					  if(top_index==-1)
					  {
						  thick_lip_up[lipid_num][1] += grid_thick[aux_ind];
						  thick_lip_up[lipid_num][4] = thick_lip_up[lipid_num][4] + 1.0;
						  time_saver_up=1;
					  }
					  if(bottom_index==-1)
					  {
						  thick_lip_down[lipid_num][1] += grid_thick[aux_ind];
						  thick_lip_down[lipid_num][4]=thick_lip_down[lipid_num][4]+1.0;
						  time_saver_down=1;
					  }
					  if(time_saver_up*time_saver_down==0)
					  {
						  for(k=0;k<lipid_num+1;k++)
						  {
							if(thick_lip_up[k][0]==top_index)
							{
								thick_lip_up[k][1] += grid_thick[aux_ind];
						                thick_lip_up[k][4]=thick_lip_up[k][4]+1.0;
								time_saver_up=1;
							}
							if(thick_lip_down[k][0]==bottom_index)
							{
								thick_lip_down[k][1] += grid_thick[aux_ind];
						                thick_lip_down[k][4]=thick_lip_down[k][4]+1.0;
								time_saver_down=1;
							}
							if(time_saver_up*time_saver_down==1)
							{
								break;
							}
						  }
					  }
				  }
				  else	//no protein
				  {
					  for(k=0;k<lipid_num;k++)
					  {
						if(thick_lip_up[k][0]==top_index)
						{

							thick_lip_up[k][1] += grid_thick[aux_ind];
						        thick_lip_up[k][4] = thick_lip_up[k][4] + 1.0;
							time_saver_up=1;
						}
						if(thick_lip_down[k][0]==bottom_index)
						{
							thick_lip_down[k][1] += grid_thick[aux_ind];
						        thick_lip_down[k][4] = thick_lip_down[k][4] + 1.0;
							time_saver_down=1;
						}
						if(time_saver_up*time_saver_down==1)
						{
							break;
						}
					  }
				  }

				  // matrices and pdb
				  if(mat || pdb)
				  {
					  if(frame_num < smooth) //first frames
					  {
						  grid_smooth_avg[aux_ind] += grid_thick[aux_ind]/smooth;
					  }
					  else //further frames
					  {
						  grid_smooth_avg[aux_ind] += grid_thick[aux_ind]/smooth
								  -grid_smooth_frames[counter_smooth-1][aux_ind]/smooth;
					  }
					  grid_smooth_frames[counter_smooth-1][aux_ind] = grid_thick[aux_ind];
				  }
				  if( pdb || (curve && mat) )
				  {
					  //deal with smoothening
					  //independent of normal
					  if(frame_num < smooth) //first frames
					  {
						  z_smooth_avg_up[aux_ind] += 10*height1/smooth;
						  z_smooth_avg_down[aux_ind] += 10*height2/smooth;
					  }
					  else //further frames
					  {
						  z_smooth_avg_up[aux_ind] += 10*height1/smooth
								  -z_smooth_frames_up[counter_smooth-1][aux_ind]/smooth;
						  z_smooth_avg_down[aux_ind] += 10*height2/smooth
								  -z_smooth_frames_down[counter_smooth-1][aux_ind]/smooth;
					  }

					  //part dependent on normal
					  if(frame_num >= smooth-1)
					  {
						  if(normal==0)
						  {
							  mov_pdb_x1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_x2 = z_smooth_avg_down[aux_ind];
							  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
							  if(swapxy)
							  {
								  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
						  if(normal==1)
						  {
							  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_y2 = z_smooth_avg_down[aux_ind];
							  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
							  if(swapxy)
							  {
								  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
						  if(normal==2)
						  {
							  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_z2 = z_smooth_avg_down[aux_ind];
							  if(swapxy)
							  {
								  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
					  }

					  z_smooth_frames_up[counter_smooth-1][aux_ind] = 10*height1;
					  z_smooth_frames_down[counter_smooth-1][aux_ind] = 10*height2;

					  if( (frame_num >= smooth-1) && pdb )
					  {
						  fprintf(fp_mov_pdb_thick,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,grid_smooth_avg[aux_ind]);
						  fprintf(fp_mov_pdb_thick,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2, 1.0,grid_smooth_avg[aux_ind]);
					  }
				  }
				  if(mat)
				  {
					  if(frame_num >= smooth-1)
					  {
						  fprintf(fp_mov_mat_thick,"%f       ",grid_smooth_avg[aux_ind]);
					  }
				  }
				  grid[aux_ind] += grid_thick[aux_ind];
				  grid_sd[aux_ind] += pow(grid_thick[aux_ind],2);
			  }
			  else //no thickness
			  {
				  if( pdb || (curve && mat) )
				  {
					  //deal with smoothening
					  //independent of normal
					  if(frame_num < smooth) //first frames
					  {
						  z_smooth_avg_up[aux_ind] += 10*height1/smooth;
						  z_smooth_avg_down[aux_ind] += 10*height2/smooth;
					  }
					  else //further frames
					  {
						  z_smooth_avg_up[aux_ind] += 10*height1/smooth
								  -z_smooth_frames_up[counter_smooth-1][aux_ind]/smooth;
						  z_smooth_avg_down[aux_ind] += 10*height2/smooth
								  -z_smooth_frames_down[counter_smooth-1][aux_ind]/smooth;
					  }

					  //normal dependent part
					  if(frame_num >= smooth-1)
					  {
						  if(normal==0)
						  {
							  mov_pdb_x1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_x2 = z_smooth_avg_down[aux_ind];
							  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
							  if(swapxy)
							  {
								  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
						  if(normal==1)
						  {
							  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_y2 = z_smooth_avg_down[aux_ind];
							  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
							  if(swapxy)
							  {
								  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
						  if(normal==2)
						  {
							  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
							  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
							  mov_pdb_z1 = z_smooth_avg_up[aux_ind];
							  mov_pdb_z2 = z_smooth_avg_down[aux_ind];
							  if(swapxy)
							  {
								  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
								  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
								  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
							  }
						  }
					  }

					  z_smooth_frames_up[counter_smooth-1][aux_ind] = 10*height1;
					  z_smooth_frames_down[counter_smooth-1][aux_ind] = 10*height2;

					  if( (frame_num >= smooth-1) && pdb )
					  {
						  fprintf(fp_mov_pdb_thick,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,0.0);
						  fprintf(fp_mov_pdb_thick,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2, 1.0,0.0);
					  }
				  }
			  }
			  grid_up_avg[aux_ind] += height1;
			  grid_down_avg[aux_ind] += height2;
/************************************ THICKNESS *******************************************/


/*****************************************APL*******************************************/
			  if(apl)
			  {
				  time_saver_up=0;
				  time_saver_down=0;
				  if(is_prot) //system with protein
				  {
					  if(top_index==-1)
					  {
						  apl_lip_up[lipid_num][1] += area_of_cell;
						  apl_grid_up[aux_ind][0]=-1;
						  time_saver_up=1;
					  }
					  if(bottom_index==-1)
					  {
						  apl_lip_down[lipid_num][1] += area_of_cell;
						  apl_grid_down[aux_ind][0]=-1;
						  time_saver_down=1;
					  }
					  if(time_saver_up*time_saver_down==0)
					  {
						  for(k=0;k<lipid_num+1;k++)
						  {
							if(apl_lip_up[k][0]==top_index)
							{
								apl_lip_up[k][1] += area_of_cell;
								apl_grid_up[aux_ind][0]=top_index;
								time_saver_up=1;
							}
							if(apl_lip_down[k][0]==bottom_index)
							{
								apl_lip_down[k][1] += area_of_cell;
								apl_grid_down[aux_ind][0]=bottom_index;
								time_saver_down=1;
							}
							if(time_saver_up*time_saver_down==1)
							{
								break;
							}
						  }
					  }
				  }
				  else	//no protein
				  {
					  for(k=0;k<lipid_num;k++)
					  {
						if(apl_lip_up[k][0]==top_index)
						{

							apl_lip_up[k][1] += area_of_cell;
							apl_grid_up[aux_ind][0]=top_index;
							time_saver_up=1;
						}
						if(apl_lip_down[k][0]==bottom_index)
						{
							apl_lip_down[k][1] += area_of_cell;
							apl_grid_down[aux_ind][0]=bottom_index;
							time_saver_down=1;
						}
						if(time_saver_up*time_saver_down==1)
						{
							break;
						}
					  }
				  }
			  }
/*****************************************APL*******************************************/

		  } //??????????? END OF THE INNER GRID ELEMENT LOOP ???????????//

		  if(mat)
		  {
			  if(thick && (frame_num >= smooth-1))
			  {
				  fprintf(fp_mov_mat_thick,"\n");
			  }
		  }
	  } //??????????? END OF THE OUTTER GRID ELEMENT LOOP ???????????//
//} //end of parallel

// output thickness per lipid
if(thick)
{
    int number=0;
    if(is_prot)
    {  number=lipid_num+1;}
    else 
    {  number=lipid_num; }
    for(k=0;k<number;k++)
    {
        if( thick_lip_up[k][4]>0 )
        { thick_lip_up[k][1] = thick_lip_up[k][1]/thick_lip_up[k][4]; }
        if( thick_lip_down[k][4]>0 )
        { thick_lip_down[k][1] = thick_lip_down[k][1]/thick_lip_down[k][4]; }

	thick_lip_up[k][2] += thick_lip_up[k][1];
	thick_lip_up[k][3] += pow(thick_lip_up[k][1],2);

	thick_lip_down[k][2] += thick_lip_down[k][1];
	thick_lip_down[k][3] += pow(thick_lip_down[k][1],2);

        // set temporary counters to zero
        thick_lip_up[k][4] = 0.0;
        thick_lip_down[k][4] = 0.0;
    }

    //output thickness per lipid over time
    fprintf(thick_fp_over_time,"FRAME %d\n",frame_num);
    real foo=0.0, bar=0.0;
    int typecast_id=0;
    for(k=0;k<number;k++)
    {
        // lip up
	typecast_id = (int) thick_lip_up[k][0];
        if(typecast_id == -1) //for protein
        {
            foo = thick_lip_up[k][1];
            fprintf(thick_fp_over_time,"%f  protein_up      %f\n",frame.time,foo);
            thick_lip_up[k][1] = 0.0;
        }
        else //for lipids
        {
            foo = thick_lip_up[k][1];
            thick_lip_up[k][1]=0.0;
            if(foo != 0.0)
            {
                bar += foo/number;
                //time index thickness
                fprintf(thick_fp_over_time,"%f	%d      %f\n",frame.time,typecast_id,foo);
            }
        }

        // lip down
        typecast_id = (int) thick_lip_down[k][0];
        if(typecast_id == -1) //for protein
        {
            foo = thick_lip_down[k][1];
            fprintf(thick_fp_over_time,"%f  protein_down      %f\n",frame.time,foo);
            thick_lip_down[k][1] = 0.0;
        }
        else //for lipids
        {
            foo = thick_lip_down[k][1];
            thick_lip_down[k][1]=0.0;
            if(foo != 0.0)
            {
                bar += foo/number;
                //time index thickness
                fprintf(thick_fp_over_time,"%f	%d      %f\n",frame.time,typecast_id,foo);
            }
        }
    }
    fprintf(thick_fp_over_time,"%f	MEAN	%f\n",frame.time,bar);
}



/************************************ ORDER PARAM ****************************************/
if(order)
{
#pragma omp parallel num_threads(nt)
{ //starting parallel
#pragma omp for private(aux_ind,j,foo)
	for(i=0; i<binx; i++)
	{
		for(j=0; j<biny; j++)
	  	{
			aux_ind = get_ind(i,j,binx);
	  		//order parameters need specific grid filling
	  			//sn1
	  			for(foo=1; foo<order_atom_num1-1; foo++)
	  			{
	  				fill_grid_order(dirx,diry,dirz,binx,nliptop,nlipbot,
	  				  pbc,i,j,bin_sizex,bin_sizey,aux_ind,
	  				  lip_ind,k,l,top_ind,
	  				  order_grid_up_sn1,order_grid_down_sn1,
	  				  is_prot,bot_ind,frame.x,
	  				  idlip,nprot_top_order1,ptop_ind_order1,pbot_ind_order1,
	  				  nprot_bot_order1,idorder1,
	  				  order_atom_num1,foo,
	  				  order_lip1,lipidCOM,z_mid,order_val,
	  				  grid_up_order1[foo-1],grid_down_order1[foo-1],scale, left_x, left_y, order_count_sn1_up, order_count_sn1_down);
	  			 }
	  			 //sn2
	  			 for(foo=1; foo<order_atom_num2-1; foo++)
	  			 {
	  				 fill_grid_order(dirx,diry,dirz,binx,nliptop,nlipbot,
	  				  pbc,i,j,bin_sizex,bin_sizey,aux_ind,
	  				  lip_ind,k,l,top_ind,
	  				  order_grid_up_sn2,order_grid_down_sn2,
	  				  is_prot,bot_ind,frame.x,
	  				  idlip,nprot_top_order2,ptop_ind_order2,pbot_ind_order2,
	  				  nprot_bot_order2,idorder2,
	  				  order_atom_num2,foo,
	  				  order_lip2,lipidCOM,z_mid,order_val,
	  				  grid_up_order2[foo-1],grid_down_order2[foo-1],scale, left_x, left_y, order_count_sn2_up, order_count_sn2_down);
	  			 }
	  	}
	}
} //ending parallel

//cannot do this in parallel
if(mat)
{
	  for(j=biny-1; j>=0; j--)
	  {
		  for(i=0; i<binx; i++)
		  {
			aux_ind = get_ind(i,j,binx);
	  		//sn1
	  		for(foo=1; foo<order_atom_num1-1; foo++)
	  		{
	  			//the trick here is the same as in apl
	            //order_grid_[up|down]_sn[1|2][0] is already an average
	            if(frame_num < smooth) //first frames
	            {
	            	order_smooth_up1 =  order_grid_up_sn1[foo-1][aux_ind][0]/smooth;
	                order_smooth_down1 =  order_grid_down_sn1[foo-1][aux_ind][0]/smooth;
	            }
	            else //further frames
	            {
	            	order_smooth_up1 =  order_grid_up_sn1[foo-1][aux_ind][0]/smooth
	            		-order_smooth_up_frames1[foo-1][counter_smooth-1][aux_ind]/smooth;
	                order_smooth_down1 =  order_grid_down_sn1[foo-1][aux_ind][0]/smooth
	                	-order_smooth_down_frames1[foo-1][counter_smooth-1][aux_ind]/smooth;
	            }
	            	order_smooth_up_frames1[foo-1][counter_smooth-1][aux_ind] =
	            			order_grid_up_sn1[foo-1][aux_ind][0];
	            	order_smooth_down_frames1[foo-1][counter_smooth-1][aux_ind] =
	            			order_grid_down_sn1[foo-1][aux_ind][0];

	            if(frame_num >= smooth-1)
	            {
					fprintf(fp_mov_mat_order_up1[foo-1],"%f       ",order_smooth_up1);
					order_smooth_down1_Xinv[foo-1][binx-1-i] = order_smooth_down1;
	            }
	  		 }
	  		 //sn2
	  		for(foo=1; foo<order_atom_num2-1; foo++)
	  		{
	  			//the trick here is the same as in apl
	            //order_grid_[up|down]_sn[1|2][0] is already an average
	            if(frame_num < smooth) //first frames
	            {
	            	order_smooth_up2 =  order_grid_up_sn2[foo-1][aux_ind][0]/smooth;
	                order_smooth_down2 =  order_grid_down_sn2[foo-1][aux_ind][0]/smooth;
	            }
	            else //further frames
	            {
	            	order_smooth_up2 =  order_grid_up_sn2[foo-1][aux_ind][0]/smooth
	            			-order_smooth_up_frames2[foo-1][counter_smooth-1][aux_ind]/smooth;
	                order_smooth_down2 =  order_grid_down_sn2[foo-1][aux_ind][0]/smooth
	                		-order_smooth_down_frames2[foo-1][counter_smooth-1][aux_ind]/smooth;
	            }
	            	order_smooth_up_frames2[foo-1][counter_smooth-1][aux_ind] =
	            			order_grid_up_sn2[foo-1][aux_ind][0];
	            	order_smooth_down_frames2[foo-1][counter_smooth-1][aux_ind] =
	            			order_grid_down_sn2[foo-1][aux_ind][0];

	            if(frame_num >= smooth-1)
	            {
					fprintf(fp_mov_mat_order_up2[foo-1],"%f       ",order_smooth_up2);
					order_smooth_down2_Xinv[foo-1][binx-1-i] = order_smooth_down2;
	            }
	  		}
	  	}
  		for(foo=1; foo<order_atom_num1-1; foo++)
  		{
			  if(frame_num >= smooth-1)
			  {
				  for(low_i=0; low_i<binx; low_i++)
				  {
					  fprintf(fp_mov_mat_order_down1[foo-1],"%f       ",order_smooth_down1_Xinv[foo-1][low_i]);
				  }
				  fprintf(fp_mov_mat_order_up1[foo-1],"\n");
				  fprintf(fp_mov_mat_order_down1[foo-1],"\n");
			  }
  		}
  		for(foo=1; foo<order_atom_num2-1; foo++)
  		{
			  if(frame_num >= smooth-1)
			  {
				  for(low_i=0; low_i<binx; low_i++)
				  {
					  fprintf(fp_mov_mat_order_down2[foo-1],"%f       ",order_smooth_down2_Xinv[foo-1][low_i]);
				  }
				  fprintf(fp_mov_mat_order_up2[foo-1],"\n");
				  fprintf(fp_mov_mat_order_down2[foo-1],"\n");
			  }
  		}
	}
}
}
/************************************ ORDER PARAM ****************************************/


/************************************ APL SECOND GO *******************************************/
	  if(apl)
	  {
		  for(j=biny-1; j>=0; j--)
		  {
			  for(i=0; i<binx; i++)
			  {
				  time_saver_up=0;
				  time_saver_down=0;
				  aux_ind = get_ind(i,j,binx);
				  if(is_prot)
				  {
					  if(apl_grid_up[aux_ind][0]==-1)
					  {
						  apl_grid_up[aux_ind][1] += apl_lip_up[lipid_num][1];
						  apl_grid_up[aux_ind][2] += pow(apl_lip_up[lipid_num][1],2);
						  time_saver_up=1;
					  }
					  if(apl_grid_down[aux_ind][0]==-1)
					  {
						  apl_grid_down[aux_ind][1] += apl_lip_down[lipid_num][1];
						  apl_grid_down[aux_ind][2] += pow(apl_lip_down[lipid_num][1],2);
						  time_saver_down=1;
					  }
					  for(k=0;k<lipid_num+1;k++)
					  {
						  if(time_saver_up*time_saver_down==1)
						  {
							  break;
						  }
						  if(apl_lip_up[k][0]==apl_grid_up[aux_ind][0])
						  {
							  apl_grid_up[aux_ind][1] += apl_lip_up[k][1];
							  apl_grid_up[aux_ind][2] += pow(apl_lip_up[k][1],2);
							  time_saver_up=1;
						  }
						  if(apl_lip_down[k][0]==apl_grid_down[aux_ind][0])
						  {
							  apl_grid_down[aux_ind][1] += apl_lip_down[k][1];
							  apl_grid_down[aux_ind][2] += pow(apl_lip_down[k][1],2);
							  time_saver_down=1;
						  }
					  }
				  }
				  else
				  {
					  for(k=0;k<lipid_num;k++)
					  {
						  if(apl_lip_up[k][0]==apl_grid_up[aux_ind][0])
						  {
							  apl_grid_up[aux_ind][1] += apl_lip_up[k][1];
							  apl_grid_up[aux_ind][2] += pow(apl_lip_up[k][1],2);
							  time_saver_up=1;
						  }
						  if(apl_lip_down[k][0]==apl_grid_down[aux_ind][0])
						  {
							  apl_grid_down[aux_ind][1] += apl_lip_down[k][1];
							  apl_grid_down[aux_ind][2] += pow(apl_lip_down[k][1],2);
							  time_saver_down=1;
						  }
						  if(time_saver_up*time_saver_down==1)
						  {
							  break;
						  }
					  }
				  }

				  if(mat)
				  {
					  //the trick here is that apl_grid_[up|down] is already an average
					  if(frame_num < smooth) //first frames
					  {
						  apl_smooth_up_avg[aux_ind] =  apl_grid_up[aux_ind][1]/smooth;
						  apl_smooth_down_avg[aux_ind] =  apl_grid_down[aux_ind][1]/smooth;
					  }
					  else //further frames
					  {
						  apl_smooth_up_avg[aux_ind] = apl_grid_up[aux_ind][1]/smooth
								  -apl_smooth_up_frames[counter_smooth-1][aux_ind]/smooth;
						  apl_smooth_down_avg[aux_ind] = apl_grid_down[aux_ind][1]/smooth
								  -apl_smooth_down_frames[counter_smooth-1][aux_ind]/smooth;
					  }
					  apl_smooth_up_frames[counter_smooth-1][aux_ind] = apl_grid_up[aux_ind][1];
					  apl_smooth_down_frames[counter_smooth-1][aux_ind] = apl_grid_down[aux_ind][1];

					  if(frame_num >= smooth-1)
					  {
						  fprintf(fp_mov_mat_apl_up,"%f       ",apl_smooth_up_avg[aux_ind]);
						  apl_smooth_down_avg_Xinv[binx-1-i] = apl_smooth_down_avg[aux_ind];
					  }
				  }
			  } //???? end of inner grid loop for APL ????//

			  if(mat)
			  {
				  if(frame_num >= smooth-1)
				  {
					  for(low_i=0; low_i<binx; low_i++)
					  {
						  fprintf(fp_mov_mat_apl_down,"%f       ",apl_smooth_down_avg_Xinv[low_i]);
					  }
					  fprintf(fp_mov_mat_apl_up,"\n");
					  fprintf(fp_mov_mat_apl_down,"\n");
				  }
			  }

		  } //???? end of outter grid loop for APL ????//

		  int number=0;
		  if(is_prot)
		  {
			  number=lipid_num+1;
		  }
		  else
		  {
			  number=lipid_num;
		  }
		  for(k=0;k<number;k++)
		  {
			  apl_lip_up[k][2] += apl_lip_up[k][1];
			  apl_lip_up[k][3] += pow(apl_lip_up[k][1],2);
			  apl_lip_down[k][2] += apl_lip_down[k][1];
			  apl_lip_down[k][3] += pow(apl_lip_down[k][1],2);
		  }


		  //output apl per lipid over time
		  fprintf(apl_fp_over_time,"FRAME %d\n",frame_num);
		  real foo=0.0, bar=0.0;
		  int typecast_id=0;
		  for(k=0;k<number;k++)
		  {
			  // lip up
			  typecast_id = (int) apl_lip_up[k][0];
			  if(typecast_id == -1) //for protein
			  {
				  foo = apl_lip_up[k][1];
				  fprintf(apl_fp_over_time,"%f  protein_up      %f\n",frame.time,foo);
				  apl_lip_up[k][1] = 0.0;
			  }
			  else //for lipids
			  {
				  foo = apl_lip_up[k][1];
		  	      apl_lip_up[k][1]=0.0;
			      if(foo != 0.0)
			      {
				      bar += foo/number;
				      //time index apl
		  		      fprintf(apl_fp_over_time,"%f	%d      %f\n",frame.time,typecast_id,foo);
			      }
			  }

			  // lip down
			  typecast_id = (int) apl_lip_down[k][0];
			  if(typecast_id == -1) //for protein
			  {
                  foo = apl_lip_down[k][1];
                  fprintf(apl_fp_over_time,"%f  protein_down      %f\n",frame.time,foo);
                  apl_lip_down[k][1] = 0.0;
			  }
			  else //for lipids
			  {
				  foo = apl_lip_down[k][1];
		  	      apl_lip_down[k][1]=0.0;
			      if(foo != 0.0)
			      {
				      bar += foo/number;
				      //time index apl
		  		      fprintf(apl_fp_over_time,"%f	%d      %f\n",frame.time,typecast_id,foo);
			      }

			  }
		  }
		  fprintf(apl_fp_over_time,"%f	MEAN	%f\n",frame.time,bar);

	  }
	  /************************************ APL SECOND GO *******************************************/



	  /***************************************** CURVATURE *******************************************/
	  if( curve && mat ) //only for movie, otherwise curvature is calculated from the averaged thickness
	  {
		  //z_smooth_avg_[up|down] is in Angstroms already (for pdb output) and already divided by smooth
		  //inside of filtering.c z_smooth_avg_[up|down] is divided by the curve_mat_frame_num
		  //therefore curve_mat_frame_num=10
		  filter_verbose = FALSE;

		  if(frame_num >= smooth-1)
		  {
			  if(q_filter_low>0.0 || q_filter_high<99999.99) //absolute radius
			  {
				  filter_curve_abs(z_smooth_avg_up,filtered_up,binx,biny,bin_sizex,bin_sizey,q_filter_low,q_filter_high,curve_mat_frame_num,filter_verbose);
				  filter_curve_abs(z_smooth_avg_down,filtered_down,binx,biny,bin_sizex,bin_sizey,q_filter_low,q_filter_high,curve_mat_frame_num,filter_verbose);
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_up,gcurve_grid_up,mcurve_grid_up,mean_curve_sign_up);
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_down,gcurve_grid_down,mcurve_grid_down,mean_curve_sign_down);
			  }
			  else if(r_filter_low>0.0 || r_filter_high<1.0) //relative radius
			  {
				  filter_curve_rel(z_smooth_avg_up,filtered_up,binx,biny,bin_sizex,bin_sizey,r_filter_low,r_filter_high,curve_mat_frame_num,filter_verbose);
				  filter_curve_rel(z_smooth_avg_down,filtered_down,binx,biny,bin_sizex,bin_sizey,r_filter_low,r_filter_high,curve_mat_frame_num,filter_verbose);
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_up,gcurve_grid_up,mcurve_grid_up,mean_curve_sign_up);
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_down,gcurve_grid_down,mcurve_grid_down,mean_curve_sign_down);
			  }
			  else
			  {
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,curve_mat_frame_num,z_smooth_avg_up,gcurve_grid_up,mcurve_grid_up,mean_curve_sign_up);
				  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,curve_mat_frame_num,z_smooth_avg_down,gcurve_grid_down,mcurve_grid_down,mean_curve_sign_down);
			  }

			  for(j=biny-1; j>=0; j--) //outter loop
			  {
				  for(i=0; i<binx; i++) //inner loop
				  {
					  aux_ind = get_ind(i,j,binx);

					  fprintf(fp_mov_mat_gcurve_up,"%10.10f	",gcurve_grid_up[aux_ind]);
					  fprintf(fp_mov_mat_mcurve_up,"%10.10f	",mcurve_grid_up[aux_ind]);
					  mat_low_gcurve[binx-1-i] = gcurve_grid_down[aux_ind];
					  mat_low_mcurve[binx-1-i] = mcurve_grid_down[aux_ind];

				  } //???? end of inner grid loop for CURVATURE ????//

				  for(low_i=0; low_i<binx; low_i++)
				  {
					  fprintf(fp_mov_mat_mcurve_down,"%f       ",mat_low_mcurve[low_i]);
					  fprintf(fp_mov_mat_gcurve_down,"%f       ",mat_low_gcurve[low_i]);
				  }
				  fprintf(fp_mov_mat_mcurve_up,"\n");
				  fprintf(fp_mov_mat_mcurve_down,"\n");
				  fprintf(fp_mov_mat_gcurve_up,"\n");
				  fprintf(fp_mov_mat_gcurve_down,"\n");

			  } //???? end of outter grid loop for CURVATURE ????//
	  		}
	  }
	  /***************************************** CURVATURE *******************************************/

	  if(pdb)
	  {
		  if(frame_num >= smooth-1)
		  {
			  fprintf(fp_mov_pdb_thick,"TER\nENDMDL\n");
		  }
	  }

	  frame_num++;

	  if(mat || pdb)
	  {
		  if(counter_smooth == smooth)
		  {
			  counter_smooth = 0;
		  }
	  }


	  sfree(lipidCOM);
	  sfree(grid_up);
	  sfree(grid_down);
	  sfree(grid_thick);
	  sfree(top_ind);
	  sfree(bot_ind);
	  sfree(ptop_ind);
	  sfree(pbot_ind);

  }while (read_next_frame(oenv,trxhandle,&frame));
  //////////////////////////////////////////////////////////////
  ////////////////// main loop ended ///////////////////////////
  //////////////////////////////////////////////////////////////

  if(pdb)
  {
	  fprintf(fp_mov_pdb_thick,"END");
	  fclose(fp_mov_pdb_thick);
  }



  /*****************************************Curvature*******************************************/
  real *gausCurveUp = NULL;
  real *gausCurveDown = NULL;
  real *meanCurveUp = NULL;
  real *meanCurveDown = NULL;

  if(curve)
  {
	  snew(gausCurveUp,grid_size);
	  snew(meanCurveUp,grid_size);
	  snew(gausCurveDown,grid_size);
	  snew(meanCurveDown,grid_size);

	  filter_verbose = TRUE;

	  //if filtering needed
          if(q_filter_low>0.0 || q_filter_high<99999.99) //absolute radius
          {
                  filter_curve_abs(grid_up_avg,filtered_up,binx,biny,bin_sizex,bin_sizey,q_filter_low,q_filter_high,frame_num,filter_verbose);
                  filter_curve_abs(grid_down_avg,filtered_down,binx,biny,bin_sizex,bin_sizey,q_filter_low,q_filter_high,frame_num,filter_verbose);
                  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_up,gausCurveUp,meanCurveUp,mean_curve_sign_up);
                  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_down,gausCurveDown,meanCurveDown,mean_curve_sign_down);
          }
	  else if(r_filter_low>0.0 || r_filter_high<1.0) //relative radius
	  {
		  filter_curve_rel(grid_up_avg,filtered_up,binx,biny,bin_sizex,bin_sizey,r_filter_low,r_filter_high,frame_num,filter_verbose);
		  filter_curve_rel(grid_down_avg,filtered_down,binx,biny,bin_sizex,bin_sizey,r_filter_low,r_filter_high,frame_num,filter_verbose);
		  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_up,gausCurveUp,meanCurveUp,mean_curve_sign_up);
		  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,1,filtered_down,gausCurveDown,meanCurveDown,mean_curve_sign_down);
	  }
	  else
	  {
		  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,frame_num,grid_up_avg,gausCurveUp,meanCurveUp,mean_curve_sign_up);
		  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,frame_num,grid_down_avg,gausCurveDown,meanCurveDown,mean_curve_sign_down);
	  }
  }
  /*****************************************Curvature*******************************************/



  /*****************************************Order Parameters*******************************************/
  if(curve==-1)
  {
	  snew(gausCurveUp,grid_size);
	  snew(meanCurveUp,grid_size);
	  curvature(dirx,diry,dirz,curve_step_x,curve_step_y,bin_sizex,bin_sizey,binx,biny,frame_num,grid_up_avg,gausCurveUp,meanCurveUp,mean_curve_sign_up);
  }
  /*****************************************Order Parameters*******************************************/



  /////////////////////////////////////////////////////
  //////////////// Outputting averages ////////////////
  /////////////////////////////////////////////////////
  /* variables for the .dat matrix of lower leaflet */
  real *mat_low_apl_avg, *mat_low_apl_sd, **mat_low_order1, **mat_low_order2;
  //snew()

  /******* THICKNESS *****/
  if(thick)
  {
	  fprintf(thick_fp_avg_pdb,"TITLE     Thickness\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
	  fprintf(thick_fp_sd_pdb,"TITLE     Thickness\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
  }

  /********** APL ********/
  if(apl)
  {
	  fprintf(apl_fp_avg_pdb,"TITLE     Area per lipid\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
	  fprintf(apl_fp_sd_pdb,"TITLE     Area per lipid\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
	  snew(mat_low_apl_avg,binx);
	  snew(mat_low_apl_sd,binx);
  }

  /********** ORDER ********/
  if(order)
  {
  	  snew(mat_low_order1,order_atom_num1-2);
  	  snew(mat_low_order2,order_atom_num2-2);
	  for(i=2; i<order_atom_num1; i++)
	  {
		  fprintf(order_fp_avg_pdb_sn1[i-2],"TITLE     Order parameters Scd\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
		  fprintf(order_fp_avg_pdb_sn1[i-2],"TITLE     Order parameters Scd\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
		  snew(mat_low_order1[i-2],binx);
	  }
	  for(i=2; i<order_atom_num2; i++)
	  {
		  fprintf(order_fp_avg_pdb_sn2[i-2],"TITLE     Order parameters Scd\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
		  fprintf(order_fp_avg_pdb_sn2[i-2],"TITLE     Order parameters Scd\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
		  snew(mat_low_order2[i-2],binx);
	  }
  }

  /********** CURVATURE ********/
  if(curve)
  {
	  fprintf(gcurve_fp_avg_pdb,"TITLE     Gaussian curvature\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
	  fprintf(mcurve_fp_avg_pdb,"TITLE     Mean curvature\nCRYST1  %.3f  %.3f   %.3f  90.00  90.00  90.00\n",10*grid_x,10*grid_y,10*grid_z);
  }


  i=0; j=0;
  real apl_val_up=0.0;
  real apl_val_down=0.0;
  real apl_sd_up=0.0;
  real apl_sd_down=0.0;
  int first_bin = 0;
  int aux_ind = 0;
  for(j=biny-1; j>=0; j--)
  {
	  for(i=0; i<binx; i++)
	  {
		  aux_ind = get_ind(i,j,binx);
		  grid_up_avg[aux_ind] /= frame_num;
		  grid_down_avg[aux_ind] /= frame_num;
		  first_bin++;

		  /**************************************THICKNESS**********************************/
		  if(thick)
		  {
			  grid[aux_ind] /= frame_num;
			  if(frame_num>1)
			  {
				  grid_sd[aux_ind] = grid_sd[aux_ind]/(frame_num-1)-pow(grid[aux_ind],2)*frame_num/(frame_num-1);

				  if(grid_sd[aux_ind]<0.0) //could happen at numerical precision
				  {
					  grid_sd[aux_ind] = 0.0;
				  }
				  else
				  {
					  grid_sd[aux_ind] = sqrt(grid_sd[aux_ind]);
				  }

				  //grid_sd[aux_ind] = sqrt(fabs(grid_sd[aux_ind]/(frame_num-1)-pow(grid[aux_ind],2)*frame_num/(frame_num-1)))/grid[aux_ind];
			  }
			  else
				  if(frame_num==1)
				  {
					  grid_sd[aux_ind] = 0.0;
				  }

			  fprintf(thick_fp_avg_dat,"%f	",grid[aux_ind]);
			  fprintf(thick_fp_sd_dat,"%f	",grid_sd[aux_ind]);

			  if(normal==0)
			  {
				  mov_pdb_x1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_x2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==1)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_y2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==2)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_z2 = 10*grid_down_avg[aux_ind];
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  fprintf(thick_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,grid[aux_ind]);
			  fprintf(thick_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,grid[aux_ind]);
			  fprintf(thick_fp_sd_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,grid_sd[aux_ind]);
			  fprintf(thick_fp_sd_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,grid_sd[aux_ind]);
		  }
		  /**************************************THICKNESS**********************************/


		  /*****************************************APL*************************************/
		  if(apl)
		  {
			  apl_val_up = apl_grid_up[aux_ind][1]/frame_num;
			  apl_val_down = apl_grid_down[aux_ind][1]/frame_num;
			  if(frame_num>1)
			  {
				  apl_sd_up = apl_grid_up[aux_ind][2]/(frame_num-1)-pow(apl_val_up,2)*frame_num/(frame_num-1);
				  apl_sd_down = apl_grid_down[aux_ind][2]/(frame_num-1)-pow(apl_val_down,2)*frame_num/(frame_num-1);

				  if(apl_sd_up<0.0) //could happen at numerical precision
				  {
					  apl_sd_up = 0.0;
				  }
				  else
				  {
					  apl_sd_up = sqrt(apl_sd_up);
				  }

				  if(apl_sd_down<0.0) //could happen at numerical precision
				  {
					  apl_sd_down = 0.0;
				  }
				  else
				  {
					  apl_sd_down = sqrt(apl_sd_down);
				  }

				//	  apl_sd_up = sqrt(apl_grid_up[aux_ind][2]/(frame_num-1)-pow(apl_val_up,2)*frame_num/(frame_num-1))/apl_val_up;
				//	  apl_sd_down = sqrt(apl_grid_down[aux_ind][2]/(frame_num-1)-pow(apl_val_down,2)*frame_num/(frame_num-1))/apl_val_down;
			  }
			  else
				  if(frame_num==1)
				  {
					  apl_sd_up = 0.0;
					  apl_sd_down = 0.0;
				  }

			  if(normal==0)
			  {
				  mov_pdb_x1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_x2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==1)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_y2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==2)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_z2 = 10*grid_down_avg[aux_ind];
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  }
			  }

			  fprintf(apl_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,apl_val_up);
			  fprintf(apl_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,apl_val_down);
			  fprintf(apl_fp_sd_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,apl_sd_up);
			  fprintf(apl_fp_sd_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,apl_sd_down);

			  fprintf(apl_up_fp_avg_dat,"%f	",apl_val_up);
			  fprintf(apl_up_fp_sd_dat,"%f	",apl_sd_up);
			  mat_low_apl_avg[binx-1-i] = apl_val_down;
			  mat_low_apl_sd[binx-1-i] = apl_sd_down;
		  }
		  /*****************************************APL*******************************************/


		  /***********************************ORDER PARAM*************************************/
		  if(order)
		  {
			  for(foo=0; foo<order_atom_num1-2; foo++)
			  {
				  if(order_count_sn1_up[foo][aux_ind]==0) //grid occupied by protein only
				  { order_grid_up_sn1[foo][aux_ind][0] = order_val; }
				  else
				  { order_grid_up_sn1[foo][aux_ind][0] /= order_count_sn1_up[foo][aux_ind]; }

				  if(order_count_sn1_down[foo][aux_ind]==0) //grid occupied by protein only
				  { order_grid_down_sn1[foo][aux_ind][0] = order_val; }
				  else
				  { order_grid_down_sn1[foo][aux_ind][0] /= order_count_sn1_down[foo][aux_ind]; }

				  grid_up_order1[foo][aux_ind] /= frame_num; //coordinate along normal
				  grid_down_order1[foo][aux_ind] /= frame_num; //coordinate along normal

				  if(normal==0)
				  {
					  mov_pdb_x1 = 10*grid_up_order1[foo][aux_ind];
					  mov_pdb_x2 = 10*grid_down_order1[foo][aux_ind];
					  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
					  if(swapxy)
					  {
						  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  if(normal==1)
				  {
					  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y1 = 10*grid_up_order1[foo][aux_ind];
					  mov_pdb_y2 = 10*grid_down_order1[foo][aux_ind];
					  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
					  if(swapxy)
					  {
						  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  if(normal==2)
				  {
					  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*grid_up_order1[foo][aux_ind];
					  mov_pdb_z2 = 10*grid_down_order1[foo][aux_ind];
					  if(swapxy)
					  {
						  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  fprintf(order_fp_avg_pdb_sn1[foo],pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,order_grid_up_sn1[foo][aux_ind][0]);
				  fprintf(order_fp_avg_pdb_sn1[foo],pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,order_grid_down_sn1[foo][aux_ind][0]);

				  fprintf(order_up_fp_avg_dat_sn1[foo],"%f	",order_grid_up_sn1[foo][aux_ind][0]);
				  //fprintf(order_down_fp_avg_dat_sn1[foo],"%f	",order_grid_down_sn1[foo][aux_ind][0]);
				  mat_low_order1[foo][binx-1-i] = order_grid_down_sn1[foo][aux_ind][0];

				  //AVG over all lipids as in g_order
				  if(first_bin == 1)
				  {
				  	real avg = order_sum1[foo]/(frame_num*lipid_num);
				        real sd = order_sum1_sd[foo]/(frame_num*lipid_num) - pow(avg,2);
                                        if(sd<0.0) //could happen due to numerical precision
                                        {
                                                sd = 0.0;
                                        }
                                        else
                                        {
                                                sd = sqrt(sd);
                                        }
				        fprintf(order_fp_AVG_sn1,"%d	%f	%f\n",foo+2,avg,sd);
				  }
			  }
			  for(foo=0; foo<order_atom_num2-2; foo++)
			  {
				  if(order_count_sn2_up[foo][aux_ind]==0) //grid occupied by protein only
				  { order_grid_up_sn2[foo][aux_ind][0] = order_val; }
				  else
				  { order_grid_up_sn2[foo][aux_ind][0] /= order_count_sn2_up[foo][aux_ind]; }

				  if(order_count_sn2_down[foo][aux_ind]==0) //grid occupied by protein only
				  { order_grid_down_sn2[foo][aux_ind][0] = order_val; }
				  else
				  { order_grid_down_sn2[foo][aux_ind][0] /= order_count_sn2_down[foo][aux_ind]; }

				  grid_up_order2[foo][aux_ind] /= frame_num;
				  grid_down_order2[foo][aux_ind] /= frame_num;

				  if(normal==0)
				  {
					  mov_pdb_x1 = 10*grid_up_order2[foo][aux_ind];
					  mov_pdb_x2 = 10*grid_down_order2[foo][aux_ind];
					  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
					  if(swapxy)
					  {
						  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  if(normal==1)
				  {
					  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y1 = 10*grid_up_order2[foo][aux_ind];
					  mov_pdb_y2 = 10*grid_down_order2[foo][aux_ind];
					  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
					  if(swapxy)
					  {
						  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  if(normal==2)
				  {
					  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*grid_up_order2[foo][aux_ind];
					  mov_pdb_z2 = 10*grid_down_order2[foo][aux_ind];
					  if(swapxy)
					  {
						  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
						  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
						  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
					  }
				  }
				  fprintf(order_fp_avg_pdb_sn2[foo],pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,order_grid_up_sn2[foo][aux_ind][0]);
				  fprintf(order_fp_avg_pdb_sn2[foo],pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,order_grid_down_sn2[foo][aux_ind][0]);

				  fprintf(order_up_fp_avg_dat_sn2[foo],"%f	",order_grid_up_sn2[foo][aux_ind][0]);
				  //fprintf(order_down_fp_avg_dat_sn2[foo],"%f	",order_grid_down_sn2[foo][aux_ind][0]);
				  mat_low_order2[foo][binx-1-i] = order_grid_down_sn2[foo][aux_ind][0];

				  //AVG over all lipids as in g_order
				  if(first_bin==1)
				  {
				  	real avg = order_sum2[foo]/(frame_num*lipid_num);
				  	real sd = order_sum2_sd[foo]/(frame_num*lipid_num) - pow(avg,2);
				  	if(sd<0.0) //could happen due to numerical precision
				  	{
				  		sd = 0.0;
				  	}
				  	else
				  	{
				  		sd = sqrt(sd);
				  	}
				  	fprintf(order_fp_AVG_sn2,"%d	%f	%f\n",foo+2,avg,sd);
				  }
			  }
		  }
		  /***********************************ORDER PARAM*******************************************/


		  /*****************************************Curvature*******************************************/
		  if(curve)
		  {
			  if(normal==0)
			  {
				  mov_pdb_x1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_x2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==1)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_y2 = 10*grid_down_avg[aux_ind];
				  mov_pdb_z1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z2 = 10*left_y+10*j*bin_sizey;
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_z1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_z2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  if(normal==2)
			  {
				  mov_pdb_x1 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_x2 = 10*left_x+10*i*bin_sizex;
				  mov_pdb_y1 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_y2 = 10*left_y+10*j*bin_sizey;
				  mov_pdb_z1 = 10*grid_up_avg[aux_ind];
				  mov_pdb_z2 = 10*grid_down_avg[aux_ind];
				  if(swapxy)
				  {
					  mov_pdb_x1 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_x2 = 10*left_y+10*j*bin_sizey;
					  mov_pdb_y1 = 10*left_x+10*i*bin_sizex;
					  mov_pdb_y2 = 10*left_x+10*i*bin_sizex;
				  }
			  }
			  fprintf(gcurve_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,gcurve_scale*gausCurveUp[aux_ind]);
			  fprintf(gcurve_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2,1.0,gcurve_scale*gausCurveDown[aux_ind]);

			  fprintf(mcurve_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chA,1, mov_pdb_x1, mov_pdb_y1, mov_pdb_z1, 1.0,mcurve_scale*meanCurveUp[aux_ind]);
			  fprintf(mcurve_fp_avg_pdb,pdbform,"ATOM",(aux_ind+1)%100000,"C","XXX",chB,1, mov_pdb_x2, mov_pdb_y2, mov_pdb_z2, 1.0,mcurve_scale*meanCurveDown[aux_ind]);

			  fprintf(gcurve_up_fp_avg_dat,"%10.10f	",gausCurveUp[aux_ind]);
			  fprintf(mcurve_up_fp_avg_dat,"%10.10f	",meanCurveUp[aux_ind]);
			  mat_low_mcurve[binx-1-i] = meanCurveDown[aux_ind];
			  mat_low_gcurve[binx-1-i] = gausCurveDown[aux_ind];
		  }
		  /*****************************************Curvature*******************************************/
	  }

	  if(thick)
	  {
		  fprintf(thick_fp_avg_dat,"\n");
		  fprintf(thick_fp_sd_dat,"\n");
	  }

	  if(apl)
	  {
		  for(low_i=0; low_i<binx; low_i++)
		  {
			  fprintf(apl_down_fp_avg_dat,"%f	",mat_low_apl_avg[low_i]);
			  fprintf(apl_down_fp_sd_dat,"%f	",mat_low_apl_sd[low_i]);
		  }
		  fprintf(apl_up_fp_avg_dat,"\n");
		  fprintf(apl_up_fp_sd_dat,"\n");
		  fprintf(apl_down_fp_avg_dat,"\n");
		  fprintf(apl_down_fp_sd_dat,"\n");
	  }

	  if(order)
	  {
		  for(foo=0; foo<order_atom_num1-2; foo++)
		  {
			  for(low_i=0; low_i<binx; low_i++)
			  {
				  fprintf(order_down_fp_avg_dat_sn1[foo],"%f	",mat_low_order1[foo][low_i]);
			  }
			  fprintf(order_up_fp_avg_dat_sn1[foo],"\n");
			  fprintf(order_down_fp_avg_dat_sn1[foo],"\n");
		  }
		  for(foo=0; foo<order_atom_num2-2; foo++)
		  {
			  for(low_i=0; low_i<binx; low_i++)
			  {
				  fprintf(order_down_fp_avg_dat_sn2[foo],"%f	",mat_low_order2[foo][low_i]);
			  }
			  fprintf(order_up_fp_avg_dat_sn2[foo],"\n");
			  fprintf(order_down_fp_avg_dat_sn2[foo],"\n");
		  }
	  }

	  if(curve)
	  {
		  for(low_i=0; low_i<binx; low_i++)
		  {
			  fprintf(gcurve_down_fp_avg_dat,"%10.10f	",mat_low_gcurve[low_i]);
			  fprintf(mcurve_down_fp_avg_dat,"%10.10f	",mat_low_mcurve[low_i]);
		  }
		  fprintf(gcurve_up_fp_avg_dat,"\n");
		  fprintf(gcurve_down_fp_avg_dat,"\n");
		  fprintf(mcurve_up_fp_avg_dat,"\n");
		  fprintf(mcurve_down_fp_avg_dat,"\n");
	  }
  }


  /********** SPECIAL THICKNESS OUTPUT: LIPID INDECES ********/
  /**************** free Thickness ************/
  if(thick)
  {
	  fprintf(thick_fp_lipids_up,"#lipid_ID	mean(THICKNESS)	stdev(THICKNESS)\n");
	  fprintf(thick_fp_lipids_down,"#lipid_ID	mean(THICKNESS)	stdev(THICKNESS)\n");
	  int isprotnlip=0;
	  if(is_prot)
	  {
		  isprotnlip = lipid_num+1;
	  }
	  else
	  {
		  isprotnlip = lipid_num;
	  }
	  for(i=0;i<isprotnlip;i++)
	  {
		  real avg_up = thick_lip_up[i][2]/frame_num;
		  real avg_down = thick_lip_down[i][2]/frame_num;
		  real sd_up=0.0;
		  real sd_down=0.0;
		  if(frame_num>1)
		  {
			  sd_up = thick_lip_up[i][3]/(frame_num-1)-pow(avg_up,2)*frame_num/(frame_num-1);
			  sd_down = thick_lip_down[i][3]/(frame_num-1)-pow(avg_down,2)*frame_num/(frame_num-1);
			  if(sd_up<0.0) //could happen due to numerical precision
			  { sd_up = 0.0; }
			  else
			  { sd_up = sqrt(sd_up); }
			  if(sd_down<0.0) //could happen due to numerical precision
			  { sd_down = 0.0; }
			  else
			  { sd_down = sqrt(sd_down); }
		  }
		  int typecast_id = (int) thick_lip_up[i][0];
		  if(avg_up>0.0)
		  {
			  if(typecast_id == -1)
			  { fprintf(thick_fp_lipids_up,"PROTEIN	%f	%f\n",avg_up,sd_up); }
			  else
			  { fprintf(thick_fp_lipids_up,"%d	%f	%f\n",typecast_id,avg_up,sd_up); }
		  }
		  if(avg_down>0.0)
		  {
			  if(typecast_id == -1)
			  { fprintf(thick_fp_lipids_down,"PROTEIN	%f	%f\n",avg_down,sd_down); }
			  else
			  { fprintf(thick_fp_lipids_down,"%d	%f	%f\n",typecast_id,avg_down,sd_down); }
		  }
	  }

	  fclose(thick_fp_lipids_up);
	  fclose(thick_fp_lipids_down);
	  fclose(thick_fp_over_time);

	  fprintf(thick_fp_avg_pdb,"END");
  	  fprintf(thick_fp_sd_pdb,"END");
  	  fclose(thick_fp_avg_pdb);
  	  fclose(thick_fp_sd_pdb);
  	  fclose(thick_fp_avg_dat);
  	  fclose(thick_fp_sd_dat);
  	  sfree(grid);
  	  sfree(grid_sd);
  	  if(mat || pdb)
  	  {
  		  for(foo=0; foo<smooth; foo++)
  		  {
  			  sfree(grid_smooth_frames[foo]);
  		  }
  		  sfree(grid_smooth_frames);
  		  sfree(grid_smooth_avg);
  		  if(mat)
  		  { fclose(fp_mov_mat_thick); }
  	  }
  }
  sfree(grid_up_avg);
  sfree(grid_down_avg);


  /********** SPECIAL APL OUTPUT: LIPID INDECES ********/
  /**************** free APL ************/
  if(apl)
  {
	  fprintf(apl_fp_lipids_up,"#lipid_ID	mean(APL)	stdev(APL)\n");
	  fprintf(apl_fp_lipids_down,"#lipid_ID	mean(APL)	stdev(APL)\n");
	  int isprotnlip=0;
	  if(is_prot)
	  {
		  isprotnlip = lipid_num+1;
	  }
	  else
	  {
		  isprotnlip = lipid_num;
	  }
	  for(i=0;i<isprotnlip;i++)
	  {
		  real avg_up = apl_lip_up[i][2]/frame_num;
		  real avg_down = apl_lip_down[i][2]/frame_num;
		  real sd_up=0.0;
		  real sd_down=0.0;
		  if(frame_num>1)
		  {
			  sd_up = apl_lip_up[i][3]/(frame_num-1)-pow(avg_up,2)*frame_num/(frame_num-1);
			  sd_down = apl_lip_down[i][3]/(frame_num-1)-pow(avg_down,2)*frame_num/(frame_num-1);
			  if(sd_up<0.0) //could happen due to numerical precision
			  { sd_up = 0.0; }
			  else
			  { sd_up = sqrt(sd_up); }
			  if(sd_down<0.0) //could happen due to numerical precision
			  { sd_down = 0.0; }
			  else
			  { sd_down = sqrt(sd_down); }
		  }
		  int typecast_id = (int) apl_lip_up[i][0];
		  if(avg_up>0.0)
		  {
			  if(typecast_id == -1)
			  { fprintf(apl_fp_lipids_up,"PROTEIN	%f	%f\n",avg_up,sd_up); }
			  else
			  { fprintf(apl_fp_lipids_up,"%d	%f	%f\n",typecast_id,avg_up,sd_up); }
		  }
		  if(avg_down>0.0)
		  {
			  if(typecast_id == -1)
			  { fprintf(apl_fp_lipids_down,"PROTEIN	%f	%f\n",avg_down,sd_down); }
			  else
			  { fprintf(apl_fp_lipids_down,"%d	%f	%f\n",typecast_id,avg_down,sd_down); }
		  }
	  }

	  fprintf(apl_fp_avg_pdb,"END");
	  fprintf(apl_fp_sd_pdb,"END");
	  fclose(apl_fp_avg_pdb);
	  fclose(apl_fp_sd_pdb);
	  fclose(apl_fp_lipids_up);
	  fclose(apl_fp_lipids_down);
	  fclose(apl_fp_over_time);

	  for(i=0;i<grid_size;i++)
	  {
		  sfree(apl_grid_up[i]);
		  sfree(apl_grid_down[i]);
	  }
	  if(is_prot)
	  {
		  for(i=0;i<lipid_num+1;i++)
		  {
			  sfree(apl_lip_up[i]);
			  sfree(apl_lip_down[i]);
		  }
	  }
	  else
	  {
		  for(i=0;i<lipid_num;i++)
		  {
			  sfree(apl_lip_up[i]);
			  sfree(apl_lip_down[i]);
		  }
	  }
	  sfree(apl_lip_up);
	  sfree(apl_lip_down);
	  sfree(apl_grid_up);
	  sfree(apl_grid_down);
	  sfree(mat_low_apl_avg);
	  sfree(mat_low_apl_sd);

	  sfree(order_sum1);
	  sfree(order_sum2);
	  sfree(order_sum1_sd);
	  sfree(order_sum2_sd);

  	  if(mat)
  	  {
  		  for(foo=0; foo<smooth; foo++)
  		  {
  			  sfree(apl_smooth_up_frames[foo]);
  			  sfree(apl_smooth_down_frames[foo]);
  		  }
		  sfree(apl_smooth_up_frames);
		  sfree(apl_smooth_down_frames);
  		  sfree(apl_smooth_up_avg);
  		  sfree(apl_smooth_down_avg);
  		  sfree(apl_smooth_down_avg_Xinv);
  		  fclose(fp_mov_mat_apl_up);
  		  fclose(fp_mov_mat_apl_down);
  	  }
  	  fclose(apl_up_fp_avg_dat);
  	  fclose(apl_up_fp_sd_dat);
  	  fclose(apl_down_fp_avg_dat);
  	  fclose(apl_down_fp_sd_dat);
  }


  /**************** free Curvature ************/
  if(curve)
  {
	  sfree(gausCurveUp);
	  sfree(gausCurveDown);
	  sfree(meanCurveUp);
	  sfree(meanCurveDown);
	  sfree(mat_low_gcurve);
	  sfree(mat_low_mcurve);
	  fprintf(gcurve_fp_avg_pdb,"END");
	  fprintf(mcurve_fp_avg_pdb,"END");

  	  fclose(mcurve_fp_avg_pdb);
  	  fclose(gcurve_fp_avg_pdb);
  	  fclose(mcurve_up_fp_avg_dat);
  	  fclose(mcurve_down_fp_avg_dat);
  	  fclose(gcurve_up_fp_avg_dat);
  	  fclose(gcurve_down_fp_avg_dat);

	  sfree(mcurve_grid_up);
	  sfree(mcurve_grid_down);
	  sfree(gcurve_grid_up);
	  sfree(gcurve_grid_down);

      if(q_filter_low>0.0 || q_filter_high<99999.99) //absolute radius
      {
                 sfree(filtered_up);
                 sfree(filtered_down);
	  }
      else if(r_filter_low>0.0 || r_filter_high<1.0) //relative radius
	  {
                 sfree(filtered_up);
                 sfree(filtered_down);
	  }

  	  if(mat)
  	  {
  		  fclose(fp_mov_mat_mcurve_up);
  		  fclose(fp_mov_mat_mcurve_down);
  		  fclose(fp_mov_mat_gcurve_up);
  		  fclose(fp_mov_mat_gcurve_down);
  	  }
  }


  /**************** free Order ************/
  if(order)
  {
	  fclose(order_fp_AVG_sn1);
	  fclose(order_fp_AVG_sn2);

//	  sfree(x_rm_pbc);
	  for(i=0;i<order_atom_num1-2;i++)
	  {
		  for(j=0;j<grid_size;j++)
		  {
			  sfree(order_grid_up_sn1[i][j]);
			  sfree(order_grid_down_sn1[i][j]);
		  }
		  sfree(order_grid_up_sn1[i]);
		  sfree(order_grid_down_sn1[i]);
		  sfree(grid_up_order1[i]);
		  sfree(grid_down_order1[i]);
		  sfree(mat_low_order1[i]);
		  sfree(order_count_sn1_up[i]);
		  sfree(order_count_sn1_down[i]);
	  }
	  for(i=0;i<order_atom_num2-2;i++)
	  {
		  for(j=0;j<grid_size;j++)
		  {
			  sfree(order_grid_up_sn2[i][j]);
			  sfree(order_grid_down_sn2[i][j]);
		  }
		  sfree(order_grid_up_sn2[i]);
		  sfree(order_grid_down_sn2[i]);
		  sfree(grid_up_order2[i]);
		  sfree(grid_down_order2[i]);
		  sfree(mat_low_order2[i]);
		  sfree(order_count_sn2_up[i]);
		  sfree(order_count_sn2_down[i]);
	  }
	  sfree(order_grid_up_sn1);
	  sfree(order_grid_down_sn1);
	  sfree(order_grid_up_sn2);
	  sfree(order_grid_down_sn2);
	  sfree(grid_up_order1);
	  sfree(grid_down_order1);
	  sfree(grid_up_order2);
	  sfree(grid_down_order2);
	  sfree(mat_low_order1);
	  sfree(mat_low_order2);
	  sfree(order_count_sn1_up);
	  sfree(order_count_sn1_down);
	  sfree(order_count_sn2_up);
	  sfree(order_count_sn2_down);

	  for(i=0;i<lipid_num;i++)
	  {
		  sfree(order_lip1[i]);
		  sfree(order_lip2[i]);
	  }
	  sfree(order_lip1);
	  sfree(order_lip2);

  	  //close files and free pointers
	  //sn1
	  for(i=2; i<order_atom_num1; i++)
	  {
		  fclose(order_up_fp_avg_dat_sn1[i-2]);
		  fclose(order_down_fp_avg_dat_sn1[i-2]);
		  fprintf(order_fp_avg_pdb_sn1[i-2],"END");
		  fclose(order_fp_avg_pdb_sn1[i-2]);
	  }
	  //sn2
	  for(i=2; i<order_atom_num2; i++)
	  {
		  fclose(order_up_fp_avg_dat_sn2[i-2]);
		  fclose(order_down_fp_avg_dat_sn2[i-2]);
		  fprintf(order_fp_avg_pdb_sn2[i-2],"END");
		  fclose(order_fp_avg_pdb_sn2[i-2]);
	  }
	  sfree(order_up_fp_avg_dat_sn1);
	  sfree(order_down_fp_avg_dat_sn1);
	  sfree(order_up_fp_avg_dat_sn2);
	  sfree(order_down_fp_avg_dat_sn2);
	  sfree(order_fp_avg_pdb_sn1);
	  sfree(order_fp_avg_pdb_sn2);

	  //protein indeces
	  if(is_prot)
	  {
		  for(i=0;i<order_atom_num1-2;i++)
		  {
			  sfree(ptop_ind_order1[i]);
			  sfree(pbot_ind_order1[i]);
		  }
		  for(i=0;i<order_atom_num2-2;i++)
		  {
			  sfree(ptop_ind_order2[i]);
			  sfree(pbot_ind_order2[i]);
		  }
		  sfree(nprot_top_order1);
		  sfree(nprot_bot_order1);
		  sfree(nprot_top_order2);
		  sfree(nprot_bot_order2);
		  sfree(ptop_ind_order1);
		  sfree(pbot_ind_order1);
		  sfree(ptop_ind_order2);
		  sfree(pbot_ind_order2);
	  }
	  //mov_mat
	  if(mat)
	  {
		  for(i=0;i<order_atom_num1-2;i++)
		  {
			  for(j=0;j<smooth;j++)
			  {
				  sfree(order_smooth_up_frames1[i][j]);
				  sfree(order_smooth_down_frames1[i][j]);
			  }
			  sfree(order_smooth_up_frames1[i]);
			  sfree(order_smooth_down_frames1[i]);
			  fclose(fp_mov_mat_order_up1[i]);
			  fclose(fp_mov_mat_order_down1[i]);
			  sfree(order_smooth_down1_Xinv[i]);
		  }
		  for(i=0;i<order_atom_num2-2;i++)
		  {
			  for(j=0;j<smooth;j++)
			  {
				  sfree(order_smooth_up_frames2[i][j]);
				  sfree(order_smooth_down_frames2[i][j]);
			  }
			  sfree(order_smooth_up_frames2[i]);
			  sfree(order_smooth_down_frames2[i]);
			  fclose(fp_mov_mat_order_up2[i]);
			  fclose(fp_mov_mat_order_down2[i]);
			  sfree(order_smooth_down2_Xinv[i]);
		  }
		  sfree(order_smooth_up_frames1);
		  sfree(order_smooth_down_frames1);
		  sfree(order_smooth_up_frames2);
		  sfree(order_smooth_down_frames2);
		  sfree(fp_mov_mat_order_up1); sfree(fp_mov_mat_order_down1);
		  sfree(fp_mov_mat_order_up2); sfree(fp_mov_mat_order_down2);
		  sfree(order_smooth_down1_Xinv);
		  sfree(order_smooth_down2_Xinv);
	  }
  }


  /**************** free Diffusion ************/
  if(diffus)
  {
	  fclose(diffus_fp_up_dat);
	  fclose(diffus_fp_down_dat);
	  fclose(diffus_fp_pdb_avg);

	  for(i=0; i<grid_size; i++)
	  {
		  sfree(diffus_grid_pos_up[i]);
		  sfree(diffus_grid_pos_down[i]);
		  sfree(diffus_grid_dist_up[i]);
		  sfree(diffus_grid_dist_down[i]);
	  }

	  sfree(diffus_grid_pos_up);
	  sfree(diffus_grid_pos_down);
	  sfree(diffus_grid_dist_up);
	  sfree(diffus_grid_dist_down);

	  sfree(diffus_grid_offset_up);
	  sfree(diffus_grid_offset_down);
  }


  /**************** free pdb movie ************/
  if( pdb || (curve && mat) )
  {
	  for(foo=0; foo<smooth; foo++)
	  {
		  sfree(z_smooth_frames_up[foo]);
		  sfree(z_smooth_frames_down[foo]);
	  }
	  sfree(z_smooth_avg_up);
	  sfree(z_smooth_avg_down);
	  sfree(z_smooth_frames_up);
	  sfree(z_smooth_frames_down);
  }

  thanx(stderr);

  return 0;
}
