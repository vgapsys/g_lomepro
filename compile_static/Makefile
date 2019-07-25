# 1. path to gromacs/src/tools (source folder)
SRCTOOLS     = /path_where_gromacs_was_downloaded/gromacs-4.6.2/gromacs/src/tools
# 2. path to gromacs/include/gromacs (folder where gromacs is compiled)
INCLUDE      = /path_where_gromacs_was_compiled/gromacs-462/include/gromacs
# 3. path to gromacs/lib64 (folder where gromacs is compiled)
GMXLIB       = /path_where_gromacs_was_compiled/gromacs-462/lib
# 4. path to fftw library
FFTW         = /path_where_fftw_was_compiled/fftw/compiled/lib64
# 5. path to fftw header files
FFTWINCL     = /path_where_fftw_was_compiled/fftw/compiled/include

# Variables set by the configuration script:
LIBS         = -lgmx -lm -lfftw3 -lfftw3f 
LDFLAGS      = -L$(GMXLIB) -L$(FFTW)
CFLAGS       = -I$(SRCTOOLS) -I$(INCLUDE) -I$(FFTWINCL)

CC           = gcc
LD           = $(CC)

g_lomepro_v1.0.2: g_lomepro.c
	$(LD) -fopenmp $(CFLAGS) $(LDFLAGS) g_lomepro.c get_index.c grid.c protein_atoms.c order.c curvature.c filtering.c dens.c -o $@ $(LIBS)

