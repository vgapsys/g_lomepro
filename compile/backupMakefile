# 1. path to gromacs/src/tools (source folder)
SRCTOOLS     = /home/vgapsys/Downloads/gromacs-4.6.2/gromacs/src/tools
# 2. path to gromacs/include/gromacs (folder where gromacs is compiled)
INCLUDE      = /home/vgapsys/gromacs-462_clean/include/gromacs
# 3. path to gromacs/lib64 (folder where gromacs is compiled)
GMXLIB       = /home/vgapsys/gromacs-462_clean/lib
# 4. path to fftw library
FFTW         = /home/vgapsys/Downloads/fftw/compiled/lib64
# 5. path to fftw header files
FFTWINCL     = /home/vgapsys/Downloads/fftw/compiled/include

# Variables set by the configuration script:
LIBS         = -lgmx -lm -lfftw3 -lfftw3f 
LDFLAGS      = -L$(GMXLIB) -L$(FFTW)
CFLAGS       = -I$(SRCTOOLS) -I$(INCLUDE) -I$(FFTWINCL)

CC           = gcc
LD           = $(CC)

g_lomepro: g_lomepro.c
	$(LD) -fopenmp $(CFLAGS) $(LDFLAGS) g_lomepro.c get_index.c grid.c protein_atoms.c order.c curvature.c filtering.c dens.c -o $@ $(LIBS) -fPIC

