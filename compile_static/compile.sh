
cp ../* .

rm g_lomepro_static
gcc-8 -static -fopenmp -I/home/vgapsys/Downloads/gromacs-4.5.5/src/tools -I/home/vgapsys/Downloads/gromacs-4.5.5/include -I/home/vgapsys/gromacs-455/include/gromacs/ -L/home/vgapsys/gromacs-455/lib64 -L/home/vgapsys/lib/lapack_blas -L/home/vgapsys/Downloads/tirpc/tirpc/lib64 -L/home/vgapsys/Downloads/scipy/scipy-0.7.2/build/temp.linux-x86_64-2.6/ -L/home/vgapsys/scripts/compilation_files/static_lib -L/home/vgapsys/Downloads/fftw/compiled/lib64 -L/home/vgapsys/Downloads/fftw/compiled_f/lib64 -L/usr/lib64/gcc/x86_64-suse-linux/5  g_lomepro.c get_index.c grid.c protein_atoms.c order.c curvature.c filtering.c dens.c  -lgmxana  -lgmx -lm -lgomp -llapack -lblas  -ltirpc -lfftw3 -lfftw3f -ldl -o g_lomepro_static

#gcc-8 -static -fopenmp -I/home/vgapsys/Downloads/gromacs-4.5.5/src/tools -I/home/vgapsys/Downloads/gromacs-4.5.5/include -I/home/vgapsys/gromacs-455/include/gromacs/ -L/home/vgapsys/owl/gromacs-462_newsc/lib -L/home/vgapsys/lib/lapack_blas -L/home/vgapsys/Downloads/tirpc/tirpc/lib64 -L/home/vgapsys/Downloads/scipy/scipy-0.7.2/build/temp.linux-x86_64-2.6/ -L/home/vgapsys/scripts/compilation_files/static_lib -L/home/vgapsys/Downloads/fftw/compiled/lib64 -L/home/vgapsys/Downloads/fftw/compiled_f/lib64 -L/usr/lib64/gcc/x86_64-suse-linux/5  g_lomepro.c get_index.c grid.c protein_atoms.c order.c curvature.c filtering.c dens.c  -lgmxana  -lgmx -lm -lgomp -llapack -lblas  -ltirpc -lfftw3 -lfftw3f -ldl -o g_lomepro_static


cp g_lomepro_static ../.


