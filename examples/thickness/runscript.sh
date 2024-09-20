
#export LOMEPRO=path_to_lomepro

mkdir thickness
cd thickness
echo 9 | $LOMEPRO/g_lomepro_static -s ../init.pdb -f ../trj_10fr.xtc -n ../index.ndx -thick -lip_num 124 -binx 100 -biny 100 
cd ../



