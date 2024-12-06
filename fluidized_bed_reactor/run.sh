module load PrgEnv-cray
source /kfs2/projects/biocon/ofoam_cray_mpich/OpenFOAM-dev/etc/bashrc
./Allclean
blockMesh -dict system/bmesh2d
rm -rf 0
cp -r 0.orig 0
setFields
decomposePar -fileHandler collated
touch soln.foam
srun -n 4 multiphaseEulerFoam -parallel -fileHandler collated
