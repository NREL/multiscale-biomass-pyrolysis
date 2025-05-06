cp -r init 0
blockMesh
setFields
decomposePar
mpirun -np 16 pyroFoam -parallel
