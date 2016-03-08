make clean
if [ -z "$1" ]
    then
    ./configure --with-problem=kinematic_ABC_flow --enable-mpi
fi
make 
