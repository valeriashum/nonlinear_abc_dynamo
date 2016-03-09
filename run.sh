make clean
if [ -z "$1" ]
    then
    ./configure --with-problem=kinematic_ABC_flow_box_1 --enable-mpi --enable-debug 
fi
make -j 
