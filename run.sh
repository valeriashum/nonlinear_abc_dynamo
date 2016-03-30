make clean
if [ -z "$1" ]
    then
    ./configure --with-problem=u_iii_hydro_dynamo_box_10x2x1 --enable-mpi --enable-debug 
fi
make -j 
