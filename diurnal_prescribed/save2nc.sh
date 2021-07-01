# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./nc_files/*.nc
rm buoyancy_npy/*.npy, dwdt_npy/*.npy
mkdir nc_files
mkdir buoyancy_npy dwdt_npy

# excute py file
echo "$BASENAME" | python save2nc.py
