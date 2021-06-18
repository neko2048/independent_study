# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./nc_files/*.nc
mkdir nc_files

# excute py file
echo "$BASENAME" | python save2nc.py
