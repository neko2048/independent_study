# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`
echo "$BASENAME" | python -W ignore th_plume.py
