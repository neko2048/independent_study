# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./th_dwdt*.jpg
rm ./th_dwdt/*.jpg
mkdir th_dwdt

# excute py file
echo "$BASENAME" | python -W ignore th_prof.py

# create video
rm th_dwdt.mkv
echo "create video..."
ffmpeg -framerate 5 -i th_dwdt%06d.jpg th_dwdt.mp4

# move jpg to folder
mv th_dwdt*.jpg ./th_dwdt/
