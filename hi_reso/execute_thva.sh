# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./th*.jpg ./vapor*.jpg
rm ./th/* ./vapor/*
mkdir th vapor

# excute py file
echo "$BASENAME" | python -W ignore thvapor.py

# create video
rm th.mkv
rm vapor.mkv
echo "create video..."
ffmpeg -framerate 5 -i th%06d.jpg th.mp4
ffmpeg -framerate 5 -i vapor%06d.jpg vapor.mp4

# move jpg to folder
mv th*.jpg ./th/
mv vapor*.jpg ./vapor/
