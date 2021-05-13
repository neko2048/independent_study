# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./th*.png #./vapor*.png
rm ./th/* #./vapor/*
mkdir th #vapor

# excute py file
echo "$BASENAME" | python -W ignore thvapor.py

# create video
rm th.mkv
#rm vapor.mkv
echo "create video..."
ffmpeg -framerate 5 -i th%06d.png th.mp4
#ffmpeg -framerate 5 -i vapor%06d.png vapor.mkv

# move png to folder
mv th*.png ./th/
#mv vapor*.png ./vapor/
