# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./core/* ./cloud/* ./
rm ./core_img/* ./cloud_img/*
#rm ./buoyancy_img/* ./buoyancy*.jpg
mkdir core cloud # save data
mkdir core_img cloud_img buoyancy_img # save image

# execute py file
echo "$BASENAME" | python -W ignore CC_buoyancy.py

# crate video
rm core.mkv cloud.mkv buoyancy.mp4
ffmpeg -framerate 5 -i core%06d.png core.mp4
ffmpeg -framerate 5 -i cloud%06d.png cloud.mp4
ffmpeg -framerate 5 -pattern_type glob -i "buoyancy*.jpg" -qscale 0 buoyancy.mp4
# move png to folder
mv core*.png core_img/
mv cloud*.png cloud_img/
mv buoyancy*.jpg buoyancy_img

# create time series figure
#echo "$BASENAME" | python CC_time.py
