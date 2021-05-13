# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm ./core/* ./cloud/* ./
rm ./core_img/* ./cloud_img/*
mkdir core cloud # save data
mkdir core_img cloud_img buoyancy_img # save image

# execute py file
echo "$BASENAME" | python -W ignore CC_buoyancy.py

# crate video
rm core.mkv cloud.mkv
ffmpeg -framerate 5 -i core%06d.png core.mkv
ffmpeg -framerate 5 -i cloud%06d.png cloud.mkv
ffmpeg -framerate 5 -i buoyancy%06d.jpg buoyancy.mkv

# move png to folder
mv core*.png core_img/
mv cloud*.png cloud_img/
mv buoyancy.jpg buoyancy_img

# create time series figure
echo "$BASENAME" | python CC_time.py
