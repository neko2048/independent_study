# retrieve currrent dir name
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

# check folder exists or not
rm Lbuoyancy*.jpg LW.jpg
rm ./LB_img/* ./LW_img/*
rm Core_L*.jpg Cloud_L*.jpg 

mkdir LB_img LW_img
# execute py file
echo "$BASENAME" | python laplace.py

# crate video
rm LB.mkv LW.mkv
rm CoreLB*.mkv CloudLB*.mkv
rm CoreLW*.mkv CoreLW*.mkv
ffmpeg -framerate 5 -i Lbuoyancy%06d.jpg LB.mkv
ffmpeg -framerate 5 -i LW%06d.jpg LW.mkv
ffmpeg -framerate 5 -i Cloud_LB%06d.jpg CloudLB.mkv
ffmpeg -framerate 5 -i Core_LB%06d.jpg CoreLB.mkv
ffmpeg -framerate 5 -i Cloud_LW%06d.jpg CloudLW.mkv
ffmpeg -framerate 5 -i Core_LW%06d.jpg CoreLW.mkv

# move jpg to folder
mv Lbuoyancy*.jpg LB_img/
mv LW*.jpg LW_img/
mv Cloud_LB*.jpg LB_img/
mv Core_LB*.jpg LB_img/
mv Cloud_LW*.jpg LW_img/
mv Core_LW*.jpg LW_img/
