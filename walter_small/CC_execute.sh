rm core*.png
rm cloud*.png
python bouyancy.py
ffmpeg -framerate 10 -i core%06d.png walter_core.mkv
ffmpeg -framerate 10 -i cloud%06d.png walter_cloud.mkv
mv core*.png core_img/
mv cloud*.png cloud_img/

