rm th*.png
rm vapor*.png
python -W ignore draw.py
ffmpeg -framerate 10 -i th%06d.png walter_th.mkv
ffmpeg -framerate 10 -i vapor%06d.png walter_vapor.mkv
mv th*.png th/
mv vapor*.png vapor/
