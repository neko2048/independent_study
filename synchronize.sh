#!bin/bash
echo 'type the folder name u want to copy file to from template:'
read name

echo 'synchronizing...'
cp ./template/*.sh ./$name/
cp ./template/*.py ./$name/

echo 'end'
