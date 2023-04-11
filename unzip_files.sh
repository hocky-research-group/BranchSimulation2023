#!/bin/bash

dirlist=$(echo full_junction morph)

for i in $dirlist;
do
	echo unzip -d ./$i/ ./$i/*.zip
done
