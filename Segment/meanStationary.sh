#!/bin/bash

mu1='median.png'		

rm -fr meanStationary
mkdir -p meanStationary
fldrs="2_Checkerboard_4_2012_11_15_13_49
4_Checkerboard_7_2012_11_19_12_11
4_Uniform_2012_11_15_8_25"

for fldr in $fldrs
do
convert data/$fldr/*.jpg -set colorspace RGB -strip -colorspace Gray -evaluate-sequence Median -auto-level meanStationary/$fldr.jpg
done