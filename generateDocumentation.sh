#!/bin/bash

#script to generate documentation

rm -rf documentation/* 
rmdir documentation 
mkdir documentation
cd matlab_implementation
rm -rf Doc/*
rmdir Doc
doxygen Doxyfile
mv Doc/html ../documentation/matlab
mkdir ../documentation/matlab/assets
mkdir -p ../documentation/matlab/assets/continuous_discrete/examples
cp continuous_discrete/examples/*.png ../documentation/matlab/assets/continuous_discrete/examples # 2>/dev/null
mkdir -p ../documentation/matlab/assets/cubature/examples
cp cubature/examples/*.png ../documentation/matlab/assets/cubature/examples #2>/dev/null
mkdir -p ../documentation/matlab/assets/discrete_discrete/examples
cp discrete_discrete/examples/*.png ../documentation/matlab/assets/discrete_discrete/examples #2>/dev/null
mkdir -p ../documentation/matlab/assets/ensemble_stochastic/examples
cp ensemble_stochastic/examples/*.png ../documentation/matlab/assets/ensemble_stochastic/examples #2>/dev/null
mkdir -p ../documentation/matlab/assets/particle/examples
cp particle/examples/*.png ../documentation/matlab/assets/particle/examples #2>/dev/null
cd ..

cd c++_implementation
rm -rf Doc/*
rmdir Doc
doxygen Doxyfile
mv html ../documentation/cpp
cp assets/*.png ../documentation/cpp
cd ..
