#!/bin/bash

cd include
file=Eigen.h
echo "#ifndef CUSTOMEIGEN" > $file
echo "#define CUSTOMEIGEN"  >> $file
echo "#ifndef EIGEN_MATRIXBASE_PLUGIN" >> $file
echo "#define EIGEN_MATRIXBASE_PLUGIN \"$(pwd)/mathWrapper/eigenSimple.h\"" >> $file
echo "#endif" >> $file
echo "#include <random>" >> $file
echo "#include <Eigen/Dense>" >> $file
echo "#include <Eigen/Cholesky>" >> $file
echo "#endif" >> $file
