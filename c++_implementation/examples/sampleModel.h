#include <iostream>
#include <math.h>
#include "../include/Eigen.h"
#include "../include/KalmanFilter.h"
#include "../include/filterModel.h"
#include "../include/model.h"
#include <fstream>
#include <string>
#include <vector>

class processNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
    Eigen::MatrixXd v;
    Eigen::MatrixXd sqrtMat;
    public:
        processNoise(){
            Eigen::MatrixXd tmp(2,2);
            tmp<< 8.47412478346666e-04,1.22566712734456e-03,1.22566712734457e-03,2.45572145988027e-03; 
            v = tmp;
            Eigen::LLT<Eigen::MatrixXd> lltOfV(v);
            sqrtMat = lltOfV.matrixL();
        }
        Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
            return v;
        }
        Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
            return  sqrtMat; 
        }
};

class sensorNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
    Eigen::MatrixXd v;
    Eigen::MatrixXd sqrtMat;
    public:
        sensorNoise(){
            Eigen::MatrixXd tmp(1,1);
            tmp<<0.25;
            v = tmp;
            Eigen::MatrixXd tmp2(1,1);
            tmp2<<0.5;
            sqrtMat = tmp2; 
        }
        Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
            return v;
        }
        Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
            return sqrtMat; 
        }
};

class transitionJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
    Eigen::MatrixXd v;
    public:
        transitionJac(){
            Eigen::MatrixXd tmp(2,2);
            tmp<<0.99523153695480593,0.93095517465814814,-0.00930955174658148,0.86357459355850885;
            v = tmp;
        }
        Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
            return v;
        }
};

class measurementJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
    Eigen::MatrixXd v;
    public:
        measurementJac(){
            Eigen::MatrixXd tmp(1,2);
            tmp<<1.0,0.0;
            v = tmp;
        }
        Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
            return v;
        }
};

class stateModel: public discreteModel<Eigen::VectorXd>{
    Eigen::MatrixXd v;
    public:
        stateModel(){
            Eigen::MatrixXd tmp(2,2);
            tmp<<0.99523153695480593,0.93095517465814814,-0.00930955174658148,0.86357459355850885;
            v = tmp;
        }
        Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
            return v*val;
        }
};

class measurementModel: public discreteModel<Eigen::VectorXd>{
    Eigen::MatrixXd v;
    public:
        measurementModel(){
            Eigen::MatrixXd tmp(1,2);tmp<<1.0,0.0;
            v = tmp;
        }
        Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
            return v*val;
        }
};
