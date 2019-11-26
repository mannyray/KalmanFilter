#include <iostream>
#include <math.h>
#include "../include/mathWrapper/boost.h"
#include "../include/KalmanFilter.h"
#include "../include/filterModel.h"
#include "../include/model.h"
#include <fstream>
#include <string>
#include <iomanip>     
#include <vector>

class processNoise: public discreteNoiseCovariance<vectorBoost,matrixBoost>{
    boost::numeric::ublas::matrix<double> v;
    boost::numeric::ublas::matrix<double> sqrtMat;
    public:
        processNoise(){
            boost::numeric::ublas::matrix<double> tmp(2,2);
            tmp(0,0) = 8.47412478346666e-04;
            tmp(0,1) = 1.22566712734456e-03;
            tmp(1,0) = 1.22566712734457e-03;
            tmp(1,1) = 2.45572145988027e-03; 
            v = tmp;
            
            boost::numeric::ublas::matrix<double> tmp2(2,2);
            tmp2(0,0) = 0.0291103500210263;
            tmp2(0,1) = 0;
            tmp2(1,0) = 0.0421041700446496;
            tmp2(1,1) = 0.0261335096137412;
            sqrtMat = tmp2;
        }
        matrixBoost function(const vectorBoost &est, int t) override{
            return matrixBoost(v);
        }
        matrixBoost sqrt(const vectorBoost &est, int t) override{
            return  matrixBoost(sqrtMat); 
        }
};
class sensorNoise: public discreteNoiseCovariance<vectorBoost,matrixBoost>{
    boost::numeric::ublas::matrix<double> v;
    boost::numeric::ublas::matrix<double> sqrtMat;
    public:
        sensorNoise(){
            boost::numeric::ublas::matrix<double> tmp(1,1);
            tmp(0,0) = 0.25;
            v = tmp;
            boost::numeric::ublas::matrix<double> tmp2(1,1);
            tmp2(0,0) = 0.5;
            sqrtMat = tmp2; 
        }
        matrixBoost function(const vectorBoost &est, int t) override{
            return v;
        }
        matrixBoost sqrt(const vectorBoost &est, int t) override{
            return sqrtMat; 
        }
};

class transitionJac: public jacobianDiscrete<vectorBoost,matrixBoost>{
    boost::numeric::ublas::matrix<double> v;
    public:
        transitionJac(){
            boost::numeric::ublas::matrix<double> tmp(2,2);
            tmp(0,0) = 0.99523153695480593;
            tmp(0,1) = 0.93095517465814814;
            tmp(1,0) = -0.00930955174658148;
            tmp(1,1) = 0.86357459355850885;
            v = tmp;
        }
       	matrixBoost function(const vectorBoost & val, int t){
            return matrixBoost(v);
        }
};

class measurementJac: public jacobianDiscrete<vectorBoost,matrixBoost>{
    boost::numeric::ublas::matrix<double> v;
    public:
        measurementJac(){
            boost::numeric::ublas::matrix<double> tmp(1,2);
            tmp(0,0) = 1;
            tmp(0,1) = 0;
            v = tmp;
        }
        matrixBoost function(const vectorBoost & val, int t){
            return matrixBoost(v);
        }
};

class stateModel: public discreteModel<vectorBoost>{
    boost::numeric::ublas::matrix<double> v;
    public:
        stateModel(){
            boost::numeric::ublas::matrix<double> tmp(2,2);
            tmp(0,0) = 0.99523153695480593;
            tmp(0,1) = 0.93095517465814814;
            tmp(1,0) = -0.00930955174658148;
            tmp(1,1) = 0.86357459355850885;
            v = tmp;
        }
        vectorBoost function(const vectorBoost & val, const int time) const override{
            return vectorBoost(prod(v,val.getSystemValue()));
        }
};


class measurementModel: public discreteModel<vectorBoost>{
    boost::numeric::ublas::matrix<double> v;
    public:
        measurementModel(){
            boost::numeric::ublas::matrix<double> tmp(1,2);
            tmp(0,0) = 1.0;
            tmp(0,1) = 0.0;
            v = tmp;
        }
        vectorBoost function(const vectorBoost & val, const int time) const override{
            return vectorBoost(prod(v,val.getSystemValue()));
        }
};
