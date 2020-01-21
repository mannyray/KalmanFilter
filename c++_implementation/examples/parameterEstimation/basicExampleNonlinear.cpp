#include <iostream>
#include <sstream>
#include "../../include/Eigen.h"
#include "../../include/KalmanFilter.h"
#include "../../include/filterModel.h"
#include "../../include/model.h"
#include <fstream>
#include <math.h>
#include <iomanip>     
#include<stdio.h>
#include<stdlib.h>


class processNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd Q;
	Eigen::MatrixXd Q_root;
	public:
		processNoise(double v){
			Eigen::MatrixXd tmp1(1,1);
			tmp1<<v;
			Q = tmp1;
			Eigen::MatrixXd tmp2(1,1);
			tmp2<<std::sqrt(v);
			Q_root = tmp2;
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
			return Q;
		}
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return  Q_root; 
		}
};

class sensorNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd R;
	Eigen::MatrixXd R_root;
	public:
		sensorNoise(double v){
			Eigen::MatrixXd tmp1(1,1);
			tmp1<<v;
			R = tmp1;
			Eigen::MatrixXd tmp2(1,1);
			tmp2<<std::sqrt(v);
			R_root = tmp2;
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
			return R;
		}
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return R_root; 
		}
};

class transitionJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
	double rate;
	double max_pop;
	Eigen::VectorXd add;
	public:
		transitionJac(double rate, double max_pop):rate(rate),max_pop(max_pop){
			Eigen::VectorXd tmp(1);
			tmp<<1+rate;
			add = tmp;
		}
		Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
			return add - (2*rate/max_pop)*val;
		}
};

class measurementJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
	
	public:
		Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
			return Eigen::MatrixXd::Identity(1,1);
		}
};

class stateModel: public discreteModel<Eigen::VectorXd>{
	double rate;
	double max_pop;
	Eigen::VectorXd one;
	public:
		stateModel(double rate, double max_pop):rate(rate),max_pop(max_pop){
			Eigen::VectorXd tmp(1);tmp<<1;
			one = tmp;
		}
		Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
			return (val + rate*val*(one-val/max_pop));
		}
};


class measurementModel: public discreteModel<Eigen::VectorXd>{
	public:
		measurementModel(){
		}
		Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
			return val;
		}
};

using namespace std;

int main(int argc, char**argv){

	if(argc != 6){
		cout<<"Incorrect arguments."<<endl; 
		return -1;
	}

	string fileName;
	double actualProcessNoise = 0;
	double actualSensorNoise = 0;
	double assumedProcessNoise = 0;
	double assumedSensorNoise = 0;
	string filePrefix = "";
	try{
		actualProcessNoise = strtod(argv[1],NULL);
		actualSensorNoise = strtod(argv[2],NULL);
		assumedProcessNoise = strtod(argv[3],NULL);
		assumedSensorNoise = strtod(argv[4],NULL);
		filePrefix = argv[5];
	}catch(...){
		cout<<"Incorrect arguments."<<endl; 
	}

	double rate = 0.01;
	double max_pop = 100;

	stateModel tm(rate,max_pop);
	measurementModel mm;
	processNoise pnActual(actualProcessNoise);
	sensorNoise snActual(actualSensorNoise);
	processNoise pnAssumed(assumedProcessNoise);
	sensorNoise snAssumed(assumedSensorNoise);
	transitionJac tj(rate,max_pop);

	measurementJac mj;	
	int stateCount = 1;
	int sensorCount = 1;
	discreteDiscreteFilterModel<Eigen::VectorXd,Eigen::MatrixXd>ddfmActual(&tm,&mm,&pnActual,&snActual,&tj,&mj,stateCount,sensorCount);
	discreteDiscreteFilterModel<Eigen::VectorXd,Eigen::MatrixXd>ddfmAssumed(&tm,&mm,&pnAssumed,&snAssumed,&tj,&mj,stateCount,sensorCount);

	int initialTime = 0;
	Eigen::VectorXd initialState(1);
	initialState<<(max_pop/2.0);
	discreteDiscreteFilterSolver<Eigen::VectorXd,Eigen::MatrixXd> ddfs(&ddfmActual,initialTime,initialState);
	
	int steps = 2000;
	ddfs.setSeed(assumedProcessNoise*assumedSensorNoise*actualProcessNoise*actualSensorNoise);
	ddfs.solve(steps);

	
	Eigen::VectorXd * states = ddfs.getSolvedStates();
	Eigen::VectorXd * measurements = ddfs.getSolvedMeasurements();

	
	Eigen::VectorXd initialEstimate(1);
	initialEstimate<<(max_pop/2.0);
	Eigen::MatrixXd initialCovariance = Eigen::MatrixXd::Identity(1,1);

	discreteDiscreteKalmanFilter<Eigen::VectorXd,Eigen::MatrixXd> filter(initialTime, initialEstimate, initialCovariance,&ddfmAssumed);


	Eigen::VectorXd* estimatesBefore = new Eigen::VectorXd[steps];
	Eigen::VectorXd* estimates = new Eigen::VectorXd[steps];
	Eigen::MatrixXd* covariances = new Eigen::MatrixXd[steps];
	for(int i = 0; i < steps; i++){
		covariances[i] = filter.getCurrentCovariance();	
		filter.predict(1);
		estimatesBefore[i] = filter.getCurrentEstimate();
		filter.update(measurements[i]);
		estimates[i] = filter.getCurrentEstimate();
	}

	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
	std::ofstream outputfileEstimates;
	stringstream ss1; ss1<<filePrefix<<"estimates.txt";
	outputfileEstimates.open(ss1.str());
	for(int i = 0; i < steps; i++){
		outputfileEstimates<<estimates[i].format(HeavyFmt)<<std::endl;
	}outputfileEstimates.close();

	std::ofstream outputfileEstimatesBefore;
	stringstream ss1_5; ss1_5<<filePrefix<<"estimatesBefore.txt";
	outputfileEstimatesBefore.open(ss1_5.str());
	for(int i = 0; i < steps; i++){
		outputfileEstimatesBefore<<estimatesBefore[i].format(HeavyFmt)<<std::endl;
	}outputfileEstimatesBefore.close();

	std::ofstream outputfileMeasurements;
	stringstream ss2; ss2<<filePrefix<<"measurements.txt";
	outputfileMeasurements.open(ss2.str());
	for(int i = 0; i < steps; i++){
		outputfileMeasurements<<measurements[i].format(HeavyFmt)<<std::endl;
	}outputfileMeasurements.close();

	std::ofstream outputfileStates;
	stringstream ss3; ss3<<filePrefix<<"states.txt";
	outputfileStates.open(ss3.str());
	for(int i = 0; i < steps; i++){
		outputfileStates<<states[i].format(HeavyFmt)<<std::endl;
	}outputfileStates.close();
	
	std::ofstream outputfileCovariances;
	stringstream ss4; ss4<<filePrefix<<"covariances.txt";
	outputfileCovariances.open(ss4.str());
	for(int i = 0; i < steps; i++){
		outputfileCovariances<<covariances[i].format(HeavyFmt)<<std::endl;
	}outputfileCovariances.close();

	delete [] estimates;
	delete [] states;
	delete [] measurements;
	delete [] covariances;
}
