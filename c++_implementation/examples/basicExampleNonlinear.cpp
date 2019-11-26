#include <iostream>
#include "../include/mathWrapper/double.h"
#include "../include/KalmanFilter.h"
#include "../include/filterModel.h"
#include "../include/model.h"
#include <fstream>
#include <math.h>
#include <iomanip>     

class processNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
	matrixDouble Q;
	public:
 	       processNoise(){
			Q = matrixDouble(1);
	        }
	matrixDouble function(const vectorDouble &est, int t) override{
		return Q;
	}
	matrixDouble sqrt(const vectorDouble &est, int t) override{
		return  Q; 
	}
};
class sensorNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
	matrixDouble R;
	matrixDouble R_root;
	public:
		sensorNoise(){
	            R = matrixDouble(3);
	            R_root = matrixDouble(std::sqrt(3)); 
	        }
		matrixDouble function(const vectorDouble &est, int t) override{
			return R;
		}
		matrixDouble sqrt(const vectorDouble &est, int t) override{
			return R_root; 
		}
};

class transitionJac: public jacobianDiscrete<vectorDouble,matrixDouble>{
	double rate;
	double max_pop;
	public:
		transitionJac(double rate, double max_pop):rate(rate),max_pop(max_pop){

		}
		matrixDouble function(const vectorDouble & val, int t){
			return matrixDouble(1 + rate - (2*rate*val.getSystemValue())/max_pop);
		}
};

class measurementJac: public jacobianDiscrete<vectorDouble,matrixDouble>{
	public:
		matrixDouble function(const vectorDouble & val, int t){
			return matrixDouble(1);
		}
};

class stateModel: public discreteModel<vectorDouble>{
	double rate;
	double max_pop;
	public:
		stateModel(double rate, double max_pop):rate(rate),max_pop(max_pop){
		}
		vectorDouble function(const vectorDouble & val, const int time) const override{
			return vectorDouble(val + rate*val*(1-val.getSystemValue()/max_pop));
		}
};


class measurementModel: public discreteModel<vectorDouble>{
	public:
		measurementModel(){
		}
		vectorDouble function(const vectorDouble & val, const int time) const override{
			return val;
		}
};



int main(){
	double rate = 0.01;
	double max_pop = 100;

	stateModel tm(rate,max_pop);
	measurementModel mm;
	processNoise pn;
	sensorNoise sn;
	transitionJac tj(rate,max_pop);
	measurementJac mj;	
	int stateCount = 1;
	int sensorCount = 1;
	discreteDiscreteFilterModel<vectorDouble,matrixDouble> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);
	
	
	
	int initialTime = 0;
	vectorDouble initialState(max_pop/2.0);
	discreteDiscreteFilterSolver<vectorDouble,matrixDouble> ddfs(&ddfm,initialTime,initialState);
	
	int steps = 1000;
	ddfs.setSeed(0);
	ddfs.solve(steps);
	
	vectorDouble * states = ddfs.getSolvedStates();
	vectorDouble * measurements = ddfs.getSolvedMeasurements();
	
	vectorDouble initialEstimate(max_pop/2.0);
	matrixDouble initialCovariance(1);

	discreteDiscreteKalmanFilter<vectorDouble,matrixDouble> filter(initialTime, initialEstimate, initialCovariance,&ddfm);

	vectorDouble* estimates = new vectorDouble[steps];
	for(int i = 0; i < steps; i++){
		filter.predict(1);
		filter.update(measurements[i]);
		
		estimates[i] = filter.getCurrentEstimate();
		//save the estimate or save the covariance, save only every 10th estimate,do nothing,
		//just run predict if you're simulating infrequent measurements..it's up to you and your context
	}

	std::ofstream outputfileEstimates;
	outputfileEstimates.open("estimatesNonlin.txt");
	for(int i = 0; i < steps; i++){
		outputfileEstimates<<std::setprecision(20)<<estimates[i]<<std::endl;
	}outputfileEstimates.close();
	

	std::ofstream outputfileMeasurements;
	outputfileMeasurements.open("measurementsNonlin.txt");
	for(int i = 0; i < steps; i++){
		outputfileMeasurements<<std::setprecision(20)<<measurements[i]<<std::endl;
	}outputfileMeasurements.close();
	
	std::ofstream outputfileStates;
	outputfileStates.open("statesNonlin.txt");
	for(int i = 0; i < steps; i++){
		outputfileStates<<std::setprecision(20)<<states[i]<<std::endl;
	}outputfileStates.close();
	
	delete [] estimates;
	delete [] states;
	delete [] measurements;
}
