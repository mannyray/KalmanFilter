#include "sampleModelBoost.h"

int main(){
	stateModel tm;
	measurementModel mm;
	processNoise pn;
	sensorNoise sn;
	transitionJac tj;
	measurementJac mj;	
	int stateCount = 2;
	int sensorCount = 1;
	discreteDiscreteFilterModel<vectorBoost,matrixBoost> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);
	
	int initialTime = 0;
	boost::numeric::ublas::vector<double> tmp1(2);
	tmp1(0) = 0;
	tmp1(1) = 0;
	vectorBoost initialState(tmp1);
	discreteDiscreteFilterSolver<vectorBoost,matrixBoost> ddfs(&ddfm,initialTime,initialState);
	 
	int steps = 10000;
	ddfs.setSeed(0);
	ddfs.solve(steps);
	
	vectorBoost * states = ddfs.getSolvedStates();
	vectorBoost * measurements = ddfs.getSolvedMeasurements();
	

	boost::numeric::ublas::vector<double> tmp2(2);
	tmp2(0) = 0;
	tmp2(1) = 0;
	vectorBoost initialEstimate(tmp2);
	
	boost::numeric::ublas::matrix<double> tmp3(2,2);
	tmp3(0,0) = 1;
	tmp3(0,1) = 0;
	tmp3(1,0) = 0;
	tmp3(1,1) = 0.01;
	matrixBoost initialCovariance(tmp3);
	discreteDiscreteKalmanFilterSquareRoot<vectorBoost,matrixBoost> filter(initialTime, initialEstimate, initialCovariance,&ddfm);

	vectorBoost * estimates = new vectorBoost[steps];
	for(int i = 0; i < steps; i++){
		filter.predict(1);
		filter.update(measurements[i]);
		
		estimates[i] = filter.getCurrentEstimate();
		//save the estimate or save the covariance, save only every 10th estimate,do nothing,
		//just run predict if you're simulating infrequent measurements..it's up to you and your context
	}

	std::ofstream outputfileEstimates;
	outputfileEstimates.open("estimatesBoostSqrt.txt");
	for(int i = 0; i < steps; i++){
		outputfileEstimates<<std::setprecision(20)<<estimates[i]<<std::endl;
	}outputfileEstimates.close();
	

	std::ofstream outputfileMeasurements;
	outputfileMeasurements.open("measurementsBoostSqrt.txt");
	for(int i = 0; i < steps; i++){
		outputfileMeasurements<<std::setprecision(20)<<measurements[i]<<std::endl;
	}outputfileMeasurements.close();
	
	std::ofstream outputfileStates;
	outputfileStates.open("statesBoostSqrt.txt");
	for(int i = 0; i < steps; i++){
		outputfileStates<<std::setprecision(20)<<states[i]<<std::endl;
	}outputfileStates.close();
	
	delete [] estimates;
	delete [] states;
	delete [] measurements;
}
