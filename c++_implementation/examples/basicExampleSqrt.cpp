#include "sampleModel.h"

int main(){
	stateModel tm;
	measurementModel mm;
	processNoise pn;
	sensorNoise sn;
	transitionJac tj;
	measurementJac mj;	
	int stateCount = 2;
	int sensorCount = 1;
	discreteDiscreteFilterModel<Eigen::VectorXd,Eigen::MatrixXd> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);
	
	int initialTime = 0;
	Eigen::VectorXd initialState(2);
	initialState<<0,0;
	discreteDiscreteFilterSolver<Eigen::VectorXd,Eigen::MatrixXd> ddfs(&ddfm,initialTime,initialState);
	
	int steps = 10000;
	ddfs.setSeed(0);
	ddfs.solve(steps);
	
	Eigen::VectorXd * states = ddfs.getSolvedStates();
	Eigen::VectorXd * measurements = ddfs.getSolvedMeasurements();
	
	Eigen::VectorXd initialEstimate(2);
	initialEstimate<<0,0;
	Eigen::MatrixXd initialCovariance(2,2);
	initialCovariance<<1,0,0,0.01;

	discreteDiscreteKalmanFilterSquareRoot<Eigen::VectorXd,Eigen::MatrixXd> 
		filter(initialTime, initialEstimate, initialCovariance,&ddfm);

	Eigen::VectorXd * estimates = new Eigen::VectorXd[steps];
	for(int i = 0; i < steps; i++){
		filter.predict(1);
		filter.update(measurements[i]);
		
		estimates[i] = filter.getCurrentEstimate();
		//save the estimate or save the covariance, save only every 10th estimate,do nothing,
		//just run predict if you're simulating infrequent measurements..it's up to you and your context
	}
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
	std::ofstream outputfileEstimates;
	outputfileEstimates.open("estimatesSQRT.txt");
	for(int i = 0; i < steps; i++){
		outputfileEstimates<<estimates[i].format(HeavyFmt)<<std::endl;
	}outputfileEstimates.close();
	

	std::ofstream outputfileMeasurements;
	outputfileMeasurements.open("measurementsSQRT.txt");
	for(int i = 0; i < steps; i++){
		outputfileMeasurements<<measurements[i].format(HeavyFmt)<<std::endl;
	}outputfileMeasurements.close();
	
	std::ofstream outputfileStates;
	outputfileStates.open("statesSQRT.txt");
	for(int i = 0; i < steps; i++){
		outputfileStates<<states[i].format(HeavyFmt)<<std::endl;
	}outputfileStates.close();
	
	delete [] estimates;
	delete [] states;
	delete [] measurements;
}
