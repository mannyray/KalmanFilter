#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "../../include/continuousDiscreteExtendedKalmanFilter.h"
#include <assert.h>

using namespace std;
using namespace Eigen;

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

template<typename M>
M load_csv (const std::string & path, bool isRow = true) {
	std::ifstream indata;
	indata.open(path);
	std::string line;
	std::vector<double> values;
	uint rows = 0;
	while (std::getline(indata, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		while (std::getline(lineStream, cell, ',')) {
			values.push_back(std::stod(cell));
		}
		++rows;
	}
	if(isRow){
		return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
	}else{
		return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, ColMajor>>(values.data(), rows, values.size()/rows);
	}
}


class sample_code: public continuousDiscreteExtendedKalmanFilter{
	MatrixXd f_mat;
	double rate = 0.5; 
	double max_pop = 100.0;
	MatrixXd ident;
	public:
		VectorXd f(VectorXd x, double t) override{
			return rate*x*(ident - (1.0/max_pop)*x);
		}
		MatrixXd f_jacobian(VectorXd x, double t) override{
			return rate*ident - (1.0/max_pop)*(2*x*rate);
		}
		void set_rk4steps(int rk4){
			dscb_measurement = rk4;
		}

		void set_f(MatrixXd f){
			f_mat = f;
		}
		MatrixXd returnEstimate(){
			return filtered_data;
		}

		MatrixXd ** returnCovariances(){
			return covariances;
		}

		sample_code(int state_count, int sensor_count, int observation_count, double dt_between_observations, Eigen::VectorXd x_0, Eigen::MatrixXd P_0, Eigen::MatrixXd measurements, Eigen::MatrixXd observation_matrix, Eigen::MatrixXd Q, Eigen::MatrixXd R):continuousDiscreteExtendedKalmanFilter(state_count,sensor_count,observation_count,dt_between_observations,x_0,P_0,measurements,observation_matrix,Q,R){ident = MatrixXd(1,1);ident(0,0) = 1;}


};


int main(){
	int state_count = 1;
	int sensor_count = 1;
	double rk4steps = 1;   

	double t_start = 0;
	double t_final = 10;
	double observation_count = 1000;
	double dt_between_observations = (t_final - t_start)/observation_count;

	MatrixXd C(1,1);C(0,0) = 1;
	MatrixXd Q_c(1,1);Q_c(0,0) = 10;
	MatrixXd R_d(1,1);R_d(0,0) = 0.01/dt_between_observations;
	MatrixXd x_0(1,1); x_0(0,0) = 50;	
	MatrixXd P_0(1,1); P_0(0,0) = 1;	
	MatrixXd measurements  = load_csv<MatrixXd>("measurements.txt");  

	sample_code t(state_count,sensor_count,observation_count,dt_between_observations,x_0,P_0,measurements,C,Q_c,R_d);

	t.set_rk4steps(rk4steps);

	t.filter();
	t.saveData("output");	
}
