#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "../../include/discreteDiscreteExtendedKalmanFilter.h"
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


class test: public discreteDiscreteExtendedKalmanFilter{
	MatrixXd f_mat;
	MatrixXd qc,lc,cc,lc_jac,cc_jac;
	MatrixXd * A;
	MatrixXd * J;
	int mode;
	public:
		VectorXd f(VectorXd x, double t) override{
			VectorXd nonlin( ((state_count+1)*state_count) /2);
			int counter = 0;
			for(int i = 0; i < state_count; i++){
				for(int j = i; j < state_count; j++){
					nonlin(counter) = x(i)*x(j);
					counter++;
				}
			}
			return qc*nonlin + lc*x + cc;
		}
		MatrixXd f_jacobian(VectorXd x, double t) override{
			MatrixXd kron_x(state_count*state_count,state_count);
			kron_x.setZero();
			for(int i = 0; i < state_count; i++){
				for(int j = 0; j < state_count; j++){
					kron_x(i*state_count + j,i) = x(j);
				}
			}
			return lc_jac*kron_x + cc_jac;
		}

		void set_f(MatrixXd Qc, MatrixXd Lc, MatrixXd Cc, MatrixXd Lc_jac, MatrixXd Cc_jac){
			qc = Qc;
			lc = Lc;
			cc = Cc;
			lc_jac = Lc_jac;
			cc_jac = Cc_jac;
		}
		MatrixXd returnEstimate(){
			return filtered_data;
		}

		MatrixXd ** returnCovariances(){
			return covariances;
		}

		test(int state_count, int sensor_count, int observation_count, Eigen::VectorXd x_0, Eigen::MatrixXd P_0, Eigen::MatrixXd measurements, Eigen::MatrixXd observation_matrix, Eigen::MatrixXd Q, Eigen::MatrixXd R):discreteDiscreteExtendedKalmanFilter(state_count,sensor_count,observation_count,x_0,P_0,measurements,observation_matrix,Q,R,true){}
};


int main(int argc, char *argv[]){
	std::ifstream nameFileout;

	MatrixXd C = load_csv<MatrixXd>("C.txt");
	
	stringstream ssn1;ssn1<<"qc_nonlin.txt";
	stringstream ssn2;ssn2<<"lc_nonlin.txt";
	stringstream ssn3;ssn3<<"cc_nonlin.txt";
	stringstream ssn4;ssn4<<"lc_jac.txt";
	stringstream ssn5;ssn5<<"cc_jac.txt";
	MatrixXd qc = load_csv<MatrixXd>(ssn1.str());
	MatrixXd lc = load_csv<MatrixXd>(ssn2.str());
	MatrixXd cc = load_csv<MatrixXd>(ssn3.str());
	MatrixXd lc_jac = load_csv<MatrixXd>(ssn4.str());
	MatrixXd cc_jac = load_csv<MatrixXd>(ssn5.str());

	MatrixXd Q_d = load_csv<MatrixXd>("Q_d.txt");
	MatrixXd R_d = load_csv<MatrixXd>("R_d.txt");
	MatrixXd x_0 = load_csv<MatrixXd>("x_0.txt");	
	MatrixXd P_0 = load_csv<MatrixXd>("P_0.txt");	

	MatrixXd measurement_count = load_csv<MatrixXd>("measurement_count.txt");
	MatrixXd state_count_mat = load_csv<MatrixXd>("state_count.txt");

	int state_count = state_count_mat(0,0);
	int sensor_count = 1;
	int observation_count = measurement_count(0,0);

	MatrixXd measurements  = load_csv<MatrixXd>("measurements.txt");  
	test t(state_count,sensor_count,observation_count,x_0,P_0,measurements,C,Q_d,R_d);
	t.set_f(qc,lc,cc,lc_jac,cc_jac);
	t.filter();

	t.saveData("output_data");	
}
