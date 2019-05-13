#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <string>


class KalmanFilter{
	protected:
		//length of state vector
		int state_count;

		//length of measurement vector
		int sensor_count;

		//the amount of total observations
		int observation_count;

		//the time difference between each observation
		double dt_between_observations;

		//current time. Assuming that first observation is at time 0.
		double current_time = 0;

		//initial state. Dimension: 'state_count'X1
		Eigen::VectorXd x_0;	

		//initial covariance. Dimension: 'state_count'X'state_count'	
		Eigen::MatrixXd P_0;
	
		//current covariance matrix P
		Eigen::MatrixXd P;
	
		//current estimate of state x
		Eigen::VectorXd x;

		//the observations themselves. There are 'observation_count' observations
		//and 'sensor_count' the amount of measurements coming in with each observation.
		//Therefore the required dimension of measurements is 
		//'sensor_count'X'observation_count'
		Eigen::MatrixXd measurements;

		//the matrix that converts state to measurement.
		//Dimensions: 'sensor_count'X'state_count'
		Eigen::MatrixXd observation_matrix;

		//the filtered data.
		//Dimensions:'state_count'X'observation_count'
		//the nth column represents the estimated state
		//at t=dt_between*(measurements - 1);
		Eigen::MatrixXd filtered_data;

		//If set to true, covariances is filled with computed covariances.
		//For large state count this takes up a lot of memory which is why
		//it is set to zero.
		bool saveCovariances = false;

		//The computed covariances matrices. (array of matrices)
		Eigen::MatrixXd ** covariances;
		
		//The traces of the covariances
		Eigen::MatrixXd traces;

		//process covariance matrix Q  
		//Dimension: 'state_count'X'state_count'
		Eigen::MatrixXd Q;

		//sensor covariance matrix R
		//Dimension: 'sensor_count'X'sensor_count'
		Eigen::MatrixXd R;

		//output format of the data saved to file - set in the constructor
		//to CSV format. This is not an efficient storage format 
		//as a double now is written with all the digits in ascii form.
		Eigen::IOFormat CSVFormat;	
	
		void writeToCSVfile(std::string name, Eigen::MatrixXd & matrix){
			std::ofstream file(name.c_str(),std::fstream::app);
			file << matrix.transpose().format(CSVFormat)<<std::endl;
		}
		
		void writeToCSVfile(std::string name, Eigen::MatrixXd * matrix){
			std::ofstream file(name.c_str(),std::fstream::app);
			file << matrix->format(CSVFormat)<<std::endl;
		}
	protected:
		//return covariance P. For square root implementatios the implementation
		//of the function will change
		virtual Eigen::MatrixXd getP(){
			return P;
		}

		//The following two methods are very specific for a general 
		//Kalman filter class such as this and can be moved to a more specific
		//sub class from which discreteDiscreteExtendedKalmanFilter
		//and continuousDiscreteExtendedKalmanFilter could inherit from
		//(instead of inheriting from this one). However, this was not done
		//in this case which technically goes against the purest forms of OOP
		//spirits.

		//Lower triangulate matrix 'in' using Householder algorithm.
		//Similar to a QR transform: A = QR (Q unitary, R is uppertriangular)
		void lower_triangulate(Eigen::MatrixXd &in){
			Eigen::VectorXd v; 
			Eigen::VectorXd w; 
			for(int i = 0; i < in.rows(); i++){
				int length = in.cols() - i;
				v = in.row(i).tail(length);
				double mew = v.norm();
				if(mew!=0){//(-2*(v(0) < 0)+1) sign of v(0)
					double beta  = v(0)  + (-2*(v(0) < 0)+1)*mew;	
					v = v/beta;
				}
				v(0) = 1;
				w = -2.0/(v.squaredNorm())*(in.block(i,i, in.rows()-i , in.cols()-i)*v);
				in.block(i,i, in.rows()-i , in.cols()-i) = in.block(i,i, in.rows()-i , in.cols()-i) + w*v.transpose();
			}
		}

		//When adding two matrices that can be expressed as square roots:
		//A = A1*(A1)^T
		//B = A2*(A2)^T
		//where ^T is the transpose, then the function add_roots(A1,A2) returns
		//A3 such that
		//A + B = A3*(A3)^T
		Eigen::MatrixXd add_roots(Eigen::MatrixXd  A1, Eigen::MatrixXd  A2){
			Eigen::MatrixXd tmp(state_count,state_count*2);
			tmp.block(0,0,state_count,state_count) = A1;
			tmp.block(0,state_count,state_count,state_count) = A2;
			lower_triangulate(tmp);
			return tmp.block(0,0,state_count,state_count);
		}

	public:
		//write computed data to file 
		void saveData(std::string file_name){
			if(saveCovariances == true){
				std::stringstream ss2; ss2<<file_name<<"_covariance_traces.txt";
				std::stringstream ss3; ss3<<file_name<<"_covariances.txt";
				Eigen::MatrixXd tmp(observation_count,1);
				for(int i = 0; i < observation_count; i++){
					double trace = 0;
					Eigen::MatrixXd tmpp = *(covariances[i]);
					for(int j = 0; j < state_count; j++){
						trace+=tmpp(j,j);
					}
					tmp(i,0) = trace; 
					writeToCSVfile(ss3.str(),tmpp);
				}	
				writeToCSVfile(ss2.str(),tmp);
			}else{
				std::stringstream ss2; ss2<<file_name<<"_covariance_traces.txt";
				writeToCSVfile(ss2.str(),traces);
			}
			std::stringstream ss1; ss1<<file_name<<"_filtered_data.txt";
			writeToCSVfile(ss1.str(),filtered_data);
		}
		//Function f - the assumed model for the Kalman filter
		//to be implemented in specific Kalman filter sub class
		//depending on context
		virtual Eigen::VectorXd f(Eigen::VectorXd x, double t) = 0;
		
		//Jacobian of f at state x at t
		virtual Eigen::MatrixXd f_jacobian(Eigen::VectorXd x, double t) = 0;

		//constructor
		KalmanFilter(int state_count, int sensor_count, int observation_count, double dt_between_observations, Eigen::VectorXd X_0, Eigen::MatrixXd p_0, Eigen::MatrixXd measur, Eigen::MatrixXd obs_mat, Eigen::MatrixXd q, Eigen::MatrixXd r,bool saveCovariances = false):state_count(state_count),sensor_count(sensor_count),observation_count(observation_count),dt_between_observations(dt_between_observations),saveCovariances(saveCovariances){ 
			//the format with which the data will be saved. It is not the most compact form as it prints
			//out the entire number. This can be changed to something more compact.
			CSVFormat = Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");	
			
			measurements = measur;
			x_0 = X_0;
			P_0 = p_0; 
			R = r;
			Q = q;
			observation_matrix = obs_mat;
			P = P_0; x = x_0;
			
			//dimensions are checked. However, the semipositive definitness of P_0,Q,R are not checked
			//this is up to the user.
			assert(P_0.rows() == P_0.cols() && P_0.rows() == state_count &&
				"covariance P has wrong dimensions");
			assert(x_0.rows() == state_count && x_0.cols() == 1 &&
				"state x has wrong dimensions");
			assert(measurements.rows() == sensor_count && measurements.cols() == observation_count &&
				"measurements has wrong dimension");
			assert(observation_matrix.rows() == sensor_count && observation_matrix.cols() == state_count &&
				"observation matrix has wrong dimensions"); 
			assert(Q.rows() == Q.cols() && Q.rows() == state_count &&
				"process noise matrix incorrect");
			assert(R.rows() == R.cols() && R.rows() == sensor_count &&
				"sensor noise matrix incorrect");
			
			//define data and space for its 
			traces = Eigen::MatrixXd(observation_count,1);
			filtered_data = Eigen::MatrixXd(state_count, observation_count);
			if(saveCovariances==true){
				covariances = new Eigen::MatrixXd*[observation_count];    
				for(int i = 0; i < observation_count; i++){
					covariances[i] = new Eigen::MatrixXd(state_count,state_count); 
				}
			}
		}

		//the predict and update phase defined in this class
		//imply that all inheriting classes have measurements
		//that are discrete in time.		

		//the predict phase of the filter
		virtual void predict() = 0;

		//the update phase of the filter
		virtual void update(int i) = 0;

		//filters through all the measurements and produces the 
		//estimates and or covariances
		virtual void filter() = 0;
		
		//destructor
		~KalmanFilter(){
			if(saveCovariances==true){
				for(int i=0; i < observation_count;i++){
					delete covariances[i];
				}delete [] covariances;
			}
		}
};

#endif
