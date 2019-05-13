#ifndef DISCRETEDISCRETEEXTENDEDKALMANFILTER_H
#define DISCRETEDISCRETEEXTENDEDKALMANFILTER_H

#include "KalmanFilter.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>

class discreteDiscreteExtendedKalmanFilter: public KalmanFilter{
	protected:
		//root of covariance P: P = P_root*P_root'.
		//It is the user's responsibility to set this correctly. In case 
		//the matrix is positive definite one case use chol(P)' or in case
		//P is identity one can take square roots of the diagonal entries.
		Eigen::MatrixXd P_root;
		
		//root of process noise matrix Q: Q = Q_root*Q_root'
		//It is the user's responsibility to set this correctly. In case 
		//the matrix is positive definite one case use chol(Q)' or in case
		//Q is identity one can take square roots of the diagonal entries.
		Eigen::MatrixXd Q_root;
		
		//root of sensor noise matrix R: R = R_root*R_root'
		//It is the user's responsibility to set this correctly. In case 
		//the matrix is positive definite one case use chol(R)' or in case
		//R is identity one can take square roots of the diagonal entries.
		Eigen::MatrixXd R_root;

	public:
		void predict() override{
			//time is just an increasing integer (k) for the discrete filter
			double time = current_time;

			Eigen::MatrixXd tmp(state_count,state_count*2);
			tmp.block(0,0,state_count,state_count) = f_jacobian(x,current_time)*P_root;
  		 	tmp.block(0,state_count,state_count,state_count) = Q_root;
			lower_triangulate(tmp);		

			x = f(x,current_time);//x_{k+1} = f(x_{k},k}
			P_root = tmp.block(0,0,state_count,state_count);
			current_time += 1;
		}

		void update(int iter) override{
			Eigen::MatrixXd tmp(state_count + sensor_count, state_count + sensor_count);
			tmp.setZero();
			tmp.block(0,0,sensor_count, sensor_count) = R_root;
			tmp.block(0,sensor_count, sensor_count, state_count) = observation_matrix*P_root;
			tmp.block(sensor_count, sensor_count,state_count, state_count) = P_root;

			lower_triangulate(tmp);
			Eigen::VectorXd e_k = (tmp.block(0,0,sensor_count,sensor_count)).colPivHouseholderQr().solve(
				(measurements.col(iter-1) - observation_matrix*x));
			
			x = x + tmp.block(sensor_count, 0,state_count, sensor_count)*e_k;
			P_root = tmp.block(sensor_count,sensor_count,state_count,state_count); 
		}	
		
		Eigen::MatrixXd getP() override{
			return P_root*P_root.transpose();
		}
		
		void filter() override{
			Eigen::VectorXd estimate_x = x_0;
			Eigen::MatrixXd P = P_0;
			Eigen::MatrixXd I = Eigen::MatrixXd::Identity(state_count, state_count);

			Eigen::MatrixXd observation_matrix_T = observation_matrix.transpose();
			Eigen::MatrixXd K(state_count,sensor_count);
			filtered_data.col(0) = x_0;
			if(saveCovariances==true){
				*(covariances[0]) = getP();
			}else{
				double tr = 0;
				Eigen::MatrixXd tmp = getP();
				for(int kk = 0; kk < state_count; kk++){
					tr+=tmp(kk,kk);
				}
				traces(0,0) = tr;
			}

			for(int iter = 1; iter < observation_count; iter++){
				predict();
				update(iter);

				if(saveCovariances==true){
					*(covariances[iter]) = getP();
				}
				else{
					double tr = 0;
					Eigen::MatrixXd tmp = getP();
					for(int kk = 0; kk < state_count; kk++){
						tr+=tmp(kk,kk);
					}
					traces(iter,0) = tr;
				}
				filtered_data.col(iter) = x;
			}
		}
		
		discreteDiscreteExtendedKalmanFilter(int state_count, int sensor_count, int observation_count, Eigen::VectorXd x_0, Eigen::MatrixXd P_0, Eigen::MatrixXd measurements, Eigen::MatrixXd observation_matrix, Eigen::MatrixXd Q, Eigen::MatrixXd R, bool saveCovariances = false):KalmanFilter(state_count,sensor_count,observation_count,0,x_0,P_0,measurements,observation_matrix,Q,R,saveCovariances){

			P_root = P;
			R_root = R;
			Q_root = Q;
			
			for(int i=0; i < state_count; i++){
				P_root(i,i) = sqrt(P(i,i));
				Q_root(i,i) = sqrt(Q(i,i));
			}		
			for(int i=0; i < sensor_count; i++){
				R_root(i,i) = sqrt(R(i,i));
			}
		}
		~discreteDiscreteExtendedKalmanFilter(){};
};
#endif
