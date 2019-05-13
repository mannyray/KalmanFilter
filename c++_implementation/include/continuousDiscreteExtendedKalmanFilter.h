#ifndef CONTINUOUSDISCRETEEXTENDEDKALMANFILTER_H
#define CONTINUOUSDISCRETEEXTENDEDKALMANFILTER_H

#include "KalmanFilter.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>

class continuousDiscreteExtendedKalmanFilter: public KalmanFilter{
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
		
		//default step count between measurements
		//dt =  dt_between_observations/dscb_measurement 
		double dscb_measurement = 5.0;

		//RK4 method to do a single step where (d/dt)x = f(x,t)
		void evolve_state(const double dt, const double time, const Eigen::VectorXd &estimate_x_n,Eigen::VectorXd &estimate_x_np1){
				Eigen::VectorXd k1,k2,k3,k4;
				k1 = dt*f(estimate_x_n,time);
				k2 = dt*f(estimate_x_n + 0.5*k1, time + 0.5*dt);
				k3 = dt*f(estimate_x_n + 0.5*k2, time + 0.5*dt);
				k4 = dt*f(estimate_x_n + k3, time + dt);
				estimate_x_np1 = estimate_x_n + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
		}
		
		//RK4 method to do a single step where (d/dt)PHI = jac PHI
		Eigen::MatrixXd evolve_phi(const double dt,const Eigen::MatrixXd  &jac, const Eigen::MatrixXd &Phi){
			Eigen::MatrixXd jac_a_2 = jac*jac;
			Eigen::MatrixXd jac_a_3 = jac*jac_a_2;
			Eigen::MatrixXd jac_a_4 = jac*jac_a_3;
			return (Eigen::MatrixXd::Identity(state_count,state_count) + dt*jac 
				+ (1/2.0)*dt*dt*jac_a_2 + (1/6.0)*dt*dt*dt*jac_a_3 
				+ (1/24.0)*dt*dt*dt*dt*jac_a_4)*Phi;
		}

	public:
		void predict() override{
			double time = current_time;
			double dt = dt_between_observations/dscb_measurement;

			Eigen::VectorXd estimate_x_n = x;
			Eigen::VectorXd estimate_x_np1 = x;
	
			Eigen::MatrixXd P_root_n = P_root;
			Eigen::MatrixXd P_root_np1; 
		
			Eigen::MatrixXd M_inv(state_count,state_count);M_inv.setZero();
			Eigen::MatrixXd jac(state_count,state_count); 
	
			Eigen::VectorXd e_np1(state_count);e_np1.setZero();

			Eigen::MatrixXd phi = Eigen::MatrixXd::Identity(state_count, state_count);
			Eigen::MatrixXd phi_sum = Eigen::MatrixXd::Zero(state_count,state_count); 
			while(time < current_time + dt_between_observations){
				if(time + dt > current_time + dt_between_observations){
					//in case dt is too big for final gap
					dt = current_time + dt_between_observations - time;
				}
			
				jac = f_jacobian(estimate_x_n, time);

				evolve_state(dt,time,estimate_x_n,estimate_x_np1);
				Eigen::MatrixXd phi_next = evolve_phi(dt,jac,phi);
				estimate_x_n = estimate_x_np1;

				//trapezoidal sum
				phi_sum = add_roots(phi_sum,sqrt(dt*0.5)*add_roots(phi*Q_root, phi_next*Q_root));
				
				phi = phi_next;
				time = time + dt;
			}
			Eigen::MatrixXd tmp(state_count,state_count*2);
			tmp.block(0,0,state_count,state_count) = phi_sum;
  		 	tmp.block(0,state_count,state_count,state_count) = phi*P_root_n;
			lower_triangulate(tmp);		

			current_time += dt_between_observations;
			x = estimate_x_n;
			P_root = tmp.block(0,0,state_count,state_count);
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
		
		continuousDiscreteExtendedKalmanFilter(int state_count, int sensor_count, int observation_count, double dt_between_observations, Eigen::VectorXd x_0, Eigen::MatrixXd P_0, Eigen::MatrixXd measurements, Eigen::MatrixXd observation_matrix, Eigen::MatrixXd Q, Eigen::MatrixXd R):KalmanFilter(state_count,sensor_count,observation_count,dt_between_observations,x_0,P_0,measurements,observation_matrix,Q,R){
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
		~continuousDiscreteExtendedKalmanFilter(){};
};
#endif
