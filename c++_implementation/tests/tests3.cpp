#include <gtest/gtest.h> 
#include "../include/model.h"
#include "../include/KalmanFilter.h"
#include "../include/filterModel.h"
#include "../include/mathWrapper/double.h"
#include "../include/mathWrapper/eigen.h"
#include "../include/mathWrapper/boost.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>



#include <math.h>
#include "../include/Eigen.h" 

 
#include <iostream>


TEST(discreteDiscreteKalmanFilter, test){
	
	int stateCount = 1;
	int sensorCount = 1;
	int initialTime = 0;
	vectorDouble initialEstimate(1);
	matrixDouble initialCovariance(1);

	Eigen::VectorXd x(1);x<<1;	
	vectorEigen initialEstimateEigen(x);
	Eigen::MatrixXd y(1,1);y<<1;
	matrixEigen initialCovarianceEigen(y);

	boost::numeric::ublas::matrix<double> initialCovarianceBoost(1,1);
	initialCovarianceBoost(0,0)=1;
	boost::numeric::ublas::vector<double> initialEstimateBoost(1);
	initialEstimateBoost(0) = 1;

	Eigen::VectorXd initialEstimateEigenSimple(1);initialEstimateEigenSimple<<1;	
	Eigen::MatrixXd initialCovarianceEigenSimple(1,1);initialCovarianceEigenSimple<<1;




	class stateModel: public discreteModel<vectorDouble>{
		public:
			vectorDouble function(const vectorDouble & val, const int time) const override{
				vectorDouble result(val.getSystemValue() + 0.1);
				return result;
			}
	};
	class stateModelEigen: public discreteModel<vectorEigen>{
		Eigen::VectorXd velocity;
		public:
			stateModelEigen(){
				Eigen::VectorXd vel(1);vel<<0.1;
				velocity = vel;
			}
			vectorEigen function(const vectorEigen & val, const int time) const override{
				vectorEigen result(val.getSystemValue() + velocity);
				return result;
			}
	};
	class stateModelBoost: public discreteModel<vectorBoost>{
		boost::numeric::ublas::vector<double> velocity;
		public:
			stateModelBoost(){
				boost::numeric::ublas::vector<double> vel(1);vel(0) = 0.1;
				velocity = vel;
			}
			vectorBoost function(const vectorBoost & val, const int time) const override{
				vectorBoost result(val.getSystemValue() + velocity);
				return result;
			}
	};
	class stateModelEigenSimple: public discreteModel<Eigen::VectorXd>{
		Eigen::VectorXd velocity;
		public:
			stateModelEigenSimple(){
				Eigen::VectorXd vel(1);vel<< 0.1;
				velocity = vel;
			}
			Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
				return val + velocity;
			}
	};




	class measurementModel: public discreteModel<vectorDouble>{
		public:
			vectorDouble function(const vectorDouble & val, const int time) const override{
				vectorDouble result(val.getSystemValue());
				return result;
			}
	};
	class measurementModelEigen: public discreteModel<vectorEigen>{
		public:
			vectorEigen function(const vectorEigen & val, const int time) const override{
				vectorEigen result(val.getSystemValue());
				return result;
			}
	};
	class measurementModelBoost: public discreteModel<vectorBoost>{
		public:
			vectorBoost function(const vectorBoost & val, const int time) const override{
				vectorBoost result(val.getSystemValue());
				return result;
			}
	};
	class measurementModelEigenSimple: public discreteModel<Eigen::VectorXd>{
		public:
			Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
				return val;
			}
	};





	class transitionJac: public jacobianDiscrete<vectorDouble,matrixDouble>{
		public:
			matrixDouble function(const vectorDouble & val, int t){
				return matrixDouble(1);
			}
	};
	class measurementJac: public jacobianDiscrete<vectorDouble,matrixDouble>{
		public:
			matrixDouble function(const vectorDouble & val, int t){
				return matrixDouble(1);
			}
	};
	class transitionJacEigen: public jacobianDiscrete<vectorEigen,matrixEigen>{
		Eigen::MatrixXd v;
		public:
			transitionJacEigen(){
				Eigen::MatrixXd tmp(1,1);tmp<<1;
				v = tmp;
			}

			matrixEigen function(const vectorEigen & val, int t){
				return matrixEigen(v); 
			}
	};
	class measurementJacEigen: public jacobianDiscrete<vectorEigen,matrixEigen>{
		Eigen::MatrixXd v;
		public:
			measurementJacEigen(){
				Eigen::MatrixXd tmp(1,1);tmp<<1;
				v = tmp;
			}

			matrixEigen function(const vectorEigen & val, int t){
				return matrixEigen(v);
			}
	};
	class transitionJacBoost: public jacobianDiscrete<vectorBoost,matrixBoost>{
		boost::numeric::ublas::matrix<double> v;
		public:
			transitionJacBoost(){
				boost::numeric::ublas::matrix<double> tmp(1,1);tmp(0,0) = 1;
				v = tmp;

			}

			matrixBoost function(const vectorBoost & val, int t){
				return matrixBoost(v); 
			}
	};
	class measurementJacBoost: public jacobianDiscrete<vectorBoost,matrixBoost>{
		boost::numeric::ublas::matrix<double> v;
		public:
			measurementJacBoost(){
				boost::numeric::ublas::matrix<double> tmp(1,1);tmp(0,0) = 1;
				v = tmp;
			}

			matrixBoost function(const vectorBoost & val, int t){
				return matrixBoost(v);
			}
	};
	class transitionJacEigenSimple: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
		Eigen::MatrixXd v;
		public:
			transitionJacEigenSimple(){
				Eigen::MatrixXd tmp(1,1);tmp << 1;
				v = tmp;

			}

			Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
				return v; 
			}
	};
	class measurementJacEigenSimple: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
		Eigen::MatrixXd v;
		public:
			measurementJacEigenSimple(){
				Eigen::MatrixXd tmp(1,1);tmp<< 1;
				v = tmp;
			}

			Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
				return v;
			}
	};




	class processNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
		matrixDouble function(const vectorDouble &est,int t) override{
			return matrixDouble(0.01);
		}	
		matrixDouble sqrt(const vectorDouble &est, int t) override{
			return matrixDouble(0.1);
		}
	};
	class sensorNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
		matrixDouble function(const vectorDouble &est,int t) override{
			return matrixDouble(0.01);
		}	
		matrixDouble sqrt(const vectorDouble &est, int t) override{
			return matrixDouble(0.1);
		}
	};
	class processNoiseEigen: public discreteNoiseCovariance<vectorEigen,matrixEigen>{
		Eigen::MatrixXd v;
		Eigen::MatrixXd v2;
		public:
		processNoiseEigen(){
			Eigen::MatrixXd tmp(1,1);tmp<<0.01;
			Eigen::MatrixXd tmp2(1,1);tmp2<<0.1;
			v = tmp;
			v2 = tmp2;
		}
		matrixEigen function(const vectorEigen &est,int t) override{
			return matrixEigen(v);
		}	
		matrixEigen sqrt(const vectorEigen &est, int t) override{
			return matrixEigen(v2);
		}
	};
	class sensorNoiseEigen: public discreteNoiseCovariance<vectorEigen,matrixEigen>{
		Eigen::MatrixXd v;
		Eigen::MatrixXd v2;
		public:
		sensorNoiseEigen(){
			Eigen::MatrixXd tmp(1,1);tmp<<0.01;
			Eigen::MatrixXd tmp2(1,1);tmp2<<0.1;
			v = tmp;
			v2 = tmp2;
		}
		matrixEigen function(const vectorEigen &est, int t) override{
			return matrixEigen(v);
		}	
		matrixEigen sqrt(const vectorEigen &est, int t) override{
			return matrixEigen(v2);
		}
	};
	class processNoiseBoost: public discreteNoiseCovariance<vectorBoost,matrixBoost>{
		boost::numeric::ublas::matrix<double> v;
		boost::numeric::ublas::matrix<double> v2;
		public:
		processNoiseBoost(){
			boost::numeric::ublas::matrix<double> tmp(1,1);tmp(0,0) = 0.01;
			boost::numeric::ublas::matrix<double> tmp2(1,1);tmp2(0,0) = 0.01;
			v = tmp;
			v2 = tmp2;
		}
		matrixBoost function(const vectorBoost &est,int t) override{
			return matrixBoost(v);
		}	
		matrixBoost sqrt(const vectorBoost &est, int t) override{
			return matrixBoost(v2);
		}
	};
	class sensorNoiseBoost: public discreteNoiseCovariance<vectorBoost,matrixBoost>{
		boost::numeric::ublas::matrix<double> v;
		boost::numeric::ublas::matrix<double> v2;
		public:
		sensorNoiseBoost(){
			boost::numeric::ublas::matrix<double> tmp(1,1);tmp(0,0) = 0.01;
			boost::numeric::ublas::matrix<double> tmp2(1,1);tmp2(0,0) = 0.01;
			v = tmp;
			v2 = tmp2;
		}
		matrixBoost function(const vectorBoost &est,int t) override{
			return matrixBoost(v);
		}	
		matrixBoost sqrt(const vectorBoost &est, int t) override{
			return matrixBoost(v2);
		}
	};
	class processNoiseEigenSimple: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
		Eigen::MatrixXd v;
		Eigen::MatrixXd v2;
		public:
		processNoiseEigenSimple(){
			Eigen::MatrixXd tmp(1,1);tmp<< 0.01;
			Eigen::MatrixXd tmp2(1,1);tmp2<<0.1;
			v = tmp;
			v2 = tmp2;
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est,int t) override{
			return v;
		}	
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return v2;
		}
	};
	class sensorNoiseEigenSimple: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
		Eigen::MatrixXd v;
		Eigen::MatrixXd v2;
		public:
		sensorNoiseEigenSimple(){
			Eigen::MatrixXd tmp(1,1);tmp<< 0.01;
			Eigen::MatrixXd tmp2(1,1);tmp2<<0.1;
			v = tmp;
			v2 = tmp2;
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est,int t) override{
			return v; 
		}	
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return v2;
		}
	};
	
	
	stateModel tm;
	measurementModel mm;
	processNoise pn;
	sensorNoise sn;
	transitionJac tj;
	measurementJac mj;	
	discreteDiscreteFilterModel<vectorDouble,matrixDouble> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);

	stateModelEigen tmEigen;
	measurementModelEigen mmEigen;
	processNoiseEigen pnEigen;
	sensorNoiseEigen snEigen;
	transitionJacEigen tjEigen;
	measurementJacEigen mjEigen;	
	discreteDiscreteFilterModel<vectorEigen,matrixEigen> ddfmEigen(&tmEigen,&mmEigen,&pnEigen,&snEigen,&tjEigen,&mjEigen,
		stateCount,sensorCount);

	stateModelBoost tmBoost;
	measurementModelBoost mmBoost;
	processNoiseBoost pnBoost;
	sensorNoiseBoost snBoost;
	transitionJacBoost tjBoost;
	measurementJacBoost mjBoost;	
	discreteDiscreteFilterModel<vectorBoost,matrixBoost> ddfmBoost(&tmBoost,&mmBoost,&pnBoost,&snBoost,&tjBoost,&mjBoost,
		stateCount,sensorCount);

	stateModelEigenSimple tmEigenSimple;
	measurementModelEigenSimple mmEigenSimple;
	processNoiseEigenSimple pnEigenSimple;
	sensorNoiseEigenSimple snEigenSimple;
	transitionJacEigenSimple tjEigenSimple;
	measurementJacEigenSimple mjEigenSimple;	
	discreteDiscreteFilterModel<Eigen::VectorXd,Eigen::MatrixXd> ddfmEigenSimple(&tmEigenSimple,&mmEigenSimple,
		&pnEigenSimple,&snEigenSimple,&tjEigenSimple,&mjEigenSimple,stateCount,sensorCount);
		

	discreteDiscreteKalmanFilter<vectorDouble,matrixDouble>
		filter(initialTime, initialEstimate, initialCovariance, &ddfm);

	discreteDiscreteKalmanFilter<vectorEigen,matrixEigen>
		filterEigen(initialTime, initialEstimateEigen, initialCovarianceEigen, &ddfmEigen);

	discreteDiscreteKalmanFilter<vectorBoost,matrixBoost>
		filterBoost(initialTime, initialEstimateBoost, initialCovarianceBoost, &ddfmBoost);

	discreteDiscreteKalmanFilter<Eigen::VectorXd,Eigen::MatrixXd>
		filterEigenSimple( initialTime, initialEstimateEigenSimple, initialCovarianceEigenSimple,
			&ddfmEigenSimple);


	for(int i = 0; i < 10; i++){
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),filterEigen.getCurrentEstimate().getSystemValue()[0],1e-9);
		ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),
			filterEigen.getCurrentCovariance().getSystemValue().coeff(0,0),1e-9);
		ASSERT_NEAR(filter.getCurrentTime(),filterEigen.getCurrentTime(),1e-9);

		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),filterBoost.getCurrentEstimate().getSystemValue()[0],1e-9);
		boost::numeric::ublas::matrix<double> tmp = filterBoost.getCurrentCovariance().getSystemValue();
		ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),tmp(0,0),1e-9);
		ASSERT_NEAR(filter.getCurrentTime(),filterBoost.getCurrentTime(),1e-9);


		filter.predict(1);
		filterEigen.predict(1);
		filterBoost.predict(1);
		filterEigenSimple.predict(1);
	}



	for(int i = 0; i < 10; i++){
		vectorDouble measurement(i);	
		Eigen::VectorXd tmp(1);tmp<<i;
		vectorEigen tmp2(tmp);
		boost::numeric::ublas::vector<double> tmpBoost(1);tmpBoost(0)=i;
		vectorBoost tmp2Boost(tmpBoost);
		
		double sk = 1*filter.getCurrentCovariance().getSystemValue()*1 + 0.01;
		double kk = filter.getCurrentCovariance().getSystemValue()*(1.0/sk);
		double curEst = filter.getCurrentEstimate().getSystemValue(); 

		filter.update(measurement);
		filterEigen.update(tmp2);
		filterBoost.update(tmp2Boost);
		filterEigenSimple.update(tmp);
		
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),filterEigen.getCurrentEstimate().getSystemValue()[0],1e-9);
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),filterEigenSimple.getCurrentEstimate().coeff(0,0),1e-9);
		ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),
			filterEigen.getCurrentCovariance().getSystemValue().coeff(0,0),1e-9);
		boost::numeric::ublas::matrix<double> t = filterBoost.getCurrentCovariance().getSystemValue();
		ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),t(0,0),1e-9);
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),filterBoost.getCurrentEstimate().getSystemValue()[0],1e-9);
	}
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
