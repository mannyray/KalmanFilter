#include <gtest/gtest.h> 
#include "../include/model.h"
#include "../include/KalmanFilter.h"
#include "../include/filterModel.h"
#include "../include/mathWrapper/double.h"
#include <math.h>
#include "../include/Eigen.h" 
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/operators.hpp>
#include <iostream>


TEST(discreteDiscreteKalmanFilter, test){
	
	int stateCount = 1;
	int sensorCount = 1;
	int initialTime = 0;
	//VECTOR_double initialEstimate(1);
	vectorDouble initialEstimate(1);
	//MATRIX_double initialCovariance(1);
	matrixDouble initialCovariance(1);

	class stateModel: public discreteModel<vectorDouble>{
		public:
			vectorDouble function(const vectorDouble & val, const int time) const override{
				vectorDouble result(val.getSystemValue() + 0.1);
				return result;
			}
	};

	class measurementModel: public discreteModel<vectorDouble>{
		public:
			vectorDouble function(const vectorDouble & val, const int time) const override{
				vectorDouble result(val.getSystemValue());
				return result;
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
	
	class processNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
		matrixDouble function(const vectorDouble &v, int t) override{
			return matrixDouble(0.01);
		}	
		matrixDouble sqrt(const vectorDouble &v, int t) override{
			return matrixDouble(0.1);
		}
	};

	class sensorNoise: public discreteNoiseCovariance<vectorDouble,matrixDouble>{
		matrixDouble function(const vectorDouble &v, int t) override{
			return matrixDouble(0.01);
		}	
		matrixDouble sqrt(const vectorDouble &v, int t) override{
			return matrixDouble(0.1);
		}
	};

	stateModel tm;
	measurementModel mm;
	processNoise pn;
	sensorNoise sn;
	transitionJac tj;
	measurementJac mj;

	discreteDiscreteFilterModel<vectorDouble,matrixDouble> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);
	discreteDiscreteKalmanFilter<vectorDouble,matrixDouble>
		filter(initialTime,initialEstimate, initialCovariance, &ddfm);
	for(int i = 0; i < 10; i++){
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),1+0.1*i,1e-9);
		ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),1+0.01*i,1e-9);
		ASSERT_NEAR(filter.getCurrentTime(),i,1e-9);
		filter.predict(1);
	}

	for(int i = 0; i < 10; i++){
		vectorDouble measurement(i);	
		
		double sk = 1*filter.getCurrentCovariance().getSystemValue()*1 + 0.01;
		double kk = filter.getCurrentCovariance().getSystemValue()*(1.0/sk);
		double curEst = filter.getCurrentEstimate().getSystemValue(); 

		filter.update(measurement);
		
		ASSERT_NEAR(filter.getCurrentEstimate().getSystemValue(),curEst+kk*(double(i) - curEst),1e-9);
		//ASSERT_NEAR(filter.getCurrentCovariance().getSystemValue(),1+0.01*i,1e-9);
		//ASSERT_NEAR(filter.getCurrentTime(),i,1e-9);
	}
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
