#include <gtest/gtest.h> 
#include "../include/mathWrapper/double.h"
#include "../include/mathWrapper/boost.h"
#include "../include/Eigen.h"
#include "../include/mathWrapper/eigen.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


#include <iostream>

TEST(discreteDiscreteKalmanFilter, test){
	int totalCount = 10000;
	
	vectorDouble * vdR_arr = new vectorDouble[totalCount];
	Eigen::VectorXd * vesR_arr = new Eigen::VectorXd[totalCount];
	vectorBoost * vbR_arr = new vectorBoost[totalCount];
	vectorEigen * veR_arr = new vectorEigen[totalCount];

	//covariance: [1.00000,0.20000;0.20000,2.00000]
	//chol(covariance): [1.0000,0.2000;0.0000,1.4000] 

	Eigen::MatrixXd cholCov(2,2);cholCov<<1.0,0.2,0.0,1.4;Eigen::VectorXd mean1(2);mean1<<2.0,1.3;
	boost::numeric::ublas::matrix<double> cholCov2(2,2);boost::numeric::ublas::vector<double> mean2(2);
	mean2(0) = 2.0; mean2(1) = 1.3;
	cholCov2(0,0)=1;cholCov2(0,1)=0.2;cholCov2(1,0)=0.0;cholCov2(1,1)=1.4;
	matrixEigen cholCov3(cholCov);vectorEigen mean3(mean1);
	for(int i = 0; i < totalCount; i++){
		vectorDouble vdR = vectorDouble::randomVector(1);
		vdR_arr[i] = vdR*matrixDouble(0.1) + vectorDouble(5.0);

		Eigen::VectorXd vesR = Eigen::VectorXd::randomVector(2); 
		vesR_arr[i] = cholCov*vesR + mean1;

		vectorBoost vbR = vectorBoost::randomVector(2);
		vbR_arr[i] = matrixBoost(cholCov2)*vbR + mean2;
	
		vectorEigen veR = vectorEigen::randomVector(2);
		veR_arr[i] = cholCov3*veR + mean3;
	}
	
	vectorDouble vdSum(0);
	Eigen::VectorXd evsSum(2);evsSum<<0.0,0.0;
	boost::numeric::ublas::vector<double> zero(2);zero(0)=0;zero(1)=0;
	vectorBoost vbSum(zero);
	vectorEigen evSum(evsSum);
	for(int i = 0; i < totalCount; i++){
		vdSum = vdSum + vdR_arr[i];	
		evsSum = evsSum + vesR_arr[i];
		vbSum = vbSum + vbR_arr[i];
		evSum = evSum + veR_arr[i];
	}
	double vdSumMean = vdSum.getSystemValue()/double(totalCount);
	ASSERT_NEAR((vdSumMean)/5.0,1.0,0.001);
	
	double vdCovariance = 0;

	for(int i = 0; i < totalCount; i++){
		vdCovariance+=(vdR_arr[i].getSystemValue()-vdSumMean)*(vdR_arr[i].getSystemValue() - vdSumMean);
	}
	vdCovariance = vdCovariance/totalCount;
	ASSERT_NEAR(vdCovariance/0.01,1,0.01);




}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
