#include <gtest/gtest.h> 
#include "../include/model.h"
#include <math.h>
#include "../include/Eigen.h" 
#include <vector>
#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/operators.hpp>


#include <iostream>

//Testing continuousSolverRK4<double> and continuousModel<double>
//---tests to make sure continuousSolverRK4<T>.solve(double t) does not
//accept values with t < 0 (can't solve back in time)
//---tests to make sure that RK4 method can solve basic heat equation PDE
//correctly according to exact solution while testing that the solved
//state and updated current time are properly stored  
TEST(continuousSolverRK4, RK4Test){

	double exponent = -0.1;
	class testModel: public continuousModel<double>{
		public:
			double function(const double & val, const double time) const override{
				double rhs = -0.1*val;
				return rhs;
			}
	};
	testModel instance;

	double initialCondition = 100;
	double initialTime = 0;
	continuousSolverRK4<double> testModelSolver(&instance,initialTime,initialCondition,0.001);

	double times[9] = {0,1,1,1,1,1,1,1,1};
	int length = 9;
	//assert time can't be negative
	double rollingSum = 0;
	for(int i=0; i<length; i++){
		rollingSum+=times[i];
		testModelSolver.solve(times[i]);
		ASSERT_NEAR(rollingSum,testModelSolver.getCurrentTime(),1e-6);
		ASSERT_NEAR(exp(exponent*rollingSum)*initialCondition,testModelSolver.getCurrentState(),1e-6);
		if(i!=length-1){continue;}
	}
	ASSERT_THROW(testModelSolver.solve(-1),std::out_of_range);
}

//CUSTOM VECTOR CLASS	
class customVector{
	private:
		std::vector<double> entries;
	public:
		int getSize() const{
			return entries.size();
		}
		
		double getValueAtIndex(int i) const{
			if(i < 0 || i >= getSize()){
				throw std::out_of_range("");
			}
			return entries[i];
		}
		void setValueAtIndex(int i, double value){
			if(i < 0 || i >= getSize()){
				throw std::out_of_range("");
			}
			entries[i] = value;
		}

		customVector(const std::vector<double> & vec){
			if(vec.size() == 0){
				throw std::out_of_range("");
			}
			entries = vec;
		}

		customVector(int size){
			if(size <= 0){
				throw std::out_of_range("");
			}
			entries = std::vector<double>(size);
		}
};


customVector operator+(const customVector &lhs,const customVector &rhs){
	if(lhs.getSize() != rhs.getSize()){
		throw std::range_error("");	
	}
	customVector res(rhs.getSize());
	for(int i = 0; i < lhs.getSize(); i++){
		res.setValueAtIndex(i,lhs.getValueAtIndex(i)+rhs.getValueAtIndex(i));
	}
	return res;
};

customVector operator*(const double &lhs, const customVector &rhs){
	customVector res(rhs.getSize());
	for(int i = 0; i < rhs.getSize(); i++){
		res.setValueAtIndex(i,rhs.getValueAtIndex(i)*lhs);
	}
	return res;
}


TEST(continuousSolverRK4,differentLibraryTest){
	class testModel: public continuousModel<customVector>{
		public:
		customVector function(const customVector & val, const double time) const{
			customVector result(2);
			result.setValueAtIndex(0,-0.1*val.getValueAtIndex(0));
			result.setValueAtIndex(1,-0.2*val.getValueAtIndex(1));
			return result;
		}
	};
	testModel instance;
	std::vector<double> t;
	t.push_back(100);t.push_back(100);
	customVector initialCondition(t);
	double initialTime = 0;
	continuousSolverRK4<customVector> testModelSolver(&instance,initialTime,initialCondition,1);

	class testModel2: public continuousModel<Eigen::VectorXd>{
		private:
			Eigen::MatrixXd modelMatrix;
		public:
			testModel2(){
				Eigen::MatrixXd tmp(2,2);
				tmp<<-0.1,0,0,-0.2;
				modelMatrix = tmp;
			}
			Eigen::VectorXd function(const Eigen::VectorXd & val, const double time) const{
				Eigen::VectorXd gg = modelMatrix*val;
				return modelMatrix*val;
			}
	};

	testModel2 instance2;;
	Eigen::VectorXd initialCondition2(2); initialCondition2<<100,100;
	instance2.function(initialCondition2,0.0);
	continuousSolverRK4<Eigen::VectorXd> testModelSolver2(&instance2,initialTime,initialCondition2,1);


/*
	class testModel3: public continuousModel<boost::numeric::ublas::vector<double>>{
		private:
			boost::numeric::ublas::matrix<double> modelMatrix;
		public:
			testModel3(){
				boost::numeric::ublas::matrix<double> tmp(2,2);
				tmp(0,0) = -0.1;
				tmp(0,1) = 0;
				tmp(1,0) = 0;
				tmp(1,1) = -0.2;
				modelMatrix = tmp;
			}
			boost::numeric::ublas::vector<double>function(boost::numeric::ublas::vector<double> & val,const double time)const{
				//modelMatrix*val;
				return prod(modelMatrix,val);
			}
	};

	testModel3 instance3 = testModel3();
	boost::numeric::ublas::vector<double> initialCondition3(2); initialCondition3(0) = 100; initialCondition3(1) = 100;
	continuousSolverRK4<boost::numeric::ublas::vector<double>> testModelSolver3(&instance3,initialTime,initialCondition3,0.001);
*/	

	double times[9] = {0,1,1,1,1,1,1,1,1};
	int length = 9;
	double rollingSum = 0;
	for(int i=0; i<length; i++){
		rollingSum+=times[i];
		testModelSolver.solve(times[i]);
		testModelSolver2.solve(times[i]);
		Eigen::VectorXd tmp2 = testModelSolver2.getCurrentState();
		double tmp3 = tmp2[0]; 
		double tmp4 = tmp2[1];
		ASSERT_NEAR(testModelSolver.getCurrentState().getValueAtIndex(0),tmp3,1e-9);
		ASSERT_NEAR(testModelSolver.getCurrentState().getValueAtIndex(1),tmp4,1e-9);
	}
	ASSERT_THROW(testModelSolver.solve(-1),std::out_of_range);
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
