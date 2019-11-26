//! [example1]
#include <Eigen/Dense>
#include <iostream>
#include "../../include/model.h"


int main(){
	//using the Eigen library to implement the exact same class (conceptually) as above	
	class odeSimpleExampleEigen: public continuousModel<Eigen::VectorXd>{
		private:
			Eigen::MatrixXd modelMatrix;
		public:
			//can add constructor to make your model more adjustable and reusable
			odeSimpleExampleEigen(double exp){
				Eigen::MatrixXd tmp(1,1);
				tmp<<exp;
				modelMatrix = tmp;
			}
			Eigen::VectorXd function(const Eigen::VectorXd & val, const double time) const override{
				return modelMatrix*val;
			}
	};

	double exponent = -0.1;
	odeSimpleExampleEigen ose2(exponent);
	
	double timeEvaluatedAt = 10;
	double valueAtTime = 1;
	Eigen::VectorXd tmp(1); tmp<<valueAtTime;

	std::cout<<ose2.function(tmp,timeEvaluatedAt) <<std::endl;
	/* produces the output:
	 * -0.1
	 */
};
//! [example1]
