//! [example1]
#include <Eigen/Dense>
#include <iostream>
#include "../../include/model.h"


int main(){
	//using the Eigen library to implement the exact same class (conceptually) as above	
	class discreteSimpleExampleEigen: public discreteModel<Eigen::VectorXd>{
		private:
			Eigen::VectorXd velocity;
		public:
			//can add constructor to make your model more adjustable and reusable
			discreteSimpleExampleEigen(double vel){
				Eigen::VectorXd tmp(1);
				tmp<<vel;
				velocity = tmp;
			}
			Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
				return val+velocity;
			}
	};

	double velocity = 0.1;
	discreteSimpleExampleEigen dse2(velocity);
	
	int timeEvaluatedAt = 10;
	double valueAtTime = 1;
	Eigen::VectorXd tmp(1); tmp<<valueAtTime;

	std::cout<<dse2.function(tmp,timeEvaluatedAt)<<std::endl;
	/* produces the output:
	 * 1.1
	 */
};
//! [example1]
