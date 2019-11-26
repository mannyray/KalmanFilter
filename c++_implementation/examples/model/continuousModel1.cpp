//! [example1]
#include <iostream>
#include "../../include/model.h"


int main(){
	//class will model one of the most simplest differential equations dy/dt = -0.1*y
	class odeSimpleExample: public continuousModel<double>{
		//implementing the pure virtual function
		public:
			double function(const double & val, const double time) const override{
				double rhs = -0.1*val;
				return rhs;
               		}
	};	

	odeSimpleExample ose;//actual instance of class
	double timeEvaluatedAt = 10;
	double valueAtTime = 1;

	std::cout<<ose.function(valueAtTime,timeEvaluatedAt)<<std::endl;
	/* produces the output:
	 * -0.1
	 */
};
//! [example1]
