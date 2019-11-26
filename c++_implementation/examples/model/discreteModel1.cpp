//! [example1]
#include <iostream>
#include "../../include/model.h"


int main(){
	//class will model an object travelling at constant speed of 0.1 per time unit    
	class discreteSimpleExample: public discreteModel<double>{
		//implementing the pure virtual function
		public:
			double function(const double & val, const int time) const override{
				return val + 0.1;
                	}
	};	

	discreteSimpleExample dse;//actual instance of class
	double timeEvaluatedAt = 10;
	double valueAtTime = 1;

	std::cout<<dse.function(valueAtTime,timeEvaluatedAt)<<std::endl;
	/* produces the output:
	 * 1.1
	 */
};
//! [example1]
