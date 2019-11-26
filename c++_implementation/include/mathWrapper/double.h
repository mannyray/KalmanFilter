#ifndef VECTORDOUBLE
#define VECTORDOUBLE

#include <stdexcept>
#include <random>


class vectorDouble{
	private:
		double value;
	public:
		vectorDouble(){
			value = 0;
		}

		vectorDouble(double value):value(value){}

		double getSystemValue() const{
			return value;
		}


		static vectorDouble randomVector(int dimension,uint32_t seed){
			std::mt19937 generator(seed);
			std::normal_distribution<double> distribution(0.0,1.0);
			double value = distribution(generator); 
			return vectorDouble(value);
		}
};




vectorDouble operator-(const vectorDouble &lhs,const vectorDouble &rhs){
	vectorDouble res(lhs.getSystemValue() - rhs.getSystemValue());
	return res;
};

vectorDouble operator+(const vectorDouble &lhs,const vectorDouble &rhs){
	vectorDouble res(lhs.getSystemValue() + rhs.getSystemValue());
	return res;
};

std::ostream& operator<<(std::ostream& stream, const vectorDouble & vec) {
        double v = vec.getSystemValue();
	stream<<v;
        return stream;
}




class matrixDouble{
        private:
                double value;

        public:
                matrixDouble(){
                        value = 0;
                }

                matrixDouble(double value):value(value){}


                static matrixDouble Identity(int dimension, int dimension2){
                        return matrixDouble(1);
                }

                matrixDouble transpose(){
                        matrixDouble t(value);
                        return t;
                }

                matrixDouble inverse(){
                        if(value == 0){
                                throw std::overflow_error("Divide by zero exception");
                        }
                        matrixDouble t(1.0/value);
                        return t;
                }

                double getSystemValue() const{
                        return value;
                }
};

matrixDouble operator*(const matrixDouble &lhs, const matrixDouble &rhs){
        matrixDouble res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}

matrixDouble operator+(const matrixDouble &lhs,const matrixDouble &rhs){
        matrixDouble res(lhs.getSystemValue() + rhs.getSystemValue());
        return res;
};

matrixDouble operator-(const matrixDouble &lhs,const matrixDouble &rhs){
        matrixDouble res(lhs.getSystemValue() - rhs.getSystemValue());
        return res;
};

std::ostream& operator<<(std::ostream& stream, const matrixDouble & vec) {
        double v = vec.getSystemValue();
	stream<<v;
        return stream;
}




vectorDouble operator*(const matrixDouble &lhs, const vectorDouble &rhs){
	vectorDouble res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}

vectorDouble operator*(const vectorDouble &lhs, const matrixDouble &rhs){
	vectorDouble res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}





#endif
