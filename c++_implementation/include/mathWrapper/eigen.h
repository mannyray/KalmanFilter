#ifndef VECTOREIGEN
#define VECTOREIGEN
#include "../Eigen.h"
#include <stdexcept>


class vectorEigen{
	private:
		Eigen::VectorXd value;
	public:
		/*
		 * Reason this empty constructor is allowed is because
		 * of the discreteSolver<K>/continuousSolver<K> class that allows
		 * a constructor with just the discreteModel<K>/continuousModel<K>
		 * and not initial conditions which means the private currentState
		 * variables within those classes end up defined with an empty constructor.
		 */
		vectorEigen(){}

		vectorEigen(Eigen::VectorXd value):value(value){}

		Eigen::VectorXd getSystemValue() const{
			return value;
		}

		static vectorEigen randomVector(int dimension,uint32_t seed){
			std::mt19937 generator(seed);
			std::normal_distribution<double> distribution(0.0,1.0);

			Eigen::VectorXd tmp(dimension);
			for(int i = 0; i < dimension; i++){
				tmp[i] = distribution(generator);
			}
			return vectorEigen(tmp);
		}


};

vectorEigen operator-(const vectorEigen &lhs,const vectorEigen &rhs){
	vectorEigen res(lhs.getSystemValue() - rhs.getSystemValue());
	return res;
};

vectorEigen operator+(const vectorEigen &lhs,const vectorEigen &rhs){
	vectorEigen res(lhs.getSystemValue() + rhs.getSystemValue());
	return res;
};

std::ostream& operator<<(std::ostream& stream, const vectorEigen & vec) {
	Eigen::VectorXd v = vec.getSystemValue();
	stream<<v;
        return stream;
}




class matrixEigen{
        private:
		Eigen::MatrixXd value;

        public:
		matrixEigen(){}
	

                matrixEigen(Eigen::MatrixXd value):value(value){}


                static matrixEigen Identity(int dimension,int dimension2){
                        return matrixEigen(Eigen::MatrixXd::Identity(dimension,dimension));
                }

                matrixEigen transpose(){
                        matrixEigen t(value.transpose());
                        return t;
                }

                matrixEigen inverse(){
                        matrixEigen t(value.inverse());
                        return t;
                }

		Eigen::MatrixXd getSystemValue() const{
                        return value;
                }
};

matrixEigen operator*(const matrixEigen &lhs, const matrixEigen &rhs){
        matrixEigen res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}

matrixEigen operator+(const matrixEigen &lhs,const matrixEigen &rhs){
        matrixEigen res(lhs.getSystemValue() + rhs.getSystemValue());
        return res;
};

matrixEigen operator-(const matrixEigen &lhs,const matrixEigen &rhs){
        matrixEigen res(lhs.getSystemValue() - rhs.getSystemValue());
        return res;
};

std::ostream& operator<<(std::ostream& stream, const matrixEigen & mat) {
	Eigen::MatrixXd v = mat.getSystemValue();
	stream<<v;
        return stream;
}




vectorEigen operator*(const matrixEigen &lhs, const vectorEigen &rhs){
	vectorEigen res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}

vectorEigen operator*(const vectorEigen &lhs, const matrixEigen &rhs){
	vectorEigen res(lhs.getSystemValue()*rhs.getSystemValue());
	return res;
}
#endif
