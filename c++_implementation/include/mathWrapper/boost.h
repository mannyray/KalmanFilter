#ifndef VECTORBOOST
#define VECTORBOOST

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <lapacke.h>

#include <random>

#include <iostream>


class matrixBoost{
        private:

                boost::numeric::ublas::matrix<double> value;
		boost::numeric::ublas::range r1;
		boost::numeric::ublas::range r2;
		boost::numeric::ublas::matrix<double>valBack;
		bool subMatrixCopy = false;

		//http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
		bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse) {
			using namespace boost::numeric::ublas;
			typedef permutation_matrix<std::size_t> pmatrix;
			// create a working copy of the input
			matrix<double> A(input);
			// create a permutation matrix for the LU-factorization
			pmatrix pm(A.size1());
			// perform LU-factorization
			int res = lu_factorize(A,pm);
			if( res != 0 ) return false;
			// create identity matrix of "inverse"
			inverse.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));
			// backsubstitute to get the inverse
			lu_substitute(A, pm, inverse);
			return true;
		}

		void setSystemValue(boost::numeric::ublas::matrix<double> v){
			value = v;
		}

        public:
                matrixBoost(){}

		matrixBoost(int rowCount, int colCount){
			value = boost::numeric::ublas::zero_matrix<double>(rowCount,colCount);
		}

                matrixBoost(boost::numeric::ublas::matrix<double> value):value(value){}

                static matrixBoost Identity(int dimension,int dimension2){
                        return matrixBoost(boost::numeric::ublas::identity_matrix<double>(dimension,dimension2));
                }

                matrixBoost transpose(){
                        matrixBoost t(boost::numeric::ublas::trans(value));
                        return t;
                }

                matrixBoost inverse(){
			boost::numeric::ublas::matrix<double> inverse(value.size1(),value.size2());
			InvertMatrix(value,inverse);
                        matrixBoost t(inverse);
                        return t;
                }

                boost::numeric::ublas::matrix<double> getSystemValue() const{
                        return value;
                }
		
		matrixBoost & block(int i1, int i2, int l1, int l2){
			if(subMatrixCopy){
				std::swap(valBack,value);
			}
			subMatrixCopy = true;
			r1 = boost::numeric::ublas::range(i1,i1+l1);
			r2 = boost::numeric::ublas::range(i2,i2+l2);
			valBack = boost::numeric::ublas::project(value,r1,r2);
			std::swap(valBack,value);
			return *this;
		}

		matrixBoost &operator=( const matrixBoost other){
			if(subMatrixCopy){
				subMatrixCopy = false;
				boost::numeric::ublas::project(valBack,r1,r2) = other.getSystemValue();
				std::swap(valBack,value);
			}else{
				value = other.getSystemValue();
			}
			return *this;//better be no chaining of equal signs
		}

		static void lowerTriangulate(matrixBoost & val){
			boost::numeric::ublas::matrix<double> in= val.getSystemValue();
			int rows = in.size1();
			int cols = in.size2();
			for(int i = 0; i < rows; i++){
				int length = cols - i;
				boost::numeric::ublas::vector<double> v = (row(in,i));
				v = boost::numeric::ublas::subrange(v,i,cols);
				double mew = boost::numeric::ublas::norm_2(v);
				if(mew!=0){
					double beta  = v(0)  + (-2*(v(0) < 0)+1)*mew;
					v = v/beta;
				}
				v(0) = 1;
				boost::numeric::ublas::range t1(i,rows);
				boost::numeric::ublas::range t2(i,cols);
				boost::numeric::ublas::vector<double> w = -2.0/
					(boost::numeric::ublas::norm_2(v)*boost::numeric::ublas::norm_2(v))*
					(prod(project(in,t1,t2),v));
				
				boost::numeric::ublas::matrix<double> vM(v.size(),1);
				for(int i = 0; i < v.size(); i++){
					vM(i,0) = v(i); 
				}
				boost::numeric::ublas::matrix<double> wM(w.size(),1);
				for(int i = 0; i < w.size(); i++){
					wM(i,0) = w(i);
				}
				project(in,t1,t2) = project(in,t1,t2) +  prod(wM,boost::numeric::ublas::trans(vM));	
			}
			val.setSystemValue(in);
		}
};

class vectorBoost{
        private:
		boost::numeric::ublas::vector<double> value;
        public:
                vectorBoost(){
                }

                vectorBoost(boost::numeric::ublas::vector<double> value):value(value){}

                boost::numeric::ublas::vector<double> getSystemValue() const{
                        return value;
                }
	
		static vectorBoost randomVector(int dimension,uint32_t seed){
			std::mt19937 generator(seed);
			std::normal_distribution<double> distribution(0.0,1.0);
			boost::numeric::ublas::vector<double> tmp(dimension);
			for(int i = 0; i < dimension; i++){
				tmp(i) = distribution(generator);
			}
			return vectorBoost(tmp);
		}

		static vectorBoost solve(matrixBoost A, vectorBoost b){
			boost::numeric::ublas::matrix<double> Alocal = A.getSystemValue();
			boost::numeric::ublas::vector<double> blocal = b.getSystemValue();
			boost::numeric::ublas::permutation_matrix<size_t> pm(Alocal.size1());

			boost::numeric::ublas::lu_factorize(Alocal,pm);
			boost::numeric::ublas::lu_substitute(Alocal, pm, blocal);
			return vectorBoost(blocal);
		}	
};


vectorBoost operator-(const vectorBoost &lhs,const vectorBoost &rhs){
        vectorBoost res(lhs.getSystemValue() - rhs.getSystemValue());
        return res;
};

vectorBoost operator+(const vectorBoost &lhs,const vectorBoost &rhs){
        vectorBoost res(lhs.getSystemValue() + rhs.getSystemValue());
        return res;
};

std::ostream& operator<<(std::ostream& stream, const vectorBoost & vec) {
	boost::numeric::ublas::vector<double> v = vec.getSystemValue(); 
	for(int i = 0; i < v.size()-1; i++){
		stream<<double(v(i))<<'\n';
	}
	stream<<v(v.size()-1);
	return stream;
}




matrixBoost operator*(const matrixBoost &lhs, const matrixBoost &rhs){
        matrixBoost res(boost::numeric::ublas::prod(lhs.getSystemValue(),rhs.getSystemValue()));
        return res;
}

matrixBoost operator+(const matrixBoost &lhs,const matrixBoost &rhs){
        matrixBoost res(lhs.getSystemValue() + rhs.getSystemValue());
        return res;
};

matrixBoost operator-(const matrixBoost &lhs,const matrixBoost &rhs){
        matrixBoost res(lhs.getSystemValue() - rhs.getSystemValue());
        return res;
};

std::ostream& operator<<(std::ostream& stream, const matrixBoost & mat) {
	boost::numeric::ublas::matrix<double> m = mat.getSystemValue(); 
	for(int i = 0; i < m.size1(); i++){
		for(int j = 0; j < m.size2(); j++){
			stream<<m(i,j)<<" ";
		}stream<<'\n';
	}	
	return stream;
}






vectorBoost operator*(const matrixBoost &lhs, const vectorBoost &rhs){
        vectorBoost res(boost::numeric::ublas::prod(lhs.getSystemValue(),rhs.getSystemValue()));
        return res;
}

vectorBoost operator*(const vectorBoost &lhs, const matrixBoost &rhs){
        vectorBoost res(boost::numeric::ublas::prod(lhs.getSystemValue(),rhs.getSystemValue()));
        return res;
}
#endif
