#ifndef EIGENSIMPLE
#define EIGENSIMPLE

EIGEN_DEVICE_FUNC static EIGEN_STRONG_INLINE
Derived randomVector(int dimension, uint32_t seed){
	std::mt19937 generator(seed);
	std::normal_distribution<double> distribution(0.0,1.0);
	Derived res(dimension);
	for(int i = 0; i < dimension; i++){
		res[i] = distribution(generator);                            
	}
	return res;
}


void static lowerTriangulate(Eigen::MatrixBase<Derived> &in){
	for(int i = 0; i < in.rows(); i++){
		int length = in.cols() - i;
		Derived v = in.row(i).tail(length).transpose();
		double mew = v.norm();
		if(mew!=0){
			double beta  = v(0)  + (-2*(v(0) < 0)+1)*mew;	
			v = v/beta;
		}
		v(0) = 1;
		Derived w = -2.0/(v.squaredNorm())*(in.block(i,i, in.rows()-i , in.cols()-i)*v);
		in.block(i,i, in.rows()-i , in.cols()-i) = in.block(i,i, in.rows()-i , in.cols()-i) + w*v.transpose();
	}
}

//Eigen::Matrix<double, -1, 1> A,Eigen::Matrix<double, -1, -1> b){
static Derived solve(Derived A, Derived b){
	return A.colPivHouseholderQr().solve(b);	
}






#endif
