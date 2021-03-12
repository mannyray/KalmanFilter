public class FunctionLinear extends Function{

	Matrix A;

	public FunctionLinear(Matrix A){
		this.A = A;
	}

        public Matrix next( Matrix x, double time){
                return Matrix.multiply(this.A,x);
        }

        public Matrix jacobian(Matrix x, double time){
                return this.A;
        }
}
