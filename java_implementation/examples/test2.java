import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

public class test2{

	static class FunctionLinear extends Function{
	
		Matrix A = Matrix.zero(2,2);

		public FunctionLinear(){
			A.set(0,0,-1);A.set(1,0,-0.1);A.set(0,1,0.2);A.set(1,1,-1);
			A = Matrix.plus(Matrix.identity(2,2),Matrix.scalarMultiply(0.002,A));
		}

		public Matrix next( Matrix x, double time){
		        return Matrix.multiply(A,x);
		}

		public Matrix jacobian(Matrix x, double time){
		        return A;
		}
	}

        public static void main(String[] args) {
        	int steps = 10000;
        
		Matrix[] measurements = new Matrix[steps];
		Matrix[] states = new Matrix[steps];
				
		Matrix Q = Matrix.zero(2,2);
		Q.set(0,0,0.000000002);Q.set(1,0,0);Q.set(0,1,0);Q.set(1,1,0.000000002);
		Matrix Q_root = Matrix.zero(2,2);//root is cholesky transpose
		Q_root.set(0,0,0.000044721359549995793928);Q_root.set(1,1,0.000044721359549995793928);


		Matrix R = Matrix.zero(1,1);
		R.set(0,0,5);
		Matrix R_root = Matrix.zero(1,1);
		R_root.set(0,0,2.2360679774997896964);
		Matrix C = Matrix.zero(1,2);
		C.set(0,0,1);C.set(0,1,0.1);
		Matrix x_0 = Matrix.zero(2,1);
		x_0.set(0,0,100);x_0.set(1,0,80);
		
		Matrix x = new Matrix(x_0);
		FunctionLinear fl = new FunctionLinear();
		
		//generate the measurements to simulate real world data
		for(int i=0; i<steps; i++){
			x = Matrix.plus(fl.next(x,i),Matrix.multiply(Q_root,Matrix.normalRandom(2,1)));
			states[i] = x;
			Matrix measurement = Matrix.plus(Matrix.multiply(C,x), Matrix.multiply(R_root,Matrix.normalRandom(1,1)));
			measurements[i] = measurement;
		}
		
		Matrix P_0 = Matrix.identity(2,2);
		Matrix P_0_root = Matrix.identity(2,2);
		
		Matrix[][] result = KalmanFilter.ddekf(new FunctionLinear(),1,0,2,1,steps,C,Q_root,R_root,P_0_root,x_0,measurements);

		//save measurements, estimates and covariances for verification with python_verify.py
		try {
			FileWriter myWriterEstimates = new FileWriter("estimates_java.txt");
			FileWriter myWriterMeasurements = new FileWriter("measurements_java.txt");
			FileWriter myWriterCovariances = new FileWriter("covariances_java.txt");
			for(int i=0; i<steps; i++){
				myWriterEstimates.write(result[0][i].get(0,0) + "," + result[0][i].get(1,0)+"\n");
				myWriterMeasurements.write(measurements[i].get(0,0)+"\n");
				myWriterCovariances.write(result[1][i].get(0,0)+","+result[1][i].get(0,1)+","+result[1][i].get(1,0)+","+result[1][i].get(1,1)+"\n");
			}
			myWriterEstimates.close();
			myWriterMeasurements.close();
			myWriterCovariances.close();
		} catch (IOException e) {}
	}
}
