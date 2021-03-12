import java.io.IOException;
import java.util.*;

public class KalmanFilter{
	
	public Matrix[] predictPhase( Function f, double time,  Matrix P_0_sqrt, Matrix x, Matrix Q_root){
		int state_count = x.rowCount();
		Matrix estimate = f.next(x,time);
		Matrix jacobian = f.jacobian(x,time);

		Matrix tmp = Matrix.zero(state_count, state_count*2);
		tmp.setSubmatrix(0,0,state_count,state_count,Matrix.multiply(jacobian,P_0_sqrt));
		tmp.setSubmatrix(0,state_count,state_count,state_count,Q_root);
		
		Matrix QR_2[] = QR.QR(tmp.transpose());
		Matrix Q = QR_2[0];
		Matrix R = QR_2[1];
		
		Matrix covariance_sqrt = R.transpose();
		covariance_sqrt = covariance_sqrt.getSubmatrix(0,0,state_count,state_count);
		
		Matrix result[] = new Matrix[2];
		result[0] = estimate;
		result[1] = covariance_sqrt;
		return result;
	}

	public Matrix[] updatePhase( Matrix R_root, Matrix P_root, Matrix C, Matrix estimate, Matrix measurement){
		int measurement_count = C.rowCount();
		int state_count = estimate.rowCount();

		Matrix tmp = Matrix.zero(state_count+measurement_count,state_count+measurement_count);
		tmp.setSubmatrix(0,0,measurement_count,measurement_count, R_root);
		tmp.setSubmatrix(0,measurement_count,measurement_count,state_count,Matrix.multiply(C,P_root));
		tmp.setSubmatrix(measurement_count,measurement_count,state_count,state_count,P_root);

		Matrix QR_2[] = QR.QR(tmp.transpose());
		Matrix Q = QR_2[0];
		Matrix R = QR_2[1];
		R = R.transpose();

		Matrix X = R.getSubmatrix(0,0,measurement_count,measurement_count);
		Matrix Y = R.getSubmatrix(measurement_count,0,state_count,measurement_count);
		Matrix Z = R.getSubmatrix(measurement_count,measurement_count,state_count,state_count);

		Matrix estimate_next = Matrix.plus(estimate,Matrix.multiply(Y,Gaussian.solve(X,Matrix.minus(measurement,Matrix.multiply(C,estimate)))));
		Matrix covariance_sqrt = Z;

		Matrix result[] = new Matrix[2];
		result[0] = estimate_next;
		result[1] = covariance_sqrt;
		return result;
	}

	public Matrix[][] ddekf( Function f, double dt_between_measurements, double start_time, int state_count, int sensor_count, int measurement_count, Matrix C, Matrix Q_root, Matrix R_root, Matrix P_0_root, Matrix x_0, Matrix[] measurements){
		Matrix x_km1_p = x_0;
		Matrix P_root_km1_p = P_0_root;

		Matrix[] estimates = new Matrix[measurement_count+1];
		Matrix[] covariances = new Matrix[measurement_count+1];

		estimates[0] = x_km1_p;//TODO copy?
		covariances[0] = P_root_km1_p;

		double current_time = start_time;

		for(int k=0; k<measurement_count; k++){
			Matrix resPredict[] = predictPhase(f,current_time,P_root_km1_p,x_km1_p,Q_root);
			Matrix x_k_m = resPredict[0];
			Matrix P_root_km = resPredict[1];

			Matrix resUpdate[] = updatePhase(R_root,P_root_km,C,x_k_m,measurements[k]);
			
			x_km1_p = resUpdate[0];
			P_root_km1_p = resUpdate[1];

			current_time = current_time + dt_between_measurements;

			estimates[k+1] = x_km1_p;
			covariances[k+1] = Matrix.multiply(P_root_km1_p,P_root_km1_p.transpose());
		}

		Matrix result [][] = new Matrix[2][];		
		result[0] = estimates;
		result[1] = covariances;
		return result;
	}
}
