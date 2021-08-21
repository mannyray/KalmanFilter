import java.io.IOException;
import java.util.*;

public class QR{

	private static double sign(double x){
		if(x < 0){
			return -1;
		}
		else{
			return 1;
		}
	}

	public static Matrix [] QR( Matrix A ){
		int rowCount = A.rowCount();
		int colCount = A.colCount();

		Matrix R = new Matrix(A);
		Matrix Q = Matrix.identity(colCount,colCount);

		for(int j=0; j<colCount; j++){
			Matrix x = R.getSubmatrix(j,j,rowCount-j,1); 
			Matrix e1 = Matrix.zero(rowCount-j,1);
			e1.set(0,0,1);
			Matrix v_k = Matrix.plus( Matrix.scalarMultiply(QR.sign(x.get(0,0))*Matrix.norm(x),e1), x);
			v_k = Matrix.scalarMultiply(1.0/Matrix.norm(v_k),v_k);

			if(rowCount < colCount ){
				for(int i=0; i<colCount; i++){
					Q.setSubmatrix(j,i,(rowCount-j),1,Matrix.minus(Q.getSubmatrix(j,i,(rowCount-j),1), Matrix.scalarMultiply(2.0,Matrix.multiply(v_k,Matrix.multiply(v_k.transpose(),Q.getSubmatrix(j,i,(rowCount-j),1))))));
				}
			}

			R.setSubmatrix(j,j,rowCount-j,colCount-j,Matrix.minus( R.getSubmatrix(j,j,rowCount-j,colCount-j), Matrix.scalarMultiply(2.0,Matrix.multiply( Matrix.multiply(v_k,v_k.transpose()), R.getSubmatrix(j,j,rowCount-j,colCount-j)))));
		}

		Matrix result[] = new Matrix[2];
		result[0] = Q.transpose();
		result[1] = R;
		return result;
	}
}
