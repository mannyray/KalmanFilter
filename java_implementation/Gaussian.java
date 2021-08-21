import java.io.IOException;
import java.util.*;
import java.lang.Math;

public class Gaussian{
	
	public static Matrix [] computeGaussian( Matrix A, boolean shifts ){
		int rowCount = A.rowCount();
		int colCount = A.colCount();
		if(rowCount != colCount){
			throw new RuntimeException("computeGaussian expects a square matrix");
		}

		Matrix U = new Matrix(A);
		Matrix L = Matrix.identity(rowCount,rowCount);
		Matrix P = Matrix.identity(rowCount,rowCount);

		for(int column=0; column<colCount; column++){
			if(shifts){
				double maxVal = Math.abs(U.get(column,column));
				int maxIndex = column;

				for(int row=column; row<rowCount; row++){
					if(Math.abs(U.get(row,column)) > maxVal){
						maxVal = Math.abs(U.get(row,column));
						maxIndex = row;
					}
				}
				if(maxIndex!=column){
					U.swapRows(column,maxIndex,column,colCount);
					L.swapRows(column,maxIndex,0,column);
					P.swapRows(column,maxIndex,0,colCount);
				}
			}

			for(int row=column+1; row<colCount; row++){
				L.set(row,column,U.get(row,column)/U.get(column,column));
				for(int column2=column; column2<rowCount; column2++){
					U.set(row,column2,U.get(row,column2)-L.get(row,column)*U.get(column,column2));
				}
			}
		}

		if(shifts){
			Matrix result [] = new Matrix[3];
			result[0] = P;
			result[1] = L;
			result[2] = U;
			return result;
		}
		else{
			Matrix result [] = new Matrix[2];
			result[0] = L;
			result[1] = U;
			return result;
		}
	}

	public static Matrix[] computeGaussian( Matrix A ){
		return computeGaussian(A,true);
	}

	public static Matrix solve(Matrix A, Matrix b){
		Matrix arr [] = computeGaussian( A );
		Matrix P = arr[0];
		Matrix L = arr[1];
		Matrix U = arr[2];

		return TriangleSolve.upperTriangleSolve(U,TriangleSolve.lowerTriangleSolve(L,Matrix.multiply(P,b)));
	}
}
