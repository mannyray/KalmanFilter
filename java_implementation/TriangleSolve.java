import java.io.IOException;
import java.util.*;

public class TriangleSolve{

	public static Matrix lowerTriangleSolve( Matrix A, Matrix b ){
		int rowCount = A.rowCount();
		int colCount = A.colCount();
		if(rowCount!=colCount){
			throw new RuntimeException("lowerTriangleSolve expecting square matrix");
		}
		if(colCount!=b.rowCount()){
			throw new RuntimeException("lowerTriangleSolve expecting matrix and vector to match");
		}
		if(b.colCount()!=1){
			A.print();
			b.print();
			throw new RuntimeException("lowerTriangleSolve expecting a vector and not matrix");
		}

		for(int i=0; i<rowCount; i++){
			if(A.get(i,i)==0){
				throw new RuntimeException("lowerTriangleSolve full rank matrix");
			}
		}

		Matrix x = Matrix.zero(rowCount,1);
		for(int row=0; row<rowCount; row++){
			double right_sum=b.get(row,0);
			double left_non_x_sum=0;
			for(int col=0; col<colCount; col++){
				left_non_x_sum = left_non_x_sum + x.get(col,0)*A.get(row,col);
			}
			x.set(row,0,(right_sum-left_non_x_sum)/A.get(row,row));
		}
		return x;	
	}	

	public static Matrix upperTriangleSolve( Matrix A, Matrix b ){
		Matrix local_A = new Matrix(A);
		Matrix local_b = new Matrix(b);

		//flip columns
		local_A.flipColumns();

		//flip rows
		local_A.flipRows();
		local_b.flipRows();

		Matrix x = lowerTriangleSolve( local_A, local_b );
		x.flipRows();

		return x;
	}
}
