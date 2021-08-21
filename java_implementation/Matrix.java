import java.io.IOException;
import java.util.*;
import java.lang.Math;
import java.util.Random;

public class Matrix{
	private int rowCount;
	private int colCount;
	private double[][] data;

	public int rowCount(){
		return rowCount;
	}

	public int colCount(){
		return colCount;
	}

	public double get(int rowIndex, int colIndex){
		if(rowIndex < 0 || rowIndex >= this.rowCount()){
			throw new RuntimeException("rowIndex out of bounds for get");
		}
		if(colIndex < 0 || colIndex >= this.colCount()){
			throw new RuntimeException("colIndex out of bounds for get");
		}
		return data[rowIndex][colIndex];
	}

	public void set(int rowIndex, int colIndex, double value){
		if(rowIndex < 0 || rowIndex >= this.rowCount()){
			throw new RuntimeException("rowIndex out of bounds for set");
		}
		if(colIndex < 0 || colIndex >= this.colCount()){
			throw new RuntimeException("colIndex out of bounds for set");
		}
		data[rowIndex][colIndex] = value; 
	}

	public Matrix getSubmatrix(int row, int column, int rowCount, int colCount){
		if(row < 0 || row >= this.rowCount()){
			throw new RuntimeException("row out of bounds for getSubmatrix");
		}
		if(column < 0 || column >= this.colCount()){
			throw new RuntimeException("column out of bounds for getSubmatrix");
		}
		if( row + rowCount > this.rowCount() ){
			throw new RuntimeException("rowCount out of bounds for getSubmatrix");
		}
		if( column + colCount > this.colCount() ){
			throw new RuntimeException("colCount out of bounds for getSubmatrix");
		}

		Matrix result = new Matrix(rowCount,colCount);
		for(int i=0;i<rowCount;i++){
			for(int j=0;j<colCount;j++){
				result.set(i,j,this.data[row+i][column+j]);
			}
		}
		return result;
	}

	public void setSubmatrix(int row, int column, int rowCount, int colCount, Matrix matrix){
		if(row < 0 || row >= this.rowCount()){
			throw new RuntimeException("row out of bounds for setSubmatrix");
		}
		if(column < 0 || column >= this.colCount()){
			throw new RuntimeException("column out of bounds for setSubmatrix");
		}
		if( row + rowCount > this.rowCount() ){
			throw new RuntimeException("rowCount out of bounds for setSubmatrix");
		}
		if( column + colCount > this.colCount() ){
			throw new RuntimeException("colCount out of bounds for setSubmatrix");
		}
		if(matrix.colCount() != colCount){
			throw new RuntimeException("matrix.colCount() not equal to colCount for setSubmatrix");
		}
		if(matrix.rowCount() != rowCount){
			throw new RuntimeException("matrix.rowCount() not equal to rowCount for setSubmatrix");
		}

		for(int i = 0; i<rowCount; i++){
			for(int j = 0; j<colCount; j++){
				this.data[i+row][j+column] = matrix.get(i,j);
			}
		}
	}

	public boolean equals( Matrix matrix ){
		if(this.rowCount() != matrix.rowCount() || this.colCount() != matrix.colCount()){
			return false;	
		}
		else{
			for(int i = 0; i < this.rowCount(); i++){
				for(int j = 0; j < this.colCount(); j++){
					if(this.get(i,j)!=matrix.get(i,j)){
						return false;
					}
				}
			}
			return true;
		}
	}

	public Matrix transpose(){
		Matrix result = new Matrix(this.colCount(),this.rowCount());
		for(int i=0; i<this.colCount(); i++){
			for(int j=0; j<this.rowCount(); j++){
				result.set(i,j,this.get(j,i));
			}
		}
		return result;
	}

	public Matrix(double[][] data){
		this.rowCount = data.length;
		this.colCount = data[0].length;
		for(int i=0; i<this.rowCount; i++){
			if(data[i].length != this.colCount){
				throw new RuntimeException("data matrix not properly formatted");
			}
		}
		this.data = new double[this.rowCount][this.colCount];

		for(int i=0; i<this.rowCount; i++){
			for(int j=0; j<this.colCount; j++){
				this.data[i][j] = data[i][j];
			}
		}
	}

	public Matrix( Matrix matrix ){
		this.colCount = matrix.colCount();
		this.rowCount = matrix.rowCount();
		this.data = new double[matrix.rowCount()][matrix.colCount()];
		for(int i = 0; i < matrix.rowCount(); i++){
			for(int j = 0; j < matrix.colCount; j++){
				this.set(i,j,matrix.get(i,j));
			}
		}
	}

	public Matrix( int rows, int columns){
		this.rowCount = rows;
		this.colCount = columns;
		this.data = new double[this.rowCount][this.colCount];
		for(int i=0; i<rowCount; i++){
			for(int j=0; j<colCount; j++){
				data[i][j] = 0;
			}
		}
	}

	public static Matrix zero(int rows, int columns){
		Matrix result = new Matrix(rows,columns);
		for(int i = 0; i < result.rowCount(); i++){
			for(int j=0; j<result.colCount();j++){
				result.set(i,j,0);
			}
		}
		return result;
	}

	public static Matrix identity(int rows, int columns){
		Matrix result = Matrix.zero(rows,columns);		
		for(int i=0; i<result.colCount(); i++){
			if(i==result.rowCount){
				break;
			}
			result.set(i,i,1);
		}

		return result;
	}

	public static Matrix random(int rows, int columns){
		Random generate = new Random();
		Matrix result = new Matrix(rows,columns);
		for(int i=0; i<result.rowCount(); i++){
			for(int j=0; j<result.colCount();j++){
				result.set(i,j,generate.nextGaussian());	
			}
		}
		return result;
	}

	public static Matrix scalarMultiply(double scalar, Matrix matrix){
		Matrix result = new Matrix(matrix.rowCount(),matrix.colCount());
		for(int i=0; i<result.rowCount(); i++){
			for(int j=0; j<result.colCount(); j++){
				result.set(i,j,scalar*matrix.get(i,j));
			}
		}
		return result;
	}

	public static Matrix plus(Matrix first,Matrix second){
		if(first.rowCount()!=second.rowCount() || first.colCount()!=second.colCount()){
			throw new RuntimeException("matrix dimensions do not match up in plus");
		}

		Matrix result = new Matrix(first.rowCount(),first.colCount());
		for(int i=0; i<first.rowCount(); i++){
			for(int j=0; j<first.colCount(); j++){
				result.set(i,j,first.get(i,j) + second.get(i,j));
			}
		}
		return result;
	}

	public static Matrix minus(Matrix first,Matrix second){
		if(first.rowCount()!=second.rowCount() || first.colCount()!=second.colCount()){
			throw new RuntimeException("matrix dimensions do not match up in minus");
		}
		Matrix result = new Matrix(first.rowCount(),first.colCount());
		for(int i=0; i<first.rowCount(); i++){
			for(int j=0; j<first.colCount(); j++){
				result.set(i,j,first.get(i,j) - second.get(i,j));
			}
		}
		return result;
	}


	public static Matrix multiply(Matrix first, Matrix second){
		if(first.colCount() != second.rowCount()){
			throw new RuntimeException( "matrix dimensions do not match up in multiply");
		}
		Matrix result = Matrix.zero(first.rowCount(),second.colCount());
		for(int i=0; i<first.rowCount(); i++){
			for(int j=0; j<second.colCount(); j++){
				for(int k=0; k<second.rowCount() ; k++){
					result.set(i,j,result.get(i,j) + first.get(i,k)*second.get(k,j));
				}
			}
		}
		return result;
	}

	public static double dot(Matrix first, Matrix second){
		if(first.rowCount()!=second.rowCount() || first.colCount()!=second.colCount()){
			throw new RuntimeException("matrix dimensions do not match up in dot");
		}
		double result = 0;
		for(int i=0; i<first.rowCount(); i++){
			for(int j=0; j<first.colCount(); j++){
				result = result + first.get(i,j)*second.get(i,j);
			}
		}
		return result;
	}

	//frobenius norm
	public static double norm(Matrix first){
		return Math.sqrt(Matrix.dot(first,first));
	}

	public void flipColumns(){
		int colCount = this.colCount();
		int rowCount = this.rowCount();
		for(int col=0; col<colCount/2; col++){
			for(int row=0; row<rowCount; row++){
				double tmp = this.get(row,col);
				this.set(row,col,this.get(row,colCount-col-1));
				this.set(row,colCount-col-1,tmp);
			}
		}
	}

	public void flipRows(){
		int rowCount = this.rowCount();
		int colCount = this.colCount();
		for(int row=0; row<rowCount/2;row++){
			for(int col=0; col<colCount; col++){
				double tmp = this.get(row,col);
				this.set(row,col,this.get(rowCount-row-1,col));
				this.set(rowCount-row-1,col,tmp);
			}
		}
	}

	public void swapRows(int row1, int row2, int startCol, int endCol){
		int rowCount = this.rowCount();
		int colCount = this.colCount();
		if(row1<0 || row2<0){
			throw new RuntimeException("swapRows row indices have to be greater than zero");
		}
		if( !( row1<rowCount && row2<rowCount )){
			throw new RuntimeException("swapRows row indeces have to be less rowCount");
		}
		if(startCol<0 || endCol<0 ){
			throw new RuntimeException("swapRows col indices have to be greater than zero");
		}
		if( !(startCol <= endCol ) ){
			throw new RuntimeException("swapRows col start has to be less or equal to end");
		}
		if( ( endCol > colCount ) ){
			throw new RuntimeException("swapRows ");
		}
		
		for(int col=startCol; col<endCol; col++){
			double tmp = this.get(row1,col);
			this.set(row1,col,this.get(row2,col));
			this.set(row2,col,tmp);
		}
	}

	public static Matrix normalRandom(int rows, int cols){
		if(rows<=0){
			throw new RuntimeException("normalRandom needs positive amount of rows");
		}
		if(cols<=0){
			throw new RuntimeException("normalRandom needs positive amount of cols");
		}

		Random r = new Random();

		Matrix result =Matrix.zero(rows,cols);
		for(int row=0; row<rows; row++){
			for(int col=0; col<cols; col++){
				result.set(row,col,r.nextGaussian());
			}
		}
		return result;
	}

	public void print(){
		for(int i = 0; i < this.rowCount(); i++){
			for(int j = 0; j < this.colCount(); j++){
				System.out.print(this.data[i][j]);	
				System.out.print(" ");
			}
			System.out.println();
		}
	}
}
