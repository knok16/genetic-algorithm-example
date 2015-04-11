package com.knok16.test;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Scanner;

import com.knok16.GeneticFramework.FitnessFunction;

public class MyFitnessFunction implements FitnessFunction {

	public long prepareTime = 0;
	public long checkTime = 0;
	public long sortingTime = 0;
	
	private static final int BIT_TO_INT = 8;
	private int pathLength = 0;
	private int[] path = null;
	private int[] seq = null;
	private int vertexCount;
	private int[][] matrix;
	
	
	public MyFitnessFunction(String filename) throws FileNotFoundException{
		super();
		Scanner in = new Scanner(new FileReader(filename));
		this.vertexCount = in.nextInt();
		this.matrix = new int[this.vertexCount][this.vertexCount];
		for(int i=0;i<this.vertexCount;i++){
			for(int j=0;j<this.vertexCount;j++){
				this.matrix[i][j] = in.nextInt();
			}
		}
		in.close();
		this.pathLength = vertexCount;
		this.path = new int[this.pathLength];
		this.seq = new int[this.pathLength];
	}
	
	@Override
	public int getArity() {
		return this.pathLength*BIT_TO_INT;
	}

	@Override
	public long run(long[] genom) {

		long old = System.currentTimeMillis(); //time
		
		int offset=0;
		int vertexNumber=0;
		int index = 0;
		long tmp = 0;
				
		for (int i=0;i<this.pathLength/BIT_TO_INT;i++){
			offset = i<<3;
			tmp = genom[i];
			for (int j=0;j<BIT_TO_INT;j++){
				vertexNumber = (int)(tmp & 255);
				index = offset+j;
				this.path[index] = vertexNumber;
				this.seq[index] = index;
				tmp >>= 8;
			}
		}
		
		this.prepareTime += (System.currentTimeMillis()-old); //time
		old = System.currentTimeMillis(); //time
		
		qsort(this.path,this.seq,0,this.pathLength-1);
		
		this.sortingTime += (System.currentTimeMillis()-old); //time
		old = System.currentTimeMillis(); //time
		
		long pathLength = this.checkPath(this.seq);
		
		this.checkTime += (System.currentTimeMillis()-old); //time
		
		return (Long.MAX_VALUE-pathLength);
	}
	
	private void qsort(int[] arrayToSort, int[] arrayToMix,int l, int r){
		int i = l;
		int j = r;
		int tmp = 0;
		int pivot = arrayToSort[(l+r)>>1];

		while (i <= j) {
			while (arrayToSort[i] < pivot) {i+=1;}
			while (arrayToSort[j] > pivot) {j-=1;}

			if (i <= j) {
				tmp = arrayToSort[i];
				arrayToSort[i] = arrayToSort[j];
				arrayToSort[j] = tmp;
				tmp = arrayToMix[i];
				arrayToMix[i] = arrayToMix[j];
				arrayToMix[j] = tmp;
				i+=1;
				j-=1;
			}
		}
		if (l < j){
			qsort(arrayToSort, arrayToMix, l, j);
		}
		if (i < r){
			qsort(arrayToSort, arrayToMix, i, r);
		}
	}
	
	public long checkPath(int[] path){
		
		long result = 0;
		int pathLength = path.length;
		int predVertex = path[0];
		int nextVertex = 0;
		
		for (int i=1;i<pathLength;i++){
			nextVertex = path[i];
			result += this.matrix[predVertex][nextVertex];			
			predVertex = nextVertex;
		}
		
		return result;
	}
	
	public static void generateRandomFile(String filename,int n) throws IOException{
		Random random = new Random();
		BufferedWriter out = new BufferedWriter(new FileWriter(filename));
		out.write(n+"\n");
		
		int[][] matrix = new int[n][n];
		
		for(int i=0;i<n-1;i++){
			for(int j=i+1;j<n;j++){
				matrix[i][j] = random.nextInt(256);
				matrix[j][i] = matrix[i][j];
			}
		}
		
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				out.write(matrix[i][j]+" ");
			}
			out.write("\n");
		}
		
		out.flush();
		out.close();
	}
}
