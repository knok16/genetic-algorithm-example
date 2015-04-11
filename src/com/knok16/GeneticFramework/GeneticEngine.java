package com.knok16.GeneticFramework;

import java.util.Random;

public class GeneticEngine{

	public long timeToSelection = 0;
	public long timeToCrossing = 0;
	public long timeToMutate = 0;
	public long timeToFF = 0;

	public static final double CHANCE_TO_FULLNESS = 0.999d;
	public static final SelectionType DEFAULT_SELECTION_TYPE = SelectionType.TOURNEY;
	public static final CrossingType DEFAULT_CROSSING_TYPE = CrossingType.ONE_POINT_RECOMBINATION;
	public static final boolean DEFAULT_USE_MUTATION = true;
	public static final long DEFAULT_GENERATION_COUNT = 10000L;

	public static final int OCTET_LENGTH = 64; // for long
	public static final int MASK_FOR_MOD = OCTET_LENGTH - 1;
	public static final int SHIFT_FOR_DIVISION;
	static {
		int shiftForDivision = 0;
		int tmp = OCTET_LENGTH;
		while (tmp > 1) {
			tmp >>= 1;
			shiftForDivision++;
		}
		SHIFT_FOR_DIVISION = shiftForDivision;
	}

	public enum SelectionType {
		TOURNEY, ROULETTE_WHEEL
	}

	public enum CrossingType {
		ONE_POINT_RECOMBINATION, TWO_POINT_RECOMBINATION, ELEMENTWISE_RECOMBINATION, ONE_ELEMENT_EXCHANGE
	}

	private FitnessFunction fitnessFunction;
	private int genomLength;
	private int sizeOfArray;
	private long generationCount;
	private int individualCount;
	private SelectionType selectionType;
	private CrossingType crossingType;
	private boolean useMutation;
	private double mutationPercent;
	private long[][] genomListParents;
	private long[][] genomListOffsprings;
	private long[] actual;
	private long[] fitnessFunctionResult;
	private long currentGeneration = 0;

	private Random random = new Random(System.currentTimeMillis());

	public GeneticEngine(FitnessFunction fitnessFunction) {
		this.fitnessFunction = fitnessFunction;
		this.genomLength = fitnessFunction.getArity();
		this.sizeOfArray = (int) Math.ceil((double) this.genomLength / OCTET_LENGTH);
		this.generationCount = DEFAULT_GENERATION_COUNT;
		this.individualCount = (int) (1 + Math.log(1 / Math.pow(1 - CHANCE_TO_FULLNESS, 1 / genomLength)) / Math.log(2));
		this.selectionType = DEFAULT_SELECTION_TYPE;
		this.crossingType = DEFAULT_CROSSING_TYPE;
		this.useMutation = DEFAULT_USE_MUTATION;
		this.mutationPercent = genomLength * (1 - Math.pow((1 - 10 * Math.pow((1 / 2), (genomLength - 1))),(1 / genomLength)));
	}

	// Main loop
	public long[] run() {
		//Preparing structuress
		this.genomListParents = new long[this.individualCount][];
		this.genomListOffsprings = new long[this.individualCount][];
		this.fitnessFunctionResult = new long[this.individualCount];
		this.actual = new long[this.individualCount];
		for (int i = 0; i < this.individualCount; i++) {
			this.actual[i] = -1;
		}
		
		//Generate 1st generation
		this.generateFirstGeneration();
		
		while (this.currentGeneration < this.generationCount) {

			this.selection();
			this.crossing();
			if (this.useMutation) {
				this.mutation();
			}
			
			long[][] tmp = this.genomListParents;
			this.genomListParents = this.genomListOffsprings;
			this.genomListOffsprings = tmp;
			
			this.currentGeneration++;
		}

		long bestFitnessFunctionResult = 0;
		long[] bestGenom = null;
		for (long[] genom : this.genomListParents) {
			long fitnessFunctionResult = this.fitnessFunction.run(genom);
			if (bestFitnessFunctionResult <= fitnessFunctionResult) {
				bestGenom = genom;
				bestFitnessFunctionResult = fitnessFunctionResult;
			}
		}

		return bestGenom;
	}

	// Generate First Generation
	private void generateFirstGeneration() {
		for (int i = 0; i < this.individualCount; i++) {
			this.genomListParents[i] = this.generateGenom();
		}
	}

	// Generate Genom - Generate 1 genom
	private long[] generateGenom() {
		long[] result = new long[this.sizeOfArray];
		for (int i = 0; i < this.sizeOfArray; i++) {
			result[i] = this.random.nextLong();
		}
		return result;
	}
	
	// Selection - Select genoms for crossing
	private void selection(){
		long old = System.currentTimeMillis(); // time
		
		switch (selectionType) {
		case ROULETTE_WHEEL:{
			
			float[] wheel = new float[this.individualCount];
			wheel[0] = this.getFitnessFunctionResult(0); 
			for (int i=1;i<this.individualCount;i++){
				wheel[i] = wheel[i-1] + this.getFitnessFunctionResult(i);
			}
			float all = wheel[this.individualCount-1];
			
			for (int i=0;i<this.individualCount;i++){
				float index = Math.abs(this.random.nextFloat())*all;
				
				int l = 0;
				int r = individualCount-1;
				int c = 0;
				while (l!=r){
					c = (l+r) >> 1;
					if (wheel[c]<index){
						l=c;
					}else{
						r=c;
					}
				}
				this.genomListOffsprings[i] = this.genomListParents[l].clone();
			}
			break;
		}
		case TOURNEY:{
			for (int i=0;i<this.individualCount;i++){
				int index1 = random.nextInt(individualCount);
				int index2 = random.nextInt(individualCount);
	
				long ffTime = System.currentTimeMillis(); // time
	
				long fr1 = this.getFitnessFunctionResult(index1);
				long fr2 = this.getFitnessFunctionResult(index2);
	
				this.timeToFF += (System.currentTimeMillis() - ffTime); // time

				this.genomListOffsprings[i] = fr1 > fr2 ? this.genomListParents[index1].clone() : this.genomListParents[index2].clone();
			}
			break;
		}
		default:
			throw new UnsupportedOperationException();
		}
		
		this.timeToSelection += (System.currentTimeMillis() - old); // time
	}

	// Crossing - Crossing all genom in generation
	private void crossing() {
		long old = System.currentTimeMillis(); // time

		for (int i = 0; i < individualCount / 2; i++) {
			int index1 = i << 1;
			int index2 = index1 | 1;
			cross(this.genomListOffsprings[index1], this.genomListOffsprings[index2]);
		}

		this.timeToCrossing += (System.currentTimeMillis() - old); // time
	}

	// Get Fitness Function Result [with cache]
	private long getFitnessFunctionResult(int genomNumber) {
		if (this.actual[genomNumber] != this.currentGeneration) {
			this.fitnessFunctionResult[genomNumber] = this.fitnessFunction.run(this.genomListParents[genomNumber]);
			this.actual[genomNumber] = this.currentGeneration;
		}
		return this.fitnessFunctionResult[genomNumber];
	}

	// Cross - Crossing 2 genom
	private void cross(long[] genom1, long[] genom2) {
		switch (crossingType) {
		case ONE_POINT_RECOMBINATION:{
			int index = this.random.nextInt(this.genomLength);
			int outerOffset = index >> SHIFT_FOR_DIVISION;
			int innerOffset = OCTET_LENGTH - (index & MASK_FOR_MOD);
			long tmp = 0;

			if (innerOffset < 63) {
				long mask = 1L << (innerOffset + 1) - 1;
				long swapMask =  (genom1[outerOffset] ^ genom2[outerOffset]) & mask;
				genom1[outerOffset] ^= swapMask;
				genom2[outerOffset] ^= swapMask;
				outerOffset++;
			}
			for (int i=outerOffset;i<this.sizeOfArray;i++){
				tmp = genom1[i];
				genom1[i] = genom2[i];
				genom2[i] = tmp;
			}
			
			break;
		}
		case TWO_POINT_RECOMBINATION:{
			int index1 = this.random.nextInt(this.genomLength);
			int index2 = this.random.nextInt(this.genomLength);
			int startIndex = Math.min(index1, index2);
			int endIndex = Math.max(index1, index2);
			int startOuterOffset = startIndex >> SHIFT_FOR_DIVISION;
			int startInnerOffset = OCTET_LENGTH - (startIndex & MASK_FOR_MOD);
			int endOuterOffset = endIndex >> SHIFT_FOR_DIVISION;
			int endInnerOffset = OCTET_LENGTH - (endIndex & MASK_FOR_MOD);
			long tmp = 0;

			if (startInnerOffset < OCTET_LENGTH-1) {
				long mask = 1L << (startInnerOffset + 1) - 1;
				long swapMask =  (genom1[startOuterOffset] ^ genom2[startOuterOffset]) & mask;
				genom1[startOuterOffset] ^= swapMask;
				genom2[startOuterOffset] ^= swapMask;
				startOuterOffset++;
			}
			for (int i=startOuterOffset;i<=endOuterOffset;i++){
				tmp = genom1[i];
				genom1[i] = genom2[i];
				genom2[i] = tmp;
			}
			if (endInnerOffset > 0) {
				long mask = 1L << endInnerOffset - 1;
				long swapMask =  (genom1[endOuterOffset] ^ genom2[endOuterOffset]) & mask;
				genom1[endOuterOffset] ^= swapMask;
				genom2[endOuterOffset] ^= swapMask;
			}
			
			break;
		}
		case ELEMENTWISE_RECOMBINATION:{
			for (int outerOffset = 0; outerOffset < this.sizeOfArray; outerOffset++) {
				long mask = this.random.nextLong();
				long swapMask = (genom1[outerOffset] ^ genom2[outerOffset]) & mask;

				genom1[outerOffset] ^= swapMask;
				genom2[outerOffset] ^= swapMask;
			}
			break;
		}
		case ONE_ELEMENT_EXCHANGE:{
			int index = this.random.nextInt(this.genomLength);
			int outerOffset = index >> SHIFT_FOR_DIVISION;
			int innerOffset = OCTET_LENGTH - (index & MASK_FOR_MOD);
			long mask = 1L << innerOffset;
			long swapMask = (genom1[outerOffset] ^ genom2[outerOffset]) & mask;

			genom1[outerOffset] ^= swapMask;
			genom2[outerOffset] ^= swapMask;
			break;
		}
		default:
			throw new UnsupportedOperationException();
		}
	}

	// Mutation - Mutate all genom in generation
	private void mutation() {
		long old = System.currentTimeMillis(); // time

		for (long[] genom : this.genomListOffsprings) {
			if (random.nextDouble() <= mutationPercent) {
				mutate(genom);
			}
		}

		this.timeToMutate += (System.currentTimeMillis() - old); // time
	}

	// Mutate - Mutate 1 genom
	private void mutate(long[] genom) {
		int index = this.random.nextInt(this.genomLength);
		int outerOffset = index >> SHIFT_FOR_DIVISION;
		int innerOffset = (index & MASK_FOR_MOD);
		long mask = 1L << innerOffset;
		genom[outerOffset] ^= mask;
	}

	// Setters and Getters
	public long getGenerationCount() {
		return generationCount;
	}

	public void setGenerationCount(long generationCount) {
		this.generationCount = generationCount;
	}

	public int getIndividualCount() {
		return individualCount;
	}

	public void setIndividualCount(int individualCount) {
		this.individualCount = individualCount;
	}

	public SelectionType getSelectionType() {
		return selectionType;
	}

	public void setSelectionType(SelectionType selectionType) {
		this.selectionType = selectionType;
	}

	public CrossingType getCrossingType() {
		return crossingType;
	}

	public void setCrossingType(CrossingType crossingType) {
		this.crossingType = crossingType;
	}

	public boolean getUseMutation() {
		return useMutation;
	}

	public void setUseMutation(boolean useMutation) {
		this.useMutation = useMutation;
	}

	public double getMutationPercent() {
		return mutationPercent;
	}

	public void setMutationPercent(double mutationPercent) {
		this.mutationPercent = mutationPercent;
	}
}