package com.knok16.GeneticFramework;

public interface FitnessFunction {
	int getArity();
	long run(long[] genom);
}
