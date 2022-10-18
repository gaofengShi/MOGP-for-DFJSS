/*
  Copyright 2010 by Sean Luke and George Mason University
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/

package ec.multiobjective.nsga2;

import ec.*;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.util.SortComparator;
import yimei.jss.gp.GPRuleEvolutionState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/* 
 * NSGA2Evaluator.java
 * 
 * Created: Sat Oct 16 00:19:57 EDT 2010
 * By: Faisal Abidi and Sean Luke
 */

/**
 * The NSGA2Evaluator is a simple generational evaluator which evaluates every
 * single member of the population (in a multithreaded fashion). Then it reduces
 * the population size to an <i>archive</i> consisting of the best front ranks.
 * When there isn't enough space to fit another front rank, individuals in that
 * final front rank vie for the remaining slots in the archive based on their
 * sparsity.
 * 
 * <p>
 * The evaluator is also responsible for calculating the rank and sparsity
 * values stored in the NSGA2MultiObjectiveFitness class and used largely for
 * statistical information.
 * 
 * <p>
 * NSGA-II has fixed archive size (the population size), and so ignores the
 * 'elites' declaration. However it will adhere to the 'reevaluate-elites'
 * parameter in SimpleBreeder to determine whether to force fitness
 * reevaluation.
 *
 */

public class NSGA2Evaluator extends SimpleEvaluator {
	/**
	 * The original population size is stored here so NSGA2 knows how large to
	 * create the archive (it's the size of the original population -- keep in mind
	 * that NSGA2Breeder had made the population larger to include the children.
	 */
	public int originalPopSize[];
	public int gen;


	public void setup(final EvolutionState state, final Parameter base) {
		super.setup(state, base);

		Parameter p = new Parameter(Initializer.P_POP);
		int subpopsLength = state.parameters.getInt(p.push(Population.P_SIZE), null, 1);
		Parameter p_subpop;
		originalPopSize = new int[subpopsLength];
		for (int i = 0; i < subpopsLength; i++) {
			p_subpop = p.push(Population.P_SUBPOP).push("" + i).push(Subpopulation.P_SUBPOPSIZE);
			originalPopSize[i] = state.parameters.getInt(p_subpop, null, 1);
		}
	}

	double currentMin;
	double currentMax;


	/**
	 * Evaluates the population, then builds the archive and reduces the population
	 * to just the archive.
	 */
	public void evaluatePopulation(final EvolutionState state) {
		super.evaluatePopulation(state);  //the same with the simpleEvalutor, during the first generation, evaluate N individuals; after that, evaluate 2N individuals

		//Gaofeng shi 21.09.2022. set the range of fitness of task.  ++++++++++++++++++++++
		ArrayList<Double> currentFitness0 = new ArrayList<>();
		ArrayList<Double> currentFitness1 = new ArrayList<>();
		GPRuleEvolutionState.rangeFitnessTask0.clear();
		for(int i = 0; i < state.population.subpops[0].individuals.length; i++){
			if(((MultiObjectiveFitness)state.population.subpops[0].individuals[i].fitness).getObjective(0) != Double.MAX_VALUE &
					((MultiObjectiveFitness)state.population.subpops[0].individuals[i].fitness).getObjective(1)!= Double.MAX_VALUE){
				currentFitness0.add(((MultiObjectiveFitness)state.population.subpops[0].individuals[i].fitness).getObjective(0));
				currentFitness1.add(((MultiObjectiveFitness)state.population.subpops[0].individuals[i].fitness).getObjective(1));
			}
		}
		GPRuleEvolutionState.rangeFitnessTask0.add(Collections.max(currentFitness0)-Collections.min(currentFitness0)); //Gaofeng 21.09.2022 for normalzation of fitness.
		GPRuleEvolutionState.rangeFitnessTask0.add(Collections.max(currentFitness1)-Collections.min(currentFitness1)); //Gaofeng 21.09.2022 for normalzation of fitness.

//		System.out.println("the range of fitness: "+ GPRuleEvolutionState.rangeFitnessTask0);
		//++++++++++++++++++++++++++++++++++++

		for (int x = 0; x < state.population.subpops.length; x++) {
			state.population.subpops[x].individuals = buildArchive(state, x); // trade the individuals in archive as the population
		}
	}


	/**
	 * Build the auxiliary fitness data and reduce the subpopulation to just the
	 * archive, which is returned.
	 */
	//achieve a archive has the same size of original population
	public Individual[] buildArchive(EvolutionState state, int subpop) {
		Individual[] dummy = new Individual[0]; //allocates an array which has 0 elements.

		// Gaofeng Shi 2022.08.10 calculate the alpha dominance
		// use the adaptived dominated alpha chose the front for the rule size colculation.
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++
		ArrayList<Individual> frontAlpha = new ArrayList<>();

		// Gaofeng shi 2022.08.24 use front as thr reulesize reference.
		MultiObjectiveFitness typicfitness = (MultiObjectiveFitness) (state.population.subpops[subpop].individuals[0].fitness);
		frontAlpha = typicfitness.partitionIntoParetoFront(state.population.subpops[subpop].individuals, null,null);


		alphaUpdate0(frontAlpha);

		System.out.println("alpha:  "+alphaMatrix[0]);
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		ArrayList ranks = assignFrontRanks(state.population.subpops[subpop]); //after this, get different several ranks
		//each rank cosists of the corresponding individuals


		ArrayList newSubpopulation = new ArrayList(); //a new one, size = 0
		int size = ranks.size();
		for (int i = 0; i < size; i++) { //do for each rank separately
			Individual[] rank = (Individual[]) ((ArrayList) (ranks.get(i))).toArray(dummy);
			assignSparsity(rank); //assign sparity value for each individual in this rank
			if (rank.length + newSubpopulation.size() >= originalPopSize[subpop]) {
				// first sort the rank by sparsity---from the small one to large one
				ec.util.QuickSort.qsort(rank, new SortComparator() {
					public boolean lt(Object a, Object b) { //Returns true if a < b, else false 
						Individual i1 = (Individual) a;
						Individual i2 = (Individual) b;
						return (((NSGA2MultiObjectiveFitness) i1.fitness).sparsity > ((NSGA2MultiObjectiveFitness) i2.fitness).sparsity);
					}

					public boolean gt(Object a, Object b) { //Returns true if a > b, else false 
						Individual i1 = (Individual) a;
						Individual i2 = (Individual) b;
						return (((NSGA2MultiObjectiveFitness) i1.fitness).sparsity < ((NSGA2MultiObjectiveFitness) i2.fitness).sparsity);
					}
				});

				// then put the m sparsest individuals in the new population
				int m = originalPopSize[subpop] - newSubpopulation.size(); //how many positions left for new individuals
				for (int j = 0; j < m; j++) 
					newSubpopulation.add(rank[j]); //add some individuals based on the sparisity

				// and bail
				break;
			} else {
				// dump in everyone
				for (int j = 0; j < rank.length; j++) //add the while rank directly
					newSubpopulation.add(rank[j]);
			}
		}

		Individual[] archive = (Individual[]) (newSubpopulation.toArray(dummy));

		// maybe force reevaluation
		NSGA2Breeder breeder = (NSGA2Breeder) (state.breeder);
		if (breeder.reevaluateElites[subpop])
			for (int i = 0; i < archive.length; i++)
				archive[i].evaluated = false;

		return archive;
	}

	/**
	 * Divides inds into ranks and assigns each individual's rank to be the rank it
	 * was placed into. Each front is an ArrayList.
	 */
	public ArrayList assignFrontRanks(Subpopulation subpop) {
		Individual[] inds = subpop.individuals; //inds includes all individuals

		ArrayList frontsByRank = MultiObjectiveFitness.partitionIntoRanks(inds);

		int numRanks = frontsByRank.size();
		for (int rank = 0; rank < numRanks; rank++) {
			ArrayList front = (ArrayList) (frontsByRank.get(rank));
			int numInds = front.size();
			for (int ind = 0; ind < numInds; ind++)
				((NSGA2MultiObjectiveFitness) (((Individual) (front.get(ind))).fitness)).rank = rank;
		}
		return frontsByRank;
	}

	/**
	 * Computes and assigns the sparsity values of a given front.
	 */
	public void assignSparsity(Individual[] front) {
		int numObjectives = ((NSGA2MultiObjectiveFitness) front[0].fitness).getObjectives().length;

		//front here means different fronts
		for (int i = 0; i < front.length; i++) //front.length means how many individuals in this ranking(front)
			((NSGA2MultiObjectiveFitness) front[i].fitness).sparsity = 0; //the first individual in this front, the sparsity 
		//is assigned to 0

		for (int i = 0; i < numObjectives; i++) {
			final int o = i;
			// 1. Sort front by each objective.
			// 2. Sum the manhattan distance of an individual's neighbours over
			// each objective.
			// NOTE: No matter which objectives objective you sort by, the
			// first and last individuals will always be the same (they maybe
			// interchanged though). This is because a Pareto front's
			// objective values are strictly increasing/decreasing.
			ec.util.QuickSort.qsort(front, new SortComparator() {
				public boolean lt(Object a, Object b) { //less than
					Individual i1 = (Individual) a;
					Individual i2 = (Individual) b;
					return (((NSGA2MultiObjectiveFitness) i1.fitness)
							.getObjective(o) < ((NSGA2MultiObjectiveFitness) i2.fitness).getObjective(o));
				}

				public boolean gt(Object a, Object b) { //great than
					Individual i1 = (Individual) a;
					Individual i2 = (Individual) b;
					return (((NSGA2MultiObjectiveFitness) i1.fitness)
							.getObjective(o) > ((NSGA2MultiObjectiveFitness) i2.fitness).getObjective(o));
				}
			});

			// Compute and assign sparsity.
			// the first and last individuals are the sparsest.
			((NSGA2MultiObjectiveFitness) front[0].fitness).sparsity = Double.POSITIVE_INFINITY;
			((NSGA2MultiObjectiveFitness) front[front.length - 1].fitness).sparsity = Double.POSITIVE_INFINITY;
			for (int j = 1; j < front.length - 1; j++) {
				NSGA2MultiObjectiveFitness f_j = (NSGA2MultiObjectiveFitness) (front[j].fitness);
				NSGA2MultiObjectiveFitness f_jplus1 = (NSGA2MultiObjectiveFitness) (front[j + 1].fitness);
				NSGA2MultiObjectiveFitness f_jminus1 = (NSGA2MultiObjectiveFitness) (front[j - 1].fitness);

//				System.out.println(f_j.maxObjective[o] - f_j.minObjective[o]);  //1
				// store the NSGA2Sparsity in sparsity
				f_j.sparsity += (f_jplus1.getObjective(o) - f_jminus1.getObjective(o))
						/ (f_j.maxObjective[o] - f_j.minObjective[o]);
			}
		}
	}

	//fzhang 2018.11.6 NSGA-II
	public void evaluatePopulationgp(final EvolutionState state) {
		super.evaluatePopulation(state);// evaluate population for GP with NSGA2
	}


	/**
	 * Gaofeng shi 2022.08.10 set the alpha dominance to unbalance objectives.
	 */
	//+++++++++++++++++++++++++++++++++++++++++++++++++COMP589++++++++++++++++++++++++
	public static double[] alphaMatrix = new double[]{0.0,0.0};
	public double learningRate = 45.0;
	public double maxFit =0;
	public double minFit = 9999999;
	public double maxSize = 0;
	public double minSize = 10;

	/** The alpha update scheme which only consider for rule size.
	 *
	 */

	public void alphaUpdate0(ArrayList<Individual> front){

		ArrayList<Double> sizeList = new ArrayList<>();
		for (int j =0; j < front.size(); j++){
			sizeList.add(((MultiObjectiveFitness)((Individual)front.get(j)).fitness).getObjective(1));
		}
		double currentMaxSize= Collections.max(sizeList);
		double currentMinSize = Collections.min(sizeList);


		if (currentMaxSize > maxSize)
			maxSize = currentMaxSize;
		if (currentMinSize < minSize)
			minSize = currentMinSize;

		double meanSize = Mean(sizeList);

		double sizeRatio = (meanSize - minSize)/(maxSize - minSize);

		if (sizeRatio<0.5){
			double degrees = Math.toDegrees(Math.atan(alphaMatrix[0])) + (learningRate * (0.5 - sizeRatio) / 0.5);

			if (degrees >= 89.9999)
				degrees = 89.9999;

			alphaMatrix[0] = Math.tan(Math.toRadians(degrees));
		}
		if (sizeRatio > 0.5){ // bias to large size
			double degrees = Math.toDegrees(Math.atan(alphaMatrix[0])) - (learningRate * (sizeRatio - 0.5) / 0.5);

			if (degrees <= 0)
				degrees = 0;
			alphaMatrix[0] = Math.tan(Math.toRadians(degrees));
		}
		System.out.println("Size range-----------"+minSize+"--"+maxSize);
        System.out.println("Mean size------------"+meanSize);
        System.out.println("SizeRatio------------"+sizeRatio);

	}

	/** The alpha update scheme which consider both eff and size,
	 * if bias to size, increase alpha,
	 * if bias to eff, decrease alpha*/
	public void alphaUpdate(ArrayList<Individual> front){

		ArrayList<Double> fitlist = new ArrayList<>();
		ArrayList<Double> sizelist = new ArrayList<>();

		for ( int i=0; i < front.size(); i++){
			fitlist.add(((MultiObjectiveFitness)((Individual)front.get(i)).fitness).getObjective(0));
			sizelist.add(((MultiObjectiveFitness)((Individual)front.get(i)).fitness).getObjective(1));
		}

		double currentMaxFit = Collections.max(fitlist);
		double currentMinFit = Collections.min(fitlist);

		double currentMaxSize = Collections.max(sizelist);
		double currentMinSize = Collections.min(sizelist);

		if (currentMaxFit > maxFit) maxFit = currentMaxFit;
		if (currentMinFit < minFit) minFit = currentMinFit;
		if (currentMaxSize > maxSize) maxSize = currentMaxSize;
		if (currentMinSize < minSize) minSize = currentMinSize;

		double meanFit = Mean(fitlist);
		double meanSize = Mean(sizelist);

		double fitRatio = (maxFit-meanFit)/(meanFit-minFit);
		double sizeRatio = (maxSize-meanSize)/(meanSize-minSize);

		if ((fitRatio > 1.0)&(sizeRatio < 1.0)){
			alphaMatrix[0] = alphaMatrix[0]-learningRate;
		}
		if ((fitRatio < 1.0)&(sizeRatio > 1.0)){
			alphaMatrix[0] = alphaMatrix[0]+learningRate;
		}

//		if ((alphaMatrix[0] <= 0) || (gen>=45)) alphaMatrix[0]=0.0;
		if (alphaMatrix[0] <= 0) alphaMatrix[0]=0.0;

		//TODO new update secheme which use the fit ratio and size ratio as the gridiant,
		// then identify how should i update alpha (alpha * ??? or alpha + 0.1*gridiant),
		// only meanfit and size is insufficient to identify bias, need another measurments

		System.out.println("Fit range------------");
		System.out.println(maxFit);
		System.out.println(minFit);
		System.out.println("Size range------------");
		System.out.println(maxSize);
		System.out.println(minSize);
		System.out.println("Fit ratio------------");
		System.out.println(fitRatio);
		System.out.println("Size ratio------------");
		System.out.println(sizeRatio);
		System.out.println("alpha------------");
		System.out.println(alphaMatrix[0]);

	}
	/** Calculates the mean value from a list of numeric values.
	 ** param values List of numeric values
	 * return Median value
	 */
	public static double Mean(ArrayList<Double> values) {
		double sum = 0;
		for (int i = 0; i < values.size();i++){
			sum += values.get(i);
		}
		return sum / values.size();
	}
}