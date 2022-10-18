/*
  Copyright 2018 by BINZI
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/

package yimei.jss.algorithm.multiobjective;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.Subpopulation;
import ec.multiobjective.MultiObjectiveFitness;
import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
import ec.simple.SimpleBreeder;
import ec.util.Parameter;
import ec.util.SortComparator;

import java.util.ArrayList;

/* 
 * SimpleBreederelite.java
 * 
 * Created: 2018
 * By: BINZI
 */

/**
 * LoadElites based on NSGA2's rank and sparity
 *
 *
 * @author BINZI
 * @version 1.0
 */

public class NSGA2BreederElite extends SimpleBreeder {

	/**
	 * A private helper function for breedPopulation which loads elites into a
	 * subpopulation.
	 */

	protected void loadElites(EvolutionState state, Population newpop) {
		// are our elites small enough?
		for (int x = 0; x < state.population.subpops.length; x++) {
			if (numElites(state, x) > state.population.subpops[x].individuals.length)
				state.output.error(
						"The number of elites for subpopulation " + x + " exceeds the actual size of the subpopulation",
						new Parameter(EvolutionState.P_BREEDER).push(P_ELITE).push("" + x));
			if (numElites(state, x) == state.population.subpops[x].individuals.length)
				state.output.warning(
						"The number of ewarninglites for subpopulation " + x + " is the actual size of the subpopulation",
						new Parameter(EvolutionState.P_BREEDER).push(P_ELITE).push("" + x));
		}
		state.output.exitIfErrors();

		// we assume that we're only grabbing a small number (say <10%), so
		// it's not being done multithreaded
		for (int sub = 0; sub < state.population.subpops.length; sub++) // ����ÿһ��subpop
		{
			if (!shouldBreedSubpop(state, sub, 0)) // don't load the elites for this one, we're not doing breeding of it
			{
				continue;
			}

			Individual[] dummy = new Individual[0];// α��Ⱥ

			ArrayList ranks = assignFrontRanks(state.population.subpops[sub]);// ����ranks

			int size = ranks.size();
			/*for (int i = 0; i < size; i++) {
				Individual[] rank = (Individual[]) ((ArrayList) (ranks.get(i))).toArray(dummy);// ȡ��i��rank
				assignSparsity(rank);// ���ŷ�ɢ�� [0,2]
			}*/
			// ���ϼ��������и����rank��sparsity
			if (numElites(state, sub) == 1) {
				int best = 0;
				Individual[] oldinds = state.population.subpops[sub].individuals;
				for (int x = 1; x < oldinds.length; x++)
					if (oldinds[x].fitness.betterThan(oldinds[best].fitness))
						best = x;
				Individual[] inds = newpop.subpops[sub].individuals;
				inds[inds.length - 1] = (Individual) (oldinds[best].clone());
			} 
			else if (numElites(state, sub) > 0) // we'll need to sort
			{
				Individual[] inds = newpop.subpops[sub].individuals;
				Individual[] oldinds = state.population.subpops[sub].individuals;
				int num = 0;// �ѷ���ľ�Ӣ����
				for (int i = 0; i < size; i++) {
					Individual[] rank = (Individual[]) ((ArrayList) (ranks.get(i))).toArray(dummy);// ȡ��i��rank
					
					if (rank.length + num >= numElites(state, sub))// ֻȡ����
					{
						// �����еĸ�������,�Ӵ�С
						ec.util.QuickSort.qsort(rank, new SortComparator() {
							public boolean lt(Object a, Object b) {
								Individual i1 = (Individual) a;
								Individual i2 = (Individual) b;
								return (((NSGA2MultiObjectiveFitness) i1.fitness).sparsity > ((NSGA2MultiObjectiveFitness) i2.fitness).sparsity);
							}

							public boolean gt(Object a, Object b) {
								Individual i1 = (Individual) a;
								Individual i2 = (Individual) b;
								return (((NSGA2MultiObjectiveFitness) i1.fitness).sparsity < ((NSGA2MultiObjectiveFitness) i2.fitness).sparsity);
							}
						});
						// end of sort

						// load the top N individuals
						   //int m = numElites(state, sub)- num;//��Ҫm��
						   int j=0;
						   for (int x = inds.length - num-1; x>= inds.length - numElites(state, sub) ; x--)
			               {
							   inds[x] = (Individual) (rank[j].clone()); 
			                    j++;
			               }
						break;
					} 
					else
					{
						   int j=0;
						for (int x = inds.length - num-1; x>= inds.length - rank.length-num; x--)
			               {
							   inds[x] = (Individual) (rank[j].clone()); 
							   j++;
			               }
						   num+=rank.length;
					}

				}
			}
		}
	}

	public ArrayList assignFrontRanks(Subpopulation subpop) {
		Individual[] inds = subpop.individuals;
		ArrayList frontsByRank = MultiObjectiveFitness.partitionIntoRanks(inds);// ����rank

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

		for (int i = 0; i < front.length; i++)
			((NSGA2MultiObjectiveFitness) front[i].fitness).sparsity = 0;

		for (int i = 0; i < numObjectives; i++) {
			final int o = i;
			// 1. Sort front by each objective.
			// 2. Sum the manhattan distance of an individual's neighbours over
			// each objective.
			// NOTE: No matter which objectives objective you sort by, the
			// first and last individuals will always be the same (they maybe
			// interchanged though). This is because a Pareto front's
			// objective values are strictly increasing/decreasing.
			ec.util.QuickSort.qsort(front, new SortComparator()// ����i��Ŀ���front�и�������
			{
				public boolean lt(Object a, Object b) {
					Individual i1 = (Individual) a;
					Individual i2 = (Individual) b;
					return (((NSGA2MultiObjectiveFitness) i1.fitness)
							.getObjective(o) < ((NSGA2MultiObjectiveFitness) i2.fitness).getObjective(o));
				}

				public boolean gt(Object a, Object b) {
					Individual i1 = (Individual) a;
					Individual i2 = (Individual) b;
					return (((NSGA2MultiObjectiveFitness) i1.fitness)
							.getObjective(o) > ((NSGA2MultiObjectiveFitness) i2.fitness).getObjective(o));
				}
			});
			final double min = ((MultiObjectiveFitness) front[0].fitness).getObjective(o);
			final double max = ((MultiObjectiveFitness) front[front.length - 1].fitness).getObjective(o);
			// Compute and assign sparsity.
			// the first and last individuals are the sparsest.
			((NSGA2MultiObjectiveFitness) front[0].fitness).sparsity = Double.POSITIVE_INFINITY;
			((NSGA2MultiObjectiveFitness) front[front.length - 1].fitness).sparsity = Double.POSITIVE_INFINITY;
			for (int j = 1; j < front.length - 1; j++)// �м�ĸ���
			{
				NSGA2MultiObjectiveFitness f_j = (NSGA2MultiObjectiveFitness) (front[j].fitness);// �ҵ���Ӧ��
				NSGA2MultiObjectiveFitness f_jplus1 = (NSGA2MultiObjectiveFitness) (front[j + 1].fitness);// ������һ������Ӧ��
				NSGA2MultiObjectiveFitness f_jminus1 = (NSGA2MultiObjectiveFitness) (front[j - 1].fitness);// ������һ������Ӧ��
				if (max == min)// �����ڸ�Ŀ�꣬�����С��һ����˵���ڸ�Ŀ���Ͼ���Ϊ0
				{
					f_j.sparsity += 0;
				} else {
					// store the NSGA2Sparsity in sparsity
					f_j.sparsity += (f_jplus1.getObjective(o) - f_jminus1.getObjective(o)) / (max - min);
				}
			}
		}
	}
}
