/*
  Copyright 2010 by Sean Luke and George Mason University
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/

package ec.multiobjective.nsga2;

import ec.EvolutionState;
import ec.Fitness;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Code;
import yimei.jss.gp.GPRuleEvolutionState;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;

/* 
 * NSGA2MultiObjectiveFitness.java
 * 
 * Created: Thu Feb 04 2010
 * By: Faisal Abidi and Sean Luke
 */

/**
 * NSGA2MultiObjectiveFitness is a subclass of MultiObjeciveFitness which
 * adds auxiliary fitness measures (sparsity, rank) largely used by MultiObjectiveStatistics.
 * It also redefines the comparison measures to compare based on rank, and break ties
 * based on sparsity. 
 *
 */

public class NSGA2MultiObjectiveFitness extends MultiObjectiveFitness
    {
    public static final String NSGA2_RANK_PREAMBLE = "Rank: ";
    public static final String NSGA2_SPARSITY_PREAMBLE = "Sparsity: ";

    public String[] getAuxilliaryFitnessNames() { return new String[] { "Rank", "Sparsity" }; }
    public double[] getAuxilliaryFitnessValues() { return new double[] { rank, sparsity }; }
        
    /** Pareto front rank measure (lower ranks are better) */
    public int rank;

    /** Sparsity along front rank measure (higher sparsity is better) */
    public double sparsity;


    public String fitnessToString()
        {
        return super.fitnessToString() + "\n" + NSGA2_RANK_PREAMBLE + Code.encode(rank) + "\n" + NSGA2_SPARSITY_PREAMBLE + Code.encode(sparsity);
        }

    public String fitnessToStringForHumans()
        {
        return super.fitnessToStringForHumans() + "\n" + NSGA2_RANK_PREAMBLE + rank + "\n" + NSGA2_SPARSITY_PREAMBLE + sparsity;
        }

    public void readFitness(final EvolutionState state, final LineNumberReader reader) throws IOException
        {
        super.readFitness(state, reader);
        rank = Code.readIntegerWithPreamble(NSGA2_RANK_PREAMBLE, state, reader);
        sparsity = Code.readDoubleWithPreamble(NSGA2_SPARSITY_PREAMBLE, state, reader);
        }

    public void writeFitness(final EvolutionState state, final DataOutput dataOutput) throws IOException
        {
        super.writeFitness(state, dataOutput);
        dataOutput.writeInt(rank);
        dataOutput.writeDouble(sparsity);
        writeTrials(state, dataOutput);
        }

    public void readFitness(final EvolutionState state, final DataInput dataInput) throws IOException
        {
        super.readFitness(state, dataInput);
        rank = dataInput.readInt();
        sparsity = dataInput.readDouble();
        readTrials(state, dataInput);
        }

    public boolean equivalentTo(Fitness _fitness)
        {
        NSGA2MultiObjectiveFitness other = (NSGA2MultiObjectiveFitness) _fitness;
        return (rank == ((NSGA2MultiObjectiveFitness) _fitness).rank) &&
            (sparsity == other.sparsity);
        }

    /**
     * We specify the tournament selection criteria, Rank (lower
     * values are better) and Sparsity (higher values are better)
     */

    public boolean betterThan(Fitness _fitness)
        {
        NSGA2MultiObjectiveFitness other = (NSGA2MultiObjectiveFitness) _fitness;
        // Rank should always be minimized.
        if (rank < ((NSGA2MultiObjectiveFitness) _fitness).rank)
            return true;
        else if (rank > ((NSGA2MultiObjectiveFitness) _fitness).rank)
            return false;
                
        // otherwise try sparsity
        return (sparsity > other.sparsity);
        }


    /**
     * Gaofeng shi 2022.08.08 comp589, paretonDominates with alpha dominance
     */
    public boolean paretoDominates1(MultiObjectiveFitness other)
    {
        // Gaofeng shi 2022.08.10 the alpha dominance.
        double[] alpha = new double[]{0,0};
        alpha[0] = NSGA2Evaluator.alphaMatrix[1]; // set the alpha dominance for the number of objectives. objective0 is the performance.
        alpha[1] = NSGA2Evaluator.alphaMatrix[0];

        boolean abeatsb = false;

        if (objectives.length != other.objectives.length)
            throw new RuntimeException("Attempt made to compare two multiobjective fitnesses; but they have different numbers of objectives.");

        for (int x = 0; x < objectives.length; x++)
        {
            if (maximize[x] != other.maximize[x])  // uh oh
                throw new RuntimeException(
                        "Attempt made to compare two multiobjective fitnesses; but for objective #" + x +
                                ", one expects higher values to be better and the other expectes lower values to be better.");

            if (maximize[x]) //fzhang 2018.11.2  maximize[0] and maximize[1] are false. It means we are looking at minimising.
            {
                if (objectives[x] + alpha[x] * objectives[0]
                        > other.objectives[x] + alpha[x] * other.objectives[0] )
                    abeatsb = true;
                else if (objectives[x] + alpha[x] * objectives[0]
                        < other.objectives[x] + alpha[x] * other.objectives[0] )
                    return false;
            }
            else
            {
                if (objectives[x] + alpha[x] * objectives[0]
                        < other.objectives[x] + alpha[x] * other.objectives[0] ) //objective[0] and objective[1]: [2514.5953504074932, 498.03951619189513]
                    abeatsb = true;
                else if (objectives[x] + alpha[x] * objectives[0]
                        > other.objectives[x] + alpha[x] * other.objectives[0] )
                    return false;
            }
        }

        return abeatsb;
    }

        /**
         * Gaofeng Shi 2022.08.11 Comp589. alpha domiance sorting.
         * This is a copy of Shaolin code. doing the sorting following the paper.
         * @param other
         * @return
         */
    public boolean paretoDominates(MultiObjectiveFitness other){
        // Gaofeng shi 2022.08.10 the alpha dominance.
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        double[] alpha = NSGA2Evaluator.alphaMatrix; // set the alpha dominance for the number of objectives. objective0 is the performance.
        boolean abeatsb = false;

        if (objectives.length != other.objectives.length)
            throw new RuntimeException("Attempt made to compare two multiobjective fitnesses; but they have different numbers of objectives.");
        // get a index list of objectives
        ArrayList<Integer> indexArray = new ArrayList<>();
        for (int e = 0; e < objectives.length; e++)
            indexArray.add(e);

        for (int x = 0; x < objectives.length; x++)
        {
            if (maximize[x] != other.maximize[x])  // uh oh
                throw new RuntimeException(
                        "Attempt made to compare two multiobjective fitnesses; but for objective #" + x +
                                ", one expects higher values to be better and the other expectes lower values to be better.");
            //identify the dominance relationship by this equation  w(x,x`) = fi(x) - fi(x') + sum{alpha_ij * (fj(x) - fj(x`))}
            double sum = 0;
            for (int i = 0; i < indexArray.size(); i++){
                if (indexArray.get(i) == x)
                    continue;
                sum += alpha[i] * ((objectives[indexArray.get(i)] - other.objectives[indexArray.get(i)])
                / GPRuleEvolutionState.rangeFitnessTask0.get(i));  //Gaofeng 21.09.2022 normalzation of fitness.

            }
            double w = (objectives[x] - other.objectives[x])/GPRuleEvolutionState.rangeFitnessTask0.get(x) + sum; //Gaofeng 21.09.2022 normalzation of fitness.
            if (maximize[x])
            {
                if (w > 0)
                    abeatsb = true;
                else if (w < 0)
                    return false;
            }
            else
            {
                if (w < 0)
                    abeatsb = true;
                else if (w > 0)
                    return false;
            }
        }
        return abeatsb;
    }

    }


