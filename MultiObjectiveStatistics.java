/*
  Copyright 2010 by Sean Luke and George Mason University
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/

package ec.multiobjective;

import java.util.ArrayList;
import java.util.Arrays;
import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.multiobjective.MultiObjectiveFitness;
import ec.multiobjective.nsga2.NSGA2Evaluator;
import ec.simple.SimpleStatistics;
import ec.util.*;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;

import java.io.*;
import java.util.Collections;
import java.util.List;

import static yimei.jss.gp.GPRuleEvolutionState.phenoCharacterisation;

/* 
 * MultiObjectiveStatistics.java
 * 
 * Created: Thu Feb 04 2010
 * By: Faisal Abidi and Sean Luke
 *
 */

/*
 * MultiObjectiveStatistics are a SimpleStatistics subclass which overrides the finalStatistics
 * method to output the current Pareto Front in various ways:
 *
 * <ul>
 * <li><p>Every individual in the Pareto Front is written to the end of the statistics log.
 * <li><p>A summary of the objective values of the Pareto Front is written to stdout.
 * <li><p>The objective values of the Pareto Front are written in tabular form to a special
 * Pareto Front file specified with the parameters below.  This file can be easily read by
 * gnuplot or Excel etc. to display the Front (if it's 2D or perhaps 3D).
 * 
 * <p>
 * <b>Parameters</b><br>
 * <table>
 * <tr>
 * <td valign=top><i>base</i>.<tt>front</tt><br>
 * <font size=-1>String (a filename)</font></td>
 * <td valign=top>(The Pareto Front file, if any)</td>
 * </tr>
 * </table>
 */

public class MultiObjectiveStatistics extends SimpleStatistics {
    /**
     * front file parameter
     */
    public static final String P_PARETO_FRONT_FILE = "front";
    public static final String P_SILENT_FRONT_FILE = "silent.front";

    public boolean silentFront;

    /**
     * The pareto front log
     */
    public int frontLog = 0;  // stdout by default

    ArrayList<Double> fitnessNondominated00;
    ArrayList<Double> fitnessNondominated01;
    ArrayList<Double> fitnessNondominated10;
    ArrayList<Double> fitnessNondominated11;

    public static ArrayList<Integer> numNondominatedInds = new ArrayList<>();

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        silentFront = state.parameters.getBoolean(base.push(P_SILENT), null, false);
        // yes, we're stating it a second time.  It's correct logic.
        silentFront = state.parameters.getBoolean(base.push(P_SILENT_FRONT_FILE), null, silentFront);

        File frontFile = state.parameters.getFile(base.push(P_PARETO_FRONT_FILE), null);

        if (silentFront) {
            frontLog = Output.NO_LOGS;
        } else if (frontFile != null) {
            try {
                frontLog = state.output.addLog(frontFile, !compress, compress);
            } catch (IOException i) {
                state.output.fatal("An IOException occurred while trying to create the log " + frontFile + ":\n" + i);
            }
        } else
            state.output.warning("No Pareto Front statistics file specified, printing to stdout at end.", base.push(P_PARETO_FRONT_FILE));
    }


    /**
     * Logs the best individual of the run.
     */
    public void finalStatistics(final EvolutionState state, final int result) {
        bypassFinalStatistics(state, result);  // just call super.super.finalStatistics(...)

        //fzhang 15.11.2018  output seed information in out.stat
    /*    Parameter p;
		// Get the job seed.
		p = new Parameter("seed").push(""+0);
        double jobSeed = state.parameters.getLongWithDefault(p, null, 0);
        state.output.println("seed: "+jobSeed,statisticslog);
        */
        if (doFinal)
            if (state.generation == 1 || state.generation == 11 || state.generation == 21 ||
                    state.generation == 31 || state.generation == 41 || state.generation == 51) {
                state.output.println("\n\n\n PARETO FRONTS", statisticslog);
            }

        //2021.10.4 define the variables for saving the fitness of non-dominated inds---rank = 0
        // =================start=============
        fitnessNondominated00 = new ArrayList<>();
        fitnessNondominated01 = new ArrayList<>();
        fitnessNondominated10 = new ArrayList<>();
        fitnessNondominated11 = new ArrayList<>();
        //================end====================


        for (int s = 0; s < state.population.subpops.length; s++) {
            MultiObjectiveFitness typicalFitness = (MultiObjectiveFitness) (state.population.subpops[s].individuals[0].fitness);
            if (doFinal)
                if (state.generation == 1 || state.generation == 11 || state.generation == 21 ||
                        state.generation == 31 || state.generation == 41 || state.generation == 51) {
                    state.output.println("\n\nPareto Front of Subpopulation " + s, statisticslog);
                }
            // build front

            // Gaofeng Shi 2022.08.25. the final front apply the traditional non domatated sorting.
            double tempAlpha = NSGA2Evaluator.alphaMatrix[0];
            NSGA2Evaluator.alphaMatrix[0] = 0.0;
            ArrayList front = typicalFitness.partitionIntoParetoFront(state.population.subpops[s].individuals, null, null);
            NSGA2Evaluator.alphaMatrix[0] = tempAlpha;

            numNondominatedInds.add(front.size());

            // sort by objective[0]
            Object[] sortedFront = front.toArray();
            QuickSort.qsort(sortedFront, new SortComparator() {
                public boolean lt(Object a, Object b) {
                    return (((MultiObjectiveFitness) (((Individual) a).fitness)).getObjective(0) <
                            (((MultiObjectiveFitness) ((Individual) b).fitness)).getObjective(0));
                }

                public boolean gt(Object a, Object b) {
                    return (((MultiObjectiveFitness) (((Individual) a).fitness)).getObjective(0) >
                            ((MultiObjectiveFitness) (((Individual) b).fitness)).getObjective(0));
                }
            });

            // print out front to statistics log
            if (doFinal)
                for (int i = 0; i < sortedFront.length; i++) {
                    if (state.generation == 1 || state.generation == 11 || state.generation == 21 ||
                            state.generation == 31 || state.generation == 41 || state.generation == 51) {
                        ((Individual) (sortedFront[i])).printIndividualForHumans(state, statisticslog);
                    }
                    //==================================start=======================================
                    if (s == 0) {
                        fitnessNondominated00.add(((MultiObjectiveFitness) ((Individual) (sortedFront[i])).fitness).objectives[0]);
                        fitnessNondominated01.add(((MultiObjectiveFitness) ((Individual) (sortedFront[i])).fitness).objectives[1]);
                    } else if (s == 1) {
                        fitnessNondominated10.add(((MultiObjectiveFitness) ((Individual) (sortedFront[i])).fitness).objectives[0]);
                        fitnessNondominated11.add(((MultiObjectiveFitness) ((Individual) (sortedFront[i])).fitness).objectives[1]);
                    } else {
                        System.out.println("there are only two populations, exceed...");
                    }
                    //===================================end====================================
                }
            // write short version of front out to disk
            if (!silentFront) {
                if (state.generation == 1 || state.generation == 11 || state.generation == 21 ||
                        state.generation == 31 || state.generation == 41 || state.generation == 51) {
                    state.output.println("Generation " + state.generation, frontLog);

                    if (state.population.subpops.length > 1)
                        state.output.println("Subpopulation " + s, frontLog);
                    for (int i = 0; i < sortedFront.length; i++) {
                        Individual ind = (Individual) (sortedFront[i]);
                        MultiObjectiveFitness mof = (MultiObjectiveFitness) (ind.fitness);
                        double[] objectives = mof.getObjectives();

                        String line = "";
                        for (int f = 0; f < objectives.length; f++)
                            line += (objectives[f] + " ");
                        state.output.println(line, frontLog);
                    }
                }
            }

            if (state.generation == 51) {
                writeNumNondominatedIndsTasks(state, numNondominatedInds);
            }
            //2021.10.5 only print the fitness of Pareto front at a numner of generations
            //=============================================start==========================
            if (state.generation == 1 || state.generation == 11 || state.generation == 21 ||
                    state.generation == 31 || state.generation == 41 || state.generation == 51) {
                writeFitnessNonDominatedInds(state, fitnessNondominated00, fitnessNondominated01, fitnessNondominated10, fitnessNondominated11);
            }

        }
    }
    //get seed
    protected long jobSeed;

    //fzhang 25.6.2018 in order to save the rulesize in each generation
    List<Long> aveSeqRulesizeTree0 = new ArrayList<>();
    List<Long> aveRouRulesizeTree1 = new ArrayList<>();
    public void preEvaluationStatistics(final EvolutionState state)
    {
        for(int x=0;x<children.length;x++)
            children[x].preEvaluationStatistics(state);

        //fzhang 17.6.2018  get the seed value
        Parameter p;
        // Get the job seed.
        p = new Parameter("seed").push(""+0);
        jobSeed = state.parameters.getLongWithDefault(p, null, 0);

        // fzhang 15.6.2018 1. save the individual size in population
        // 2. calculate the average size of individuals in population
        // check the average size of sequencing and routing rules in population

        //int indSizePop = 0; // in order to check whether SeqSizePop1 and RouSizePop2 are calculated correctly
        // should be the sum of SeqSizePop1 and RouSizePop2
        long aveSeqSizeTree0 = 0;
        long aveRouSizeTree1 = 0;
/*		for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {
			GPIndividual indi = (GPIndividual) state.population.subpops[0].individuals[ind];
			SeqSizeTree0 += indi.trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
			RouSizeTree1 += indi.trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
		}
		aveSeqSizeTree0 = SeqSizeTree0 / state.population.subpops[0].individuals.length;
		aveRouSizeTree1 = RouSizeTree1 / state.population.subpops[0].individuals.length;

		aveSeqRulesizeTree0.add(aveSeqSizeTree0);
		aveRouRulesizeTree1.add(aveRouSizeTree1);*/

        //2020.2.12 save the average routing and sequencing rule size for all subpops, respectively.
        for (int subpop = 0; subpop < state.population.subpops.length; subpop++){
            //fzhang 15.6.2018  in order to check the average size of sequencing and routing rules in population
            int SeqSizeTree0 = 0;
            int RouSizeTree1 = 0;

            for(int inds = 0; inds < state.population.subpops[subpop].individuals.length; inds++){
                GPIndividual indi = (GPIndividual) state.population.subpops[subpop].individuals[inds];
                SeqSizeTree0 += indi.trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
                RouSizeTree1 += indi.trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
            }
            aveSeqSizeTree0 = SeqSizeTree0 / state.population.subpops[subpop].individuals.length;
            aveRouSizeTree1 = RouSizeTree1 / state.population.subpops[subpop].individuals.length;

            aveSeqRulesizeTree0.add(aveSeqSizeTree0);
            aveRouRulesizeTree1.add(aveRouSizeTree1);
        }

        if(state.generation == state.numGenerations-1) {
            //fzhang  15.6.2018  save the size of rules in each generation
//                File rulesizeFile = new File(out_dir + "/job." + jobSeed + ".aveGenRulesize.csv"); // jobSeed = 0
            File rulesizeFile = new File("job." + jobSeed + ".aveGenRulesize.csv"); // jobSeed = 0


			/*try {
				BufferedWriter writer = new BufferedWriter(new FileWriter(rulesizeFile));
				writer.write("Gen,aveSeqRuleSize,aveRouRuleSize,avePairSize");
				writer.newLine();
				for (int gen = 0; gen < aveSeqRulesizeTree0.size(); gen++) {
					writer.write(gen + "," + aveSeqRulesizeTree0.get(gen) + "," + aveRouRulesizeTree1.get(gen) + "," +
							(aveSeqRulesizeTree0.get(gen) + aveRouRulesizeTree1.get(gen))/2);
					writer.newLine();
				}
				writer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}*/

/*                //2020.2.12 save the information of all the 3 subpops
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter(rulesizeFile));
                    writer.write("Gen,aveSeqRuleSize0,aveRouRuleSize0,aveSeqRuleSize1,aveRouRuleSize1,aveSeqRuleSize2,aveRouRuleSize2");
                    writer.newLine();
                    for (int gen = 0; gen < aveSeqRulesizeTree0.size(); gen+=state.population.subpops.length) {
                        writer.write(gen/state.population.subpops.length + "," + aveSeqRulesizeTree0.get(gen) + "," + aveRouRulesizeTree1.get(gen) + "," +
                                aveSeqRulesizeTree0.get(gen+1) + "," + aveRouRulesizeTree1.get(gen+1) + "," +
                                aveSeqRulesizeTree0.get(gen+2) + "," + aveRouRulesizeTree1.get(gen+2));
                        writer.newLine();
                    }
                    writer.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }*/


            //2020.7.31 save the information of all the 2 subpops
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter(rulesizeFile));
                writer.write("Gen,aveSeqRuleSize0,aveRouRuleSize0,aveSeqRuleSize1,aveRouRuleSize1");
                writer.newLine();
                for (int gen = 0; gen < aveSeqRulesizeTree0.size(); gen+=state.population.subpops.length) {
                    writer.write(gen/state.population.subpops.length + "," + aveSeqRulesizeTree0.get(gen) + "," + aveRouRulesizeTree1.get(gen) );
//                                + "," +
//                                aveSeqRulesizeTree0.get(gen+1) + "," + aveRouRulesizeTree1.get(gen+1));  //original Gaofeng shi. Because it is for multi tasks.
                    writer.newLine();
                }
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

			/*System.out.println(SeqSizeTree0);
			System.out.println(RouSizeTree1);
			System.out.println(aveSeqRulesizeTree0.get(state.generation));
			System.out.println(aveRouRulesizeTree1.get(state.generation));*/

            //fzhang 15.6.2018 in order to check whether SeqSizePop1 and RouSizePop2 are calculated correctly (YES)
	 	/*	for (int pop = 0; pop < state.population.subpops.length; pop++) {
	 			for (int ind = 0; ind < state.population.subpops[pop].individuals.length; ind++) {
	 				indSizePop += state.population.subpops[pop].individuals[ind].size();
	 			}
	 		}
	 		System.out.println(indSizePop);*/
        }
    }



    //2021.10.4 write the information of fitness of non-dominated individuals
    public void writeFitnessNonDominatedInds(final EvolutionState state, ArrayList<Double> fitnessNondominated00, ArrayList<Double> fitnessNondominated01, ArrayList<Double> fitnessNondominated10, ArrayList<Double> fitnessNondominated11) {

        Parameter p;
        p = new Parameter("seed").push("" + 0);
        jobSeed = state.parameters.getIntWithDefault(p, null, 0);
        File distanceFitnessGapFile = new File("job." + jobSeed + "." + (state.generation - 1) + ".fitnessNondominatedInds.csv"); // jobSeed = 0
        ArrayList<Integer> sizeArray = new ArrayList<>();
        sizeArray.add(fitnessNondominated00.size());
        sizeArray.add(fitnessNondominated01.size());
        sizeArray.add(fitnessNondominated10.size());
        sizeArray.add(fitnessNondominated11.size());

        int maxLength = Collections.max(sizeArray);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(distanceFitnessGapFile));
//                writer.write("fitnessNondominated00,fitnessNondominated01,fitnessNondominated10,fitnessNondominated11"); //for multi tasks
            writer.write("fitnessNondominated00,fitnessNondominated01");
            writer.newLine();
            for (int num = 0; num < maxLength; num++) {
                String toWrite = "";
                //write max-flowtime and mean-flowtime
                if (fitnessNondominated00.size() > num) {
                    toWrite += fitnessNondominated00.get(num) + "," + fitnessNondominated01.get(num);
                }
                //Gaofeng shi 2022.09.03 only one task.
//                    else{
//                        toWrite += ",";
//                    }
//
//                    //write max-tardiness and mean-tardiness
//                    if(fitnessNondominated10.size() > num){
//                        toWrite += "," + fitnessNondominated10.get(num) + "," + fitnessNondominated11.get(num);
//                    }
//                    else{
//                        toWrite += ",";
//                    }

                writer.write(toWrite);
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //2021.10.5 the number of non-dominated inds at all generations
    public void writeNumNondominatedIndsTasks(EvolutionState state, ArrayList<Integer> numNondominatedInds) {
        Parameter p;
        // Get the job seed.
        p = new Parameter("seed").push("" + 0);
        jobSeed = state.parameters.getLongWithDefault(p, null, 0);
//		File indSizeForTask = new File(out_dir + "/job." + jobSeed + ".parentPopIndSizeTask.csv");
        File numNondominatedIndsFile = new File("job." + jobSeed + ".numNondominatedInds.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(numNondominatedIndsFile));
//                writer.write("generation,numNondominatedIndsT0,numNondominatedIndsT1"); //for multi tasks
            writer.write("generation,numNondominatedIndsT0");
            writer.newLine();
            for (int cutPoint = 0; cutPoint < numNondominatedInds.size(); cutPoint += state.population.subpops.length) {
//                    writer.write(cutPoint / state.population.subpops.length + "," + numNondominatedInds.get(cutPoint) + "," + numNondominatedInds.get(cutPoint + 1));//for multi tasks
                writer.write(cutPoint / state.population.subpops.length + "," + numNondominatedInds.get(cutPoint));
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Gaofeng Shi 2022.09.06 remove the duplicated individual from front.
    public static ArrayList<Individual> removeDuplicates(ArrayList<Individual> front){

        ArrayList<Individual> newList = new ArrayList<>();
        int count = 0;
        for(int i=0; i < front.size(); i++){
            int[][] frontCRT = phenotypicForSurrogate.muchBetterPhenotypicPopulation(front,phenoCharacterisation); //calculate characterise of front.
            Boolean duplicate = false;
            for (int j=0; j< newList.size(); j++){
                int[][] tempArcCRT = phenotypicForSurrogate.muchBetterPhenotypicPopulation(newList,phenoCharacterisation); //calculate characterise archive temp.
                double  distances = PhenoCharacterisation.distance(frontCRT[i],tempArcCRT[j]); // calculate inds distance.
                if (distances==0) {
                    duplicate = true;
                    count += 1;
                    break;
                }
            }
            if(duplicate == false){
                Individual ind = (Individual)front.get(i).clone();
                newList.add(ind);
            }
        }
//            System.out.println("The duplicated inds: "+count);
        return newList;
    }
}
