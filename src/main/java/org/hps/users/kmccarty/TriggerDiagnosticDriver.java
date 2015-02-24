package org.hps.users.kmccarty;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.hps.readout.ecal.TriggerModule;
import org.hps.readout.ecal.triggerbank.AbstractIntData;
import org.hps.readout.ecal.triggerbank.SSPCluster;
import org.hps.readout.ecal.triggerbank.SSPData;
import org.hps.readout.ecal.triggerbank.SSPPairTrigger;
import org.hps.readout.ecal.triggerbank.SSPSinglesTrigger;
import org.hps.readout.ecal.triggerbank.SSPTrigger;
import org.hps.readout.ecal.triggerbank.TIData;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

public class TriggerDiagnosticDriver extends Driver {
	// Store the LCIO collection names for the needed objects.
	private String bankCollectionName = "TriggerBank";
	private String clusterCollectionName = "EcalClusters";
	
	// Store the lists of parsed objects.
	private TIData tiBank;
	private SSPData sspBank;
	private List<Cluster> reconClusters = new ArrayList<Cluster>();
	private List<SSPCluster> sspClusters;
	private List<List<PairTrigger<Cluster[]>>> reconPairsTriggers = new ArrayList<List<PairTrigger<Cluster[]>>>(2);
	private List<List<PairTrigger<SSPCluster[]>>> sspPairsTriggers = new ArrayList<List<PairTrigger<SSPCluster[]>>>(2);
	private List<List<SinglesTrigger<Cluster>>> reconSinglesTriggers = new ArrayList<List<SinglesTrigger<Cluster>>>(2);
	private List<List<SinglesTrigger<SSPCluster>>> sspSinglesTriggers = new ArrayList<List<SinglesTrigger<SSPCluster>>>(2);
	
	// Trigger modules for performing trigger analysis.
	private TriggerModule[] singlesTrigger = new TriggerModule[2];
	private TriggerModule[] pairsTrigger = new TriggerModule[2];
	
	// Verification settings.
	private int nsa = 100;
	private int nsb = 20;
	private int windowWidth = 200;
	private int hitAcceptance = 1;
	private double energyAcceptance = 0.03;
	private boolean performClusterVerification = true;
	private boolean performSinglesTriggerVerification = true;
	private boolean performPairTriggerVerification = true;
	
	// Efficiency tracking variables.
	private int globalFound = 0;
	private int globalMatched = 0;
	private int globalPosition = 0;
	private int globalEnergy = 0;
	private int globalHitCount = 0;
	
	private int singlesSSPTriggers = 0;
	private int singlesReconMatched = 0;
	private int singlesReconTriggers = 0;
	private int singlesInternalMatched = 0;
	private int singlesReportedTriggers = 0;
	
	private int pairSSPTriggers = 0;
	private int pairReconMatched = 0;
	private int pairReconTriggers = 0;
	private int pairInternalMatched = 0;
	private int pairReportedTriggers = 0;
	
	private int[] globalEnergyMinCut = new int[2];
	private int[] globalEnergyMaxCut = new int[2];
	private int[] globalHitCountCut = new int[2];
	private int[] globalSinglesTimeCut = new int[2];
	
	private int[] globalEnergySumCut = new int[2];
	private int[] globalEnergyDiffCut = new int[2];
	private int[] globalEnergySlopeCut = new int[2];
	private int[] globalCoplanarityCut = new int[2];
	private int[] globalPairTimeCut = new int[2];
	
	// Diagnostic plots.
    private AIDA aida = AIDA.defaultInstance();
    IHistogram1D clusterTimePlotF;
    IHistogram1D clusterEnergyDiffPlotF;
    IHistogram1D clusterHitDiffPlotF;
    IHistogram2D energyHitDiffPercentPlotF;
    
    IHistogram1D clusterTimePlotA;
    IHistogram1D clusterEnergyDiffPlotA;
    IHistogram1D clusterHitDiffPlotA;
    IHistogram2D energyHitDiffPercentPlotA;
    
    IHistogram1D reconClusterTotalEnergyPlot;
    IHistogram1D reconClusterHitCountPlot;
    IHistogram1D reconClusterTimePlot;
    IHistogram2D reconClusterDistribution;
    IHistogram1D sspClusterTotalEnergyPlot;
    IHistogram1D sspClusterHitCountPlot;
    IHistogram2D sspClusterDistribution;
	
    IHistogram1D cyclesRemoved;
    IHistogram2D cyclesRemovedEnergy;
    IHistogram2D cyclesRemovedSeed;
    IHistogram2D[] hitsInCycle = new IHistogram2D[9];
    IHistogram1D repeatedHitsPlot;
    IHistogram2D cyclesRemovedEnergyPercent;
    
    // Internal state variables.
    private static final int STATE_CLUSTER_UNDEFINED = -1;
    private static final int STATE_CLUSTER_FAIL_ENERGY = 1;
    private static final int STATE_CLUSTER_FAIL_HIT_COUNT = 2;
    private static final int STATE_CLUSTER_SUCCESS_MATCH = 3;
    
    // Verbose settings.
    private boolean clusterFail = false;
    private boolean singlesFail = false;
    private boolean pairFail = false;
    private boolean verbose = false;
    private boolean printClusterFail = true;
    private boolean printSinglesTriggerFail = true;
    private boolean printPairTriggerFail = true;
    private StringBuffer outputBuffer = new StringBuffer();
    
	/**
	 * Define the trigger modules. This should be replaced by parsing
	 * the DAQ configuration at some point.
	 */
	@Override
	public void startOfData() {
		// Instantiate the diagnostic plots.
	    energyHitDiffPercentPlotF = aida.histogram2D("Trigger Diagnostics :: Failed Cluster Energy-Hit Difference (Percentage)", 11, -5, 6, 25, 0.75, 1.25);
	    clusterTimePlotF = aida.histogram1D("Trigger Diagnostics :: Failed Cluster Time Distribution", 2 * windowWidth / 5, 0, windowWidth);
	    clusterHitDiffPlotF = aida.histogram1D("Trigger Diagnostics :: Failed Cluster Hit Count Difference Distribution", 11, -5, 6);
	    clusterEnergyDiffPlotF = aida.histogram1D("Trigger Diagnostics :: Failed Cluster Energy Difference Distribution", 25, 0.75, 1.25);
	    
	    energyHitDiffPercentPlotA = aida.histogram2D("Trigger Diagnostics :: All Cluster Energy-Hit Difference (Percentage)", 11, -5, 6, 25, 0.75, 1.25);
	    clusterTimePlotA = aida.histogram1D("Trigger Diagnostics :: All Cluster Time Distribution", 2 * windowWidth / 5, 0, windowWidth);
	    clusterHitDiffPlotA = aida.histogram1D("Trigger Diagnostics :: All Cluster Hit Count Difference Distribution", 11, -5, 6);
	    clusterEnergyDiffPlotA = aida.histogram1D("Trigger Diagnostics :: All Cluster Energy Difference Distribution", 25, 0.75, 1.25);
	    
	    reconClusterTotalEnergyPlot = aida.histogram1D("Trigger Diagnostics :: Recon Cluster Total Energy Distribution", 55, 0.0, 2.2);
	    reconClusterHitCountPlot= aida.histogram1D("Trigger Diagnostics :: Recon Cluster Hit Count Distribution", 10, -0.5, 9.5);
	    reconClusterTimePlot = aida.histogram1D("Trigger Diagnostics :: Recon Cluster Time Distribution", 2 * windowWidth / 5, 0, windowWidth);
	    reconClusterDistribution = aida.histogram2D("Trigger Diagnostics :: Recon Cluster Distribution", 46, -23, 23, 11, -5.5, 5.5);
	    sspClusterTotalEnergyPlot = aida.histogram1D("Trigger Diagnostics :: SSP Cluster Total Energy Distribution", 55, 0.0, 2.2);
	    sspClusterHitCountPlot= aida.histogram1D("Trigger Diagnostics :: SSP Cluster Hit Count Distribution", 10, -0.5, 9.5);
	    sspClusterDistribution = aida.histogram2D("Trigger Diagnostics :: SSP Cluster Distribution", 46, -23, 23, 11, -5.5, 5.5);
		
	    cyclesRemoved = aida.histogram1D("Trigger Diagnostics :: Hit Cycles Removed From Seed", 10, -0.5, 9.5);
	    cyclesRemovedEnergy = aida.histogram2D("Trigger Diagnostics :: Energy Percent vs. Clock Cycles Removed", 10, -0.5, 9.5, 100, 0.0, 1.0);
	    cyclesRemovedSeed = aida.histogram2D("Trigger Diagnostics :: Seed Percent vs. Clock Cycles Removed", 10, -0.5, 9.5, 100, 0.0, 1.5);
	    
	    repeatedHitsPlot = aida.histogram1D("Trigger Diagnostics :: Repeated Crystal Positions in Cluster", 16, -0.5, 15.5);
	    cyclesRemovedEnergyPercent = aida.histogram2D("Trigger Diagnostics :: Cycle Energy Percent vs. Clock Cycles Removed", 10, -0.5, 9.5, 100, 0.0, 1.0);
	    for(int i = 0; i < hitsInCycle.length; i++) {
	    	hitsInCycle[i] = aida.histogram2D("Trigger Diagnostics :: Hits in Cycle vs. Cycles Removed (Cluster Size " + i + ")", 5, -0.5, 4.5, 10, -0.5, 9.5);
	    }
	    
		// Print the cluster verification header.
		System.out.println();
		System.out.println();
		System.out.println("======================================================================");
		System.out.println("=== Cluster/Trigger Verification Settings ============================");
		System.out.println("======================================================================");
		
		// Define the first singles trigger.
		singlesTrigger[0] = new TriggerModule();
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.500);
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 8.191);
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 0);
		
		// Define the second singles trigger.
		singlesTrigger[1] = new TriggerModule();
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.000);
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 8.191);
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 0);
		
		// Define the first pairs trigger.
		pairsTrigger[0] = new TriggerModule();
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.000);
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 8.191);
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 0);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SUM_LOW, 0.000);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SUM_HIGH, 8.191);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_DIFFERENCE_HIGH, 8.191);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_LOW, 0.000);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_F, 0.001);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_COPLANARITY_HIGH, 180);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_TIME_COINCIDENCE, 8);
		
		// Define the second pairs trigger.
		pairsTrigger[1] = new TriggerModule();
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.000);
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 8.191);
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 0);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SUM_LOW, 0.000);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SUM_HIGH, 8.191);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_DIFFERENCE_HIGH, 8.191);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_LOW, 0.000);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_F, 0.001);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_COPLANARITY_HIGH, 180);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_TIME_COINCIDENCE, 8);
		
		// Instantiate the triggers lists.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			reconPairsTriggers.add(new ArrayList<PairTrigger<Cluster[]>>());
			sspPairsTriggers.add(new ArrayList<PairTrigger<SSPCluster[]>>());
			reconSinglesTriggers.add(new ArrayList<SinglesTrigger<Cluster>>());
			sspSinglesTriggers.add(new ArrayList<SinglesTrigger<SSPCluster>>());
		}
		
		// Print the initial settings.
		logSettings();
	}
	
	/**
	 * Prints the total run statistics.
	 */
	@Override
	public void endOfData() {
		// Print the cluster/trigger verification header.
		System.out.println();
		System.out.println();
		System.out.println("======================================================================");
		System.out.println("=== Cluster/Trigger Verification Results =============================");
		System.out.println("======================================================================");
		
		// Print the cluster verification data.
		System.out.println("Cluster Verification:");
		System.out.printf("\tRecon Clusters     :: %d%n", globalFound);
		System.out.printf("\tClusters Matched   :: %d%n", globalMatched);
		System.out.printf("\tFailed (Position)  :: %d%n", globalPosition);
		System.out.printf("\tFailed (Energy)    :: %d%n", globalEnergy);
		System.out.printf("\tFailed (Hit Count) :: %d%n", globalHitCount);
		if(globalFound == 0) { System.out.printf("\tCluster Efficiency :: N/A%n"); }
		else { System.out.printf("\tCluster Efficiency :: %7.3f%%%n", 100.0 * globalMatched / globalFound); }
		
		// Print the singles trigger verification data.
		System.out.println();
		System.out.println("Singles Trigger Verification:");
		System.out.printf("\tSSP Cluster Sim Triggers   :: %d%n", singlesSSPTriggers);
		System.out.printf("\tRecon Cluster Sim Triggers :: %d%n", singlesReconTriggers);
		System.out.printf("\tSSP Reported Triggers      :: %d%n", singlesReportedTriggers);

		if(singlesSSPTriggers == 0) {
			System.out.printf("\tInternal Efficiency        :: %d / %d (N/A)%n",
					singlesInternalMatched, singlesSSPTriggers);
		} else {
			System.out.printf("\tInternal Efficiency        :: %d / %d (%7.3f%%)%n",
					singlesInternalMatched, singlesSSPTriggers, (100.0 * singlesInternalMatched / singlesSSPTriggers));
		}
		if(singlesReconTriggers == 0) {
			System.out.printf("\tTrigger Efficiency         :: %d / %d (N/A)%n",
					singlesReconMatched, singlesReconTriggers);
		} else {
			System.out.printf("\tTrigger Efficiency         :: %d / %d (%7.3f%%)%n",
					singlesReconMatched, singlesReconTriggers, (100.0 * singlesReconMatched / singlesReconTriggers));
		}
		
		// Print the individual cut performances.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			System.out.println();
			System.out.printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(singlesSSPTriggers == 0) {
				System.out.printf("\tCluster Energy Lower Bound :: %d / %d%n", globalEnergyMinCut[triggerNum], singlesSSPTriggers);
				System.out.printf("\tCluster Energy Upper Bound :: %d / %d%n", globalEnergyMaxCut[triggerNum], singlesSSPTriggers);
				System.out.printf("\tCluster Hit Count          :: %d / %d%n", globalHitCountCut[triggerNum], singlesSSPTriggers);
			} else {
				System.out.printf("\tCluster Energy Lower Bound :: %d / %d (%7.3f%%)%n",
						globalEnergyMinCut[triggerNum], singlesSSPTriggers, (100.0 * globalEnergyMinCut[triggerNum] / singlesSSPTriggers));
				System.out.printf("\tCluster Energy Upper Bound :: %d / %d (%7.3f%%)%n",
						globalEnergyMaxCut[triggerNum], singlesSSPTriggers, (100.0 * globalEnergyMaxCut[triggerNum] / singlesSSPTriggers));
				System.out.printf("\tCluster Hit Count          :: %d / %d (%7.3f%%)%n",
						globalHitCountCut[triggerNum], singlesSSPTriggers, (100.0 * globalHitCountCut[triggerNum] / singlesSSPTriggers));
			}
		}
		
		// Print the pair trigger verification data.
		System.out.println();
		System.out.println("Pair Trigger Verification:");
		System.out.printf("\tSSP Cluster Sim Triggers   :: %d%n", pairSSPTriggers);
		System.out.printf("\tRecon Cluster Sim Triggers :: %d%n", pairReconTriggers);
		System.out.printf("\tSSP Reported Triggers      :: %d%n", pairReportedTriggers);

		if(pairSSPTriggers == 0) {
			System.out.printf("\tInternal Efficiency        :: %d / %d (N/A)%n",
					pairInternalMatched, pairSSPTriggers);
		} else {
			System.out.printf("\tInternal Efficiency        :: %d / %d (%7.3f%%)%n",
					pairInternalMatched, pairSSPTriggers, (100.0 * pairInternalMatched / pairSSPTriggers));
		}
		if(pairReconTriggers == 0) {
			System.out.printf("\tTrigger Efficiency         :: %d / %d (N/A)%n",
					pairReconMatched, pairReconTriggers);
		} else {
			System.out.printf("\tTrigger Efficiency         :: %d / %d (%7.3f%%)%n",
					pairReconMatched, pairReconTriggers, (100.0 * pairReconMatched / pairReconTriggers));
		}
		
		// Print the individual cut performances.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			System.out.println();
			System.out.printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(pairSSPTriggers == 0) {
				System.out.printf("\tPair Energy Sum            :: %d / %d%n", globalEnergySumCut[triggerNum], pairSSPTriggers);
				System.out.printf("\tPair Energy Difference     :: %d / %d%n", globalEnergyDiffCut[triggerNum], pairSSPTriggers);
				System.out.printf("\tPair Energy Slope          :: %d / %d%n", globalEnergySlopeCut[triggerNum], pairSSPTriggers);
				System.out.printf("\tPair Coplanarity           :: %d / %d%n", globalCoplanarityCut[triggerNum], pairSSPTriggers);
				System.out.printf("\tPair Trigger Time          :: %d / %d%n", globalPairTimeCut[triggerNum], pairSSPTriggers);
			} else {
				System.out.printf("\tPair Energy Sum            :: %d / %d (%7.3f%%)%n",
						globalEnergySumCut[triggerNum], pairSSPTriggers, (100.0 * globalEnergySumCut[triggerNum] / pairSSPTriggers));
				System.out.printf("\tPair Energy Difference     :: %d / %d (%7.3f%%)%n",
						globalEnergyDiffCut[triggerNum], pairSSPTriggers, (100.0 * globalEnergyDiffCut[triggerNum] / pairSSPTriggers));
				System.out.printf("\tPair Energy Slope          :: %d / %d (%7.3f%%)%n",
						globalEnergySlopeCut[triggerNum], pairSSPTriggers, (100.0 * globalEnergySlopeCut[triggerNum] / pairSSPTriggers));
				System.out.printf("\tPair Coplanarity           :: %d / %d (%7.3f%%)%n",
						globalCoplanarityCut[triggerNum], pairSSPTriggers, (100.0 * globalCoplanarityCut[triggerNum] / pairSSPTriggers));
				System.out.printf("\tPair Trigger Time          :: %d / %d (%7.3f%%)%n",
						globalPairTimeCut[triggerNum], pairSSPTriggers, (100.0 * globalPairTimeCut[triggerNum] / pairSSPTriggers));
			}
		}
	}
	
	/**
	 * Gets the banks and clusters from the event.
	 */
	@Override
	public void process(EventHeader event) {
		// ==========================================================
		// ==== Initialize the Event ================================
		// ==========================================================
		
		// Reset the output buffer and print flags.
		outputBuffer = new StringBuffer();
		clusterFail = false;
		singlesFail = false;
		pairFail = false;
		
		println("======================================================================");
		println("==== Cluster/Trigger Verification ====================================");
		println("======================================================================");
		
		// Clear the list of triggers from previous events.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			sspSinglesTriggers.get(triggerNum).clear();
			reconSinglesTriggers.get(triggerNum).clear();
			sspPairsTriggers.get(triggerNum).clear();
			reconPairsTriggers.get(triggerNum).clear();
		}
		
		
		
		// ==========================================================
		// ==== Obtain Reconstructed Clusters =======================
		// ==========================================================
		
		// Get the reconstructed clusters.
		if(event.hasCollection(Cluster.class, clusterCollectionName)) {
			// Get the reconstructed clusters.
			List<Cluster> allClusters = event.get(Cluster.class, clusterCollectionName);
			
			// Keep only the clusters that can be verified.
			println();
			println("Process cluster for verifiability:");
			reconClusters.clear();
			for(Cluster reconCluster : allClusters) {
				// Check that the cluster is within the safe region of the
				// FADC readout window. If it is not, it will likely have
				// inaccurate energy or hit values and may not produce the
				// expected results.
				printf("\t%s", clusterToString(reconCluster));
				if(isVerifiable(reconCluster)) {
					reconClusters.add(reconCluster);
					println(" [  verifiable  ]");
					
					int repeatedHits = 0;
					double seedTime = getReconTime(reconCluster);
					Map<Integer, Double> cycleMap = new HashMap<Integer, Double>();
					Map<Integer, Integer> hitCountMap = new HashMap<Integer, Integer>();
					Set<Point> hasHitSet = new HashSet<Point>();
					for(CalorimeterHit hit : reconCluster.getCalorimeterHits()) {
						double timeDiff = (hit.getTime() - seedTime) / 4.0;
						cyclesRemoved.fill(timeDiff);
						cyclesRemovedEnergy.fill(timeDiff, hit.getCorrectedEnergy() / reconCluster.getEnergy());
						cyclesRemovedSeed.fill(timeDiff, hit.getCorrectedEnergy() / reconCluster.getCalorimeterHits().get(0).getCorrectedEnergy());
						Point hitLocation = new Point(hit.getIdentifierFieldValue("ix"), hit.getIdentifierFieldValue("iy"));
						if(hasHitSet.contains(hitLocation)) { repeatedHits++; }
						else { hasHitSet.add(hitLocation); }
						Integer positionHits = hitCountMap.get((int) timeDiff);
						if(positionHits != null) { hitCountMap.put((int) timeDiff, positionHits + 1); }
						else { hitCountMap.put((int) timeDiff, 1); }
						Double cycleEnergy = cycleMap.get((int) timeDiff);
						if(cycleEnergy != null) { cycleMap.put((int) timeDiff, cycleEnergy + hit.getCorrectedEnergy()); }
						else { cycleMap.put((int) timeDiff, hit.getCorrectedEnergy()); }
					}
					repeatedHitsPlot.fill(repeatedHits);
					for(Entry<Integer, Double> entry : cycleMap.entrySet()) {
						cyclesRemovedEnergyPercent.fill(entry.getKey(), entry.getValue());
					}

					int hitCount = reconCluster.getCalorimeterHits().size() - 1;
					if(hitCount >= 0 && hitCount < 9) {
						for(Entry<Integer, Integer> entry : hitCountMap.entrySet()) {
							hitsInCycle[hitCount].fill(entry.getKey(), entry.getValue());
						}
					}
					
				} else { println(" [ unverifiable ]"); }
				
			}
			
			// Output the number of verifiable clusters found.
			if(reconClusters.size() == 1) { println("1 verifiable reconstructed cluster found."); }
			else { printf("%d verifiable reconstructed clusters found.%n", reconClusters.size()); }
			
			// Output the number of unverifiable clusters found.
			int unverifiableClusters = allClusters.size() - reconClusters.size();
			if(unverifiableClusters == 1) { println("1 unverifiable reconstructed cluster found."); }
			else { printf("%d unverifiable reconstructed clusters found.%n", unverifiableClusters); }
		} else {
			reconClusters = new ArrayList<Cluster>(0);
			printf("No reconstructed clusters were found for collection \"%s\" in this event.%n", clusterCollectionName);
		}
		
		
		
		// ==========================================================
		// ==== Obtain SSP and TI Banks =============================
		// ==========================================================
		
		// Get the SSP clusters.
		if(event.hasCollection(GenericObject.class, bankCollectionName)) {
			// Get the bank list.
			List<GenericObject> bankList = event.get(GenericObject.class, bankCollectionName);
			
			// Search through the banks and get the SSP and TI banks.
			for(GenericObject obj : bankList) {
				// If this is an SSP bank, parse it.
				if(AbstractIntData.getTag(obj) == SSPData.BANK_TAG) {
					sspBank = new SSPData(obj);
				}
				
				// Otherwise, if this is a TI bank, parse it.
				else if(AbstractIntData.getTag(obj) == TIData.BANK_TAG) {
					tiBank = new TIData(obj);
					
					if(tiBank.isPulserTrigger()) { println("Trigger type :: Pulser"); }
					else if(tiBank.isSingle0Trigger() || tiBank.isSingle1Trigger()) { println("Trigger type :: Singles"); }
					else if(tiBank.isPair0Trigger() || tiBank.isPair1Trigger()) { println("Trigger type :: Pair"); }
					else if(tiBank.isCalibTrigger()) { println("Trigger type :: Cosmic"); }
				}
			}
			
			// If there is an SSP bank, get the list of SSP clusters.
			if(sspBank != null) {
				sspClusters = sspBank.getClusters();
				if(sspClusters.size() == 1) {
					println("1 SSP cluster found.");
				} else {
					printf("%d SSP clusters found.%n", sspClusters.size());
				}
			}
		}
		
		
		
		// ==========================================================
		// ==== Establish Event Integrity ===========================
		// ==========================================================
		
		// Check that all of the required objects are present.
		if(sspBank == null) {
			println("No SSP bank found for this event. No verification will be performed.");
			return;
		} if(tiBank == null) {
			println("No TI bank found for this event. No verification will be performed.");
			return;
		}
		
		
		
		// ==========================================================
		// ==== Perform Event Verification ==========================
		// ==========================================================
		
		// Perform the cluster verification step.
		if(performClusterVerification) { clusterVerification(); }
		
		// Construct lists of triggers for the SSP clusters and the
		// reconstructed clusters.
		if(performSinglesTriggerVerification) {
			constructSinglesTriggers();
			singlesTriggerVerification();
		}
		if(performPairTriggerVerification) {
			constructPairTriggers();
			pairTriggerVerification();
		}
		
		
		
		// ==========================================================
		// ==== Perform Event Write-Out =============================
		// ==========================================================
		
		if(verbose ||(clusterFail && printClusterFail) ||
				(singlesFail && printSinglesTriggerFail) ||
				(pairFail && printPairTriggerFail)) {
			System.out.println(outputBuffer.toString());
		}
	}
	
	public void setPrintOnClusterFailure(boolean state) {
		printClusterFail = state;
	}
	
	public void setPrintOnSinglesFailure(boolean state) {
		printSinglesTriggerFail = state;
	}
	
	public void setPrintOnPairFailure(boolean state) {
		printPairTriggerFail = state;
	}
	
	public void setVerbose(boolean state) {
		verbose = state;
	}
	
	/**
	 * Attempts to match all reconstructed clusters that are safely
	 * within the integration window with clusters reported by the SSP.
	 * Method also tracks the ratio of valid reconstructed clusters to
	 * matches found.<br/>
	 * <br/>
	 * Note that unmatched SSP clusters are ignored. Since these may
	 * or may not correspond to reconstructed clusters that occur in
	 * the forbidden time region, it is impossible to say whether or
	 * not these legitimately failed to match or not.
	 */
	private void clusterVerification() {
		// ==========================================================
		// ==== Initialize Cluster Verification =====================
		// ==========================================================
		
		// Print the cluster verification header.
		println();
		println();
		println("======================================================================");
		println("=== Cluster Verification =============================================");
		println("======================================================================");
		
		// If there are no reconstructed clusters, than there is nothing
		// that can be verified.
		if(reconClusters.isEmpty()) {
			println("No reconstructed clusters are present. Skipping event...");
			return;
		}
		
		// Track the number of cluster pairs that were matched and that
		// failed by failure type.
		int eventMatched = 0;
		int eventPosition = 0;
		int eventEnergy = 0;
		int eventHitCount = 0;
		
		// Track the matched clusters.
		StringBuffer eventMatchedText = new StringBuffer();
		
		
		
		// ==========================================================
		// ==== Produce the Cluster Position Mappings ===============
		// ==========================================================
		
		// Create maps to link cluster position to the list of clusters
		// that were found at that location.
		Map<Point, List<Cluster>> reconClusterMap = new HashMap<Point, List<Cluster>>(reconClusters.size());
		Map<Point, List<SSPCluster>> sspClusterMap = new HashMap<Point, List<SSPCluster>>(reconClusters.size());
		
		// Populate the reconstructed cluster map.
		for(Cluster reconCluster : reconClusters) {
			// Get the cluster position.
			Point position = new Point(getReconXIndex(reconCluster), getReconYIndex(reconCluster));
			
			// Get the list for this cluster position.
			List<Cluster> reconList = reconClusterMap.get(position);
			if(reconList == null) {
				reconList = new ArrayList<Cluster>();
				reconClusterMap.put(position, reconList);
			}
			
			// Add the cluster to the list.
			reconList.add(reconCluster);
		}
		
		// Populate the SSP cluster map.
		for(SSPCluster sspCluster : sspClusters) {
			// Get the cluster position.
			Point position = new Point(sspCluster.getXIndex(), sspCluster.getYIndex());
			
			// Get the list for this cluster position.
			List<SSPCluster> sspList = sspClusterMap.get(position);
			if(sspList == null) {
				sspList = new ArrayList<SSPCluster>();
				sspClusterMap.put(position, sspList);
			}
			
			// Add the cluster to the list.
			sspList.add(sspCluster);
		}
		
		
		
		// ==========================================================
		// ==== Perform Cluster Matching ============================
		// ==========================================================
		
		// For each reconstructed cluster, attempt to match the clusters
		// with SSP clusters at the same position.
		positionLoop:
		for(Entry<Point, List<Cluster>> clusterSet : reconClusterMap.entrySet()) {
			// Get the reconstructed and SSP clusters at this position.
			List<Cluster> reconList = clusterSet.getValue();
			List<SSPCluster> sspList = sspClusterMap.get(clusterSet.getKey());
			
			// Print the crystal position header.
			println();
			printf("Considering clusters at (%3d, %3d)%n", clusterSet.getKey().x, clusterSet.getKey().y);
			
			// If there are no SSP clusters, then matching fails by
			// reason of position. The remainder of the loop may be
			// skipped, since there is nothing to check.
			if(sspList == null || sspList.isEmpty()) {
				eventPosition += reconList.size();
				clusterFail = true;
				continue positionLoop;
			}
			
			// If there are more reconstructed clusters than there are
			// SSP clusters, than a number equal to the difference must
			// fail by means of positions.
			if(sspList.size() < reconList.size()) {
				clusterFail = true;
				eventPosition += (sspList.size() - reconList.size());
			}
			
			// Get all possible permutations of SSP clusters.
			List<List<Pair>> permutations = getPermutations(reconList, sspList);
			
			// Print the information for this crystal position.
			printf("\tRecon Clusters :: %d%n", reconList.size());
			printf("\tSSP Clusters   :: %d%n", sspList.size());
			printf("\tPermutations   :: %d%n", permutations.size());
			
			// Track the best results found over all permutations.
			int positionMatched = -1;
			int postionEnergy = -1;
			int postionHitCount = -1;
			StringBuffer positionMatchedText = new StringBuffer();
			
			// Track the plotted values for the current best permutation.
			List<Double> positionAllTimes = new ArrayList<Double>();
			List<Integer> positionAllHitDifference = new ArrayList<Integer>();
			List<Double> positionAllEnergyDiffPercent = new ArrayList<Double>();
			List<Double> positionFailTimes = new ArrayList<Double>();
			List<Integer> positionFailHitDifference = new ArrayList<Integer>();
			List<Double> positionFailEnergyDiffPercent = new ArrayList<Double>();
			
			// Iterate over the permutations and find the permutation
			// that produces the best possible result when compared to
			// the reconstructed clusters.
			int permIndex = 0;
			for(List<Pair> pairs : permutations) {
				// Update the current permutation number.
				permIndex++;
				
				// Track the results of this permutation.
				int permutationMatched = 0;
				int permutationEnergy = 0;
				int permutationHitCount = 0;
				StringBuffer permutationMatchedText = new StringBuffer();
				
				// Track the plot values for this permutation.
				List<Double> permutationAllTimes = new ArrayList<Double>();
				List<Integer> permutationAllHitDifference = new ArrayList<Integer>();
				List<Double> permutationAllEnergyDiffPercent = new ArrayList<Double>();
				List<Double> permutationFailTimes = new ArrayList<Double>();
				List<Integer> permutationFailHitDifference = new ArrayList<Integer>();
				List<Double> permutationFailEnergyDiffPercent = new ArrayList<Double>();
				
				// Try to match each pair.
				pairLoop:
				for(Pair pair : pairs) {
					// Track the state of the current pair.
					int pairState = STATE_CLUSTER_UNDEFINED;
					
					// Print the current reconstructed/SSP cluster pair.
					printf("\tP%d :: %s --> %s", permIndex,
							pair.reconCluster == null ? "None" : clusterToString(pair.reconCluster),
							pair.sspCluster == null ? "None" : clusterToString(pair.sspCluster));
					
					// If either cluster in the pair is null, there
					// are not enough clusters to perform this match.
					if(pair.reconCluster == null || pair.sspCluster == null) {
						printf(" [ %18s ]%n", "failure: unpaired cluster");
						continue pairLoop;
					}
					
					// Check if the reconstructed cluster has an energy
					// within the allotted threshold of the SSP cluster.
					if(pair.sspCluster.getEnergy() >= pair.reconCluster.getEnergy() * (1 - energyAcceptance) &&
							pair.sspCluster.getEnergy() <= pair.reconCluster.getEnergy() * (1 + energyAcceptance)) {
						// Check that the hit count of the reconstructed
						// is within the allotted threshold of the SSP
						// cluster.
						if(pair.sspCluster.getHitCount() >= pair.reconCluster.getCalorimeterHits().size() - hitAcceptance &&
								pair.sspCluster.getHitCount() <= pair.reconCluster.getCalorimeterHits().size() + hitAcceptance) {
							// Designate the pair as a match.
							pairState = STATE_CLUSTER_SUCCESS_MATCH;
						} else { pairState = STATE_CLUSTER_FAIL_HIT_COUNT; }
					} else { pairState = STATE_CLUSTER_FAIL_ENERGY; }
					
					// Update the permutation variables based on the
					// pair state.
					if(pairState == STATE_CLUSTER_SUCCESS_MATCH) {
						// Having passed the position, energy, and hit
						// count tests, this cluster pair is designated
						// a match.
						permutationMatched++;
						
						// Output the matched cluster to a string to be
						// printed in the event summary.
						permutationMatchedText.append(String.format("\t%s --> %s%n",
								clusterToString(pair.reconCluster),
								clusterToString(pair.sspCluster)));
						printf(" [ %18s ]%n", "success: matched");
					} else {
						// Update the tracker variable and print the
						// outcome appropriate to the failure type.
						if(pairState == STATE_CLUSTER_FAIL_ENERGY) {
							permutationEnergy++;
							printf(" [ %18s ]%n", "failure: energy");
						} else if(pairState == STATE_CLUSTER_FAIL_HIT_COUNT) {
							permutationHitCount++;
							printf(" [ %18s ]%n", "failure: hit count");
						} else {
							throw new IllegalArgumentException("Unknown cluster match state!");
						}
						
						// Track the plotted values for failed pairs.
						permutationFailTimes.add(getReconTime(pair.reconCluster));
						permutationFailHitDifference.add(getHitCountDifference(pair.sspCluster, pair.reconCluster));
						permutationFailEnergyDiffPercent.add(getEnergyPercentDifference(pair.sspCluster, pair.reconCluster));
					}
					
					// Track the plotted values for all pairs.
					permutationAllTimes.add(getReconTime(pair.reconCluster));
					permutationAllHitDifference.add(getHitCountDifference(pair.sspCluster, pair.reconCluster));
					permutationAllEnergyDiffPercent.add(getEnergyPercentDifference(pair.sspCluster, pair.reconCluster));
				} // End Pair Loop
				
				// Print the results of the permutation.
				printf("\t\tPermutation Matched   :: %d%n", permutationMatched);
				printf("\t\tPermutation Energy    :: %d%n", permutationEnergy);
				printf("\t\tPermutation Hit Count :: %d%n", permutationHitCount);
				
				// Check whether the results from this permutation
				// exceed the quality of the last best results. A
				// greater number of matches is always better. If the
				// matches are the same, select the one with fewer
				// failures due to energy.
				if((permutationMatched > positionMatched) || (permutationMatched == positionMatched && permutationHitCount < postionEnergy)) {
					// Update the statistics.
					positionMatched = permutationMatched;
					postionEnergy = permutationEnergy;
					postionHitCount = permutationHitCount;
					positionMatchedText = permutationMatchedText;
					
					// Set the plot values.
					positionAllTimes = permutationAllTimes;
					positionAllHitDifference = permutationAllHitDifference;
					positionAllEnergyDiffPercent = permutationAllEnergyDiffPercent;
					positionFailTimes = permutationFailTimes;
					positionFailHitDifference = permutationFailHitDifference;
					positionFailEnergyDiffPercent = permutationFailEnergyDiffPercent;
				}
			} // End Permutation Loop
			
			// Print the final results for the position.
			printf("\tPosition Matched   :: %d%n", positionMatched);
			printf("\tPosition Energy    :: %d%n", postionEnergy);
			printf("\tPosition Hit Count :: %d%n", postionHitCount);
			
			// Add the results from the best-matched permutation
			// to the event efficiency results.
			eventMatched += positionMatched;
			eventEnergy += postionEnergy;
			eventHitCount += postionHitCount;
			eventMatchedText.append(positionMatchedText.toString());
			
			// Update the plots.
			for(int index = 0; index < positionFailHitDifference.size(); index++) {
				// Update the fail plots.
				clusterTimePlotF.fill(positionFailTimes.get(index));
				clusterHitDiffPlotF.fill(positionFailHitDifference.get(index));
				clusterEnergyDiffPlotF.fill(positionFailEnergyDiffPercent.get(index));
				energyHitDiffPercentPlotF.fill(positionFailHitDifference.get(index), positionFailEnergyDiffPercent.get(index));
				
				// Update the all plots.
				clusterTimePlotA.fill(positionAllTimes.get(index));
				clusterHitDiffPlotA.fill(positionAllHitDifference.get(index));
				clusterEnergyDiffPlotA.fill(positionAllEnergyDiffPercent.get(index));
				energyHitDiffPercentPlotA.fill(positionAllHitDifference.get(index), positionAllEnergyDiffPercent.get(index));
			}
		} // End Crystal Position Loop
		
		// Add the event results to the global results.
		globalMatched += eventMatched;
		globalPosition += eventPosition;
		globalEnergy += eventEnergy;
		globalHitCount += eventHitCount;
		globalFound += reconClusters.size();
		
		
		
		// ==========================================================
		// ==== Output Event Summary ================================
		// ==========================================================
		
		// Print the valid reconstructed clusters and populate their
		// distribution graphs.
		println();
		println("Verified Reconstructed Clusters:");
		if(!reconClusters.isEmpty()) {
			for(Cluster reconCluster : reconClusters) {
				printf("\t%s%n", clusterToString(reconCluster));
			    reconClusterTotalEnergyPlot.fill(reconCluster.getEnergy());
			    reconClusterHitCountPlot.fill(getReconHitCount(reconCluster));
			    reconClusterTimePlot.fill(getReconTime(reconCluster));
			    int ix = getReconXIndex(reconCluster);
			    reconClusterDistribution.fill(ix > 0 ? ix - 1 : ix, getReconYIndex(reconCluster));
			}
		} else { println("\tNone"); }
		
		// Print the SSP clusters and populate their distribution graphs.
		println("SSP Clusters:");
		if(!sspClusters.isEmpty()) {
			for(SSPCluster sspCluster : sspClusters) {
				printf("\t%s%n", clusterToString(sspCluster));
			    sspClusterTotalEnergyPlot.fill(sspCluster.getEnergy());
			    sspClusterHitCountPlot.fill(sspCluster.getHitCount());
			    int ix = sspCluster.getXIndex();
			    sspClusterDistribution.fill(ix > 0 ? ix - 1 : ix, sspCluster.getYIndex());
			}
		} else { println("\tNone"); }
		
		// Print the matched clusters.
		println("Matched Clusters:");
		if(eventMatchedText.length() != 0) {
			print(eventMatchedText.toString());
		} else { println("\tNone"); }
		
		// Print event statistics.
		println();
		println("Event Statistics:");
		printf("\tRecon Clusters     :: %d%n", reconClusters.size());
		printf("\tClusters Matched   :: %d%n", eventMatched);
		printf("\tFailed (Position)  :: %d%n", eventPosition);
		printf("\tFailed (Energy)    :: %d%n", eventEnergy);
		printf("\tFailed (Hit Count) :: %d%n", eventHitCount);
		printf("\tCluster Efficiency :: %3.0f%%%n", 100.0 * eventMatched / reconClusters.size());
		
		// Note whether there was a cluster match failure.
		if(eventMatched - reconClusters.size() != 0) {
			clusterFail = true;
		}
	}
	
	/**
	 * Checks triggers simulated on SSP clusters against the SSP bank's
	 * reported triggers to verify that the trigger is correctly applying
	 * cuts to the clusters it sees. Additionally compares triggers
	 * simulated on reconstructed clusters to measure trigger efficiency.
	 */
	private void singlesTriggerVerification() {
		// ==========================================================
		// ==== Initialize Singles Trigger Verification =============
		// ==========================================================
		
		// Print the cluster verification header.
		println();
		println();
		println("======================================================================");
		println("=== Singles Trigger Verification =====================================");
		println("======================================================================");
		
		// Track the number of triggers seen and the number found.
		int sspReportedTriggers = 0;
		int sspInternalMatched = 0;
		int reconTriggersMatched = 0;
		
		// Track the number of times a given cut caused a trigger to
		// fail to match.
		int[] eventEnergyMin = new int[2];
		int[] eventEnergyMax = new int[2];
		int[] eventHitCount = new int[2];
		int[] eventTime = new int[2];
		
		
		
		// ==========================================================
		// ==== Output Event Summary ================================
		// ==========================================================
		
		// Get the list of triggers reported by the SSP.
		List<SSPTrigger> sspTriggers = sspBank.getTriggers();
		
		// Output the SSP cluster singles triggers.
		println();
		println("SSP Cluster Singles Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(SinglesTrigger<SSPCluster> simTrigger : sspSinglesTriggers.get(triggerNum)) {
				printf("\tTrigger %d :: %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d%n",
						(triggerNum + 1), sspClusterPositionString(simTrigger.getTriggerSource()),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0);
			}
		}
		if(sspSinglesTriggers.get(0).size() + sspSinglesTriggers.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the reconstructed cluster singles triggers.
		println("Reconstructed Cluster Singles Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(SinglesTrigger<Cluster> simTrigger : reconSinglesTriggers.get(triggerNum)) {
				printf("\tTrigger %d :: %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d%n",
						(triggerNum + 1), reconClusterPositionString(simTrigger.getTriggerSource()),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0);
			}
		}
		if(reconSinglesTriggers.get(0).size() + reconSinglesTriggers.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the SSP reported triggers.
		println("SSP Reported Singles Triggers");
		for(SSPTrigger sspTrigger : sspTriggers) {
			if(sspTrigger instanceof SSPSinglesTrigger) {
				// Cast the trigger to a singles trigger.
				SSPSinglesTrigger sspSingles = (SSPSinglesTrigger) sspTrigger;
				
				// Increment the number of SSP cluster singles triggers.
				sspReportedTriggers++;
				
				// Get the trigger properties.
				int triggerNum = sspSingles.isFirstTrigger() ? 1 : 2;
				
				// Print the trigger.
				printf("\tTrigger %d :: %3d ns :: EClusterLow: %d; EClusterHigh %d; HitCount: %d%n",
						triggerNum, sspSingles.getTime(), sspSingles.passCutEnergyMin() ? 1 : 0,
						sspSingles.passCutEnergyMax() ? 1 : 0, sspSingles.passCutHitCount() ? 1 : 0);
			}
		}
		if(sspReportedTriggers == 0) { println("\tNone"); }
		
		
		
		// ==========================================================
		// ==== SSP Internal Logic Verification =====================
		// ==========================================================
		
		// Track which SSP triggers have been matched to avoid matching
		// multiple reconstructed SSP cluster triggers to the same SSP
		// trigger.
		Set<SSPSinglesTrigger> sspTriggerSet = new HashSet<SSPSinglesTrigger>();
		Set<SinglesTrigger<SSPCluster>> simTriggerSet = new HashSet<SinglesTrigger<SSPCluster>>();
		
		// Iterate over the triggers.
		println();
		println("SSP Reported Trigger --> SSP Cluster Trigger Match Status");
		for(SSPTrigger sspTrigger : sspTriggers) {
			// If the trigger is a singles trigger, convert it.
			if(sspTrigger instanceof SSPSinglesTrigger) {
				// Cast the trigger to a singles trigger.
				SSPSinglesTrigger sspSingles = (SSPSinglesTrigger) sspTrigger;
				int triggerNum = sspSingles.isFirstTrigger() ? 0 : 1;
				boolean matchedTrigger = false;
				
				// Iterate over the SSP cluster simulated triggers and
				// look for a trigger that matches.
				matchLoop:
				for(SinglesTrigger<SSPCluster> simTrigger : sspSinglesTriggers.get(triggerNum)) {
					// If the current SSP trigger has already been
					// matched, skip it.
					if(sspTriggerSet.contains(sspSingles)) { continue matchLoop; }
					
					// Otherwise, check whether the reconstructed SSP
					// cluster trigger matches the SSP trigger.
					if(compareSSPSinglesTriggers(sspSingles, simTrigger)) {
						matchedTrigger = true;
						sspTriggerSet.add(sspSingles);
						simTriggerSet.add(simTrigger);
						sspInternalMatched++;
						break matchLoop;
					}
				}
				
				printf("\tTrigger %d :: %3d :: EClusterLow: %d; EClusterHigh %d; HitCount: %d :: Matched: %5b%n",
						(triggerNum + 1), sspSingles.getTime(), sspSingles.passCutEnergyMin() ? 1 : 0,
						sspSingles.passCutEnergyMax() ? 1 : 0, sspSingles.passCutHitCount() ? 1 : 0,
						matchedTrigger);
			}
		}
		
		// Iterate over the unmatched simulated triggers again and the
		// unmatched SSP reported trigger that most closely matches it.
		reportedTriggerLoop:
		for(SSPTrigger sspTrigger : sspTriggers) {
			// If the trigger is a singles trigger, convert it.
			if(sspTrigger instanceof SSPSinglesTrigger) {
				// Cast the trigger to a singles trigger.
				SSPSinglesTrigger sspSingles = (SSPSinglesTrigger) sspTrigger;
				
				// If this reported trigger has already been matched,
				// ignore it and continue to the next.
				if(sspTriggerSet.contains(sspSingles)) { continue reportedTriggerLoop; }
				
				// Otherwise, obtain information about the trigger.
				int triggerNum = sspSingles.isFirstTrigger() ? 0 : 1;
				int numMatched = -1;
				boolean foundBest = false;
				boolean[] matchedCut = new boolean[3];
				
				// Iterate over the simulated SSP triggers and find the
				// trigger which best matches this trigger, but does not
				// have a match already.
				matchLoop:
				for(SinglesTrigger<SSPCluster> simTrigger : sspSinglesTriggers.get(triggerNum)) {
					// If this trigger is at a different time, skip it.
					if(sspSingles.getTime() != simTrigger.getTriggerSource().getTime()) {
						continue matchLoop;
					}
					
					// If this trigger has been matched, skip it.
					if(simTriggerSet.contains(simTrigger)) { continue matchLoop; }
					
					// Check each of the cuts.
					boolean[] tempMatchedCut = new boolean[3];
					tempMatchedCut[0] = (simTrigger.getStateClusterEnergyLow()  == sspSingles.passCutEnergyMin());
					tempMatchedCut[1] = (simTrigger.getStateClusterEnergyHigh() == sspSingles.passCutEnergyMax());
					tempMatchedCut[2] = (simTrigger.getStateHitCount()          == sspSingles.passCutHitCount());
					
					// Check each cut and see if this is a closer match
					// than the previous best match.
					int tempNumMatched = 0;
					for(boolean passed : tempMatchedCut) { if(passed) { tempNumMatched++; } }
					
					// If the number of matched cuts exceeds the old
					// best result, this becomes the new best result.
					if(tempNumMatched > numMatched) {
						foundBest = true;
						numMatched = tempNumMatched;
						matchedCut = tempMatchedCut;
					}
				}
				
				// If some match was found, note what caused it to not
				// qualify as a complete match.
				if(foundBest) {
					if(!matchedCut[0]) { eventEnergyMin[triggerNum]++; }
					if(!matchedCut[1]) { eventEnergyMax[triggerNum]++; }
					if(!matchedCut[2]) { eventHitCount[triggerNum]++; }
				}
				
				// If there was no match found, it means that there were
				// no triggers that were both unmatched and at the same
				// time as this simulated trigger.
				else { eventTime[triggerNum]++; }
			}
		}
		
		
		
		// ==========================================================
		// ==== SSP Singles Trigger Efficiency ======================
		// ==========================================================
		
		// Reset the SSP matched trigger set.
		sspTriggerSet.clear();
		
		// Iterate over the reconstructed cluster singles triggers.
		println();
		println("Recon Cluster Trigger --> SSP Reported Trigger Match Status");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(SinglesTrigger<Cluster> simTrigger : reconSinglesTriggers.get(triggerNum)) {
				
				printf("\tTrigger %d :: %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d%n",
						(triggerNum + 1), reconClusterPositionString(simTrigger.getTriggerSource()),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0);
				printf("\t\tCluster Energy :: %.3f GeV; Hit Count :: %d%n", simTrigger.getTriggerSource().getEnergy(),
						simTrigger.getTriggerSource().getCalorimeterHits().size());
				printf("\t\tCluster Energy Cut :: [ %.3f GeV, %.3f GeV ]; Hit Count :: [ %.0f, INF )%n",
						singlesTrigger[triggerNum].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW),
						singlesTrigger[triggerNum].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH),
						singlesTrigger[triggerNum].getCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW));
				
				// Iterate over the SSP reported triggers and compare
				// them to the reconstructed cluster simulated trigger.
				matchLoop:
				for(SSPTrigger sspTrigger : sspTriggers) {
					// Only compare singles triggers.
					if(sspTrigger instanceof SSPSinglesTrigger) {
						// Cast the trigger.
						SSPSinglesTrigger sspSingles = (SSPSinglesTrigger) sspTrigger;
						
						printf("\t\t\tTrigger %d :: %3d ns :: EClusterLow: %d; EClusterHigh %d; HitCount: %d",
								sspSingles.isFirstTrigger() ? 1 : 2,
								sspSingles.getTime(),
								sspSingles.passCutEnergyMin() ? 1 : 0,
								sspSingles.passCutEnergyMax() ? 1 : 0,
								sspSingles.passCutHitCount() ? 1 : 0);
						
						// Only compare triggers if they are from the
						// same trigger source.
						if((triggerNum == 0 && sspSingles.isSecondTrigger())
								|| (triggerNum == 1 && sspSingles.isFirstTrigger())) {
							print(" [ fail; source    ]%n");
							continue matchLoop;
						}
						
						// Only compare the singles trigger if it was
						// not already matched to another trigger.
						if(sspTriggerSet.contains(sspSingles)) {
							print(" [ fail; matched   ]%n");
							continue matchLoop;
						}
						
						// Test each cut.
						if(sspSingles.passCutEnergyMin() != simTrigger.getStateClusterEnergyLow()) {
							print(" [ fail; E_min     ]%n");
							continue matchLoop;
						} if(sspSingles.passCutEnergyMax() != simTrigger.getStateClusterEnergyHigh()) {
							print(" [ fail; E_max     ]%n");
							continue matchLoop;
						} if(sspSingles.passCutHitCount() != simTrigger.getStateHitCount()) {
							print(" [ fail; hit count ]%n");
							continue matchLoop;
						}
						
						// If all the trigger flags match, then the
						// triggers are a match.
						reconTriggersMatched++;
						sspTriggerSet.add(sspSingles);
						print(" [ success         ]%n");
						break matchLoop;
					}
				}
			}
		}
		
		
		
		// ==========================================================
		// ==== Output Event Results ================================
		// ==========================================================
		
		// Get the number of SSP and reconstructed cluster simulated
		// triggers.
		int sspSimTriggers = sspSinglesTriggers.get(0).size() + sspSinglesTriggers.get(1).size();
		int reconSimTriggers = reconSinglesTriggers.get(0).size() + reconSinglesTriggers.get(1).size();
		
		// Print event statistics.
		println();
		println("Event Statistics:");
		printf("\tSSP Cluster Sim Triggers   :: %d%n", sspSimTriggers);
		printf("\tRecon Cluster Sim Triggers :: %d%n", reconSimTriggers);
		printf("\tSSP Reported Triggers      :: %d%n", sspReportedTriggers);
		if(sspSimTriggers == 0) {
			printf("\tInternal Efficiency        :: %d / %d (N/A)%n",
					sspInternalMatched, sspSimTriggers);
		} else {
			printf("\tInternal Efficiency        :: %d / %d (%3.0f%%)%n",
					sspInternalMatched, sspSimTriggers, (100.0 * sspInternalMatched / sspSimTriggers));
		}
		if(reconSimTriggers == 0) {
			printf("\tTrigger Efficiency         :: %d / %d (N/A)%n",
					reconTriggersMatched, reconSimTriggers);
		} else {
			printf("\tTrigger Efficiency         :: %d / %d (%3.0f%%)%n",
					reconTriggersMatched, reconSimTriggers, (100.0 * reconTriggersMatched / reconSimTriggers));
		}
		
		// Print the individual cut performances.
		println();
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
		printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(sspSimTriggers == 0) {
				printf("\tCluster Energy Lower Bound :: %d / %d%n", eventEnergyMin[triggerNum], sspSimTriggers);
				printf("\tCluster Energy Upper Bound :: %d / %d%n", eventEnergyMax[triggerNum], sspSimTriggers);
				printf("\tCluster Hit Count          :: %d / %d%n", eventHitCount[triggerNum], sspSimTriggers);
			} else {
				printf("\tCluster Energy Lower Bound :: %d / %d (%3.0f%%)%n",
						eventEnergyMin[triggerNum], sspSimTriggers, (100.0 * eventEnergyMin[triggerNum] / sspSimTriggers));
				printf("\tCluster Energy Upper Bound :: %d / %d (%3.0f%%)%n",
						eventEnergyMax[triggerNum], sspSimTriggers, (100.0 * eventEnergyMax[triggerNum] / sspSimTriggers));
				printf("\tCluster Hit Count          :: %d / %d (%3.0f%%)%n",
						eventHitCount[triggerNum], sspSimTriggers, (100.0 * eventHitCount[triggerNum] / sspSimTriggers));
			}
		}
		
		// Update the global trigger tracking variables.
		singlesSSPTriggers += sspSimTriggers;
		singlesReconMatched += reconTriggersMatched;
		singlesReconTriggers += reconSimTriggers;
		singlesInternalMatched += sspInternalMatched;
		singlesReportedTriggers += sspReportedTriggers;
		
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			globalEnergyMinCut[triggerNum] += eventEnergyMin[triggerNum];
			globalEnergyMaxCut[triggerNum] += eventEnergyMax[triggerNum];
			globalHitCountCut[triggerNum] += eventHitCount[triggerNum];
			globalSinglesTimeCut[triggerNum] += eventTime[triggerNum];
		}
		
		// Note whether the was a singles trigger match failure.
		if((reconTriggersMatched - reconSimTriggers != 0) || (sspInternalMatched - sspSimTriggers != 0)) {
			singlesFail = true;
		}
	}
	
	/**
	 * Checks triggers simulated on SSP clusters against the SSP bank's
	 * reported triggers to verify that the trigger is correctly applying
	 * cuts to the clusters it sees. Additionally compares triggers
	 * simulated on reconstructed clusters to measure trigger efficiency.
	 */
	private void pairTriggerVerification() {
		// ==========================================================
		// ==== Initialize Pair Trigger Verification ===============
		// ==========================================================
		
		// Print the cluster verification header.
		println();
		println();
		println("======================================================================");
		println("=== Pair Trigger Verification ========================================");
		println("======================================================================");
		
		// Track the number of triggers seen and the number found.
		int sspReportedTriggers = 0;
		int sspInternalMatched = 0;
		int reconTriggersMatched = 0;
		
		int[] eventEnergySum = new int[2];
		int[] eventEnergyDiff = new int[2];
		int[] eventEnergySlope = new int[2];
		int[] eventCoplanarity = new int[2];
		int[] eventTime = new int[2];
		
		
		
		// ==========================================================
		// ==== Output Event Summary ================================
		// ==========================================================
		
		// Get the list of triggers reported by the SSP.
		List<SSPTrigger> sspTriggers = sspBank.getTriggers();
		
		// Output the SSP cluster pair triggers.
		println();
		println("SSP Cluster Pair Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(PairTrigger<SSPCluster[]> simTrigger : sspPairsTriggers.get(triggerNum)) {
				printf("\tTrigger %d :: %s, %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d; ESumLow: %d, ESumHigh: %d, EDiff: %d, ESlope: %d, Coplanarity: %d%n",
						(triggerNum + 1), sspClusterPositionString(simTrigger.getTriggerSource()[0]),
						sspClusterPositionString(simTrigger.getTriggerSource()[1]),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0,
						simTrigger.getStateEnergySumLow() ? 1 : 0,
						simTrigger.getStateEnergySumHigh() ? 1 : 0,
						simTrigger.getStateEnergyDifference() ? 1 : 0,
						simTrigger.getStateEnergySlope() ? 1 : 0,
						simTrigger.getStateCoplanarity() ? 1 : 0);
			}
		}
		if(sspPairsTriggers.get(0).size() + sspPairsTriggers.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the reconstructed cluster singles triggers.
		println("Reconstructed Cluster Pair Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(PairTrigger<Cluster[]> simTrigger : reconPairsTriggers.get(triggerNum)) {
				printf("\tTrigger %d :: %s, %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d; ESumLow: %d, ESumHigh: %d, EDiff: %d, ESlope: %d, Coplanarity: %d%n",
						(triggerNum + 1), reconClusterPositionString(simTrigger.getTriggerSource()[0]),
						reconClusterPositionString(simTrigger.getTriggerSource()[1]),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0,
						simTrigger.getStateEnergySumLow() ? 1 : 0,
						simTrigger.getStateEnergySumHigh() ? 1 : 0,
						simTrigger.getStateEnergyDifference() ? 1 : 0,
						simTrigger.getStateEnergySlope() ? 1 : 0,
						simTrigger.getStateCoplanarity() ? 1 : 0);
			}
		}
		if(reconPairsTriggers.get(0).size() + reconPairsTriggers.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the SSP reported triggers.
		println("SSP Reported Pair Triggers");
		for(SSPTrigger sspTrigger : sspTriggers) {
			if(sspTrigger instanceof SSPPairTrigger) {
				// Cast the trigger to a singles trigger.
				SSPPairTrigger sspPair = (SSPPairTrigger) sspTrigger;
				
				// Increment the number of SSP cluster singles triggers.
				sspReportedTriggers++;
				
				// Get the trigger properties.
				int triggerNum = sspPair.isFirstTrigger() ? 1 : 2;
				
				// Print the trigger.
				printf("\tTrigger %d :: %3d ns :: ESum: %d, EDiff: %d, ESlope: %d, Coplanarity: %d%n",
						triggerNum, sspPair.getTime(),
						sspPair.passCutEnergySum() ? 1 : 0, sspPair.passCutEnergyDifference() ? 1 : 0,
						sspPair.passCutEnergySlope() ? 1 : 0, sspPair.passCutCoplanarity() ? 1 : 0);
			}
		}
		if(sspReportedTriggers == 0) { println("\tNone"); }
		
		
		
		// ==========================================================
		// ==== SSP Internal Logic Verification =====================
		// ==========================================================
		
		// Track which SSP triggers have been matched to avoid matching
		// multiple reconstructed SSP cluster triggers to the same SSP
		// trigger.
		Set<SSPPairTrigger> sspTriggerSet = new HashSet<SSPPairTrigger>();
		Set<PairTrigger<SSPCluster[]>> simTriggerSet = new HashSet<PairTrigger<SSPCluster[]>>();
		
		// Iterate over the triggers.
		println();
		println("SSP Reported Trigger --> SSP Cluster Trigger Match Status");
		for(SSPTrigger sspTrigger : sspTriggers) {
			// If the trigger is a pair trigger, convert it.
			if(sspTrigger instanceof SSPPairTrigger) {
				// Cast the trigger to a pair trigger.
				SSPPairTrigger sspPair = (SSPPairTrigger) sspTrigger;
				int triggerNum = sspPair.isFirstTrigger() ? 0 : 1;
				boolean matchedTrigger = false;
				
				// Iterate over the SSP cluster simulated triggers and
				// look for a trigger that matches.
				matchLoop:
				for(PairTrigger<SSPCluster[]> simTrigger : sspPairsTriggers.get(triggerNum)) {
					// If the current SSP trigger has already been
					// matched, skip it.
					if(sspTriggerSet.contains(sspPair)) { continue matchLoop; }
					
					// Otherwise, check whether the reconstructed SSP
					// cluster trigger matches the SSP trigger.
					if(compareSSPPairTriggers(sspPair, simTrigger)) {
						matchedTrigger = true;
						sspTriggerSet.add(sspPair);
						simTriggerSet.add(simTrigger);
						sspInternalMatched++;
						break matchLoop;
					}
				}
				
				printf("\tTrigger %d :: %3d ns :: ESum: %d, EDiff: %d, ESlope: %d, Coplanarity: %d :: Matched: %5b%n",
						triggerNum, sspPair.getTime(), sspPair.passCutEnergySum() ? 1 : 0,
						sspPair.passCutEnergyDifference() ? 1 : 0, sspPair.passCutEnergySlope() ? 1 : 0,
						sspPair.passCutCoplanarity() ? 1 : 0, matchedTrigger);
			}
		}
		
		// Iterate over the unmatched simulated triggers again and the
		// unmatched SSP reported trigger that most closely matches it.
		reportedTriggerLoop:
		for(SSPTrigger sspTrigger : sspTriggers) {
			// If the trigger is a singles trigger, convert it.
			if(sspTrigger instanceof SSPPairTrigger) {
				// Cast the trigger to a singles trigger.
				SSPPairTrigger sspPair = (SSPPairTrigger) sspTrigger;
				
				// If this reported trigger has already been matched,
				// ignore it and continue to the next.
				if(sspTriggerSet.contains(sspPair)) { continue reportedTriggerLoop; }
				
				// Otherwise, obtain information about the trigger.
				int triggerNum = sspPair.isFirstTrigger() ? 0 : 1;
				int numMatched = -1;
				boolean foundBest = false;
				boolean[] matchedCut = new boolean[3];
				
				// Iterate over the simulated SSP triggers and find the
				// trigger which best matches this trigger, but does not
				// have a match already.
				matchLoop:
				for(PairTrigger<SSPCluster[]> simTrigger : sspPairsTriggers.get(triggerNum)) {
					// Get the time of the simulated pair.
					int simTime = 0;
					if(simTrigger.getTriggerSource()[0].getYIndex() < 0) {
						simTime = simTrigger.getTriggerSource()[0].getTime();
					} else {
						simTime = simTrigger.getTriggerSource()[1].getTime();
					}
					
					// If this trigger is at a different time, skip it.
					if(sspPair.getTime() != simTime) {
						continue matchLoop;
					}
					
					// If this trigger has been matched, skip it.
					if(simTriggerSet.contains(simTrigger)) { continue matchLoop; }
					
					// Check each of the cuts.
					boolean[] tempMatchedCut = new boolean[3];
					tempMatchedCut[0] = (simTrigger.getStateEnergySum()        == sspPair.passCutEnergySum());
					tempMatchedCut[1] = (simTrigger.getStateEnergyDifference() == sspPair.passCutEnergyDifference());
					tempMatchedCut[2] = (simTrigger.getStateEnergySlope()      == sspPair.passCutEnergySlope());
					tempMatchedCut[3] = (simTrigger.getStateCoplanarity()      == sspPair.passCutCoplanarity());
					
					// Check each cut and see if this is a closer match
					// than the previous best match.
					int tempNumMatched = 0;
					for(boolean passed : tempMatchedCut) { if(passed) { tempNumMatched++; } }
					
					// If the number of matched cuts exceeds the old
					// best result, this becomes the new best result.
					if(tempNumMatched > numMatched) {
						foundBest = true;
						numMatched = tempNumMatched;
						matchedCut = tempMatchedCut;
					}
				}
				
				// If some match was found, note what caused it to not
				// qualify as a complete match.
				if(foundBest) {
					if(!matchedCut[0]) { eventEnergySum[triggerNum]++; }
					if(!matchedCut[1]) { eventEnergyDiff[triggerNum]++; }
					if(!matchedCut[2]) { eventEnergySlope[triggerNum]++; }
					if(!matchedCut[3]) { eventCoplanarity[triggerNum]++; }
				}
				
				// If there was no match found, it means that there were
				// no triggers that were both unmatched and at the same
				// time as this simulated trigger.
				else { eventTime[triggerNum]++; pairFail = true; }
			}
		}
		
		
		
		// ==========================================================
		// ==== SSP Pair Trigger Efficiency =========================
		// ==========================================================
		
		// Reset the SSP matched trigger set.
		sspTriggerSet.clear();
		
		// Iterate over the reconstructed cluster pair triggers.
		println();
		println("Recon Cluster Trigger --> SSP Reported Trigger Match Status");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(PairTrigger<Cluster[]> simTrigger : reconPairsTriggers.get(triggerNum)) {
				// Track whether the trigger was matched.
				boolean matchedTrigger = false;
				
				// Iterate over the SSP reported triggers and compare
				// them to the reconstructed cluster simulated trigger.
				matchLoop:
				for(SSPTrigger sspTrigger : sspTriggers) {
					// Only compare pair triggers.
					if(sspTrigger instanceof SSPPairTrigger) {
						// Cast the trigger.
						SSPPairTrigger sspPair = (SSPPairTrigger) sspTrigger;
						
						// Only compare the pair trigger if it was
						// not already matched to another trigger.
						if(sspTriggerSet.contains(sspPair)) { continue matchLoop; }
						
						// Compare the triggers.
						if(compareReconPairTriggers(sspPair, simTrigger)) {
							reconTriggersMatched++;
							matchedTrigger = true;
							sspTriggerSet.add(sspPair);
							break matchLoop;
						}
					}
				}
				
				// Print the trigger matching status.
				printf("\tTrigger %d :: %s, %s :: EClusterLow: %d; EClusterHigh %d; HitCount: %d; ESumLow: %d, ESumHigh: %d, EDiff: %d, ESlope: %d, Coplanarity: %d :: Matched: %5b%n",
						(triggerNum + 1), reconClusterPositionString(simTrigger.getTriggerSource()[0]),
						reconClusterPositionString(simTrigger.getTriggerSource()[1]),
						simTrigger.getStateClusterEnergyLow() ? 1 : 0,
						simTrigger.getStateClusterEnergyHigh() ? 1 : 0,
						simTrigger.getStateHitCount() ? 1 : 0,
						simTrigger.getStateEnergySumLow() ? 1 : 0,
						simTrigger.getStateEnergySumHigh() ? 1 : 0,
						simTrigger.getStateEnergyDifference() ? 1 : 0,
						simTrigger.getStateEnergySlope() ? 1 : 0,
						simTrigger.getStateCoplanarity() ? 1 : 0, matchedTrigger);
			}
		}
		
		
		
		// ==========================================================
		// ==== Output Event Results ================================
		// ==========================================================
		
		// Get the number of SSP and reconstructed cluster simulated
		// triggers.
		int sspSimTriggers = sspPairsTriggers.get(0).size() + sspPairsTriggers.get(1).size();
		int reconSimTriggers = reconPairsTriggers.get(0).size() + reconPairsTriggers.get(1).size();
		
		// Print event statistics.
		println();
		println("Event Statistics:");
		printf("\tSSP Cluster Sim Triggers   :: %d%n", sspSimTriggers);
		printf("\tRecon Cluster Sim Triggers :: %d%n", reconSimTriggers);
		printf("\tSSP Reported Triggers      :: %d%n", sspReportedTriggers);
		if(sspSimTriggers == 0) {
			printf("\tInternal Efficiency        :: %d / %d (N/A)%n",
					sspInternalMatched, sspSimTriggers);
		} else {
			printf("\tInternal Efficiency        :: %d / %d (%3.0f%%)%n",
					sspInternalMatched, sspSimTriggers, (100.0 * sspInternalMatched / sspSimTriggers));
		}
		if(reconSimTriggers == 0) {
			printf("\tTrigger Efficiency         :: %d / %d (N/A)%n",
					reconTriggersMatched, reconSimTriggers);
		} else {
			printf("\tTrigger Efficiency         :: %d / %d (%3.0f%%)%n",
					reconTriggersMatched, reconSimTriggers, (100.0 * reconTriggersMatched / reconSimTriggers));
		}
		
		// Print the individual cut performances.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			println();
			printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(sspSimTriggers == 0) {
				printf("\tPair Energy Sum            :: %d / %d%n", eventEnergySum[triggerNum], sspSimTriggers);
				printf("\tPair Energy Difference     :: %d / %d%n", eventEnergyDiff[triggerNum], sspSimTriggers);
				printf("\tPair Energy Slope          :: %d / %d%n", eventEnergySlope[triggerNum], sspSimTriggers);
				printf("\tPair Coplanarity           :: %d / %d%n", eventCoplanarity[triggerNum], sspSimTriggers);
				printf("\tPair Trigger Time          :: %d / %d%n", eventTime[triggerNum], sspSimTriggers);
			} else {
				printf("\tPair Energy Sum            :: %d / %d (%3.0f%%)%n",
						eventEnergySum[triggerNum], sspSimTriggers, (100.0 * eventEnergySum[triggerNum] / sspSimTriggers));
				printf("\tPair Energy Difference     :: %d / %d (%3.0f%%)%n",
						eventEnergyDiff[triggerNum], sspSimTriggers, (100.0 * eventEnergyDiff[triggerNum] / sspSimTriggers));
				printf("\tPair Energy Slope          :: %d / %d (%3.0f%%)%n",
						eventEnergySlope[triggerNum], sspSimTriggers, (100.0 * eventEnergySlope[triggerNum] / sspSimTriggers));
				printf("\tPair Coplanarity           :: %d / %d (%3.0f%%)%n",
						eventCoplanarity[triggerNum], sspSimTriggers, (100.0 * eventCoplanarity[triggerNum] / sspSimTriggers));
				printf("\tPair Trigger Time          :: %d / %d (%3.0f%%)%n",
						eventTime[triggerNum], sspSimTriggers, (100.0 * eventTime[triggerNum] / sspSimTriggers));
			}
		}
		
		// Update the global trigger tracking variables.
		pairSSPTriggers += sspSimTriggers;
		pairReconMatched += reconTriggersMatched;
		pairReconTriggers += reconSimTriggers;
		pairInternalMatched += sspInternalMatched;
		pairReportedTriggers += sspReportedTriggers;
		
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			globalEnergySumCut[triggerNum] += eventEnergySum[triggerNum];
			globalEnergyDiffCut[triggerNum] += eventEnergyDiff[triggerNum];
			globalEnergySlopeCut[triggerNum] += eventEnergySlope[triggerNum];
			globalCoplanarityCut[triggerNum] += eventCoplanarity[triggerNum];
			globalPairTimeCut[triggerNum] += eventTime[triggerNum];
		}
		
		// Note whether the was a singles trigger match failure.
		if((reconTriggersMatched - reconSimTriggers != 0) || (sspInternalMatched - sspSimTriggers != 0)) {
			pairFail = true;
		}
	}
	
	/**
	 * Generates and stores the singles triggers for both reconstructed
	 * and SSP clusters.
	 */
	private void constructSinglesTriggers() {
		// Run the SSP clusters through the singles trigger to determine
		// whether they pass it or not.
		for(SSPCluster cluster : sspClusters) {
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
				// For a cluster to have formed it is assumed to have passed
				// the cluster seed energy cuts. This can not be verified
				// since the SSP bank does not report individual hit. 
				boolean passSeedLow = true;
				boolean passSeedHigh = true;
				
				// The remaining cuts may be acquired from trigger module.
				boolean passClusterLow = singlesTrigger[triggerNum].clusterTotalEnergyCutLow(cluster);
				boolean passClusterHigh = singlesTrigger[triggerNum].clusterTotalEnergyCutHigh(cluster);
				boolean passHitCount = singlesTrigger[triggerNum].clusterHitCountCut(cluster);
				
				// Make a trigger to store the results.
				SinglesTrigger<SSPCluster> trigger = new SinglesTrigger<SSPCluster>(cluster);
				trigger.setStateSeedEnergyLow(passSeedLow);
				trigger.setStateSeedEnergyHigh(passSeedHigh);
				trigger.setStateClusterEnergyLow(passClusterLow);
				trigger.setStateClusterEnergyHigh(passClusterHigh);
				trigger.setStateHitCount(passHitCount);
				
				// Store the trigger.
				sspSinglesTriggers.get(triggerNum).add(trigger);
			}
		}
		
		// Run the reconstructed clusters through the singles trigger
		// to determine whether they pass it or not.
		for(Cluster cluster : reconClusters) {
			// Simulate each of the cluster singles triggers.
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
				// For a cluster to have formed it is assumed to have passed
				// the cluster seed energy cuts. This can not be verified
				// since the SSP bank does not report individual hit. 
				boolean passSeedLow = true;
				boolean passSeedHigh = true;
				
				// The remaining cuts may be acquired from trigger module.
				boolean passClusterLow = singlesTrigger[triggerNum].clusterTotalEnergyCutLow(cluster);
				boolean passClusterHigh = singlesTrigger[triggerNum].clusterTotalEnergyCutHigh(cluster);
				boolean passHitCount = singlesTrigger[triggerNum].clusterHitCountCut(cluster);
				
				// Make a trigger to store the results.
				SinglesTrigger<Cluster> trigger = new SinglesTrigger<Cluster>(cluster);
				trigger.setStateSeedEnergyLow(passSeedLow);
				trigger.setStateSeedEnergyHigh(passSeedHigh);
				trigger.setStateClusterEnergyLow(passClusterLow);
				trigger.setStateClusterEnergyHigh(passClusterHigh);
				trigger.setStateHitCount(passHitCount);
				
				// Store the trigger.
				reconSinglesTriggers.get(triggerNum).add(trigger);
			}
		}
	}
	
	/**
	 * Generates and stores the pair triggers for both reconstructed
	 * and SSP clusters.
	 */
	private void constructPairTriggers() {
		// Store cluster pairs.
		List<Cluster> topReconClusters = new ArrayList<Cluster>();
		List<Cluster> bottomReconClusters = new ArrayList<Cluster>();
		List<Cluster[]> reconPairs = new ArrayList<Cluster[]>();
		List<SSPCluster> topSSPClusters = new ArrayList<SSPCluster>();
		List<SSPCluster> bottomSSPClusters = new ArrayList<SSPCluster>();
		List<SSPCluster[]> sspPairs = new ArrayList<SSPCluster[]>();
		
		// Split the clusters into lists of top and bottom clusters.
		for(Cluster reconCluster : reconClusters) {
			if(reconCluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy") > 0) {
				topReconClusters.add(reconCluster);
			} else {
				bottomReconClusters.add(reconCluster);
			}
		}
		for(SSPCluster sspCluster : sspClusters) {
			if(sspCluster.getYIndex() > 0) {
				topSSPClusters.add(sspCluster);
			} else {
				bottomSSPClusters.add(sspCluster);
			}
		}
		
		// Form all possible top/bottom cluster pairs.
		for(Cluster topReconCluster : topReconClusters) {
			for(Cluster bottomReconCluster : bottomReconClusters) {
				Cluster[] reconPair = new Cluster[2];
				reconPair[0] = topReconCluster;
				reconPair[1] = bottomReconCluster;
				reconPairs.add(reconPair);
			}
		}
		for(SSPCluster topSSPCluster : topSSPClusters) {
			for(SSPCluster bottomSSPCluster : bottomSSPClusters) {
				SSPCluster[] sspPair = new SSPCluster[2];
				sspPair[0] = topSSPCluster;
				sspPair[1] = bottomSSPCluster;
				sspPairs.add(sspPair);
			}
		}
		
		// Simulate the pair triggers and record the results.
		for(Cluster[] reconPair : reconPairs) {
			// Simulate each of the cluster pair triggers.
			reconTriggerLoop:
			for(int triggerIndex = 0; triggerIndex < 2; triggerIndex++) {
				// Check that the pair passes the time coincidence cut.
				// If it does not, it is not a valid pair and should be
				// destroyed.
				if(!pairsTrigger[triggerIndex].pairTimeCoincidenceCut(reconPair)) {
					continue reconTriggerLoop;
				}
				
				// For a cluster to have formed it is assumed to have passed
				// the cluster seed energy cuts. This can not be verified
				// since the SSP bank does not report individual hit. 
				boolean passSeedLow = true;
				boolean passSeedHigh = true;
				
				// The remaining cuts may be acquired from trigger module.
				boolean passClusterLow = pairsTrigger[triggerIndex].clusterTotalEnergyCutLow(reconPair[0])
						&& pairsTrigger[triggerIndex].clusterTotalEnergyCutLow(reconPair[1]);
				boolean passClusterHigh = pairsTrigger[triggerIndex].clusterTotalEnergyCutHigh(reconPair[0])
						&& pairsTrigger[triggerIndex].clusterTotalEnergyCutHigh(reconPair[1]);
				boolean passHitCount = pairsTrigger[triggerIndex].clusterHitCountCut(reconPair[0])
						&& pairsTrigger[triggerIndex].clusterHitCountCut(reconPair[1]);
				boolean passPairEnergySumLow = pairsTrigger[triggerIndex].pairEnergySumCutLow(reconPair);
				boolean passPairEnergySumHigh = pairsTrigger[triggerIndex].pairEnergySumCutHigh(reconPair);
				boolean passPairEnergyDifference = pairsTrigger[triggerIndex].pairEnergyDifferenceCut(reconPair);
				boolean passPairEnergySlope = pairsTrigger[triggerIndex].pairEnergySlopeCut(reconPair);
				boolean passPairCoplanarity = pairsTrigger[triggerIndex].pairCoplanarityCut(reconPair);
				boolean passTimeCoincidence = pairsTrigger[triggerIndex].pairTimeCoincidenceCut(reconPair);
				
				// Create a trigger from the results.
				PairTrigger<Cluster[]> trigger = new PairTrigger<Cluster[]>(reconPair);
				trigger.setStateSeedEnergyLow(passSeedLow);
				trigger.setStateSeedEnergyHigh(passSeedHigh);
				trigger.setStateClusterEnergyLow(passClusterLow);
				trigger.setStateClusterEnergyHigh(passClusterHigh);
				trigger.setStateHitCount(passHitCount);
				trigger.setStateEnergySumLow(passPairEnergySumLow);
				trigger.setStateEnergySumHigh(passPairEnergySumHigh);
				trigger.setStateEnergyDifference(passPairEnergyDifference);
				trigger.setStateEnergySlope(passPairEnergySlope);
				trigger.setStateCoplanarity(passPairCoplanarity);
				trigger.setStateTimeCoincidence(passTimeCoincidence);
				
				// Add the trigger to the list.
				reconPairsTriggers.get(triggerIndex).add(trigger);
			}
		}
		
		for(SSPCluster[] sspPair : sspPairs) {
			pairTriggerLoop:
			for(int triggerIndex = 0; triggerIndex < 2; triggerIndex++) {
				// Check that the pair passes the time coincidence cut.
				// If it does not, it is not a valid pair and should be
				// destroyed.
				if(!pairsTrigger[triggerIndex].pairTimeCoincidenceCut(sspPair)) {
					continue pairTriggerLoop;
				}
				
				// For a cluster to have formed it is assumed to have passed
				// the cluster seed energy cuts. This can not be verified
				// since the SSP bank does not report individual hit. 
				boolean passSeedLow = true;
				boolean passSeedHigh = true;
				
				// The remaining cuts may be acquired from trigger module.
				boolean passClusterLow = pairsTrigger[triggerIndex].clusterTotalEnergyCutLow(sspPair[0])
						&& pairsTrigger[triggerIndex].clusterTotalEnergyCutLow(sspPair[1]);
				boolean passClusterHigh = pairsTrigger[triggerIndex].clusterTotalEnergyCutHigh(sspPair[0])
						&& pairsTrigger[triggerIndex].clusterTotalEnergyCutHigh(sspPair[1]);
				boolean passHitCount = pairsTrigger[triggerIndex].clusterHitCountCut(sspPair[0])
						&& pairsTrigger[triggerIndex].clusterHitCountCut(sspPair[1]);
				boolean passPairEnergySumLow = pairsTrigger[triggerIndex].pairEnergySumCutLow(sspPair);
				boolean passPairEnergySumHigh = pairsTrigger[triggerIndex].pairEnergySumCutHigh(sspPair);
				boolean passPairEnergyDifference = pairsTrigger[triggerIndex].pairEnergyDifferenceCut(sspPair);
				boolean passPairEnergySlope = pairsTrigger[triggerIndex].pairEnergySlopeCut(sspPair);
				boolean passPairCoplanarity = pairsTrigger[triggerIndex].pairCoplanarityCut(sspPair);
				boolean passTimeCoincidence = pairsTrigger[triggerIndex].pairTimeCoincidenceCut(sspPair);
				
				// Create a trigger from the results.
				PairTrigger<SSPCluster[]> trigger = new PairTrigger<SSPCluster[]>(sspPair);
				trigger.setStateSeedEnergyLow(passSeedLow);
				trigger.setStateSeedEnergyHigh(passSeedHigh);
				trigger.setStateClusterEnergyLow(passClusterLow);
				trigger.setStateClusterEnergyHigh(passClusterHigh);
				trigger.setStateHitCount(passHitCount);
				trigger.setStateEnergySumLow(passPairEnergySumLow);
				trigger.setStateEnergySumHigh(passPairEnergySumHigh);
				trigger.setStateEnergyDifference(passPairEnergyDifference);
				trigger.setStateEnergySlope(passPairEnergySlope);
				trigger.setStateCoplanarity(passPairCoplanarity);
				trigger.setStateTimeCoincidence(passTimeCoincidence);
				
				// Add the trigger to the list.
				sspPairsTriggers.get(triggerIndex).add(trigger);
			}
		}
	}
	
	/**
	 * Outputs all of the verification parameters currently in use by
	 * the software. A warning will be issued if the values for NSA and
	 * NSB, along with the FADC window, preclude clusters from being
	 * verified.
	 */
	private void logSettings() {
		// Output general settings.
		System.out.println("Cluster Verification Settings");
		System.out.printf("\tEnergy Threshold       :: %1.2f%n", energyAcceptance);
		System.out.printf("\tHit Threshold          :: %1d%n", hitAcceptance);
		
		// Output window settings.
		System.out.println("FADC Timing Window Settings");
		System.out.printf("\tNSB                    :: %3d ns%n", nsb);
		System.out.printf("\tNSA                    :: %3d ns%n", nsa);
		System.out.printf("\tFADC Window            :: %3d ns%n", windowWidth);
		
		// Calculate the valid clustering window.
		int start = nsb;
		int end = windowWidth - nsa;
		if(start < end) {
			System.out.printf("\tValid Cluster Window   :: [ %3d ns, %3d ns ]%n", start, end);
			performClusterVerification = true;
		} else {
			System.out.println("\tNSB, NSA, and FADC window preclude a valid cluster verification window.");
			System.out.println("\tCluster verification will not be performed!");
			performClusterVerification = false;
		}
		
		// Output the singles trigger settings.
		for(int i = 0; i < 2; i++) {
			System.out.printf("Singles Trigger %d Settings%n", (i + 1));
			System.out.printf("\tCluster Energy Low     :: %.3f GeV%n", singlesTrigger[i].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW));
			System.out.printf("\tCluster Energy High    :: %.3f GeV%n", singlesTrigger[i].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH));
			System.out.printf("\tCluster Hit Count      :: %.0f hits%n", singlesTrigger[i].getCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW));
		}
		
		// Output the pair trigger settings.
		for(int i = 0; i < 2; i++) {
			System.out.printf("Pairs Trigger %d Settings%n", (i + 1));
			System.out.printf("\tCluster Energy Low     :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW));
			System.out.printf("\tCluster Energy High    :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH));
			System.out.printf("\tCluster Hit Count      :: %.0f hits%n", pairsTrigger[i].getCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW));
			System.out.printf("\tPair Energy Sum Low    :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_ENERGY_SUM_LOW));
			System.out.printf("\tPair Energy Sum Low    :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_ENERGY_SUM_HIGH));
			System.out.printf("\tPair Energy Difference :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_ENERGY_DIFFERENCE_HIGH));
			System.out.printf("\tPair Energy Slope      :: %.3f GeV%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_ENERGY_SLOPE_LOW));
			System.out.printf("\tPair Energy Slope F    :: %.3f GeV / mm%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_ENERGY_SLOPE_F));
			System.out.printf("\tPair Coplanarity       :: %.0f Degrees%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_COPLANARITY_HIGH));
			System.out.printf("\tPair Time Coincidence  :: %.0f ns%n", pairsTrigger[i].getCutValue(TriggerModule.PAIR_TIME_COINCIDENCE));
		}
	}
	
	/**
	 * Checks whether all of the hits in a cluster are within the safe
	 * region of the FADC output window.
	 * @param reconCluster - The cluster to check.
	 * @return Returns <code>true</code> if the cluster is safe and
	 * returns <code>false</code> otherwise.
	 */
	private final boolean isVerifiable(Cluster reconCluster) {
		// Iterate over the hits in the cluster.
		for(CalorimeterHit hit : reconCluster.getCalorimeterHits()) {
			// Check that none of the hits are within the disallowed
			// region of the FADC readout window.
			if(hit.getTime() <= nsb || hit.getTime() >= (windowWidth - nsa)) {
				return false;
			}
			
			// Also check to make sure that the cluster does not have
			// any negative energy hits. These are, obviously, wrong.
			if(hit.getCorrectedEnergy() < 0.0) {
				return false;
			}
		}
		
		// If all of the cluster hits pass the time cut, the cluster
		// is valid.
		return true;
	}
	
	/**
	 * Convenience method that writes the information in a cluster to
	 * a <code>String</code>.
	 * @param cluster - The cluster.
	 * @return Returns the cluster information as a <code>String</code>.
	 */
	private static final String clusterToString(Cluster cluster) {
		return String.format("Cluster at (%3d, %3d) with %.3f GeV and %d hits at %4.0f ns.",
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix"),
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"),
				cluster.getEnergy(), cluster.getCalorimeterHits().size(),
				cluster.getCalorimeterHits().get(0).getTime());
	}
	
	/**
	 * Convenience method that writes the information in a cluster to
	 * a <code>String</code>.
	 * @param cluster - The cluster.
	 * @return Returns the cluster information as a <code>String</code>.
	 */
	private static final String clusterToString(SSPCluster cluster) {
		return String.format("Cluster at (%3d, %3d) with %.3f GeV and %d hits at %4d ns.",
				cluster.getXIndex(), cluster.getYIndex(), cluster.getEnergy(),
				cluster.getHitCount(), cluster.getTime());
	}
	
	/**
	 * Convenience method that writes the position of a reconstructed
	 * cluster in the form (ix, iy).
	 * @param cluster - The cluster.
	 * @return Returns the cluster position as a <code>String</code>.
	 */
	private static final String reconClusterPositionString(Cluster cluster) {
		return String.format("(%3d, %3d)",
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix"),
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"));
	}
	
	/**
	 * Convenience method that writes the position of an SSP bank
	 * cluster in the form (ix, iy).
	 * @param cluster - The cluster.
	 * @return Returns the cluster position as a <code>String</code>.
	 */
	private static final String sspClusterPositionString(SSPCluster cluster) {
		return String.format("(%3d, %3d)", cluster.getXIndex(), cluster.getYIndex());
	}
	
	/**
	 * Compares a trigger from the SSP bank to a trigger simulated on
	 * an SSP cluster.
	 * @param bankTrigger - The trigger from the SSP bank.
	 * @param simTrigger - The trigger from the simulation.
	 * @return Returns <code>true</code> if the triggers match and
	 * <code>false</code> if they do not.
	 */
	private static final boolean compareSSPSinglesTriggers(SSPSinglesTrigger bankTrigger, SinglesTrigger<SSPCluster> simTrigger) {
		// The bank trigger and simulated trigger must have the same
		// time. This is equivalent to the time of the triggering cluster.
		if(bankTrigger.getTime() != simTrigger.getTriggerSource().getTime()) {
			return false;
		}
		
		// If the time stamp is the same, check that the trigger flags
		// are all the same. Start with cluster energy low.
		if(bankTrigger.passCutEnergyMin() != simTrigger.getStateClusterEnergyLow()) {
			return false;
		}
		
		// Check cluster energy high.
		if(bankTrigger.passCutEnergyMax() != simTrigger.getStateClusterEnergyHigh()) {
			return false;
		}
		
		// Check cluster hit count.
		if(bankTrigger.passCutHitCount() != simTrigger.getStateHitCount()) {
			return false;
		}
		
		// If all of the tests are successful, the triggers match.
		return true;
	}
	
	/**
	 * Compares a trigger from the SSP bank to a trigger simulated on
	 * an SSP cluster.
	 * @param bankTrigger - The trigger from the SSP bank.
	 * @param simTrigger - The trigger from the simulation.
	 * @return Returns <code>true</code> if the triggers match and
	 * <code>false</code> if they do not.
	 */
	private static final boolean compareSSPPairTriggers(SSPPairTrigger bankTrigger, PairTrigger<SSPCluster[]> simTrigger) {
		// Get the time of the bottom cluster in the pair.
		int simTime = 0;
		if(simTrigger.getTriggerSource()[0].getYIndex() < 0) {
			simTime = simTrigger.getTriggerSource()[0].getTime();
		} else {
			simTime = simTrigger.getTriggerSource()[1].getTime();
		}
		
		// The bank trigger and simulated trigger must have the same
		// time. This is equivalent to the time of the triggering cluster.
		if(bankTrigger.getTime() != simTime) { return false; }
		
		// If the time stamp is the same, check that the trigger flags
		// are all the same. Start with energy sum.
		if(bankTrigger.passCutEnergySum() != simTrigger.getStateEnergySum()) {
			return false;
		}
		
		// Check pair energy difference.
		if(bankTrigger.passCutEnergyDifference() != simTrigger.getStateEnergyDifference()) {
			return false;
		}
		
		// Check pair energy slope.
		if(bankTrigger.passCutEnergySlope() != simTrigger.getStateEnergySlope()) {
			return false;
		}
		
		// Check pair coplanarity.
		if(bankTrigger.passCutCoplanarity() != simTrigger.getStateCoplanarity()) {
			return false;
		}
		
		// If all of the tests are successful, the triggers match.
		return true;
	}
	
	/**
	 * Compares a trigger from the SSP bank to a trigger simulated on
	 * an reconstructed cluster.
	 * @param bankTrigger - The trigger from the SSP bank.
	 * @param simTrigger - The trigger from the simulation.
	 * @return Returns <code>true</code> if the triggers match and
	 * <code>false</code> if they do not.
	 */
	private static final boolean compareReconPairTriggers(SSPPairTrigger bankTrigger, PairTrigger<Cluster[]> simTrigger) {
		// Check that the trigger flags are all the same. Start with
		// energy sum.
		if(bankTrigger.passCutEnergySum() != simTrigger.getStateEnergySum()) {
			return false;
		}
		
		// Check pair energy difference.
		if(bankTrigger.passCutEnergyDifference() != simTrigger.getStateEnergyDifference()) {
			return false;
		}
		
		// Check pair energy slope.
		if(bankTrigger.passCutEnergySlope() != simTrigger.getStateEnergySlope()) {
			return false;
		}
		
		// Check pair coplanarity.
		if(bankTrigger.passCutCoplanarity() != simTrigger.getStateCoplanarity()) {
			return false;
		}
		
		// If all of the tests are successful, the triggers match.
		return true;
	}
	
	/**
	 * Generates a <code>List</code> collection that contains a set
	 * of <code>ArrayList</code> collections representing a unique
	 * permutation of the entries in the argument.
	 * @param values - A collection of the entries to be permuted.
	 * @return Returns a list of lists representing the permutations.
	 */
	private static final List<List<Pair>> getPermutations(List<Cluster> reconClusters, List<SSPCluster> sspClusters) {
		// Store the SSP cluster permutations.
		List<List<SSPCluster>> permList = new ArrayList<List<SSPCluster>>();
		
		// Make sure that the two lists are the same size.
		int reconSize = reconClusters.size();
		int sspSize = sspClusters.size();
		while(sspClusters.size() < reconClusters.size()) {
			sspClusters.add(null);
		}
		while(reconClusters.size() < sspClusters.size()) {
			reconClusters.add(null);
		}
		
		// Get the SSP cluster permutations.
		permute(new ArrayList<SSPCluster>(0), sspClusters, permList);
		
		// Create pairs from the permutations.
		List<List<Pair>> pairList = new ArrayList<List<Pair>>();
		for(List<SSPCluster> permutation : permList) {
			List<Pair> pairs = new ArrayList<Pair>(reconClusters.size());
			
			for(int clusterIndex = 0; (clusterIndex < reconClusters.size() && clusterIndex < permutation.size()); clusterIndex++) {
				pairs.add(new Pair(reconClusters.get(clusterIndex), permutation.get(clusterIndex)));
			}
			
			pairList.add(pairs);
		}
		
		// Remove the extra values.
		for(int i = sspClusters.size() - 1; i >= sspSize; i--) { sspClusters.remove(i); }
		for(int i = reconClusters.size() - 1; i >= reconSize; i--) { reconClusters.remove(i); }
		
		// Return the pairs.
		return pairList;
	}
	
	/**
	 * Recursive method for permuting all entries in the argument
	 * collection <code>remainingValues</code> into the argument
	 * <code>permutedValues</code> values. Completed permutations are
	 * placed in the argument <code>permList</code>.
	 * @param permutedValues - List to store entries that have already
	 * been permuted.
	 * @param remainingValues - List to store  entries that need to be
	 * permuted.
	 * @param permList - List to store completed permutations.
	 */
	private static final void permute(List<SSPCluster> permutedValues, List<SSPCluster> remainingValues, List<List<SSPCluster>> permList) {
		// If the list of entries that still need to be sorted is empty,
		// then there is nothing to sort. Just return and empty list.
		if(remainingValues.isEmpty()) { return; }
		
		// If there is only one value left in the list of entries that
		// still need to be sorted, then just add it to the permutation
		// list and return it.
		else if(remainingValues.size() <= 1) {
			// Add the last entry.
			permutedValues.add(remainingValues.get(0));
			
			// Add the permutation to the list of completed permutations.
			permList.add(permutedValues);
		}
		
		// Otherwise, continue to get all possible permutations.
		else {
			// Iterate over the entries that have not been permuted.
			for(int i = 0; i < remainingValues.size(); i++) {
				// Make new lists to contain the permutations.
				List<SSPCluster> newPermList = new ArrayList<SSPCluster>(permutedValues.size() + 1);
				List<SSPCluster> newRemainList = new ArrayList<SSPCluster>(remainingValues.size());
				
				// Copy the current permuted entries to the new list
				// and one value from the list of entries that have
				// not been permuted yet.
				newPermList.addAll(permutedValues);
				newPermList.add(remainingValues.get(i));
				
				// The new list of entries that have not been permuted
				// should be identical, except it should now be missing
				// the entry that was moved.
				for(int index = 0; index < remainingValues.size(); index++) {
					if(index != i) { newRemainList.add(remainingValues.get(index)); }
				}
				
				// Repeat the process with the new lists.
				permute(newPermList, newRemainList, permList);
			}
		}
	}
	
	private void printf(String text, Object... args) {
		outputBuffer.append(String.format(text, args));
	}
	
	private void println() { printf(String.format("%n")); }
	
	private void println(String text) { printf(String.format("%s%n", text)); }
	
	private void print(String text) { printf(text); }
	
	/**
	 * Gets the time stamp of the cluster in nanoseconds.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the time-stamp.
	 */
	private static final double getReconTime(Cluster reconCluster) {
		return reconCluster.getCalorimeterHits().get(0).getTime();
	}
	
	/**
	 * Gets the x-index of the cluster's seed hit.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the x-index.
	 */
	private static final int getReconXIndex(Cluster reconCluster) {
		return reconCluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
	}
	
	/**
	 * Gets the y-index of the cluster's seed hit.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the y-index.
	 */
	private static final int getReconYIndex(Cluster reconCluster) {
		return reconCluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
	}
	
	/**
	 * Gets the number of hits in a reconstructed cluster.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the number of hits in the cluster.
	 */
	private static final int getReconHitCount(Cluster reconCluster) {
		return reconCluster.getCalorimeterHits().size();
	}
	
	/**
	 * Gets the difference in hit count between an SSP cluster and a
	 * reconstructed cluster.
	 * @param sspCluster - The SSP cluster.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the difference as an <code>int</code> primitive.
	 */
	private static final int getHitCountDifference(SSPCluster sspCluster, Cluster reconCluster) {
		return sspCluster.getHitCount() - getReconHitCount(reconCluster);
	}
	
	/**
	 * Solves the equation <code>|E_ssp / E_recon|</code>.
	 * @param sspCluster - The SSP cluster.
	 * @param reconCluster - The reconstructed cluster.
	 * @return Returns the solution to the equation as a <code>double
	 * </code> primitive.
	 */
	private static final double getEnergyPercentDifference(SSPCluster sspCluster, Cluster reconCluster) {
		return Math.abs((sspCluster.getEnergy() / reconCluster.getEnergy()));
	}
	
	/**
	 * Class <code>Pair</code> provides a convenient means of putting
	 * a reconstructed cluster and an SSP cluster in the same object
	 * for cluster matching.
	 * 
	 * @author Kyle McCarty <mccarty@jlab.org>
	 */
	private static class Pair {
		public final Cluster reconCluster;
		public final SSPCluster sspCluster;
		
		/**
		 * Instantiates a <code>Pair</code> consisting of the two
		 * cluster objects specified.
		 * @param reconCluster - A reconstructed cluster.
		 * @param sspCluster - An SSP bank cluster.
		 */
		public Pair(Cluster reconCluster, SSPCluster sspCluster) {
			this.reconCluster = reconCluster;
			this.sspCluster = sspCluster;
		}
	}
}