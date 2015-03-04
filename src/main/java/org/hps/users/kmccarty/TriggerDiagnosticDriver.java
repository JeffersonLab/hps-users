package org.hps.users.kmccarty;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Calendar;
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
import org.hps.readout.ecal.triggerbank.SSPNumberedTrigger;
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
	private String hitCollectionName = "EcalCalHits";
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
	private int noiseThreshold = 50;
	private long localWindowStart = 0;
	private double energyAcceptance = 0.03;
	private int localWindowThreshold = 5000;
	private boolean performClusterVerification = true;
	private boolean performSinglesTriggerVerification = true;
	private boolean performPairTriggerVerification = true;
	
	// Efficiency tracking variables.
	private ClusterMatchStatus clusterRunStats = new ClusterMatchStatus();
	private ClusterMatchStatus clusterLocalStats = new ClusterMatchStatus();
	
	private int singlesSSPTriggers = 0;
	private int singlesReconMatched = 0;
	private int singlesReconTriggers = 0;
	private int singlesInternalMatched = 0;
	private int singlesReportedTriggers = 0;
	private int singlesReportedExtras = 0;
	
	private int pairSSPTriggers = 0;
	private int pairReconMatched = 0;
	private int pairReconTriggers = 0;
	private int pairInternalMatched = 0;
	private int pairReportedTriggers = 0;
	private int pairReportedExtras = 0;
	
	private int[] globalEnergyMinCut = new int[2];
	private int[] globalEnergyMaxCut = new int[2];
	private int[] globalHitCountCut = new int[2];
	private int[] globalSinglesTimeCut = new int[2];
	
	private int[] globalEnergySumCut = new int[2];
	private int[] globalEnergyDiffCut = new int[2];
	private int[] globalEnergySlopeCut = new int[2];
	private int[] globalCoplanarityCut = new int[2];
	private int[] globalPairTimeCut = new int[2];
	
	private int failedClusterEvents = 0;
	private int failedSinglesEvents = 0;
	private int failedPairEvents = 0;
	private int totalEvents = 0;
	private int noiseEvents = 0;
	
	// Diagnostic plots.
    private AIDA aida = AIDA.defaultInstance();
    IHistogram1D cyclesRemoved;
    IHistogram2D cyclesRemovedEnergy;
    IHistogram2D cyclesRemovedSeed;
    IHistogram2D[] hitsInCycle = new IHistogram2D[9];
    IHistogram1D repeatedHitsPlot;
    IHistogram2D cyclesRemovedEnergyPercent;
    
    // Verbose settings.
    private boolean clusterFail = false;
    private boolean singlesEfficiencyFail = false;
    private boolean singlesInternalFail = false;
    private boolean pairEfficiencyFail = false;
    private boolean pairInternalFail = false;
    private boolean verbose = false;
    private boolean printClusterFail = true;
    private boolean printSinglesTriggerEfficiencyFail = true;
    private boolean printSinglesTriggerInternalFail = true;
    private boolean printPairTriggerEfficiencyFail = true;
    private boolean printPairTriggerInternalFail = true;
    private StringBuffer outputBuffer = new StringBuffer();
    
    // Cut index arrays for trigger verification.
	private static final int ENERGY_MIN = 0;
	private static final int ENERGY_MAX = 1;
	private static final int HIT_COUNT = 2;
	private static final int ENERGY_SUM = 0;
	private static final int ENERGY_DIFF = 1;
	private static final int ENERGY_SLOPE = 2;
	private static final int COPLANARITY = 3;
	private static final int TIME = 4;
    
	/**
	 * Define the trigger modules. This should be replaced by parsing
	 * the DAQ configuration at some point.
	 */
	@Override
	public void startOfData() {
		// Instantiate the diagnostic plots.
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
		
		// Set the FADC settings.
		nsa = 100;
		nsb = 20;
		windowWidth = 400;
		
		/*
		// Define the first singles trigger.
		singlesTrigger[0] = new TriggerModule();
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.010);
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 8.191);
		singlesTrigger[0].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 2);
		
		// Define the second singles trigger.
		singlesTrigger[1] = new TriggerModule();
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.010);
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 0.050);
		singlesTrigger[1].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 2);
		
		// Define the first pairs trigger.
		pairsTrigger[0] = new TriggerModule();
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.020);
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 0.055);
		pairsTrigger[0].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 1);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SUM_LOW, 0.010);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SUM_HIGH, 2.000);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_DIFFERENCE_HIGH, 1.200);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_LOW, 0.400);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_F, 0.0055);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_COPLANARITY_HIGH, 40);
		pairsTrigger[0].setCutValue(TriggerModule.PAIR_TIME_COINCIDENCE, 16);
		
		// Define the second pairs trigger.
		pairsTrigger[1] = new TriggerModule();
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_LOW, 0.010);
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_TOTAL_ENERGY_HIGH, 1.800);
		pairsTrigger[1].setCutValue(TriggerModule.CLUSTER_HIT_COUNT_LOW, 2);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SUM_LOW, 0.020);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SUM_HIGH, 2.000);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_DIFFERENCE_HIGH, 1.200);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_LOW, 0.400);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_ENERGY_SLOPE_F, 0.0055);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_COPLANARITY_HIGH, 40);
		pairsTrigger[1].setCutValue(TriggerModule.PAIR_TIME_COINCIDENCE, 16);
		*/
		
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
		
		// Print the general event failure rate.
		System.out.println("Event Failure Rate:");
		System.out.printf("\tNoise Events          :: %d / %d (%7.3f%%)%n",
				noiseEvents, totalEvents, (100.0 * noiseEvents / totalEvents));
		System.out.printf("\tCluster Events Failed :: %d / %d (%7.3f%%)%n",
				failedClusterEvents, totalEvents, (100.0 * failedClusterEvents / totalEvents));
		System.out.printf("\tSingles Events Failed :: %d / %d (%7.3f%%)%n",
				failedSinglesEvents, totalEvents, (100.0 * failedSinglesEvents / totalEvents));
		System.out.printf("\tPair Events Failed    :: %d / %d (%7.3f%%)%n",
				failedPairEvents, totalEvents, (100.0 * failedPairEvents / totalEvents));
		
		// Print the cluster verification data.
		System.out.println("Cluster Verification:");
	
		System.out.printf("\tRecon Clusters        :: %d%n", clusterRunStats.getReconClusterCount());
		System.out.printf("\tSSP Clusters          :: %d%n", clusterRunStats.getSSPClusterCount());
		System.out.printf("\tClusters Matched      :: %d%n", clusterRunStats.getMatches());
		System.out.printf("\tFailed (Position)     :: %d%n", clusterRunStats.getPositionFailures());
		System.out.printf("\tFailed (Energy)       :: %d%n", clusterRunStats.getEnergyFailures());
		System.out.printf("\tFailed (Hit Count)    :: %d%n", clusterRunStats.getHitCountFailures());
		if(clusterRunStats.getReconClusterCount() == 0) { System.out.printf("\tCluster Efficiency    :: N/A%n"); }
		else { System.out.printf("\tCluster Efficiency :: %7.3f%%%n", 100.0 * clusterRunStats.getMatches() / clusterRunStats.getReconClusterCount()); }
		
		// Print the singles trigger verification data.
		int spaces = getPrintSpaces(singlesSSPTriggers, singlesReconTriggers,
				singlesReportedTriggers, singlesInternalMatched, singlesReconMatched);
		System.out.println();
		System.out.println("Singles Trigger Verification:");
		System.out.printf("\tSSP Cluster Sim Triggers   :: %" + spaces + "d%n", singlesSSPTriggers);
		System.out.printf("\tRecon Cluster Sim Triggers :: %" + spaces + "d%n", singlesReconTriggers);
		System.out.printf("\tSSP Reported Triggers      :: %" + spaces + "d%n", singlesReportedTriggers);
		System.out.printf("\tExtra Reported Triggers    :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
				singlesReportedExtras, singlesReportedTriggers, (100.0 * singlesReportedExtras/ singlesReportedTriggers));
		
		System.out.printf("\tInternal Efficiency        :: %" + spaces + "d / %" + spaces + "d ",
				singlesInternalMatched, singlesSSPTriggers);
		if(singlesSSPTriggers == 0) { System.out.printf("(N/A)%n"); }
		else { System.out.printf("(%7.3f%%)%n", (100.0 * singlesInternalMatched / singlesSSPTriggers)); }
		
		System.out.printf("\tTrigger Efficiency         :: %" + spaces + "d / %" + spaces + "d ",
				singlesReconMatched, singlesReconTriggers);
		if(singlesReconTriggers == 0) { System.out.printf("(N/A)%n"); }
		else { System.out.printf("(%7.3f%%)%n" , (100.0 * singlesReconMatched / singlesReconTriggers)); }
		
		// Print the individual cut performances.
		int halfsinglesSSPTriggers = singlesSSPTriggers / 2;
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			System.out.println();
			System.out.printf("\tTrigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(singlesSSPTriggers == 0) {
				System.out.printf("\t\tCluster Energy Lower Bound :: %" + spaces + "d / %" + spaces + "d%n",
						globalEnergyMinCut[triggerNum], halfsinglesSSPTriggers);
				System.out.printf("\t\tCluster Energy Upper Bound :: %" + spaces + "d / %" + spaces + "d%n",
						globalEnergyMaxCut[triggerNum], halfsinglesSSPTriggers);
				System.out.printf("\t\tCluster Hit Count          :: %" + spaces + "d / %" + spaces + "d%n",
						globalHitCountCut[triggerNum], halfsinglesSSPTriggers);
			} else {
				System.out.printf("\t\tCluster Energy Lower Bound :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalEnergyMinCut[triggerNum], halfsinglesSSPTriggers, (100.0 * globalEnergyMinCut[triggerNum] / halfsinglesSSPTriggers));
				System.out.printf("\t\tCluster Energy Upper Bound :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalEnergyMaxCut[triggerNum], halfsinglesSSPTriggers, (100.0 * globalEnergyMaxCut[triggerNum] / halfsinglesSSPTriggers));
				System.out.printf("\t\tCluster Hit Count          :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalHitCountCut[triggerNum], halfsinglesSSPTriggers, (100.0 * globalHitCountCut[triggerNum] / halfsinglesSSPTriggers));
			}
		}
		
		// Print the pair trigger verification data.
		spaces = getPrintSpaces(pairReportedExtras, pairReportedTriggers, pairInternalMatched, pairSSPTriggers);
		System.out.println();
		System.out.println("Pair Trigger Verification:");
		System.out.printf("\tSSP Cluster Sim Triggers   :: %" + spaces + "d%n", pairSSPTriggers);
		System.out.printf("\tRecon Cluster Sim Triggers :: %" + spaces + "d%n", pairReconTriggers);
		System.out.printf("\tSSP Reported Triggers      :: %" + spaces + "d%n", pairReportedTriggers);
		System.out.printf("\tExtra Reported Triggers    :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
				pairReportedExtras, pairReportedTriggers, (100.0 * pairReportedExtras/ pairReportedTriggers));
		
		System.out.printf("\tInternal Efficiency        :: %" + spaces + "d / %" + spaces + "d ",
				pairInternalMatched, pairSSPTriggers);
		if(pairSSPTriggers == 0) { System.out.printf("(N/A)%n"); }
		else { System.out.printf("(%7.3f%%)%n", (100.0 * pairInternalMatched / pairSSPTriggers)); }
		
		System.out.printf("\tTrigger Efficiency         :: %" + spaces + "d / %" + spaces + "d ",
				pairReconMatched, pairReconTriggers);
		if(pairReconTriggers == 0) { System.out.printf("(N/A)%n"); }
		else { System.out.printf("(%7.3f%%)%n", (100.0 * pairReconMatched / pairReconTriggers)); }
		
		// Print the individual cut performances.
		int halfPairSSPTriggers = pairSSPTriggers / 2;
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			System.out.println();
			System.out.printf("\tTrigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
			if(pairSSPTriggers == 0) {
				System.out.printf("\t\tPair Energy Sum            :: %" + spaces + "d / %" + spaces + "d%n",
						globalEnergySumCut[triggerNum], halfPairSSPTriggers);
				System.out.printf("\t\tPair Energy Difference     :: %" + spaces + "d / %" + spaces + "d%n",
						globalEnergyDiffCut[triggerNum], halfPairSSPTriggers);
				System.out.printf("\t\tPair Energy Slope          :: %" + spaces + "d / %" + spaces + "d%n",
						globalEnergySlopeCut[triggerNum], halfPairSSPTriggers);
				System.out.printf("\t\tPair Coplanarity           :: %" + spaces + "d / %" + spaces + "d%n",
						globalCoplanarityCut[triggerNum], halfPairSSPTriggers);
				System.out.printf("\t\tPair Trigger Time          :: %" + spaces + "d / %" + spaces + "d%n",
						globalPairTimeCut[triggerNum], halfPairSSPTriggers);
			} else {
				System.out.printf("\t\tPair Energy Sum            :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalEnergySumCut[triggerNum], halfPairSSPTriggers, (100.0 * globalEnergySumCut[triggerNum] / halfPairSSPTriggers));
				System.out.printf("\t\tPair Energy Difference     :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalEnergyDiffCut[triggerNum], halfPairSSPTriggers, (100.0 * globalEnergyDiffCut[triggerNum] / halfPairSSPTriggers));
				System.out.printf("\t\tPair Energy Slope          :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalEnergySlopeCut[triggerNum], halfPairSSPTriggers, (100.0 * globalEnergySlopeCut[triggerNum] / halfPairSSPTriggers));
				System.out.printf("\t\tPair Coplanarity           :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalCoplanarityCut[triggerNum], halfPairSSPTriggers, (100.0 * globalCoplanarityCut[triggerNum] / halfPairSSPTriggers));
				System.out.printf("\t\tPair Trigger Time          :: %" + spaces + "d / %" + spaces + "d (%7.3f%%)%n",
						globalPairTimeCut[triggerNum], halfPairSSPTriggers, (100.0 * globalPairTimeCut[triggerNum] / halfPairSSPTriggers));
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
        
        // Print the verification header.
		println();
		println();
		println("======================================================================");
		println("==== Cluster/Trigger Verification ====================================");
		println("======================================================================");
		
		// Increment the total event count.
		totalEvents++;
		
		// Reset the output buffer and print flags.
		outputBuffer = new StringBuffer();
		clusterFail = false;
		singlesInternalFail = false;
		singlesEfficiencyFail = false;
		pairInternalFail = false;
		pairEfficiencyFail = false;
		
		
		
		// ==========================================================
		// ==== Output the Event Cluster List =======================
		// ==========================================================
		
		/*
		// Print the clusterer header.
		println();
		println();
		println("======================================================================");
		println("=== GTP Readout Clusterer ============================================");
		println("======================================================================");
		
        // Output all clusters in full detail.
		List<Cluster> allClusters = event.get(Cluster.class, clusterCollectionName);
    	println("Event Cluster Collection:");
        for(Cluster cluster : allClusters) {
            CalorimeterHit seedHit = cluster.getCalorimeterHits().get(0);
            int ix = seedHit.getIdentifierFieldValue("ix");
            int iy = seedHit.getIdentifierFieldValue("iy");
            double energy = cluster.getEnergy();
            double time = seedHit.getTime();
            
            printf("\tCluster --> %6.3f GeV at (%3d, %3d) and at t = %.2f%n", energy, ix, iy, time);
            
            for(CalorimeterHit hit : cluster.getCalorimeterHits()) {
                int hix = hit.getIdentifierFieldValue("ix");
                int hiy = hit.getIdentifierFieldValue("iy");
                double henergy = hit.getCorrectedEnergy();
                double htime = hit.getTime();
                printf("\t\tCompHit --> %.3f GeV at (%3d, %3d) and at t = %.2f%n", henergy, hix, hiy, htime);
            }
        }
        */
		
		
		
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
			if(verbose) { System.out.println(outputBuffer.toString()); }
			return;
		} if(tiBank == null) {
			println("No TI bank found for this event. No verification will be performed.");
			if(verbose) { System.out.println(outputBuffer.toString()); }
			return;
		}
		
		
		
		// ==========================================================
		// ==== Check the Noise Level ===============================
		// ==========================================================
		
		// Check if there are hits.
		if(event.hasCollection(CalorimeterHit.class, hitCollectionName)) {
			// Check if there are more hits than the noise threshold.
			if(event.get(CalorimeterHit.class, hitCollectionName).size() >= noiseThreshold) {
				noiseEvents++;
				println("Noise event detected. Skipping event...");
				if(verbose) { System.out.println(outputBuffer.toString()); }
				return;
			}
		}
        
        
        
		// ==========================================================
		// ==== Obtain Reconstructed Clusters =======================
		// ==========================================================
		
		// Clear the list of triggers from previous events.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			sspSinglesTriggers.get(triggerNum).clear();
			reconSinglesTriggers.get(triggerNum).clear();
			sspPairsTriggers.get(triggerNum).clear();
			reconPairsTriggers.get(triggerNum).clear();
		}
		
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
					double seedTime = TriggerDiagnosticUtil.getClusterTime(reconCluster);
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
		
		// Track how many events failed due to each type of verification.
		if(clusterFail) { failedClusterEvents++; }
		if(pairInternalFail || pairEfficiencyFail) { failedPairEvents++; }
		if(singlesInternalFail || singlesEfficiencyFail) { failedSinglesEvents++; }
		
		
		
		// ==========================================================
		// ==== Perform Event Write-Out =============================
		// ==========================================================
		
		if(verbose ||(clusterFail && printClusterFail) ||
				(singlesInternalFail && printSinglesTriggerInternalFail) ||
				(singlesEfficiencyFail && printSinglesTriggerEfficiencyFail) ||
				(pairInternalFail && printPairTriggerInternalFail) ||
				(pairEfficiencyFail && printPairTriggerEfficiencyFail)) {
			System.out.println(outputBuffer.toString());
		}
		
		
		
		// ==========================================================
		// ==== Process Local Tracked Variables =====================
		// ==========================================================
		if(Calendar.getInstance().getTimeInMillis() - localWindowStart > localWindowThreshold) {
			// Write a snapshot of the driver to the event stream.
			// TODO: Readout the snapshot!!
			
			// Clear the local statistical data.
			clusterLocalStats.clear();
		}
	}
	
	public void setPrintOnClusterFailure(boolean state) {
		printClusterFail = state;
	}
	
	public void setPrintOnSinglesEfficiencyFailure(boolean state) {
		printSinglesTriggerEfficiencyFail = state;
	}
	
	public void setPrintOnSinglesSSPFailure(boolean state) {
		printSinglesTriggerInternalFail = state;
	}
	
	public void setPrintOnPairEfficiencyFailure(boolean state) {
		printPairTriggerEfficiencyFail = state;
	}
	
	public void setPrintOnPairSSPFailure(boolean state) {
		printPairTriggerInternalFail = state;
	}
	
	public void setVerbose(boolean state) {
		verbose = state;
	}
	
	public void setHitCollectionName(String hitCollectionName) {
		this.hitCollectionName = hitCollectionName;
	}
	
	public void setClusterCollectionName(String clusterCollectionName) {
		this.clusterCollectionName = clusterCollectionName;
	}
	
	public void setBankCollectionName(String bankCollectionName) {
		this.bankCollectionName = bankCollectionName;
	}
	
	public void setNoiseThresholdCount(int noiseHits) {
		noiseThreshold = noiseHits;
	}
	
	public void setHitAcceptanceWindow(int window) {
		hitAcceptance = window;
	}
	
	public void setEnergyAcceptanceWindow(double window) {
		energyAcceptance = window;
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
		
		// Track the number of cluster pairs that were matched and that
		// failed by failure type.
		ClusterMatchEvent event = new ClusterMatchEvent();
		
		
		
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
			Point position = new Point(TriggerDiagnosticUtil.getXIndex(reconCluster),
					TriggerDiagnosticUtil.getYIndex(reconCluster));
			
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
				clusterFail = true;
				for(Cluster cluster : reconList) {
					event.pairFailPosition(cluster, null);
				}
				continue positionLoop;
			}
			
			// If there are more reconstructed clusters than there are
			// SSP clusters, than a number equal to the difference must
			// fail by means of positions.
			if(sspList.size() < reconList.size()) { clusterFail = true; }
			
			// Get all possible permutations of SSP clusters.
			List<List<Pair>> permutations = getPermutations(reconList, sspList);
			
			// Print the information for this crystal position.
			printf("\tRecon Clusters :: %d%n", reconList.size());
			printf("\tSSP Clusters   :: %d%n", sspList.size());
			printf("\tPermutations   :: %d%n", permutations.size());
			
			// Track the plotted values for the current best permutation.
			ClusterMatchEvent bestPerm = null;
			
			// Iterate over the permutations and find the permutation
			// that produces the best possible result when compared to
			// the reconstructed clusters.
			int permIndex = 0;
			for(List<Pair> pairs : permutations) {
				// Update the current permutation number.
				permIndex++;
				
				// Track the plot values for this permutation.
				ClusterMatchEvent perm = new ClusterMatchEvent();
				
				// Try to match each pair.
				pairLoop:
				for(Pair pair : pairs) {
					// Print the current reconstructed/SSP cluster pair.
					printf("\tP%d :: %s --> %s", permIndex,
							pair.reconCluster == null ? "None" : clusterToString(pair.reconCluster),
							pair.sspCluster == null ? "None" : clusterToString(pair.sspCluster));
					
					// If either cluster in the pair is null, there
					// are not enough clusters to perform this match.
					if(pair.reconCluster == null || pair.sspCluster == null) {
						printf(" [ %18s ]%n", "failure: unpaired");
						perm.pairFailPosition(pair.reconCluster, pair.sspCluster);
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
							perm.pairMatch(pair.reconCluster, pair.sspCluster);
							printf(" [ %18s ]%n", "success: matched");
						} else {
							perm.pairFailHitCount(pair.reconCluster, pair.sspCluster);
							printf(" [ %18s ]%n", "failure: hit count");
						} // End hit count check.
					} else {
						perm.pairFailEnergy(pair.reconCluster, pair.sspCluster);
						printf(" [ %18s ]%n", "failure: energy");
					} // End energy check.
				} // End Pair Loop
				
				// Print the results of the permutation.
				printf("\t\tPermutation Matched   :: %d%n", perm.getMatches());
				printf("\t\tPermutation Energy    :: %d%n", perm.getEnergyFailures());
				printf("\t\tPermutation Hit Count :: %d%n", perm.getHitCountFailures());
				
				// Check whether the results from this permutation
				// exceed the quality of the last best results. A
				// greater number of matches is always better. If the
				// matches are the same, select the one with fewer
				// failures due to energy.
				bestPerm = getBestPermutation(bestPerm, perm);
			} // End Permutation Loop
			
			// Print the final results for the position.
			printf("\tPosition Matched   :: %d%n", bestPerm.getMatches());
			printf("\tPosition Energy    :: %d%n", bestPerm.getEnergyFailures());
			printf("\tPosition Hit Count :: %d%n", bestPerm.getHitCountFailures());
			
			// Add the results from the best-matched permutation
			// to the event efficiency results.
			event.addEvent(bestPerm);
		} // End Crystal Position Loop
		
		// Add the event results to the global results.
		clusterRunStats.addEvent(event, reconClusters, sspClusters);
		clusterLocalStats.addEvent(event, reconClusters, sspClusters);
		
		
		
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
			}
		} else { println("\tNone"); }
		
		// Print the SSP clusters and populate their distribution graphs.
		println("SSP Clusters:");
		if(!sspClusters.isEmpty()) {
			for(SSPCluster sspCluster : sspClusters) {
				printf("\t%s%n", clusterToString(sspCluster));
			}
		} else { println("\tNone"); }
		
		// Print the matched clusters.
		println("Matched Clusters:");
		if(event.getMatchedPairs().size() != 0) {
			// Iterate over the matched pairs.
			for(ClusterMatchedPair pair : event.getMatchedPairs()) {
				// If the pair is a match, print it out.
				if(pair.isMatch()) {
					printf("\t%s --> %s%n",
							clusterToString(pair.getReconstructedCluster()),
							clusterToString(pair.getSSPCluster()));
				}
			}
		}
		 else { println("\tNone"); }
		
		// Get the number of position failures.
		int failPosition = event.getPositionFailures();
		if(sspClusters == null || sspClusters.isEmpty()) {
			failPosition = (reconClusters == null ? 0 : reconClusters.size());
		}
		
		// Print event statistics.
		println();
		println("Event Statistics:");
		printf("\tRecon Clusters     :: %d%n", reconClusters.size());
		printf("\tClusters Matched   :: %d%n", event.getMatches());
		printf("\tFailed (Position)  :: %d%n", failPosition);
		printf("\tFailed (Energy)    :: %d%n", event.getEnergyFailures());
		printf("\tFailed (Hit Count) :: %d%n", event.getHitCountFailures());
		printf("\tCluster Efficiency :: %3.0f%%%n", 100.0 * event.getMatches() / reconClusters.size());
		
		// Note whether there was a cluster match failure.
		if(event.getMatches() - reconClusters.size() != 0) {
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
		// Create lists of generic triggers.
		List<List<? extends Trigger<?>>> sspTriggerList = new ArrayList<List<? extends Trigger<?>>>(2);
		List<List<? extends Trigger<?>>> reconTriggerList = new ArrayList<List<? extends Trigger<?>>>(2);
		
		// Convert the simulated triggers to generic versions and add
		// them to the generic list.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			// Get the generic trigger list.
			List<? extends Trigger<?>> sspTriggers = sspSinglesTriggers.get(triggerNum);
			List<? extends Trigger<?>> reconTriggers = reconSinglesTriggers.get(triggerNum);
			
			// Add it to the generic list.
			sspTriggerList.add(sspTriggers);
			reconTriggerList.add(reconTriggers);
		}
		
		// Run generic trigger verification.
		triggerVerification(sspTriggerList, reconTriggerList, true);
	}
	
	/**
	 * Checks triggers simulated on SSP clusters against the SSP bank's
	 * reported triggers to verify that the trigger is correctly applying
	 * cuts to the clusters it sees. Additionally compares triggers
	 * simulated on reconstructed clusters to measure trigger efficiency.
	 */
	private void pairTriggerVerification() {
		// Create lists of generic triggers.
		List<List<? extends Trigger<?>>> sspTriggerList = new ArrayList<List<? extends Trigger<?>>>(2);
		List<List<? extends Trigger<?>>> reconTriggerList = new ArrayList<List<? extends Trigger<?>>>(2);
		
		// Convert the simulated triggers to generic versions and add
		// them to the generic list.
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			// Get the generic trigger list.
			List<? extends Trigger<?>> sspTriggers = sspPairsTriggers.get(triggerNum);
			List<? extends Trigger<?>> reconTriggers = reconPairsTriggers.get(triggerNum);
			
			// Add it to the generic list.
			sspTriggerList.add(sspTriggers);
			reconTriggerList.add(reconTriggers);
		}
		
		// Run generic trigger verification.
		triggerVerification(sspTriggerList, reconTriggerList, false);
	}
	
	private void triggerVerification(List<List<? extends Trigger<?>>> sspTriggerList, 
			List<List<? extends Trigger<?>>> reconTriggerList, boolean isSingles) {
		
		// ==========================================================
		// ==== Initialize Trigger Verification =====================
		// ==========================================================
		
		// Print the cluster verification header.
		println();
		println();
		println("======================================================================");
		if(isSingles) { println("=== Singles Trigger Verification ====================================="); }
		else { println("=== Pair Trigger Verification ========================================"); }
		println("======================================================================");
		
		// Track the number of triggers seen and the number found.
		int sspReportedTriggers = 0;
		int sspInternalMatched = 0;
		int reconTriggersMatched = 0;
		int sspReportedExtras = 0;
		
		// Track the number of times a given cut caused a trigger to
		// fail to match.
		int[][] triggerComp = new int[5][2];
		
		
		
		// ==========================================================
		// ==== Output Event Summary ================================
		// ==========================================================
		
		// Get the list of triggers reported by the SSP.
		List<? extends SSPNumberedTrigger> sspTriggers;
		if(isSingles) { sspTriggers = sspBank.getSinglesTriggers(); }
		else { sspTriggers = sspBank.getPairTriggers(); }
		
		// Output the SSP cluster singles triggers.
		println();
		println("SSP Cluster " + (isSingles ? "Singles" : "Pair") + " Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(Trigger<?> simTrigger : sspTriggerList.get(triggerNum)) {
				printf("\tTrigger %d :: %s :: %s%n",
						(triggerNum + 1), triggerPositionString(simTrigger),
						simTrigger.toString());
			}
		}
		if(sspTriggerList.get(0).size() + sspTriggerList.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the reconstructed cluster singles triggers.
		println("Reconstructed Cluster " + (isSingles ? "Singles" : "Pair") + " Triggers");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(Trigger<?> simTrigger : reconTriggerList.get(triggerNum)) {
				printf("\tTrigger %d :: %s :: %s%n",
						(triggerNum + 1), triggerPositionString(simTrigger),
						simTrigger.toString());
			}
		}
		if(reconTriggerList.get(0).size() + reconTriggerList.get(1).size() == 0) {
			println("\tNone");
		}
		
		// Output the SSP reported triggers.
		println("SSP Reported " + (isSingles ? "Singles" : "Pair") + " Triggers");
		for(SSPTrigger sspTrigger : sspTriggers) {
			// Increment the number of SSP cluster singles triggers.
			sspReportedTriggers++;
			
			// Print the trigger.
			printf("\t%s%n", sspTrigger.toString());
		}
		if(sspReportedTriggers == 0) { println("\tNone"); }
		
		
		
		// ==========================================================
		// ==== SSP Internal Logic Verification =====================
		// ==========================================================
		
		// Track which SSP triggers have been matched to avoid matching
		// multiple reconstructed SSP cluster triggers to the same SSP
		// trigger.
		Set<SSPNumberedTrigger> sspTriggerSet = new HashSet<SSPNumberedTrigger>();
		Set<Trigger<?>> simTriggerSet = new HashSet<Trigger<?>>();
		
		// Track the number of SSP reported triggers that are found in
		// excess of the SSP simulated triggers.
		int extraTriggers = sspTriggers.size() - (sspTriggerList.get(0).size() + sspTriggerList.get(1).size());
		if(extraTriggers > 0) {
			sspReportedExtras += extraTriggers;
			if(isSingles) { singlesInternalFail = true; }
			else { pairInternalFail = true; }
		}
		
		// Iterate over the triggers.
		println();
		println("SSP Reported Trigger --> SSP Cluster Trigger Match Status");
		for(SSPNumberedTrigger sspTrigger : sspTriggers) {
			// Get the trigger information.
			int triggerNum = sspTrigger.isFirstTrigger() ? 0 : 1;
			boolean matchedTrigger = false;
			
			// Iterate over the SSP cluster simulated triggers and
			// look for a trigger that matches.
			matchLoop:
			for(Trigger<?> simTrigger : sspTriggerList.get(triggerNum)) {
				// If the current SSP trigger has already been
				// matched, skip it.
				if(sspTriggerSet.contains(sspTrigger)) { continue matchLoop; }
				
				// Otherwise, check whether the reconstructed SSP
				// cluster trigger matches the SSP trigger.
				if(compareTriggers(sspTrigger, simTrigger)) {
					matchedTrigger = true;
					sspTriggerSet.add(sspTrigger);
					simTriggerSet.add(simTrigger);
					sspInternalMatched++;
					break matchLoop;
				}
				
				printf("\t%s :: Matched: %5b%n", sspTrigger.toString(), matchedTrigger);
			}
		}
		
		// Iterate over the unmatched simulated triggers again and the
		// unmatched SSP reported trigger that most closely matches it.
		simLoop:
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(Trigger<?> simTrigger : sspTriggerList.get(triggerNum)) {
				// Check whether this trigger has already been matched
				// or not. If it has been matched, skip it.
				if(simTriggerSet.contains(simTrigger)) { continue simLoop; }
				
				// Get the trigger time for the simulated trigger.
				double simTime = getTriggerTime(simTrigger);
				
				// Track the match statistics for each reported trigger
				// so that the closest match may be found.
				int numMatched = -1;
				boolean[] matchedCut = null;
				
				// Iterate over the reported triggers to find a match.
				reportedLoop:
				for(SSPNumberedTrigger sspTrigger : sspTriggers) {
					// If the two triggers have different times, this
					// trigger should be skipped.
					if(sspTrigger.getTime() != simTime) {
						continue reportedLoop;
					}
					
					// If this reported trigger has been matched then
					// it should be skipped.
					if(sspTriggerSet.contains(sspTrigger)) { continue reportedLoop; }
					
					// Check each of the cuts.
					boolean[] tempMatchedCut = triggerCutMatch(simTrigger, sspTrigger);
					
					// Check each cut and see if this is a closer match
					// than the previous best match.
					int tempNumMatched = 0;
					for(boolean passed : tempMatchedCut) { if(passed) { tempNumMatched++; } }
					
					// If the number of matched cuts exceeds the old
					// best result, this becomes the new best result.
					if(tempNumMatched > numMatched) {
						numMatched = tempNumMatched;
						matchedCut = tempMatchedCut;
					}
				}
				
				// If some match was found, note what caused it to not
				// qualify as a complete match.
				if(matchedCut != null) {
					for(int cutIndex = 0; cutIndex < matchedCut.length; cutIndex++) {
						if(!matchedCut[cutIndex]) { triggerComp[cutIndex][triggerNum]++; }
					}
				}
				
				// If there was no match found, it means that there were
				// no triggers that were both unmatched and at the same
				// time as this simulated trigger.
				else {
					triggerComp[TIME][triggerNum]++;
					if(isSingles) { singlesInternalFail = true; }
					else { pairInternalFail = true; }
				}
			}
		}
		
		
		
		// ==========================================================
		// ==== Trigger Efficiency ==================================
		// ==========================================================
		
		// Reset the SSP matched trigger set.
		sspTriggerSet.clear();
		
		// Iterate over the reconstructed cluster singles triggers.
		println();
		println("Recon Cluster Trigger --> SSP Reported Trigger Match Status");
		for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			for(Trigger<?> simTrigger : reconTriggerList.get(triggerNum)) {
				
				printf("\tTrigger %d :: %s :: %s%n", (triggerNum + 1),
						triggerPositionString(simTrigger), simTrigger.toString());
				
				// Iterate over the SSP reported triggers and compare
				// them to the reconstructed cluster simulated trigger.
				matchLoop:
				for(SSPNumberedTrigger sspTrigger : sspTriggers) {
					printf("\t\t\t%s", sspTrigger.toString());
					
					// Only compare triggers if they are from the
					// same trigger source.
					if((triggerNum == 0 && sspTrigger.isSecondTrigger())
							|| (triggerNum == 1 && sspTrigger.isFirstTrigger())) {
						print(" [ fail; source    ]%n");
						continue matchLoop;
					}
					
					// Only compare the singles trigger if it was
					// not already matched to another trigger.
					if(sspTriggerSet.contains(sspTrigger)) {
						print(" [ fail; matched   ]%n");
						continue matchLoop;
					}
					
					// Test each cut.
					String[][] cutNames = {
							{ "E_min", "E_max", "hit count", "null" },
							{ "E_sum", "E_diff", "E_slope", "coplanar" }
					};
					int typeIndex = isSingles ? 0 : 1;
					boolean[] matchedCuts = triggerCutMatch(simTrigger, sspTrigger);
					for(int cutIndex = 0; cutIndex < matchedCuts.length; cutIndex++) {
						if(!matchedCuts[cutIndex]) {
							printf(" [ fail; %-9s ]%n", cutNames[typeIndex][cutIndex]);
							continue matchLoop;
						}
					}
					
					// If all the trigger flags match, then the
					// triggers are a match.
					reconTriggersMatched++;
					sspTriggerSet.add(sspTrigger);
					print(" [ success         ]%n");
					break matchLoop;
				}
			}
		}
		
		
		
		// ==========================================================
		// ==== Output Event Results ================================
		// ==========================================================
		
		// Get the number of SSP and reconstructed cluster simulated
		// triggers.
		int sspSimTriggers = sspTriggerList.get(0).size() + sspTriggerList.get(1).size();
		int reconSimTriggers = reconTriggerList.get(0).size() + reconTriggerList.get(1).size();
		int halfSimTriggers = sspSimTriggers / 2;
		
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
		if(isSingles) {
			println();
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
			printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
				if(sspSimTriggers == 0) {
					printf("\tCluster Energy Lower Bound :: %d / %d%n", triggerComp[ENERGY_MIN][triggerNum], halfSimTriggers);
					printf("\tCluster Energy Upper Bound :: %d / %d%n", triggerComp[ENERGY_MAX][triggerNum], halfSimTriggers);
					printf("\tCluster Hit Count          :: %d / %d%n", triggerComp[HIT_COUNT][triggerNum], halfSimTriggers);
				} else {
					printf("\tCluster Energy Lower Bound :: %d / %d (%3.0f%%)%n",
							triggerComp[ENERGY_MIN][triggerNum], halfSimTriggers, (100.0 * triggerComp[ENERGY_MIN][triggerNum] / halfSimTriggers));
					printf("\tCluster Energy Upper Bound :: %d / %d (%3.0f%%)%n",
							triggerComp[ENERGY_MAX][triggerNum], halfSimTriggers, (100.0 * triggerComp[ENERGY_MAX][triggerNum] / halfSimTriggers));
					printf("\tCluster Hit Count          :: %d / %d (%3.0f%%)%n",
							triggerComp[HIT_COUNT][triggerNum], halfSimTriggers, (100.0 * triggerComp[HIT_COUNT][triggerNum] / halfSimTriggers));
				}
				printf("\tExcess Reported Triggers   :: %d%n", sspReportedExtras / 2);
			}
			
			// Update the global trigger tracking variables.
			singlesSSPTriggers += sspSimTriggers;
			singlesReconMatched += reconTriggersMatched;
			singlesReconTriggers += reconSimTriggers;
			singlesInternalMatched += sspInternalMatched;
			singlesReportedTriggers += sspReportedTriggers;
			singlesReportedExtras += sspReportedExtras;
			
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
				globalEnergyMinCut[triggerNum] += triggerComp[ENERGY_MIN][triggerNum];
				globalEnergyMaxCut[triggerNum] += triggerComp[ENERGY_MAX][triggerNum];
				globalHitCountCut[triggerNum] += triggerComp[HIT_COUNT][triggerNum];
				globalSinglesTimeCut[triggerNum] += triggerComp[TIME][triggerNum];
			}
		} else {
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
				println();
				printf("Trigger %d Individual Cut Failure Rate:%n", (triggerNum + 1));
				if(sspSimTriggers == 0) {
					printf("\tPair Energy Sum            :: %d / %d%n", triggerComp[ENERGY_SUM][triggerNum], halfSimTriggers);
					printf("\tPair Energy Difference     :: %d / %d%n", triggerComp[ENERGY_DIFF][triggerNum], halfSimTriggers);
					printf("\tPair Energy Slope          :: %d / %d%n", triggerComp[ENERGY_SLOPE][triggerNum], halfSimTriggers);
					printf("\tPair Coplanarity           :: %d / %d%n", triggerComp[COPLANARITY][triggerNum], halfSimTriggers);
					printf("\tPair Trigger Time          :: %d / %d%n", triggerComp[TIME][triggerNum], halfSimTriggers);
				} else {
					printf("\tPair Energy Sum            :: %d / %d (%3.0f%%)%n",
							triggerComp[ENERGY_SUM][triggerNum], halfSimTriggers, (100.0 * triggerComp[ENERGY_SUM][triggerNum] / halfSimTriggers));
					printf("\tPair Energy Difference     :: %d / %d (%3.0f%%)%n",
							triggerComp[ENERGY_DIFF][triggerNum], halfSimTriggers, (100.0 * triggerComp[ENERGY_DIFF][triggerNum] / halfSimTriggers));
					printf("\tPair Energy Slope          :: %d / %d (%3.0f%%)%n",
							triggerComp[ENERGY_SLOPE][triggerNum], halfSimTriggers, (100.0 * triggerComp[ENERGY_SLOPE][triggerNum] / halfSimTriggers));
					printf("\tPair Coplanarity           :: %d / %d (%3.0f%%)%n",
							triggerComp[COPLANARITY][triggerNum], halfSimTriggers, (100.0 * triggerComp[COPLANARITY][triggerNum] / halfSimTriggers));
					printf("\tPair Trigger Time          :: %d / %d (%3.0f%%)%n",
							triggerComp[TIME][triggerNum], halfSimTriggers, (100.0 * triggerComp[TIME][triggerNum] / halfSimTriggers));
				}
				printf("\tExcess Reported Triggers   :: %d%n", sspReportedExtras / 2);
			}
			
			// Update the global trigger tracking variables.
			pairSSPTriggers += sspSimTriggers;
			pairReconMatched += reconTriggersMatched;
			pairReconTriggers += reconSimTriggers;
			pairInternalMatched += sspInternalMatched;
			pairReportedTriggers += sspReportedTriggers;
			pairReportedExtras += sspReportedExtras;
			
			for(int triggerNum = 0; triggerNum < 2; triggerNum++) {
				globalEnergySumCut[triggerNum] += triggerComp[ENERGY_SUM][triggerNum];
				globalEnergyDiffCut[triggerNum] += triggerComp[ENERGY_DIFF][triggerNum];
				globalEnergySlopeCut[triggerNum] += triggerComp[ENERGY_SLOPE][triggerNum];
				globalCoplanarityCut[triggerNum] += triggerComp[COPLANARITY][triggerNum];
				globalPairTimeCut[triggerNum] += triggerComp[TIME][triggerNum];
			}
		}
		
		// Note whether the was a trigger match failure.
		if((reconTriggersMatched - reconSimTriggers != 0) || (sspInternalMatched - sspSimTriggers != 0)) {
			if(isSingles) { singlesEfficiencyFail = true; }
			else { pairEfficiencyFail = true; }
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
		System.out.printf("\tEnergy Threshold       :: %1.2f%%%n", energyAcceptance);
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
	 * Convenience method that writes the position of a cluster in the
	 * form (ix, iy).
	 * @param cluster - The cluster.
	 * @return Returns the cluster position as a <code>String</code>.
	 */
	private static final String clusterPositionString(Cluster cluster) {
		return String.format("(%3d, %3d)",
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix"),
				cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"));
	}
	
	/**
	 * Convenience method that writes the position of a cluster in the
	 * form (ix, iy).
	 * @param cluster - The cluster.
	 * @return Returns the cluster position as a <code>String</code>.
	 */
	private static final String clusterPositionString(SSPCluster cluster) {
		return String.format("(%3d, %3d)", cluster.getXIndex(), cluster.getYIndex());
	}
	
	/**
	 * Compares an SSP trigger with a simulated trigger. Note that only
	 * certain class combinations are supported. Triggers of the type
	 * <code>SSPSinglesTrigger</code> may be compared with triggers of
	 * the type <code>SinglesTrigger<SSPCluster></code> and triggers of
	 * the type <code>SSPPairTrigger</code> may be compared to either
	 * <code>PairTrigger<SSPCluster[]></code> triggers objects.
	 * @param bankTrigger - The SSP bank trigger.
	 * @param simTrigger - The simulated trigger.
	 * @return Returns <code>true</code> if the triggers are valid
	 * matches and <code>false</code> if they are not.
	 * @throws IllegalArgumentException Occurs if the trigger types
	 * are not of a supported type.
	 */
	@SuppressWarnings("unchecked")
	private static final boolean compareTriggers(SSPNumberedTrigger bankTrigger, Trigger<?> simTrigger) throws IllegalArgumentException {
		// Get the classes of the arguments. This is used to check the
		// generic type of the Trigger<?> object, and means that the
		// "unchecked" warnings can be safely ignored.
		Object source = simTrigger.getTriggerSource();
		
		// If the combination of classes is supported, pass the triggers
		// to the appropriate handler.
		if(bankTrigger instanceof SSPSinglesTrigger && simTrigger instanceof SinglesTrigger && source instanceof SSPCluster) {
			return compareSSPSinglesTriggers((SSPSinglesTrigger) bankTrigger, (SinglesTrigger<SSPCluster>) simTrigger);
		} else if(bankTrigger instanceof SSPPairTrigger && simTrigger instanceof PairTrigger && source instanceof SSPCluster[]) {
			return compareSSPPairTriggers((SSPPairTrigger) bankTrigger, (PairTrigger<SSPCluster[]>) simTrigger);
		}
		
		// Otherwise, the trigger combination is not supported. Produce
		// and exception.
		throw new IllegalArgumentException(String.format("Trigger type \"%s\" can not be compared to trigger type \"%s\" with source type \"%s\".",
				bankTrigger.getClass().getSimpleName(), simTrigger.getClass().getSimpleName(), source.getClass().getSimpleName()));
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
	
	private static final ClusterMatchEvent getBestPermutation(ClusterMatchEvent firstEvent, ClusterMatchEvent secondEvent) {
		// If both permutations are null, return that.
		if(firstEvent == null && secondEvent == null) {
			return null;
		}
		
		// If one permutation is null, it is not the best.
		if(firstEvent == null) { return secondEvent; }
		else if(secondEvent == null) { return firstEvent; }
		
		// A permutation is better if it has more matches.
		if(firstEvent.getMatches() > secondEvent.getMatches()) { return firstEvent; }
		else if(secondEvent.getMatches() > firstEvent.getMatches()) { return secondEvent; }
		
		// Otherwise, the permutation with the least energy failures is
		// the better permutation.
		if(firstEvent.getEnergyFailures() < secondEvent.getEnergyFailures()) { return firstEvent; }
		else if(secondEvent.getEnergyFailures() < firstEvent.getEnergyFailures()) { return secondEvent; }
		
		// If both these values are the same, then the events are identical.
		return firstEvent;
	}
	
	/**
	 * Determines the number of spaces needed to render the longest of
	 * a series of integers as a string.
	 * @param vals - The series of integers.
	 * @return Returns the number of spaces needed to render the longest
	 * integer as a base-10 string.
	 */
	private static final int getPrintSpaces(int... vals) {
		// Track the largest value.
		int largest = 0;
		
		// Iterate over the arguments and find the largest.
		for(int val : vals) {
			// Get the length of the string.
			int length = TriggerDiagnosticUtil.getDigits(val);
			
			// If it is larger, track it.
			if(length > largest) { largest = length; }
		}
		
		// Return the longer one.
		return largest;
	}
	
	/**
	 * Gets the position of the source of a <code>Trigger</code> object
	 * as text. This method only supports trigger sources of the types
	 * <code>SSPCluster</code>, <code>Cluster</code>, and arrays of size
	 * two of either type.
	 * @param trigger - The trigger from which to obtain the source.
	 * @return Returns the source of the trigger as a <code>String</code>
	 * object.
	 * @throws IllegalArgumentException Occurs if the source of the
	 * trigger is not any of the supported types.
	 */
	private static final String triggerPositionString(Trigger<?> trigger) throws IllegalArgumentException {
		// Get the trigger source.
		Object source = trigger.getTriggerSource();
		
		// Handle valid trigger sources.
		if(source instanceof SSPCluster) {
			return clusterPositionString((SSPCluster) source);
		} else if(source instanceof Cluster) {
			return clusterPositionString((Cluster) source);
		} else if(source instanceof SSPCluster[]) {
			SSPCluster[] sourcePair = (SSPCluster[]) source;
			if(sourcePair.length == 2) {
				return String.format("%s, %s", clusterPositionString(sourcePair[0]),
						clusterPositionString(sourcePair[1]));
			}
		} else if(source instanceof Cluster[]) {
			Cluster[] sourcePair = (Cluster[]) source;
			if(sourcePair.length == 2) {
				return String.format("%s, %s", clusterPositionString(sourcePair[0]),
						clusterPositionString(sourcePair[1]));
			}
		}
		
		// Otherwise, the source type is unrecognized. Throw an error.
		throw new IllegalArgumentException(String.format("Trigger source type \"%s\" is not supported.",
				trigger.getTriggerSource().getClass().getSimpleName()));
	}
	
	/**
	 * Gets the time of a simulated trigger object. Method supports
	 * triggers with source objects of type <code>SSPCluster</code>,
	 * <code>Cluster</code>, and arrays of size two composed of either
	 * object type.
	 * @param trigger - The trigger.
	 * @return Returns the time at which the trigger occurred.
	 * @throws IllegalArgumentException Occurs if the trigger source
	 * is not a supported type.
	 */
	private static final double getTriggerTime(Trigger<?> trigger) throws IllegalArgumentException {
		// Get the trigger source.
		Object source = trigger.getTriggerSource();
		
		// Get the trigger time for supported trigger types.
		if(source instanceof SSPCluster) {
			return ((SSPCluster) source).getTime();
		} else if(source instanceof Cluster) {
			return TriggerDiagnosticUtil.getClusterTime((Cluster) source);
		} else if(source instanceof SSPCluster[]) {
			// Get the pair.
			SSPCluster[] sourcePair = (SSPCluster[]) source;
			
			// Get the time of the bottom cluster.
			if(sourcePair.length == 2) {
				if(sourcePair[0].getYIndex() < 0) { return sourcePair[0].getTime(); }
				else if(sourcePair[1].getYIndex() < 0) { return sourcePair[1].getTime(); }
				else { throw new IllegalArgumentException("Cluster pairs must be formed of a top/bottom pair."); }
			}
			else { throw new IllegalArgumentException("Cluster pairs must be of size 2."); }
		} else if(source instanceof Cluster[]) {
			// Get the pair.
			Cluster[] sourcePair = (Cluster[]) source;
			int[] iy = {
				TriggerDiagnosticUtil.getYIndex(sourcePair[0]),
				TriggerDiagnosticUtil.getYIndex(sourcePair[1])
			};
			
			// Get the time of the bottom cluster.
			if(sourcePair.length == 2) {
				if(iy[0] < 0) { return TriggerDiagnosticUtil.getClusterTime(sourcePair[0]); }
				else if(iy[1] < 0) { return TriggerDiagnosticUtil.getClusterTime(sourcePair[1]); }
				else { throw new IllegalArgumentException("Cluster pairs must be formed of a top/bottom pair."); }
			}
			else { throw new IllegalArgumentException("Cluster pairs must be of size 2."); }
		}
		
		// If the source type is unrecognized, throw an exception.
		throw new IllegalArgumentException(String.format("Trigger source type \"%\" is not supported.",
				source.getClass().getSimpleName()));
	}
	
	private static final boolean[] triggerCutMatch(Trigger<?> simTrigger, SSPTrigger sspTrigger) {
		// Check that the cuts match for supported trigger types.
		if(simTrigger instanceof SinglesTrigger && sspTrigger instanceof SSPSinglesTrigger) {
			// Create an array to store the cut checks.
			boolean[] cutMatch = new boolean[3];
			
			// Cast the triggers.
			SinglesTrigger<?> simSingles = (SinglesTrigger<?>) simTrigger;
			SSPSinglesTrigger sspSingles = (SSPSinglesTrigger) sspTrigger;
			
			// Perform the check.
			cutMatch[ENERGY_MIN] = (simSingles.getStateClusterEnergyLow()  == sspSingles.passCutEnergyMin());
			cutMatch[ENERGY_MAX] = (simSingles.getStateClusterEnergyHigh() == sspSingles.passCutEnergyMax());
			cutMatch[HIT_COUNT] = (simSingles.getStateHitCount()          == sspSingles.passCutHitCount());
			
			// Return the match array.
			return cutMatch;
		} else if(simTrigger instanceof PairTrigger && sspTrigger instanceof SSPPairTrigger) {
			// Create an array to store the cut checks.
			boolean[] cutMatch = new boolean[4];
			
			// Cast the triggers.
			PairTrigger<?> simPair = (PairTrigger<?>) simTrigger;
			SSPPairTrigger sspPair = (SSPPairTrigger) sspTrigger;
			
			// Perform the check.
			cutMatch[ENERGY_SUM] = (simPair.getStateEnergySum()        == sspPair.passCutEnergySum());
			cutMatch[ENERGY_DIFF] = (simPair.getStateEnergyDifference() == sspPair.passCutEnergyDifference());
			cutMatch[ENERGY_SLOPE] = (simPair.getStateEnergySlope()      == sspPair.passCutEnergySlope());
			cutMatch[COPLANARITY] = (simPair.getStateCoplanarity()      == sspPair.passCutCoplanarity());
			
			// Return the match array.
			return cutMatch;
		}
		
		// If this point is reached, the triggers are not of a supported
		// type for cut comparison. Produce an exception.
		throw new IllegalArgumentException(String.format("Triggers of type \"%s\" can not be cut-matched with triggers of type \"%s\".",
				simTrigger.getClass().getSimpleName(), sspTrigger.getClass().getSimpleName()));
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