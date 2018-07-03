package org.hps.users.mgraham;

import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IPlotter;
import hep.aida.IPlotterStyle;
import hep.physics.matrix.MutableMatrix;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.hps.analysis.dataquality.DataQualityMonitor;

import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackerHitUtils;
import org.hps.recon.utils.TrackClusterMatcher;
import org.lcsim.constants.Constants;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.converter.compact.subdetector.HpsTracker2;
import org.lcsim.detector.converter.compact.subdetector.SvtStereoLayer;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.solids.Point3D;
import org.lcsim.detector.solids.Polygon3D;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.fit.helicaltrack.HelicalTrackCross;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.fit.helicaltrack.HelicalTrackStrip;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.geometry.subdetector.BarrelEndcapFlag;
import org.lcsim.lcio.LCIOWriter;

/**
 *
 */
public class TrackDebugger extends DataQualityMonitor {

    private static Logger LOGGER = Logger.getLogger(TrackDebugger.class.getPackage().getName());

    String finalStateParticlesColName = "FinalStateParticles";
    String readoutHitCollectionName = "EcalReadoutHits";//these are in ADC counts
    String calibratedHitCollectionName = "EcalCalHits";//these are in energy
    String clusterCollectionName = "EcalClustersCorr";
    private String notrackFile;
    private String helicalTrackHitCollectionName = "HelicalTrackHits";
    private String rotatedTrackHitCollectionName = "RotatedHelicalTrackHits";
    String[] fpQuantNames = {"nEle_per_Event", "nPos_per_Event", "nPhoton_per_Event", "nUnAssociatedTracks_per_Event", "avg_delX_at_ECal", "avg_delY_at_ECal", "avg_E_Over_P", "avg_mom_beam_elec", "sig_mom_beam_elec"};
    private String outputFile;
    private LCIOWriter writer;
    private LCIOWriter notrackwriter;
    private TrackClusterMatcher matcher = new TrackClusterMatcher();
    //private Map<SiSensor, Map<Integer, Hep3Vector>> stripPositions = new HashMap<SiSensor, Map<Integer, Hep3Vector>>(); 
    private List<HpsSiSensor> sensors = null;
    private Map<Integer, List<SvtStereoLayer>> topStereoLayers = new HashMap<Integer, List<SvtStereoLayer>>();
    private Map<Integer, List<SvtStereoLayer>> bottomStereoLayers = new HashMap<Integer, List<SvtStereoLayer>>();
    // Constants
    public static final double SENSOR_LENGTH = 98.33; // mm
    public static final double SENSOR_WIDTH = 38.3399; // mm
    private static final String SUBDETECTOR_NAME = "Tracker";
    boolean doSkim = false;
    TrackerHitUtils trackerHitUtils = new TrackerHitUtils();
    //some counters
    int nRecoEvents = 0;
    int nTotEle = 0;
    int nTotPos = 0;
    int nTotPhotons = 0;
    int nTotUnAss = 0;
    int nTotAss = 0;
    //some summers
    double sumdelX = 0.0;
    double sumdelY = 0.0;
    double sumEoverP = 0.0;
    private final String plotDir = "FinalStateParticles/";
    // double beamEnergy = 1.05; //GeV
    double clTimeMin = 30;//ns
    double clTimeMax = 50;//ns
    double deltaTimeMax = 4;//ns
    double coplanMean = 180;//degrees
    double coplanWidth = 10;//degrees
    double esumMin = 0.4;
    double esumMax = 1.2;
    double phot_nom_x = 42.52;//nominal photon position (px=0)
    boolean requirePositron = false;
    double maxPairs = 5;
    boolean requireSuperFiducial = false;
    int nbins = 50;
    double B_FIELD = 0.23;//Tesla

    double minPhi = -0.25;
    double maxPhi = 0.25;

    IHistogram1D elePx;
    IHistogram1D elePy;
    IHistogram1D elePz;
    IHistogram1D elePzBeam;
    IHistogram1D elePzBeamTop;
    IHistogram1D elePzBeamBottom;
    IHistogram1D elePTop;
    IHistogram1D elePBottom;
    IHistogram1D eleClEne;

    IHistogram1D posPx;
    IHistogram1D posPy;
    IHistogram1D posPz;
    IHistogram1D posPTop;
    IHistogram1D posPBottom;
    IHistogram1D posClEne;

    /*  photon quanties (...right now, just unassociated clusters) */
    IHistogram1D nPhotonsHisto;
    IHistogram1D enePhoton;
    IHistogram1D xPhoton;
    IHistogram1D yPhoton;

    IHistogram1D phoClEnePosSide;
    IHistogram1D phoClEneEleSide;
    IHistogram1D clEnergySumAll;
    IHistogram1D clEnergySumNoEle;
    IHistogram1D clEnergySumNoPos;
    IHistogram1D clEnergySumNoTracks;
    IHistogram1D clEnergySumBothTracks;

    /*  tracks with associated clusters */
    IHistogram1D eneOverp;
    IHistogram1D deltaXAtCal;
    IHistogram1D deltaYAtCal;
//    IHistogram2D trackXvsECalX;
//    IHistogram2D trackYvsECalY;
    IHistogram2D trackPvsECalE;
    IHistogram2D trackTvsECalT;
    IHistogram1D timeMatchDeltaT;
    /* number of unassocaited tracks/event */
    IHistogram1D nUnAssTracksHisto;
    IHistogram2D clXvsYEle;
    IHistogram2D clXvsYNoEle;
    IHistogram2D clXvsYPos;
    IHistogram2D clXvsYNoPos;
    IHistogram2D elePhiTrackVsCluster;
    IHistogram2D posPhiTrackVsCluster;
    IHistogram1D eleSlope;
    IHistogram1D elePhi;
    IHistogram1D posSlope;
    IHistogram1D posPhi;

    IHistogram1D fromClusterEleSlopeFound;
    IHistogram1D fromClusterElePhiFound;
    IHistogram1D fromClusterEleSlopeMissed;
    IHistogram1D fromClusterElePhiMissed;

    IHistogram1D fromClusterPosSlopeFound;
    IHistogram1D fromClusterPosPhiFound;
    IHistogram1D fromClusterPosSlopeMissed;
    IHistogram1D fromClusterPosPhiMissed;

    IHistogram1D fromClusterEleClusterDX;
    IHistogram1D fromClusterEleClusterDY;

    IHistogram1D eleMissedNModulesCrossed;
    IHistogram1D posMissedNModulesCrossed;
    IHistogram1D eleMissedHitsFound;
    IHistogram1D posMissedHitsFound;
    int nmodules = 6;
    IHistogram1D[] eleMissedClosestResidX = new IHistogram1D[nmodules];
    ;
    IHistogram1D[] eleMissedClosestResidY = new IHistogram1D[nmodules];
    ;
    IHistogram1D[] posMissedClosestResidX = new IHistogram1D[nmodules];
    ;
    IHistogram1D[] posMissedClosestResidY = new IHistogram1D[nmodules];
    ;
    IHistogram1D eleMissedTotalDX;
    IHistogram1D eleMissedTotalDY;
    IHistogram1D posMissedTotalDX;
    IHistogram1D posMissedTotalDY;
    IHistogram2D siClusterTimevsZ;

    String detName;

    /**
     * The B field map
     */
    FieldMap bFieldMap = null;
    /**
     * Position of the Ecal face
     */
    private double ecalPosition = 0; // mm

    /**
     * Z position to start extrapolation from
     */
    double extStartPos = 700; // mm

    /**
     * The extrapolation step size
     */
    double stepSize = 5.0; // mm
    /**
     * Name of the constant denoting the position of the Ecal face in the
     * compact description.
     */
    private static final String ECAL_POSITION_CONSTANT_NAME = "ecal_dface";

    public void setFinalStateParticlesColName(String fsp) {
        this.finalStateParticlesColName = fsp;
    }

    public void setOutputFilePath(String output) {
        this.outputFile = output;
    }

    public void setNoTracksFilePath(String output) {
        this.notrackFile = output;
    }

    public void setDoSkim(boolean doit) {
        this.doSkim = doit;
    }

    protected void detectorChanged(Detector detector) {

        super.detectorChanged(detector);
        detName = detector.getName();
        double maxFactor = 1.5;
        double feeMomentumCut = 0.75; //this number, multiplied by the beam energy, is the actual cut
        B_FIELD = detector.getFieldMap().getField(new BasicHep3Vector(0, 0, 500)).y();
        // Get the field map from the detector object
        bFieldMap = detector.getFieldMap();
        // Get the position of the Ecal from the compact description
        ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
        // Get the position of the Ecal from the compact description
        ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
        LOGGER.setLevel(Level.ALL);
        LOGGER.info("B_FIELD=" + B_FIELD);
        LOGGER.info("Setting up the plotter");
        aida.tree().cd("/");
        String trkType = "SeedTrack/";
        if (isGBL)
            trkType = "GBLTrack/";
        // Get the HpsSiSensor objects from the tracker detector element
        sensors = detector.getSubdetector(SUBDETECTOR_NAME)
                .getDetectorElement().findDescendants(HpsSiSensor.class);

        // If the detector element had no sensors associated with it, throw
        // an exception
        if (sensors.size() == 0)
            throw new RuntimeException("No sensors were found in this detector.");

        // Get the stereo layers from the geometry and build the stereo
        // layer maps
        List<SvtStereoLayer> stereoLayers
                = ((HpsTracker2) detector.getSubdetector(SUBDETECTOR_NAME).getDetectorElement()).getStereoPairs();
        for (SvtStereoLayer stereoLayer : stereoLayers)
            if (stereoLayer.getAxialSensor().isTopLayer()) {
                //System.out.println("Adding stereo layer " + stereoLayer.getLayerNumber());
                if (!topStereoLayers.containsKey(stereoLayer.getLayerNumber()))
                    topStereoLayers.put(stereoLayer.getLayerNumber(), new ArrayList<SvtStereoLayer>());
                topStereoLayers.get(stereoLayer.getLayerNumber()).add(stereoLayer);
            } else {
                if (!bottomStereoLayers.containsKey(stereoLayer.getLayerNumber()))
                    bottomStereoLayers.put(stereoLayer.getLayerNumber(), new ArrayList<SvtStereoLayer>());
                bottomStereoLayers.get(stereoLayer.getLayerNumber()).add(stereoLayer);
            }
        /*  Final State Particle Quantities   */
        /*  plot electron & positron momentum separately  */
        elePx = aida.histogram1D("Electron Px (GeV)", nbins, -0.1 * beamEnergy, 0.200 * beamEnergy);
        elePy = aida.histogram1D("Electron Py (GeV)", nbins, -0.1 * beamEnergy, 0.1 * beamEnergy);
        elePz = aida.histogram1D("Electron Pz (GeV)", nbins, 0, beamEnergy * maxFactor);
//        elePTop = aida.histogram1D( "Electron Total P (GeV):  Top", nbins, 0, beamEnergy * maxFactor);
//        elePBottom = aida.histogram1D( "Electron Total P (GeV):  Bottom", nbins, 0, beamEnergy * maxFactor);
        eleClEne = aida.histogram1D("Electron Cluster Energy", nbins, 0, beamEnergy * maxFactor);

        posPx = aida.histogram1D("Positron Px (GeV)", nbins, -0.1 * beamEnergy, 0.200 * beamEnergy);
        posPy = aida.histogram1D("Positron Py (GeV)", nbins, -0.1 * beamEnergy, 0.1 * beamEnergy);
        posPz = aida.histogram1D("Positron Pz (GeV)", nbins, 0, beamEnergy * maxFactor);
//        posPTop = aida.histogram1D( "Positron Total P (GeV):  Top", nbins, 0, beamEnergy * maxFactor);
//        posPBottom = aida.histogram1D( "Positron Total P (GeV):  Bottom", nbins, 0, beamEnergy * maxFactor);
        posClEne = aida.histogram1D("Positron Cluster Energy", nbins, 0, beamEnergy * maxFactor);
        /*  No track quantities
        
         /*  photon quanties (...right now, just unassociated clusters) */
//        nPhotonsHisto = aida.histogram1D( "Number of photons per event", 15, 0, 15);
//        enePhoton = aida.histogram1D( "Photon Energy (GeV)", nbins, 0, 2.4 * beamEnergy);
//        xPhoton = aida.histogram1D( "Photon X position (mm)", nbins, -200, 200);
//        yPhoton = aida.histogram1D( "Photon Y position (mm)", nbins, -100, 100);
        phoClEneEleSide = aida.histogram1D("Photon Cluster Energy: Electron Side", nbins, 0, beamEnergy * maxFactor);
        phoClEnePosSide = aida.histogram1D("Photon Cluster Energy: Positron Side", nbins, 0, beamEnergy * maxFactor);
        clEnergySumAll = aida.histogram1D("Cluster Energy Sum: All", nbins, 0, beamEnergy * maxFactor);
        clEnergySumNoPos = aida.histogram1D("Cluster Energy Sum: No Positron", nbins, 0, beamEnergy * maxFactor);
        clEnergySumNoEle = aida.histogram1D("Cluster Energy Sum: No Electron", nbins, 0, beamEnergy * maxFactor);
        clEnergySumNoTracks = aida.histogram1D("Cluster Energy Sum: No Tracks", nbins, 0, beamEnergy * maxFactor);
        clEnergySumBothTracks = aida.histogram1D("Cluster Energy Sum: Both Tracks", nbins, 0, beamEnergy * maxFactor);

        fromClusterEleSlopeFound = aida.histogram1D("From Cluster Found Electron Slope", nbins, -0.1, 0.1);
        fromClusterElePhiFound = aida.histogram1D("From Cluster Found Electron Phi", nbins, minPhi, maxPhi);

        fromClusterEleSlopeMissed = aida.histogram1D("From Cluster Missed Electron Slope", nbins, -0.1, 0.1);
        fromClusterElePhiMissed = aida.histogram1D("From Cluster Missed Electron Phi", nbins, minPhi, maxPhi);

        fromClusterPosSlopeFound = aida.histogram1D("From Cluster Found Positron Slope", nbins, -0.1, 0.1);
        fromClusterPosPhiFound = aida.histogram1D("From Cluster Found Positron Phi", nbins, minPhi, maxPhi);

        fromClusterPosSlopeMissed = aida.histogram1D("From Cluster Missed Positron Slope", nbins, -0.1, 0.1);
        fromClusterPosPhiMissed = aida.histogram1D("From Cluster Missed Positron Phi", nbins, minPhi, maxPhi);

        eleSlope = aida.histogram1D(" Electron Slope", nbins, -0.1, 0.1);
        elePhi = aida.histogram1D(" Electron Phi", nbins, minPhi, maxPhi);

        posSlope = aida.histogram1D(" Positron Slope", nbins, -0.1, 0.1);
        posPhi = aida.histogram1D(" Positron Phi", nbins, minPhi, maxPhi);

        elePhiTrackVsCluster = aida.histogram2D("Electron Phi From Track vs Cluster", nbins, minPhi, maxPhi, nbins, minPhi, maxPhi);
        posPhiTrackVsCluster = aida.histogram2D("Positron Phi From Track vs Cluster", nbins, minPhi, maxPhi, nbins, minPhi, maxPhi);

        eleMissedNModulesCrossed = aida.histogram1D(" N Modules Crossed: Missed Electron", 7, 0, 7);
        posMissedNModulesCrossed = aida.histogram1D(" N Modules Crossed: Missed Positron", 7, 0, 7);
        eleMissedHitsFound = aida.histogram1D(" N Modules Hits Found: Missed Electron", 7, 0, 7);
        posMissedHitsFound = aida.histogram1D(" N Modules Hits Found: Missed Positron", 7, 0, 7);
        for (int i = 0; i < 6; i++) {
            eleMissedClosestResidX[i] = aida.histogram1D(plotDir + trkType  + "/Residuals/MissedElectron" + "X Resid Layer " + (i + 1), 50, -50, 50);
            eleMissedClosestResidY[i] = aida.histogram1D(plotDir + trkType  + "/Residuals/MissedElectron" + "Y Resid Layer " + (i + 1), 50, -50, 50);
            posMissedClosestResidX[i] = aida.histogram1D(plotDir + trkType  + "/Residuals/MissedPositron" + "X Resid Layer " + (i + 1), 50, -50, 50);
            posMissedClosestResidY[i] = aida.histogram1D(plotDir + trkType  + "/Residuals/MissedPositron" + "Y Resid Layer " + (i + 1), 50, -50, 50);

        }
        eleMissedTotalDX = aida.histogram1D(plotDir + trkType  + "/Residuals/MissedElectron" + "X Resid Total", 50, 0, 100);
        eleMissedTotalDY = aida.histogram1D(plotDir + trkType + "/Residuals/MissedElectron" + "Y Resid Total", 50, 0, 100);
        posMissedTotalDX = aida.histogram1D(plotDir + trkType + "/Residuals/MissedPositron" + "X Resid Total", 50, 0, 100);
        posMissedTotalDY = aida.histogram1D(plotDir + trkType + "/Residuals/MissedPositron" + "Y Resid Total", 50, 0, 100);

        fromClusterEleClusterDX = aida.histogram1D("Fake Electron Cluster DX", nbins, -20, 20);
        fromClusterEleClusterDY = aida.histogram1D("Fake Electron Cluster DY", nbins, -2, 2);
        /*  tracks with associated clusters */
//        eneOverp = aida.histogram1D( "Cluster Energy Over TrackMomentum", nbins, 0, 2.0);
//        deltaXAtCal = aida.histogram1D( "delta X @ ECal (mm)", nbins, -50, 50.0);
//        deltaYAtCal = aida.histogram1D( "delta Y @ ECal (mm)", nbins, -50, 50.0);
//        trackXvsECalX = aida.histogram2D(plotDir +trkType+ triggerType + "/" + "track X vs ECal X", nbins, -300, 300.0, nbins, -300, 300.0);
//        trackYvsECalY = aida.histogram2D(plotDir +trkType+ triggerType + "/" + "track Y vs ECal Y", nbins, -100, 100.0, nbins, -100, 100.0);
//        trackPvsECalE = aida.histogram2D( "track mom vs ECal E", nbins, 0.1, beamEnergy * maxFactor, nbins, 0.1, beamEnergy * maxFactor);
//        trackTvsECalT = aida.histogram2D( "track T vs ECal T", 200, 0.0, 200.0, nbins, -25.0, 25.0);
//        timeMatchDeltaT = aida.histogram1D( "ECal T minus track T", 200, -25, 175);
        /* number of unassocaited tracks/event */
//        nUnAssTracksHisto = aida.histogram1D( "Number of unassociated tracks per event", 5, 0, 5);
        clXvsYEle = aida.histogram2D("Cluster X vs Y: Found Electron", nbins, -300, 100, nbins, -100, 100.0);
        clXvsYNoEle = aida.histogram2D("Cluster X vs Y: Missed Electron", nbins, -300, 100, nbins, -100, 100.0);
        clXvsYPos = aida.histogram2D("Cluster X vs Y: Found Positron", nbins, -100, 300, nbins, -100, 100.0);
        clXvsYNoPos = aida.histogram2D("Cluster X vs Y: Missed Positron", nbins, -100, 300, nbins, -100, 100.0);

        siClusterTimevsZ = aida.histogram2D("SiCluster time vs Z", 100, 0, 1000.0, nbins, -20, 20);

    }

    @Override
    public void process(EventHeader event) {
        /*  make sure everything is there */

//        List<CalorimeterHit> hits;
//        if (event.hasCollection(CalorimeterHit.class, calibratedHitCollectionName))
//            hits = event.get(CalorimeterHit.class, calibratedHitCollectionName);
//        else
//            return; //this might be a non-data event        
        if (!event.hasCollection(ReconstructedParticle.class, finalStateParticlesColName)) {
            if (debug)
                LOGGER.info(finalStateParticlesColName + " collection not found???");
            return;
        }
        if (!event.hasCollection(TrackerHit.class, rotatedTrackHitCollectionName)) {
            LOGGER.info("No " + rotatedTrackHitCollectionName + "????");
            return;
        }

        if (!event.hasCollection(Cluster.class, clusterCollectionName)) {
            LOGGER.info("No " + clusterCollectionName + "????");
            return;
        }

        List<TrackerHit> hthList = event.get(TrackerHit.class, rotatedTrackHitCollectionName);
        List<Cluster> calHits = event.get(Cluster.class, clusterCollectionName);
//        RelationalTable hitToStrips = TrackUtils.getHitToStripsTable(event);
//        RelationalTable hitToRotated = TrackUtils.getHitToRotatedTable(event);

        //check to see if this event is from the correct trigger (or "all");
        if (!matchTrigger(event))
            return;

        nRecoEvents++;
        int nPhotons = 0;  //number of photons 
        int nUnAssTracks = 0; //number of tracks w/o clusters
        List<List> pairList = new ArrayList<>();
        List<ReconstructedParticle> finalStateParticles = event.get(ReconstructedParticle.class, finalStateParticlesColName);
        if (debug)
            LOGGER.info("This events has " + finalStateParticles.size() + " final state particles");
//        for (ReconstructedParticle fsPart1 : finalStateParticles) {
        for (Cluster cl1 : calHits) {
//            if (debug)
            //               LOGGER.info("PDGID = " + fsPart1.getParticleIDUsed() + "; charge = " + fsPart1.getCharge() + "; pz = " + fsPart1.getMomentum().x());
//            if (isGBL != TrackType.isGBL(fsPart1.getType()))
//                continue;
//            if (fsPart1.getClusters().size() == 0)
//                continue;
//            Cluster cl1 = fsPart1.getClusters().get(0);
            double[] cl1Position = cl1.getPosition();
            double cl_xi = cl1Position[0];
            double cl_yi = cl1Position[1];
            double cl_zi = cl1Position[2];
            double cl_ti;
            cl_ti = ClusterUtilities.getSeedHitTime(cl1);
            double cl_Ei = cl1.getEnergy();

            if (requireSuperFiducial && !inSuperFiducialRegion(cl_xi, cl_yi))
                continue;
            if (cl_ti < clTimeMin || cl_ti > clTimeMax)
                continue;

            //          for (ReconstructedParticle fsPart2 : finalStateParticles.subList(finalStateParticles.indexOf(fsPart1), finalStateParticles.size())) {
            for (Cluster cl2 : calHits.subList(calHits.indexOf(cl1), calHits.size())) {
//                if (fsPart2.getClusters().size() == 0)
//                    continue;
//                if (fsPart1 == fsPart2)
//                    continue;
//                if (debug)
//                    LOGGER.info("PDGID = " + fsPart2.getParticleIDUsed() + "; charge = " + fsPart2.getCharge() + "; pz = " + fsPart1.getMomentum().x());
//                if (isGBL != TrackType.isGBL(fsPart2.getType()))
//                    continue;
//                Cluster cl2 = fsPart2.getClusters().get(0);
                double[] cl2Position = cl2.getPosition();
                double cl_xj = cl2Position[0];
                double cl_yj = cl2Position[1];
                double cl_zj = cl2Position[2];
                double cl_tj;
                cl_tj = ClusterUtilities.getSeedHitTime(cl2);
                double cl_Ej = cl2.getEnergy();
                if (requireSuperFiducial && !inSuperFiducialRegion(cl_xj, cl_yj))
                    continue;
                //timing cuts
                if (cl_tj < clTimeMin || cl_tj > clTimeMax)
                    continue;
                double deltaT = cl_tj - cl_ti;
                //time difference cut
                if (Math.abs(deltaT) > deltaTimeMax)
                    continue;
                //require top-bottom
                if (cl_yi * cl_yj > 0)
                    continue;
                //require left-right
                if (cl_xi * cl_xj > 0)
                    continue;

                double coplan = getECalCoplanarity(cl1, cl2);
                if (Math.abs(coplan - coplanMean) > coplanWidth)
                    continue;

                double esum = cl_Ei + cl_Ej;
                if (esum > esumMax || esum < esumMin)
                    continue;
//                //require a matched positron track!
//                if (requirePositron)
//                    if (fsPart1.getCharge() != 1 && fsPart2.getCharge() != 1)
//                        continue;

                /////////  passed all cuts....make a new pair!
                List<Cluster> pair = new ArrayList<Cluster>();
                //add electron first (.get(0)
                if (cl_xi < 0) {
                    pair.add(cl1);
                    pair.add(cl2);
                } else {
                    pair.add(cl2);
                    pair.add(cl1);
                }
                pairList.add(pair);
            }
        }

        boolean saveEvent = false;
        boolean saveNoTrack = false;
        //LOGGER.info("Found " + pairList.size() + " pairs in this event");
//        if (pairList.size() > 0)
//            LOGGER.info("Found " + pairList.size() + " pairs in this event");
        if (pairList.size() > maxPairs)
            return;
        if (doSkim && pairList.size() > 0)
            saveEvent = true;

        for (List pair : pairList) {
            ReconstructedParticle fsEleSide = (ReconstructedParticle) findFSParticleFromCluster(finalStateParticles, (Cluster) pair.get(0));
            ReconstructedParticle fsPosSide = (ReconstructedParticle) findFSParticleFromCluster(finalStateParticles, (Cluster) pair.get(1));
            double esum = fsEleSide.getEnergy() + fsPosSide.getEnergy();
            Cluster eleSide = fsEleSide.getClusters().get(0);
            Cluster posSide = fsPosSide.getClusters().get(0);
            boolean posTrack = false;
            boolean eleTrack = false;
            boolean posTrackIsPos = false;
            boolean eleTrackIsEle = false;
            clEnergySumAll.fill(esum);
            if (fsPosSide.getTracks().size() > 0) {
                posTrack = true;
                if (fsPosSide.getTracks().get(0).getCharge() < 0) //remember the sign is flipped!!!
                    posTrackIsPos = true;
//                LOGGER.info("Found Positron Tracks! ");
                Hep3Vector posMom = fsPosSide.getMomentum();
                posPx.fill(posMom.x());
                posPy.fill(posMom.y());
                posPz.fill(posMom.z());
                posClEne.fill(fsPosSide.getEnergy());
                posPhi.fill(fsPosSide.getTracks().get(0).getTrackStates().get(0).getPhi());
                posSlope.fill(fsPosSide.getTracks().get(0).getTrackStates().get(0).getTanLambda());
                clXvsYPos.fill(posSide.getPosition()[0], posSide.getPosition()[1]);

            } else {
                clXvsYNoPos.fill(posSide.getPosition()[0], posSide.getPosition()[1]);
                phoClEnePosSide.fill(fsPosSide.getEnergy());
            }
            if (fsEleSide.getTracks().size() > 0) {
                eleTrack = true;
                if (fsEleSide.getTracks().get(0).getCharge() > 0)//remember the sign is flipped!!!
                    eleTrackIsEle = true;
//                LOGGER.info("Found Electron Tracks! ");
                Hep3Vector eleMom = fsPosSide.getMomentum();
                elePx.fill(eleMom.x());
                elePy.fill(eleMom.y());
                elePz.fill(eleMom.z());
                eleClEne.fill(fsEleSide.getEnergy());
                clXvsYEle.fill(eleSide.getPosition()[0], eleSide.getPosition()[1]);
                elePhi.fill(fsEleSide.getTracks().get(0).getTrackStates().get(0).getPhi());
                eleSlope.fill(fsEleSide.getTracks().get(0).getTrackStates().get(0).getTanLambda());
            } else {
                phoClEneEleSide.fill(fsEleSide.getEnergy());
                clXvsYNoEle.fill(eleSide.getPosition()[0], eleSide.getPosition()[1]);
            }
            if (posTrack && eleTrack)
                clEnergySumBothTracks.fill(esum);
            else if (!posTrack && !eleTrack) {
                clEnergySumNoTracks.fill(esum);
                saveNoTrack = true;
            } else if ((!eleTrack && posTrackIsPos) || (!posTrack && !eleTrackIsEle))
                clEnergySumNoEle.fill(esum);
            else
                clEnergySumNoPos.fill(esum);
            //require positron and no electron...look for electron-side SVT hits

            //the slope of the track in the non-bend assuming starts at (0,0,0)
            double slpFromCluster = eleSide.getPosition()[1] / eleSide.getPosition()[2];
            //the curvature of the track using the measured cluster energy
            double curvFromCluster = TrackUtils.calculateCurvature(eleSide.getEnergy(), -1, Constants.fieldConversion * B_FIELD);
//            LOGGER.info("Curvature from Cluster = " + curvFromCluster);
            //make a fromCluster track with these 
            double phi = findBestPhi(eleSide, curvFromCluster, slpFromCluster, minPhi, maxPhi, 100);
//            LOGGER.info("Phi from Cluster = " + phi);
            double[] trackParameters = new double[5];
            Track trkFromEleCluster = new BaseTrack();
            trackParameters[BaseTrack.D0] = 0;
            trackParameters[BaseTrack.OMEGA] = curvFromCluster;
            trackParameters[BaseTrack.PHI] = phi;
            trackParameters[BaseTrack.TANLAMBDA] = slpFromCluster;
            trackParameters[BaseTrack.Z0] = 0;
            ((BaseTrack) trkFromEleCluster).setTrackParameters(trackParameters, B_FIELD);
            HelicalTrackFit htfFromEleCluster = TrackUtils.getHTF(trkFromEleCluster);
            TrackState state = TrackUtils.extrapolateTrackUsingFieldMap(trkFromEleCluster, extStartPos, ecalPosition, stepSize, bFieldMap);
            trkFromEleCluster.getTrackStates().add(state);
            if (!eleTrack && posTrackIsPos) {
                fromClusterEleSlopeMissed.fill(slpFromCluster);
                fromClusterElePhiMissed.fill(phi);
//                fromClusterEleClusterDX.fill(matcher.getXDistance(eleSide, trkFromEleCluster));
//                fromClusterEleClusterDY.fill(matcher.getYDistance(eleSide, trkFromEleCluster));
//                LOGGER.info("Delta X = " + matcher.getXDistance(eleSide, trkFromEleCluster) + "; Delta Y = " + matcher.getYDistance(eleSide, trkFromEleCluster));
            } else if (eleTrackIsEle) {
                fromClusterEleSlopeFound.fill(slpFromCluster);
                fromClusterElePhiFound.fill(phi);
                elePhiTrackVsCluster.fill(fsEleSide.getTracks().get(0).getTrackStates().get(0).getPhi(), phi);
                List<TrackerHit> hits = fsEleSide.getTracks().get(0).getTrackerHits();
                for (TrackerHit hit : hits) {
                    HelicalTrackHit hth = (HelicalTrackHit) hit;
                    HelicalTrackCross cross = (HelicalTrackCross) hth;
                    List<HelicalTrackStrip> clusterlist = cross.getStrips();
                    for (HelicalTrackStrip cl : clusterlist)
                        siClusterTimevsZ.fill(cl.origin().x(), cl.time());
                }
            }

            //the slope of the track in the non-bend assuming starts at (0,0,0)
            slpFromCluster = posSide.getPosition()[1] / posSide.getPosition()[2];
            //the curvature of the track using the measured cluster energy
            curvFromCluster = TrackUtils.calculateCurvature(posSide.getEnergy(), 1, Constants.fieldConversion * B_FIELD);
//            LOGGER.info("Curvature from Cluster = " + curvFromCluster);
            //make a fromCluster track with these 
            phi = findBestPhi(posSide, curvFromCluster, slpFromCluster, minPhi, maxPhi, 100);
//            LOGGER.info("Phi from Cluster = " + phi);
            trackParameters = new double[5];
            Track trkFromPosCluster = new BaseTrack();
            trackParameters[BaseTrack.D0] = 0;
            trackParameters[BaseTrack.OMEGA] = curvFromCluster;
            trackParameters[BaseTrack.PHI] = phi;
            trackParameters[BaseTrack.TANLAMBDA] = slpFromCluster;
            trackParameters[BaseTrack.Z0] = 0;
            ((BaseTrack) trkFromPosCluster).setTrackParameters(trackParameters, B_FIELD);
            HelicalTrackFit htfFromPosCluster = TrackUtils.getHTF(trkFromPosCluster);
            state = TrackUtils.extrapolateTrackUsingFieldMap(trkFromPosCluster, extStartPos, ecalPosition, stepSize, bFieldMap);
            trkFromPosCluster.getTrackStates().add(state);
            if (!posTrack && eleTrackIsEle) {
                fromClusterPosSlopeMissed.fill(slpFromCluster);
                fromClusterPosPhiMissed.fill(phi);
                //   fromClusterPosClusterDX.fill(matcher.getXDistance(posSide, trkFromPosCluster));
                //   fromClusterPosClusterDY.fill(matcher.getYDistance(posSide, trkFromPosCluster));
//                LOGGER.info("Delta X = " + matcher.getXDistance(posSide, trkFromPosCluster) + "; Delta Y = " + matcher.getYDistance(posSide, trkFromPosCluster));
            } else if (posTrackIsPos) {
                fromClusterPosSlopeFound.fill(slpFromCluster);
                fromClusterPosPhiFound.fill(phi);
                posPhiTrackVsCluster.fill(fsPosSide.getTracks().get(0).getTrackStates().get(0).getPhi(), phi);
                List<TrackerHit> hits = fsPosSide.getTracks().get(0).getTrackerHits();
                for (TrackerHit hit : hits) {
                    HelicalTrackHit hth = (HelicalTrackHit) hit;
                    HelicalTrackCross cross = (HelicalTrackCross) hth;
                    List<HelicalTrackStrip> clusterlist = cross.getStrips();
                    for (HelicalTrackStrip cl : clusterlist)
                        siClusterTimevsZ.fill(cl.origin().x(), cl.time());
                }
            }

            if (posTrack && !eleTrack && posTrackIsPos) {
                LOGGER.info("Found a suspicious event...WILL SAVE");
                saveEvent = true;
                //see if the trajectory taken from cluster intercepts the sensors:   
                int nSensorsCrossed = 0;
                double totalDX = 0;
                double totalDY = 0;
                List<HelicalTrackHit> newTrackHits = new ArrayList<HelicalTrackHit>();
                for (int l = 1; l < 7; l++) {
                    boolean inAcceptance = isWithinAcceptance(trkFromEleCluster, l);
                    if (inAcceptance)
                        nSensorsCrossed++;
                    double minDelta = 99999999;
                    double minX = 99999;
                    double minY = 99999;
                    HelicalTrackHit bestHit = null;
                    LOGGER.info("Is the electron trajectory within the acceptance of Layer" + l + "?  " + inAcceptance);
                    for (TrackerHit hthTH : hthList)
                        if (getLayerFromPosition(hthTH.getPosition()[0]) == l) {
                            Hep3Vector pos = new BasicHep3Vector(hthTH.getPosition());
                            SymmetricMatrix cvm = new SymmetricMatrix(3, hthTH.getCovMatrix(), true);
                            HelicalTrackHit hth;
                            hth = new HelicalTrackHit(pos, cvm, hthTH.getdEdx(), hthTH.getTime(), hthTH.getType(), null, detName, l, BarrelEndcapFlag.BARREL);
                            Map<String, Double> resid = TrackUtils.calculateTrackHitResidual(hth, htfFromEleCluster, false);
                            double resX = resid.get("resy") / resid.get("erry");//X in detector == Y in tracking
                            double resY = resid.get("resz") / resid.get("errz");//Y in detector == Z in tracking
                            LOGGER.info("hit resiudal =  " + resid);
                            double delta = Math.sqrt(resX * resX + resY * resY);
                            if (delta < minDelta) {
                                bestHit = hth;
                                minDelta = delta;
                                minX = resid.get("resy");
                                minY = resid.get("resz");
                            }
                        }
                    if (bestHit != null && minDelta < 10)
                        newTrackHits.add(bestHit);
                    totalDX += minX * minX;
                    totalDY += minY * minY;
                    eleMissedClosestResidX[l - 1].fill(minX);
                    eleMissedClosestResidY[l - 1].fill(minY);
                }
                eleMissedTotalDX.fill(Math.sqrt(totalDX));
                eleMissedTotalDY.fill(Math.sqrt(totalDY));
                eleMissedNModulesCrossed.fill(nSensorsCrossed);
                eleMissedHitsFound.fill(newTrackHits.size());
            } else if (posTrack && !eleTrack)
                LOGGER.info("Won't save this...have an electron track on positron side");
        }
        //if I got here I want to save the event.
        if (saveEvent)
            try {
                LOGGER.info("Writing Event!");
                writer.write(event);
            } catch (IOException x) {
                throw new RuntimeException("Error writing LCIO file", x);
            }
        if (saveNoTrack)
            try {
                LOGGER.info("Writing No Track Event!");
                notrackwriter.write(event);
            } catch (IOException x) {
                throw new RuntimeException("Error writing LCIO file", x);
            }
//        LOGGER.info("Done with event");

    }
//MG -- 6/26/17...this doens't work to get just best phi anymore...kludge dist definitly so it doesn't break build, but it get the best dist in X and Y
    public double findBestPhi(Cluster cl, double curv, double slp, double minPhi, double maxPhi, double nSteps) {
        double bestPhi = 999999999;
        double minDist = 9999999;
        double phiStep = (maxPhi - minPhi) / nSteps;
        for (int i = 0; i < nSteps; i++) {
            double thisPhi = phiStep * i + minPhi;;
            double[] trackParameters = new double[5];
            BaseTrack trkFromCluster = new BaseTrack();
            trackParameters[BaseTrack.D0] = 0;
            trackParameters[BaseTrack.OMEGA] = curv;
            trackParameters[BaseTrack.PHI] = thisPhi;
            trackParameters[BaseTrack.TANLAMBDA] = slp;
            trackParameters[BaseTrack.Z0] = 0;
            ((BaseTrack) trkFromCluster).setTrackParameters(trackParameters, B_FIELD);
            TrackState state = TrackUtils.extrapolateTrackUsingFieldMap(trkFromCluster, extStartPos, ecalPosition, stepSize, bFieldMap);
            trkFromCluster.getTrackStates().add(state);

//            double dist = Math.abs(matcher.getXDistance(cl, trkFromCluster));
            double dist=Math.abs(matcher.getDistance(cl, trkFromCluster));
//            LOGGER.info("Finding Best Phi:  phi = " + trackParameters[BaseTrack.PHI] + "; distance = " + dist);
            if (dist < minDist) {
//                LOGGER.info("Finding Best Phi:  Found Better phi = " + trackParameters[BaseTrack.PHI] + "; distance = " + dist + "  old minDist = " + minDist);
                minDist = dist;
                bestPhi = thisPhi;
            }
        }
        return bestPhi;
    }

    public void printDQMData() {
        LOGGER.info("FinalStateMonitoring::printDQMData");
        for (Entry<String, Double> entry : monitoredQuantityMap.entrySet())
            LOGGER.info(entry.getKey() + " = " + entry.getValue());
        LOGGER.info("*******************************");
    }

    IFitResult fitBeamEnergyPeak(IHistogram1D h1d, IFitter fitter, String range) {
//        return fitter.fit(h1d, "g", range);

//        return fitter.fit(h1d, "g+p1", init, range);
        double[] init = {20.0, 2.2, 0.12, 10, 0.0};
//        double[] init = {20.0, 2.2, 0.1};
        IFitResult ifr = null;
        try {
            ifr = fitter.fit(h1d, "g+p1", init);
        } catch (RuntimeException ex) {
            LOGGER.info(this.getClass().getSimpleName() + ":  caught exception in fitGaussian");
        }

        return ifr;
    }

    double getECalCoplanarity(Cluster cl1, Cluster cl2) {
        Cluster clTop;
        Cluster clBottom;
        if (cl1.getPosition()[1] > 0) {
            clTop = cl1;
            clBottom = cl2;
        } else {
            clTop = cl2;
            clBottom = cl1;
        }
        double[] clTopPosition = clTop.getPosition();
        double[] clBottomPosition = clBottom.getPosition();

        double topX = clTopPosition[0];
        double topY = clTopPosition[1];
        double botX = clBottomPosition[0];
        double botY = clBottomPosition[1];
        double topE = clTop.getEnergy();
        double botE = clBottom.getEnergy();
        double Esum = topE + botE;
        double cl_impact_angleTop = Math.toDegrees(Math.atan2(topY, topX - phot_nom_x));
        double cl_impact_angleBottom = Math.toDegrees(Math.atan2(botY, botX - phot_nom_x));
        if (cl_impact_angleTop < 0.)
            cl_impact_angleTop = cl_impact_angleTop + 360.;
        if (cl_impact_angleBottom < 0.)
            cl_impact_angleBottom = cl_impact_angleBottom + 360.;
        return cl_impact_angleBottom - cl_impact_angleTop;

    }

    boolean inSuperFiducialRegion(double x, double y) {
        boolean in_fid = false;
        double x_edge_low = -262.74;
        double x_edge_high = 347.7;
        double y_edge_low = 33.54;
        double y_edge_high = 75.18;

        double x_gap_low = -106.66;
        double x_gap_high = 42.17;
        double y_gap_high = 47.18;

//        #set y_edge_low to the y of the electron gap!
        y_edge_low = y_gap_high;

        y = Math.abs(y);

        if (x > x_edge_low && x < x_edge_high && y > y_edge_low && y < y_edge_high)
            if ((x > x_gap_low && x < x_gap_high && y > y_edge_low && y < y_gap_high) != true)
                in_fid = true;
        return in_fid;
    }

    boolean inFiducialRegion(double x, double y) {
        boolean in_fid = false;
        double x_edge_low = -262.74;
        double x_edge_high = 347.7;
        double y_edge_low = 33.54;
        double y_edge_high = 75.18;

        double x_gap_low = -106.66;
        double x_gap_high = 42.17;
        double y_gap_high = 47.18;

        y = Math.abs(y);

        if (x > x_edge_low && x < x_edge_high && y > y_edge_low && y < y_edge_high)
            if ((x > x_gap_low && x < x_gap_high && y > y_edge_low && y < y_gap_high) != true)
                in_fid = true;
        return in_fid;
    }

    private void setupWriter() {
        // Cleanup existing writer.
        if (writer != null)
            try {
                writer.flush();
                writer.close();
                writer = null;
            } catch (IOException x) {
                System.err.println(x.getMessage());
            }

        // Setup new writer.
        try {
            writer = new LCIOWriter(outputFile);
        } catch (IOException x) {
            throw new RuntimeException("Error creating writer", x);
        }

        try {
            writer.reOpen();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    private void setupNoTrackWriter() {
        // Cleanup existing writer.
        if (notrackwriter != null)
            try {
                notrackwriter.flush();
                notrackwriter.close();
                notrackwriter = null;
            } catch (IOException x) {
                System.err.println(x.getMessage());
            }

        // Setup new writer.
        try {
            notrackwriter = new LCIOWriter(notrackFile);
        } catch (IOException x) {
            throw new RuntimeException("Error creating writer", x);
        }

        try {
            notrackwriter.reOpen();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    protected void startOfData() {
        setupWriter();
        setupNoTrackWriter();
    }

    @Override
    public void endOfData() {
        try {
            writer.close();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
         try {
            notrackwriter.close();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    /**
     * Extrapolate a track to a layer and check that it lies within its
     * acceptance.
     *
     * @param track The track that will be extrapolated to the layer of interest
     * @param layer The layer number to extrapolate to
     * @return true if the track lies within the sensor acceptance, false
     * otherwise
     */
    private boolean isWithinAcceptance(Track track, int layer) {

        // TODO: Move this to a utility class 
        //System.out.println("Retrieving sensors for layer: " + layer);
        // Since TrackUtils.isTop/BottomTrack does not work when running off 
        // a recon file, get the detector volume that a track is associated 
        // with by using the sensor.  This assumes that a track is always
        // composed by stereo hits that lie within a single detector volume
        //HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit)track.getTrackerHits().get(0).getRawHits().get(0)).getDetectorElement();
        boolean isTop = true;
        if (track.getTrackStates().get(0).getTanLambda() < 0)
            isTop = false;
        // Get the sensors associated with the layer that the track
        // will be extrapolated to
        List<SvtStereoLayer> stereoLayers = null;

        // if (TrackUtils.isTopTrack(track, track.getTrackerHits().size())) {
        if (isTop)
            //System.out.println("Top track found.");
            stereoLayers = this.topStereoLayers.get(layer); //} else if (TrackUtils.isBottomTrack(track, track.getTrackerHits().size())) {
        else
            //System.out.println("Bottom track found.");
            stereoLayers = this.bottomStereoLayers.get(layer);

        for (SvtStereoLayer stereoLayer : stereoLayers) {
            Hep3Vector axialSensorPosition = stereoLayer.getAxialSensor().getGeometry().getPosition();
            Hep3Vector stereoSensorPosition = stereoLayer.getStereoSensor().getGeometry().getPosition();

            //System.out.println("Axial sensor position: " + axialSensorPosition.toString());
            //System.out.println("Stereo sensor position: " + stereoSensorPosition.toString());
            Hep3Vector axialTrackPos = TrackUtils.extrapolateTrack(track, axialSensorPosition.z());
            Hep3Vector stereoTrackPos = TrackUtils.extrapolateTrack(track, stereoSensorPosition.z());
//            LOGGER.info("track position on axial sensor:  " + axialTrackPos.x() + ", " + axialTrackPos.y() + ", " + axialTrackPos.z());
//            LOGGER.info("track position on stero sensor:  " + stereoTrackPos.x() + ", " + stereoTrackPos.y() + ", " + stereoTrackPos.z());
            //System.out.println("Track position at axial sensor: " + axialTrackPos.toString());
            //System.out.println("Track position at stereo sensor: " + stereoTrackPos.toString());
            boolean inAxial = this.sensorContainsTrack(axialTrackPos, stereoLayer.getAxialSensor());
            boolean inStereo = this.sensorContainsTrack(stereoTrackPos, stereoLayer.getStereoSensor());
//            LOGGER.info("in Axial = " + inAxial + "; in Stereo = " + inStereo);
            if (inAxial
                    && inStereo)
                //System.out.println("Track lies within layer acceptance.");
                return true;
        }

        return false;

        /*int layerNumber = (layer - 1)/2 + 1;
         String title = "Track Position - Layer " + layerNumber + " - Tracks Within Acceptance";
         //aida.histogram2D(title).fill(trackPos.y(), trackPos.z());
         //aida.cloud2D(title).fill(frontTrackPos.y(), frontTrackPos.z()); */
    }

    public boolean sensorContainsTrack(Hep3Vector trackPosition, HpsSiSensor sensor) {

//        if(maskBadChannels){
//            int intersectingChannel = this.findIntersectingChannel(trackPosition, sensor);
//            if(intersectingChannel == 0 || intersectingChannel == 638) return false;
//
//            if(sensor.isBadChannel(intersectingChannel) 
//                    || sensor.isBadChannel(intersectingChannel+1) 
//                    || sensor.isBadChannel(intersectingChannel-1)){
//                //this.printDebug("Track intersects a bad channel!");
//                return false;
//                    }
//        }
        ITransform3D localToGlobal = sensor.getGeometry().getLocalToGlobal();

        Hep3Vector sensorPos = sensor.getGeometry().getPosition();
        Box sensorSolid = (Box) sensor.getGeometry().getLogicalVolume().getSolid();
        Polygon3D sensorFace = sensorSolid.getFacesNormalTo(new BasicHep3Vector(0, 0, 1)).get(0);

        List<Point3D> vertices = new ArrayList<Point3D>();
        for (int index = 0; index < 4; index++)
            vertices.add(new Point3D());
        for (Point3D vertex : sensorFace.getVertices())
            if (vertex.y() < 0 && vertex.x() > 0) {
                localToGlobal.transform(vertex);
                //vertices.set(0, new Point3D(vertex.y() + sensorPos.x(), vertex.x() + sensorPos.y(), vertex.z() + sensorPos.z()));
                vertices.set(0, new Point3D(vertex.x(), vertex.y(), vertex.z()));
                //System.out.println(this.getClass().getSimpleName() + ": Vertex 1 Position: " + vertices.get(0).toString());
                //System.out.println(this.getClass().getSimpleName() + ": Transformed Vertex 1 Position: " + localToGlobal.transformed(vertex).toString());
            } else if (vertex.y() > 0 && vertex.x() > 0) {
                localToGlobal.transform(vertex);
                //vertices.set(1, new Point3D(vertex.y() + sensorPos.x(), vertex.x() + sensorPos.y(), vertex.z() + sensorPos.z()));
                vertices.set(1, new Point3D(vertex.x(), vertex.y(), vertex.z()));
                //System.out.println(this.getClass().getSimpleName() + ": Vertex 2 Position: " + vertices.get(1).toString());
                //System.out.println(this.getClass().getSimpleName() + ": Transformed Vertex 2 Position: " + localToGlobal.transformed(vertex).toString());
            } else if (vertex.y() > 0 && vertex.x() < 0) {
                localToGlobal.transform(vertex);
                //vertices.set(2, new Point3D(vertex.y() + sensorPos.x(), vertex.x() + sensorPos.y(), vertex.z() + sensorPos.z()));
                vertices.set(2, new Point3D(vertex.x(), vertex.y(), vertex.z()));
                //System.out.println(this.getClass().getSimpleName() + ": Vertex 3 Position: " + vertices.get(2).toString());
                //System.out.println(this.getClass().getSimpleName() + ": Transformed Vertex 3 Position: " + localToGlobal.transformed(vertex).toString());
            } else if (vertex.y() < 0 && vertex.x() < 0) {
                localToGlobal.transform(vertex);
                //vertices.set(3, new Point3D(vertex.y() + sensorPos.x(), vertex.x() + sensorPos.y(), vertex.z() + sensorPos.z()));
                vertices.set(3, new Point3D(vertex.x(), vertex.y(), vertex.z()));
                //System.out.println(this.getClass().getSimpleName() + ": Vertex 4 Position: " + vertices.get(3).toString());
                //System.out.println(this.getClass().getSimpleName() + ": Transformed Vertex 4 Position: " + localToGlobal.transformed(vertex).toString());
            }
        /*
         double area1 = this.findTriangleArea(vertices.get(0).x(), vertices.get(0).y(), vertices.get(1).x(), vertices.get(1).y(), trackPosition.y(), trackPosition.z()); 
         double area2 = this.findTriangleArea(vertices.get(1).x(), vertices.get(1).y(), vertices.get(2).x(), vertices.get(2).y(), trackPosition.y(), trackPosition.z()); 
         double area3 = this.findTriangleArea(vertices.get(2).x(), vertices.get(2).y(), vertices.get(3).x(), vertices.get(3).y(), trackPosition.y(), trackPosition.z()); 
         double area4 = this.findTriangleArea(vertices.get(3).x(), vertices.get(3).y(), vertices.get(0).x(), vertices.get(0).y(), trackPosition.y(), trackPosition.z()); 
         */

        double area1 = this.findTriangleArea(vertices.get(0).x(), vertices.get(0).y(), vertices.get(1).x(), vertices.get(1).y(), trackPosition.x(), trackPosition.y());
        double area2 = this.findTriangleArea(vertices.get(1).x(), vertices.get(1).y(), vertices.get(2).x(), vertices.get(2).y(), trackPosition.x(), trackPosition.y());
        double area3 = this.findTriangleArea(vertices.get(2).x(), vertices.get(2).y(), vertices.get(3).x(), vertices.get(3).y(), trackPosition.x(), trackPosition.y());
        double area4 = this.findTriangleArea(vertices.get(3).x(), vertices.get(3).y(), vertices.get(0).x(), vertices.get(0).y(), trackPosition.x(), trackPosition.y());

        if ((area1 > 0 && area2 > 0 && area3 > 0 && area4 > 0) || (area1 < 0 && area2 < 0 && area3 < 0 && area4 < 0))
            return true;
        return false;

    }

    /**
     *
     */
    public double findTriangleArea(double x0, double y0, double x1, double y1, double x2, double y2) {
        return .5 * (x1 * y2 - y1 * x2 - x0 * y2 + y0 * x2 + x0 * y1 - y0 * x1);
    }

    public int getLayerFromPosition(double posZ) {
        if (posZ < 150)
            return 1;
        if (posZ < 250)
            return 2;
        if (posZ < 350)
            return 3;
        if (posZ < 550)
            return 4;
        if (posZ < 750)
            return 5;
        if (posZ < 950)
            return 6;
        return 666;
    }

    public ReconstructedParticle findFSParticleFromCluster(List<ReconstructedParticle> fsList, Cluster cl) {
        for (ReconstructedParticle fsPart : fsList)
            if (fsPart.getClusters().size() > 0 && fsPart.getClusters().get(0) == cl)
                return fsPart;
        LOGGER.info("Couldn't find an fs particle...should never get here??????");
        return null;
    }

    public int findIntersectingChannel(Hep3Vector trackPosition, SiSensor sensor) {

        //--- Check that the track doesn't pass through a region of bad channels ---//
        //--------------------------------------------------------------------------//
        //Rotate the track position to the JLab coordinate system
        //this.printDebug("Track position in tracking frame: " + trackPosition.toString());
        Hep3Vector trackPositionDet = VecOp.mult(VecOp.inverse(this.trackerHitUtils.detToTrackRotationMatrix()), trackPosition);
        //this.printDebug("Track position in JLab frame " + trackPositionDet.toString());
        // Rotate the track to the sensor coordinates
        ITransform3D globalToLocal = sensor.getReadoutElectrodes(ChargeCarrier.HOLE).getGlobalToLocal();
        globalToLocal.transform(trackPositionDet);
        //this.printDebug("Track position in sensor electrodes frame " + trackPositionDet.toString());

        // Find the closest channel to the track position
        double deltaY = Double.MAX_VALUE;
        int intersectingChannel = 0;
        for (int physicalChannel = 0; physicalChannel < 639; physicalChannel++) {
            /*this.printDebug(SvtUtils.getInstance().getDescription(sensor) + " : Channel " + physicalChannel 
             + " : Strip Position " + stripPositions.get(sensor).get(physicalChannel));
             this.printDebug(SvtUtils.getInstance().getDescription(sensor) + ": Channel " + physicalChannel
             + " : Delta Y: " + Math.abs(trackPositionDet.x() - stripPositions.get(sensor).get(physicalChannel).x()));*/
            /*if(Math.abs(trackPositionDet.x() - stripPositions.get(sensor).get(physicalChannel).x()) < deltaY ){
             deltaY = Math.abs(trackPositionDet.x() - stripPositions.get(sensor).get(physicalChannel).x()); 
             intersectingChannel = physicalChannel;
             }*/
        }

        return intersectingChannel;
    }

}
