/**
 * Driver used to compute SVT hit efficiencies at each sensor
 * as a function of strip and y
 */
/**
 * @author mrsolt
 *
 */
package org.hps.users.mrsolt;

import static java.lang.Math.abs;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogramFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.ITree;
import hep.physics.matrix.BasicMatrix;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import org.lcsim.constants.Constants;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.solids.LineSegment3D;
import org.lcsim.detector.solids.Polygon3D;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.apache.commons.math3.util.Pair;
import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.svt.AbstractSvtDaqMapping;
import org.hps.conditions.svt.SvtChannel;
import org.hps.conditions.svt.SvtConditions;
import org.hps.conditions.svt.SvtDaqMapping;
import org.hps.conditions.svt.SvtDaqMapping.SvtDaqMappingCollection;
import org.hps.conditions.svt.SvtChannel.SvtChannelCollection;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.HpsHelicalTrackFit;
import org.hps.recon.tracking.TrackStateUtils;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackerHitUtils;
import org.hps.recon.tracking.gbl.FittedGblTrajectory;
import org.hps.recon.tracking.gbl.GBLStripClusterData;
import org.hps.recon.tracking.gbl.GblUtils;
import org.hps.recon.tracking.gbl.matrix.Matrix;

public class SvtHitEfficiency extends Driver {


    // Use JFreeChart as the default plotting backend
    static { 
        hep.aida.jfree.AnalysisFactory.register();
    }

    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree; 
    IHistogramFactory histogramFactory; 

    //List of Sensors
    private List<HpsSiSensor> sensors = null;
    
    Map<String, IHistogram1D> numberOfTracksChannel = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerChannel = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannel = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksP = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerP = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyP = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> numberOfTracksChannelEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerChannelEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannelEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksYEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerYEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksPEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerPEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPEle = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> numberOfTracksChannelPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerChannelPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannelPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksYPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerYPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksPPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksWithHitOnMissingLayerPPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPPos = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> numberOfTracksChannelCorrected = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannelCorrected = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksYCorrected = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrected = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksPCorrected = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrected = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> numberOfTracksChannelCorrectedEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannelCorrectedEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksYCorrectedEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrectedEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksPCorrectedEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrectedEle = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> numberOfTracksChannelCorrectedPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyChannelCorrectedPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksYCorrectedPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrectedPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> numberOfTracksPCorrectedPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrectedPos = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> TotalEff = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalEffEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalEffPos = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> TotalCorrectedEff = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalCorrectedEffEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalCorrectedEffPos = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPerr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelCorrectederr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrectederr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrectederr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPEleerr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelCorrectedEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrectedEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrectedEleerr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelPoserr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYPoserr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPPoserr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> hitEfficiencyChannelCorrectedPoserr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyYCorrectedPoserr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> hitEfficiencyPCorrectedPoserr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> TotalEfferr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalEffEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalEffPoserr = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> TotalCorrectedEfferr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalCorrectedEffEleerr = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> TotalCorrectedEffPoserr = new HashMap<String,IHistogram1D>();
    
    
    Map<String, IHistogram1D> residualY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> errorY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> pullY = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> residualYerrorYDiff = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram2D> residualYvsV = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsV = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsV = new HashMap<String,IHistogram2D>();
    
    Map<String, IHistogram2D> residualYvsU = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsU = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsU = new HashMap<String,IHistogram2D>();
    
    Map<String, IHistogram1D> residualYEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> errorYEle = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> pullYEle = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram2D> residualYvsVEle = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsVEle = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsVEle = new HashMap<String,IHistogram2D>();
    
    Map<String, IHistogram2D> residualYvsUEle = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsUEle = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsUEle = new HashMap<String,IHistogram2D>();
    
    Map<String, IHistogram1D> residualYPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> errorYPos = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> pullYPos = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram2D> residualYvsVPos = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsVPos = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsVPos = new HashMap<String,IHistogram2D>();
    
    Map<String, IHistogram2D> residualYvsUPos = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> errorYvsUPos = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> pullYvsUPos = new HashMap<String,IHistogram2D>();
    
    
    
    Map<String, IHistogram1D> D0 = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Z0 = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Tanlambda = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Phi0 = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Omega = new HashMap<String,IHistogram1D>();
    
    Map<String, IHistogram1D> D0_err = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Z0_err = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Tanlambda_err = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Phi0_err = new HashMap<String,IHistogram1D>();
    Map<String, IHistogram1D> Omega_err = new HashMap<String,IHistogram1D>();
    
    //Histogram Settings
    int nBins = 50;
    double maxPull = 7;
    double minPull = -maxPull;
    double maxRes = 2;
    double minRes = -maxRes;
    double maxYerror = 1;
    double maxD0 = 5;
    double minD0 = -maxD0;
    double maxZ0 = 10;
    double minZ0 = -maxZ0;
    double maxTLambda = 0.1;
    double minTLambda = -maxTLambda;
    double maxPhi0 = 0.2;
    double minPhi0 = -maxPhi0;
    double maxOmega = 0.001;
    double minOmega = -maxOmega;
    double maxD0Err = 1.5;
    double minD0Err = 0;
    double maxZ0Err = 1;
    double minZ0Err = 0;
    double maxTLambdaErr = 0.005;
    double minTLambdaErr = 0;
    double maxPhi0Err = 0.01;
    double minPhi0Err = 0;
    double maxOmegaErr = 0.0001;
    double minOmegaErr = 0;
    
    int chanExtd = 0;
    
    String atIP = "IP";
    
    //Collection Strings
    private String stripHitOutputCollectionName = "StripClusterer_SiTrackerHitStrip1D";
    private String GBLTrackCollectionName = "GBLTracks";
    
    TrackerHitUtils trackerHitUtils = new TrackerHitUtils();
    
    boolean maskBadChannels = false;
    
    protected static double bfield;
    FieldMap bFieldMap = null;
   
    //Constants
    private static final String SUBDETECTOR_NAME = "Tracker";
    //Configurable Variables
    boolean debug = false;
    
    String outputFileName = "channelEff.txt";
    boolean cleanTridents = false;
    int nLay = 6;
    double nSig = 5;
    
    SvtChannelCollection channelMap;
    SvtDaqMappingCollection daqMap;
    
    
    
    public void setOutputFileName(String outputFileName) {
        this.outputFileName = outputFileName;
    }
    
    public void setCleanTridents(boolean cleanTridents) { 
        this.cleanTridents = cleanTridents;
    }
    
    public void setNLay(int nLay) { 
        this.nLay = nLay;
    }
    
    public void setSig(double nSig) { 
        this.nSig = nSig;
    }
    
    public void setChanExtd(int chanExtd) { 
        this.chanExtd = chanExtd;
    }
    
    public void setMaskBadChannels(boolean maskBadChannels){
        this.maskBadChannels = maskBadChannels; 
    }
    
    //Beam Energy
    double ebeam;
    
    public void detectorChanged(Detector detector){
    	
    	aida.tree().cd("/");
    	tree = aida.tree();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);
        
        DatabaseConditionsManager mgr = DatabaseConditionsManager.getInstance();
        SvtConditions svtConditions = mgr.getCachedConditions(SvtConditions.class, "svt_conditions").getCachedData();
        
        channelMap = svtConditions.getChannelMap();
        daqMap = svtConditions.getDaqMap();
    
        //Set Beam Energy
        BeamEnergyCollection beamEnergyCollection = 
                this.getConditionsManager().getCachedConditions(BeamEnergyCollection.class, "beam_energies").getCachedData();        
        ebeam = beamEnergyCollection.get(0).getBeamEnergy();
        
        bfield = TrackUtils.getBField(detector).magnitude();
        bFieldMap = detector.getFieldMap();
        
        // Get the HpsSiSensor objects from the tracker detector element
        sensors = detector.getSubdetector(SUBDETECTOR_NAME)
                          .getDetectorElement().findDescendants(HpsSiSensor.class);
   
        // If the detector element had no sensors associated with it, throw
        // an exception
        if (sensors.size() == 0) {
            throw new RuntimeException("No sensors were found in this detector.");
        }
        
        D0.put(atIP,histogramFactory.createHistogram1D("D0 " + atIP, nBins, minD0, maxD0));
    	Z0.put(atIP,histogramFactory.createHistogram1D("Z0 " + atIP, nBins, minZ0, maxZ0));
    	Tanlambda.put(atIP,histogramFactory.createHistogram1D("TanLambda " + atIP, nBins, minTLambda, maxTLambda));
    	Phi0.put(atIP,histogramFactory.createHistogram1D("Phi0 " + atIP, nBins, minPhi0, maxPhi0));
    	Omega.put(atIP,histogramFactory.createHistogram1D("Omega " + atIP, nBins, minOmega, maxOmega)); 
    	
    	D0_err.put(atIP,histogramFactory.createHistogram1D("D0 Error " + atIP, nBins, minD0Err, maxD0Err));
    	Z0_err.put(atIP,histogramFactory.createHistogram1D("Z0 Error " + atIP, nBins, minZ0Err, maxZ0Err));
    	Tanlambda_err.put(atIP,histogramFactory.createHistogram1D("TanLambda Error " + atIP, nBins, minTLambdaErr, maxTLambdaErr));
    	Phi0_err.put(atIP,histogramFactory.createHistogram1D("Phi0 Error " + atIP, nBins, minPhi0Err, maxPhi0Err));
    	Omega_err.put(atIP,histogramFactory.createHistogram1D("Omega Error " + atIP, nBins, minOmegaErr, maxOmegaErr));

        for(HpsSiSensor sensor:sensors){
        	String sensorName = sensor.getName();
        	int nChan = sensor.getNumberOfChannels()+chanExtd;
        	double readoutPitch = sensor.getReadoutStripPitch();
        	double maxU = nChan * readoutPitch / 2;
        	double width = getSensorLength(sensor);
        	double maxV = width/2;
            double minV = -maxV;
        	numberOfTracksChannel.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksWithHitOnMissingLayerChannel.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannel.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksY.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksWithHitOnMissingLayerY.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyY.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksP.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerP.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyP.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Ele " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksWithHitOnMissingLayerChannelEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel Ele " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannelEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Ele " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksYEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Ele " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksWithHitOnMissingLayerYEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y Ele " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyYEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Ele " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksPEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Ele " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerPEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyPEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Pos " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksWithHitOnMissingLayerChannelPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel Pos " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannelPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Pos " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksYPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Pos " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksWithHitOnMissingLayerYPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y Pos " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyYPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Pos " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksPPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Pos " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerPPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyPPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannelCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksYCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyYCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksPCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected Ele " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannelCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Ele " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksYCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected Ele " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyYCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Ele " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksPCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected Ele " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected Pos " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyChannelCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Pos " + sensorName, nChan, -chanExtd, nChan));
        	numberOfTracksYCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected Pos " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyYCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Pos " + sensorName, nBins, -maxU, maxU));
        	numberOfTracksPCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected Pos " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	TotalEff.put(sensorName,histogramFactory.createHistogram1D("Total Eff " + sensorName, 1, 0, 1));
        	TotalEffEle.put(sensorName,histogramFactory.createHistogram1D("Total Eff Ele " + sensorName, 1, 0, 1));
        	TotalEffPos.put(sensorName,histogramFactory.createHistogram1D("Total Eff Pos " + sensorName, 1, 0, 1));
        	
        	TotalCorrectedEff.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff " + sensorName, 1, 0, 1));
        	TotalCorrectedEffEle.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Ele " + sensorName, 1, 0, 1));
        	TotalCorrectedEffPos.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Pos " + sensorName, 1, 0, 1));
        	
        	hitEfficiencyChannelerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Ele Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Ele Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Ele Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Ele Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Ele Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Ele Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Pos Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Pos Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Pos Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Pos Error " + sensorName, nChan, -chanExtd, nChan));
        	hitEfficiencyYCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Pos Error " + sensorName, nBins, -maxU, maxU));
        	hitEfficiencyPCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Pos Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	TotalEfferr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Error " + sensorName, 1, 0, 1));
        	TotalEffEleerr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Ele Error " + sensorName, 1, 0, 1));
        	TotalEffPoserr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Pos Error " + sensorName, 1, 0, 1));
        	
        	TotalCorrectedEfferr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Error " + sensorName, 1, 0, 1));
        	TotalCorrectedEffEleerr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Ele Error " + sensorName, 1, 0, 1));
        	TotalCorrectedEffPoserr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Pos Error " + sensorName, 1, 0, 1));
        	
        	residualY.put(sensorName,histogramFactory.createHistogram1D("Residual U " + sensorName, nBins, minRes, maxRes));
        	errorY.put(sensorName,histogramFactory.createHistogram1D("Error U " + sensorName, nBins, 0, maxYerror));
        	pullY.put(sensorName,histogramFactory.createHistogram1D("U Pulls " + sensorName, nBins, minPull, maxPull));
        	residualYerrorYDiff.put(sensorName,histogramFactory.createHistogram1D("Residual U - Error U " + sensorName, nBins, -10, 10));
        	
        	residualYvsV.put(sensorName,histogramFactory.createHistogram2D("Residual U vs V " + sensorName, 2*nBins, minV, maxV, nBins, minRes, maxRes));
        	errorYvsV.put(sensorName,histogramFactory.createHistogram2D("Error U vs V " + sensorName, 2*nBins, minV, maxV, nBins, 0, maxYerror));
        	pullYvsV.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs V " + sensorName, 2*nBins, minV, maxV, nBins, minPull, maxPull));
        	
        	residualYvsU.put(sensorName,histogramFactory.createHistogram2D("Residual U vs U " + sensorName, 2*nBins, -maxU, maxU, nBins, minRes, maxRes));
        	errorYvsU.put(sensorName,histogramFactory.createHistogram2D("Error U vs U " + sensorName, 2*nBins, -maxU, maxU, nBins, 0, maxYerror));
        	pullYvsU.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs U " + sensorName, 2*nBins, -maxU, maxU, nBins, minPull, maxPull));
        	
        	residualYEle.put(sensorName,histogramFactory.createHistogram1D("Residual U Electron " + sensorName, nBins, minRes, maxRes));
        	errorYEle.put(sensorName,histogramFactory.createHistogram1D("Error U Electron " + sensorName, nBins, 0, maxYerror));
        	pullYEle.put(sensorName,histogramFactory.createHistogram1D("U Pulls Electron " + sensorName, nBins, minPull, maxPull));
        	
        	residualYvsVEle.put(sensorName,histogramFactory.createHistogram2D("Residual U vs V Electron " + sensorName, 2*nBins, minV, maxV, nBins, minRes, maxRes));
        	errorYvsVEle.put(sensorName,histogramFactory.createHistogram2D("Error U vs V Electron " + sensorName, 2*nBins, minV, maxV, nBins, 0, maxYerror));
        	pullYvsVEle.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs V Electron " + sensorName, 2*nBins, minV, maxV, nBins, minPull, maxPull));
        	
        	residualYvsUEle.put(sensorName,histogramFactory.createHistogram2D("Residual U vs U Electron " + sensorName, 2*nBins, -maxU, maxU, nBins, minRes, maxRes));
        	errorYvsUEle.put(sensorName,histogramFactory.createHistogram2D("Error U vs U Electron " + sensorName, 2*nBins, -maxU, maxU, nBins, 0, maxYerror));
        	pullYvsUEle.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs U Electron " + sensorName, 2*nBins, -maxU, maxU, nBins, minPull, maxPull));
        	
        	residualYPos.put(sensorName,histogramFactory.createHistogram1D("Residual U Positron " + sensorName, nBins, minRes, maxRes));
        	errorYPos.put(sensorName,histogramFactory.createHistogram1D("Error U Positron " + sensorName, nBins, 0, maxYerror));
        	pullYPos.put(sensorName,histogramFactory.createHistogram1D("U Pulls Positron " + sensorName, nBins, minPull, maxPull));
        	
        	residualYvsVPos.put(sensorName,histogramFactory.createHistogram2D("Residual U vs V Positron " + sensorName, 2*nBins, minV, maxV, nBins, minRes, maxRes));
        	errorYvsVPos.put(sensorName,histogramFactory.createHistogram2D("Error U vs V Positron " + sensorName, 2*nBins, minV, maxV, nBins, 0, maxYerror));
        	pullYvsVPos.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs V Positron " + sensorName, 2*nBins, minV, maxV, nBins, minPull, maxPull));
        	
        	residualYvsUPos.put(sensorName,histogramFactory.createHistogram2D("Residual U vs U Positron " + sensorName, 2*nBins, -maxU, maxU, nBins, minRes, maxRes));
        	errorYvsUPos.put(sensorName,histogramFactory.createHistogram2D("Error U vs U Positron " + sensorName, 2*nBins, -maxU, maxU, nBins, 0, maxYerror));
        	pullYvsUPos.put(sensorName,histogramFactory.createHistogram2D("U Pulls vs U Positron " + sensorName, 2*nBins, -maxU, maxU, nBins, minPull, maxPull));
        	
        	D0.put(sensorName,histogramFactory.createHistogram1D("D0 " + sensorName, nBins, minD0, maxD0));
        	Z0.put(sensorName,histogramFactory.createHistogram1D("Z0 " + sensorName, nBins, minZ0, maxZ0));
        	Tanlambda.put(sensorName,histogramFactory.createHistogram1D("TanLambda " + sensorName, nBins, minTLambda, maxTLambda));
        	Phi0.put(sensorName,histogramFactory.createHistogram1D("Phi0 " + sensorName, nBins, minPhi0, maxPhi0));
        	Omega.put(sensorName,histogramFactory.createHistogram1D("Omega " + sensorName, nBins, minOmega, maxOmega));
        	
        	D0_err.put(sensorName,histogramFactory.createHistogram1D("D0 Error " + sensorName, nBins, minD0Err, maxD0Err));
        	Z0_err.put(sensorName,histogramFactory.createHistogram1D("Z0 Error " + sensorName, nBins, minZ0Err, maxZ0Err));
        	Tanlambda_err.put(sensorName,histogramFactory.createHistogram1D("TanLambda Error " + sensorName, nBins, minTLambdaErr, maxTLambdaErr));
        	Phi0_err.put(sensorName,histogramFactory.createHistogram1D("Phi0 Error " + sensorName, nBins, minPhi0Err, maxPhi0Err));
        	Omega_err.put(sensorName,histogramFactory.createHistogram1D("Omega Error " + sensorName, nBins, minOmegaErr, maxOmegaErr));
        }
    }

    public void process(EventHeader event){
		aida.tree().cd("/");	
        
		//Grab all the lists of tracks in each track collection
        List<List<Track>> trackCollections = event.get(Track.class);
        
        List<Track> tracks = event.get(Track.class,GBLTrackCollectionName);
        
        //Grab all the clusters in the event
        List<SiTrackerHitStrip1D> stripHits = event.get(SiTrackerHitStrip1D.class, stripHitOutputCollectionName);
        
        //for(List<Track> tracks:trackCollections){
        	
    		if(cleanTridents){
    			// Require an event to have exactly two tracks
    			if (tracks.size() != 2) return;

    			// Require the two tracks to be in opposite volumes
    			if (tracks.get(0).getTrackStates().get(0).getTanLambda()*tracks.get(1).getTrackStates().get(0).getTanLambda() >= 0) return;

    			// Require the two tracks to be oppositely charged
    			if (tracks.get(0).getTrackStates().get(0).getOmega()*tracks.get(1).getTrackStates().get(0).getOmega() >= 0) return;
    		}
    		
            for(Track track:tracks){
            	//boolean gbl = isGBL(track);//fix this function
            	//if(!gbl) continue;
            	//Grab the unused layer on the track
        	    int unusedLay = getUnusedSvtLayer(track.getTrackerHits());   
        	    if(unusedLay == -1) continue;
        	    
        	    List<TrackState> TStates = track.getTrackStates();
        	    TrackState tState = getTrackState(track,unusedLay);
        	    if(tState == null){
        	    	continue;
        	    }
        	    
        	    double[] covAtIP = TStates.get(0).getCovMatrix();
            	SymmetricMatrix LocCovAtIP = new SymmetricMatrix(5,covAtIP,true);
        	    
                D0.get(atIP).fill(TStates.get(0).getD0());
        		Z0.get(atIP).fill(TStates.get(0).getZ0());
        		Tanlambda.get(atIP).fill(TStates.get(0).getTanLambda());
        		Phi0.get(atIP).fill(TStates.get(0).getPhi());
        		Omega.get(atIP).fill(TStates.get(0).getOmega());
        		
        		D0_err.get(atIP).fill(Math.sqrt(LocCovAtIP.e(HelicalTrackFit.dcaIndex,HelicalTrackFit.dcaIndex)));
        		Z0_err.get(atIP).fill(Math.sqrt(LocCovAtIP.e(HelicalTrackFit.z0Index,HelicalTrackFit.z0Index)));
        		Tanlambda_err.get(atIP).fill(Math.sqrt(LocCovAtIP.e(HelicalTrackFit.slopeIndex,HelicalTrackFit.slopeIndex)));
        		Phi0_err.get(atIP).fill(Math.sqrt(LocCovAtIP.e(HelicalTrackFit.phi0Index,HelicalTrackFit.phi0Index)));
        		Omega_err.get(atIP).fill(Math.sqrt(LocCovAtIP.e(HelicalTrackFit.curvatureIndex,HelicalTrackFit.curvatureIndex)));
        	    
        	    //tState = track.getTrackStates().get(0);
    	    	
        	    Hep3Vector p = toHep3(tState.getMomentum());
        	    
        	    //See if track is within acceptance of both the axial and stereo sensors of the unused layer
        	    Pair<HpsSiSensor,Integer> axialSensorPair = isWithinSensorAcceptance(track,tState,unusedLay,true, p,bFieldMap);
        	    Pair<HpsSiSensor,Integer> stereoSensorPair = isWithinSensorAcceptance(track,tState,unusedLay,false,p,bFieldMap);        	    
        	    
        	    if(axialSensorPair == null || stereoSensorPair == null) continue;
        	    
        	    HpsSiSensor axialSensor = axialSensorPair.getFirst();
        	    HpsSiSensor stereoSensor = stereoSensorPair.getFirst();
        	    
        	    String sensorAxialName = axialSensor.getName();
    	    	String sensorStereoName = stereoSensor.getName();
        	    
        	    //Grab the sensor positions
        	    Hep3Vector axialPos = axialSensor.getGeometry().getPosition();
    	    	Hep3Vector stereoPos = stereoSensor.getGeometry().getPosition();
    	    	
    	    	//Get the track state at the axial sensor at the unused layer + 1
    	    	//If unused layer is 6, then grab the track state at L5

    	    	//double[] cov = tState.getCovMatrix();
    	    	double q = track.getCharge();

    	    	//Grab the hit position of the sensor where the track state is found above
    	    	Hep3Vector endHitPos = getHitPos(track,unusedLay);
	    		
    	    	//Extrapolate from the unused layer + 1 to the unused layer for both axial and stereo sensors
    	    	//If the unused layer is 6, then extraploate from L5 to L6
    	    	HelicalTrackFit htf = TrackUtils.getHTF(tState);
    	    	//Hep3Matrix matrix = getPerToClPrj(htf);
    	    	
    	    	//Hep3Vector axialExtrapPos = extrapolateTrackPosition(tState,endHitPos.x(),axialPos.z(),5,bFieldMap);
    	    	//Hep3Vector stereoExtrapPos = extrapolateTrackPosition(tState,endHitPos.x(),stereoPos.z(),5,bFieldMap);
    	    	//System.out.println(event.getEventNumber());
    	    	//Hep3Vector axialExtrapPos = TrackUtils.extrapolateTrackPositionToSensor(track.getTrackStates().get(0), axialSensor,bfield);
    	    	//Hep3Vector stereoExtrapPos = TrackUtils.extrapolateTrackPositionToSensor(track.getTrackStates().get(0), stereoSensor,bfield);
    	    	
    	    	//Hep3Vector axialExtrapPos = TrackUtils.extrapolateTrackPositionToSensor(track, axialSensor,sensors,bfield);
    	    	//Hep3Vector stereoExtrapPos = TrackUtils.extrapolateTrackPositionToSensor(track, stereoSensor,sensors,bfield);
    	    	
    	    	Hep3Vector axialExtrapPos = getLocationAtSensor(tState, axialSensor, bfield);
    	    	Hep3Vector stereoExtrapPos = getLocationAtSensor(tState, stereoSensor, bfield);
    	    	
    	    	if(axialExtrapPos == null || stereoExtrapPos == null){
    	    		System.out.println("Extrapolation is NULL");
    	    		continue;
    	    	}
    	    	
    	    	Hep3Vector axialExtrapPosSensor = globalToSensor(axialExtrapPos,axialSensor);
    	    	Hep3Vector stereoExtrapPosSensor = globalToSensor(stereoExtrapPos,stereoSensor);
    	    	
    	    	//Compute the extrapolation errors
    	    	//compute error correctly
    	    	//System.out.println("Layer " + unusedLay);
    	    	double yErrorAxial = computeExtrapErrorY(track,tState,axialSensor,endHitPos)[0];
    	    	double yErrorStereo = computeExtrapErrorY(track,tState,stereoSensor,endHitPos)[0];
    	    	
    	    	//Compute the channel where the track extrapolates to in each sensor
    	    	int chanAxial = axialSensorPair.getSecond();
    	    	int chanStereo = stereoSensorPair.getSecond();
    	    	
    	    	double trackP = toHep3(track.getTrackStates().get(0).getMomentum()).magnitude();
    	    	
    	    	double weightAxial = findWeight(axialExtrapPosSensor.x(),yErrorAxial,axialSensor);
    	    	double weightStereo = findWeight(stereoExtrapPosSensor.x(),yErrorStereo,stereoSensor);
    	    	
    	    	//Fill the denominator of the efficiency histos
    	    	numberOfTracksChannel.get(sensorAxialName).fill(chanAxial);
        	    numberOfTracksChannel.get(sensorStereoName).fill(chanStereo);
        	    numberOfTracksY.get(sensorAxialName).fill(axialExtrapPosSensor.x());
        	    numberOfTracksY.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
        	    numberOfTracksP.get(sensorAxialName).fill(trackP);
        	    numberOfTracksP.get(sensorStereoName).fill(trackP);
        	    
        	    numberOfTracksChannelCorrected.get(sensorAxialName).fill(chanAxial,weightAxial);
        	    numberOfTracksChannelCorrected.get(sensorStereoName).fill(chanStereo,weightStereo);
        	    numberOfTracksYCorrected.get(sensorAxialName).fill(axialExtrapPosSensor.x(),weightAxial);
        	    numberOfTracksYCorrected.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),weightStereo);
        	    numberOfTracksPCorrected.get(sensorAxialName).fill(trackP,weightAxial);
        	    numberOfTracksPCorrected.get(sensorStereoName).fill(trackP,weightStereo);
        	    
        	    if(q < 0){
        	    	numberOfTracksChannelEle.get(sensorAxialName).fill(chanAxial);
            	    numberOfTracksChannelEle.get(sensorStereoName).fill(chanStereo);
            	    numberOfTracksYEle.get(sensorAxialName).fill(axialExtrapPosSensor.x());
            	    numberOfTracksYEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
            	    numberOfTracksPEle.get(sensorAxialName).fill(trackP);
            	    numberOfTracksPEle.get(sensorStereoName).fill(trackP);
            	    
            	    numberOfTracksChannelCorrectedEle.get(sensorAxialName).fill(chanAxial,weightAxial);
            	    numberOfTracksChannelCorrectedEle.get(sensorStereoName).fill(chanStereo,weightStereo);
            	    numberOfTracksYCorrectedEle.get(sensorAxialName).fill(axialExtrapPosSensor.x(),weightAxial);
            	    numberOfTracksYCorrectedEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),weightStereo);
            	    numberOfTracksPCorrectedEle.get(sensorAxialName).fill(trackP,weightAxial);
            	    numberOfTracksPCorrectedEle.get(sensorStereoName).fill(trackP,weightStereo);
        	    }
        	    else{
        	    	numberOfTracksChannelPos.get(sensorAxialName).fill(chanAxial);
            	    numberOfTracksChannelPos.get(sensorStereoName).fill(chanStereo);
            	    numberOfTracksYPos.get(sensorAxialName).fill(axialExtrapPosSensor.x());
            	    numberOfTracksYPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
            	    numberOfTracksPPos.get(sensorAxialName).fill(trackP);
            	    numberOfTracksPPos.get(sensorStereoName).fill(trackP);
            	    
            	    numberOfTracksChannelCorrectedPos.get(sensorAxialName).fill(chanAxial,weightAxial);
            	    numberOfTracksChannelCorrectedPos.get(sensorStereoName).fill(chanStereo,weightStereo);
            	    numberOfTracksYCorrectedPos.get(sensorAxialName).fill(axialExtrapPosSensor.x(),weightAxial);
            	    numberOfTracksYCorrectedPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),weightStereo);
            	    numberOfTracksPCorrectedPos.get(sensorAxialName).fill(trackP,weightAxial);
            	    numberOfTracksPCorrectedPos.get(sensorStereoName).fill(trackP,weightStereo);
        	    }
    	    	
        	    //Fill the error histos
    	    	errorY.get(sensorAxialName).fill(yErrorAxial);
    	    	errorY.get(sensorStereoName).fill(yErrorStereo);
    	    	errorYvsV.get(sensorAxialName).fill(axialExtrapPosSensor.y(),yErrorAxial);
    	    	errorYvsV.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),yErrorStereo);
    	    	errorYvsU.get(sensorAxialName).fill(axialExtrapPosSensor.x(),yErrorAxial);
    	    	errorYvsU.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),yErrorStereo);
    	    	if(q < 0){
    	    		errorYEle.get(sensorAxialName).fill(yErrorAxial);
        	    	errorYEle.get(sensorStereoName).fill(yErrorStereo);
        	    	errorYvsVEle.get(sensorAxialName).fill(axialExtrapPosSensor.y(),yErrorAxial);
        	    	errorYvsVEle.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),yErrorStereo);
        	    	errorYvsUEle.get(sensorAxialName).fill(axialExtrapPosSensor.x(),yErrorAxial);
        	    	errorYvsUEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),yErrorStereo);
    	    	}
    	    	else{
    	    		errorYPos.get(sensorAxialName).fill(yErrorAxial);
        	    	errorYPos.get(sensorStereoName).fill(yErrorStereo);
        	    	errorYvsVPos.get(sensorAxialName).fill(axialExtrapPosSensor.y(),yErrorAxial);
        	    	errorYvsVPos.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),yErrorStereo);
        	    	errorYvsUPos.get(sensorAxialName).fill(axialExtrapPosSensor.x(),yErrorAxial);
        	    	errorYvsUPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),yErrorStereo);
    	    	}
    	    	double residualAxial = 9999;
    	    	double residualStereo = 9999;
    	    	
    	    	for(SiTrackerHitStrip1D hit:stripHits){
    	    		//Get the sensor and position of the hit
    	    		HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();
            	    double[] hitPos = hit.getPosition();
            	    Hep3Vector hitPosSensor = globalToSensor(toHep3(hitPos),sensor);
            	    //System.out.println(sensor.getName() + " hit position = " + hitPosSensor);
            	    //Check to see if the sensor of this hit is the same sensor you expect to see an axial hit
            	    if(sensorAxialName == sensor.getName()){
            	    	//Compute the residual between extrapolated track and hit position and fill histo
            	    	//double residual = axialExtrapPos.y() - hitPos[1];
            	    	double residual = axialExtrapPosSensor.x() - hitPosSensor.x();
            	    	if(Math.abs(residual) < Math.abs(residualAxial)){
            	    		residualAxial = residual;
            	    	}
            	    }
            	  //Check to see if the sensor of this hit is the same sensor you expect to see a stereo hit
            	    if(sensorStereoName == sensor.getName()){
            	    	//Compute the residual between extrapolated track and hit position and fill histo
            	    	//double residual = stereoExtrapPos.y() - hitPos[1];
            	    	double residual = stereoExtrapPosSensor.x() - hitPosSensor.x();
            	    	if(Math.abs(residual) < Math.abs(residualStereo)){
            	    		residualStereo = residual;
            	    	}
            	    }
    	    	}
    	    	
    	    	
    	    	
    	    	
        	    //Compute the residual between extrapolated track and hit position and fill histo
        	    //double residual = axialExtrapPos.y() - hitPos[1];
        	    residualY.get(sensorAxialName).fill(residualAxial);
        	    pullY.get(sensorAxialName).fill(residualAxial/yErrorAxial);
        	    residualYerrorYDiff.get(sensorAxialName).fill(Math.abs(residualAxial) - yErrorAxial);
        	    residualY.get(sensorStereoName).fill(residualStereo);
        	    pullY.get(sensorStereoName).fill(residualStereo/yErrorStereo);
        	    residualYerrorYDiff.get(sensorStereoName).fill(Math.abs(residualStereo) - yErrorStereo);
        	    
        	    residualYvsV.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial);
        	    pullYvsV.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial/yErrorAxial);
        	    residualYvsV.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo);
        	    pullYvsV.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo/yErrorStereo);
        	    
        	    residualYvsU.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial);
        	    pullYvsU.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial/yErrorAxial);
        	    residualYvsU.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo);
        	    pullYvsU.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo/yErrorStereo);
        	    
        	    if(q < 0){
        	    	residualYEle.get(sensorAxialName).fill(residualAxial);
            	    pullYEle.get(sensorAxialName).fill(residualAxial/yErrorAxial);
            	    residualYEle.get(sensorStereoName).fill(residualStereo);
            	    pullYEle.get(sensorStereoName).fill(residualStereo/yErrorStereo);
            	    
            	    residualYvsVEle.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial);
            	    pullYvsVEle.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial/yErrorAxial);
            	    residualYvsVEle.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo);
            	    pullYvsVEle.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo/yErrorStereo);
            	    
            	    residualYvsUEle.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial);
            	    pullYvsUEle.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial/yErrorAxial);
            	    residualYvsUEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo);
            	    pullYvsUEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo/yErrorStereo);
        	    }
        	    else{
        	    	residualYPos.get(sensorAxialName).fill(residualAxial);
            	    pullYPos.get(sensorAxialName).fill(residualAxial/yErrorAxial);
            	    residualYPos.get(sensorStereoName).fill(residualStereo);
            	    pullYPos.get(sensorStereoName).fill(residualStereo/yErrorStereo);
            	    
            	    residualYvsVPos.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial);
            	    pullYvsVPos.get(sensorAxialName).fill(axialExtrapPosSensor.y(),residualAxial/yErrorAxial);
            	    residualYvsVPos.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo);
            	    pullYvsVPos.get(sensorStereoName).fill(stereoExtrapPosSensor.y(),residualStereo/yErrorStereo);
            	    
            	    residualYvsUPos.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial);
            	    pullYvsUPos.get(sensorAxialName).fill(axialExtrapPosSensor.x(),residualAxial/yErrorAxial);
            	    residualYvsUPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo);
            	    pullYvsUPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x(),residualStereo/yErrorStereo);
        	    }
        	    
        	    if((Math.abs(residualAxial) < this.nSig * yErrorAxial)){
        	    	numberOfTracksWithHitOnMissingLayerChannel.get(sensorAxialName).fill(chanAxial);
        	    	numberOfTracksWithHitOnMissingLayerY.get(sensorAxialName).fill(axialExtrapPosSensor.x());
        	    	numberOfTracksWithHitOnMissingLayerP.get(sensorAxialName).fill(trackP);
        	    	if(q < 0){
        	    		numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorAxialName).fill(chanAxial);
            	    	numberOfTracksWithHitOnMissingLayerYEle.get(sensorAxialName).fill(axialExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerPEle.get(sensorAxialName).fill(trackP);
        	    	}
        	    	else{
        	    		numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorAxialName).fill(chanAxial);
            	    	numberOfTracksWithHitOnMissingLayerYPos.get(sensorAxialName).fill(axialExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerPPos.get(sensorAxialName).fill(trackP);
        	    	}
        	    }
        	  //Check to see if the sensor of this hit is the same sensor you expect to see a stereo hit
        	    //Compute the residual between extrapolated track and hit position and fill histo
        	    //double residual = stereoExtrapPos.y() - hitPos[1];
        	    if((Math.abs(residualStereo) < this.nSig * yErrorStereo)){
        	    	numberOfTracksWithHitOnMissingLayerChannel.get(sensorStereoName).fill(chanStereo);
        	    	numberOfTracksWithHitOnMissingLayerY.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
        	    	numberOfTracksWithHitOnMissingLayerP.get(sensorStereoName).fill(trackP);
        	    	if(q < 0){
        	    		numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorStereoName).fill(chanStereo);
            	    	numberOfTracksWithHitOnMissingLayerYEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerPEle.get(sensorStereoName).fill(trackP);
        	    	}
        	    	else{
        	    		numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorStereoName).fill(chanStereo);
            	    	numberOfTracksWithHitOnMissingLayerYPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerPPos.get(sensorStereoName).fill(trackP);
        	    	}
        	    }
        	    
        	    
        	    
        	    
	    		
    	    	//Loop over all tracker hits in the event
    	    	/*boolean fillAxial = false;
    	    	boolean fillStereo = false;
    	    	for(SiTrackerHitStrip1D hit:stripHits){
    	    		if(fillAxial && fillStereo) break;
    	    		//Get the sensor and position of the hit
    	    		HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();
            	    double[] hitPos = hit.getPosition();
            	    Hep3Vector hitPosSensor = globalToSensor(toHep3(hitPos),sensor);
            	    //System.out.println(sensor.getName() + " hit position = " + hitPosSensor);
            	    //Check to see if the sensor of this hit is the same sensor you expect to see an axial hit
            	    if(sensorAxialName == sensor.getName() && !fillAxial){
            	    	//Compute the residual between extrapolated track and hit position and fill histo
            	    	//double residual = axialExtrapPos.y() - hitPos[1];
            	    	double residual = axialExtrapPosSensor.x() - hitPosSensor.x();
            	    	residualY.get(sensorAxialName).fill(residual);
            	    	pullY.get(sensorAxialName).fill(residual/yErrorAxial);
            	    	//System.out.println(sensor.getName() + " hit position = " + hitPosSensor + "  Extrap Pos " + axialExtrapPosSensor);
            	    	//System.out.println("Axial residual = " + residual);
            	    	//Check to see if the residual is within desired error values
            	    	//If it is, then fill numerator of efficiency histos
            	    	residualYerrorYDiff.get(sensorAxialName).fill(Math.abs(residual) - yErrorAxial);
            	    	if((Math.abs(residual) > this.nSig * yErrorAxial)) continue;
            	    	numberOfTracksWithHitOnMissingLayerChannel.get(sensorAxialName).fill(chanAxial);
            	    	numberOfTracksWithHitOnMissingLayerY.get(sensorAxialName).fill(axialExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerP.get(sensorAxialName).fill(trackP);
            	    	if(q < 0){
            	    		numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorAxialName).fill(chanAxial);
                	    	numberOfTracksWithHitOnMissingLayerYEle.get(sensorAxialName).fill(axialExtrapPosSensor.x());
                	    	numberOfTracksWithHitOnMissingLayerPEle.get(sensorAxialName).fill(trackP);
            	    	}
            	    	else{
            	    		numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorAxialName).fill(chanAxial);
                	    	numberOfTracksWithHitOnMissingLayerYPos.get(sensorAxialName).fill(axialExtrapPosSensor.x());
                	    	numberOfTracksWithHitOnMissingLayerPPos.get(sensorAxialName).fill(trackP);
            	    	}
            	    	fillAxial = true;
            	    }
            	  //Check to see if the sensor of this hit is the same sensor you expect to see a stereo hit
            	    if(sensorStereoName == sensor.getName() && !fillStereo){
            	    	//Compute the residual between extrapolated track and hit position and fill histo
            	    	//double residual = stereoExtrapPos.y() - hitPos[1];
            	    	double residual = stereoExtrapPosSensor.x() - hitPosSensor.x();
            	    	residualY.get(sensorStereoName).fill(residual);
            	    	//System.out.println(sensor.getName() + " hit position = " + hitPosSensor + "  Extrap Pos " + stereoExtrapPosSensor);
            	    	pullY.get(sensorStereoName).fill(residual/yErrorStereo);
            	    	//Check to see if the residual is within desired error values
            	    	//If it is, then fill numerator of efficiency histos
            	    	residualYerrorYDiff.get(sensorStereoName).fill(Math.abs(residual) - yErrorStereo);
            	    	if((Math.abs(residual) > this.nSig * yErrorStereo)) continue;
            	    	numberOfTracksWithHitOnMissingLayerChannel.get(sensorStereoName).fill(chanStereo);
            	    	numberOfTracksWithHitOnMissingLayerY.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
            	    	numberOfTracksWithHitOnMissingLayerP.get(sensorStereoName).fill(trackP);
            	    	if(q < 0){
            	    		numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorStereoName).fill(chanStereo);
                	    	numberOfTracksWithHitOnMissingLayerYEle.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
                	    	numberOfTracksWithHitOnMissingLayerPEle.get(sensorStereoName).fill(trackP);
            	    	}
            	    	else{
            	    		numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorStereoName).fill(chanStereo);
                	    	numberOfTracksWithHitOnMissingLayerYPos.get(sensorStereoName).fill(stereoExtrapPosSensor.x());
                	    	numberOfTracksWithHitOnMissingLayerPPos.get(sensorStereoName).fill(trackP);
            	    	}
            	    	fillStereo = true;
            	    }
    	    	}*/
            }
        //}
    }
    
    private double findWeight(double y, double yErr, HpsSiSensor sensor){
    	double readoutPitch = sensor.getReadoutStripPitch();
    	int nChan = sensor.getNumberOfChannels();
    	boolean firstChan = firstChanIsEdge(sensor);
    	double height = readoutPitch * nChan;
    	double distanceToEdge = 0;
    	if(firstChan){
    		distanceToEdge = height/2 - y;
    	}
    	else{
    		distanceToEdge = height/2 + y;
    	}
    	double nSig = distanceToEdge/yErr;
    	//System.out.println(sensor.getName() + "  " + firstChan + "  Sigma " + nSig + "  Gauss " + computeGaussInt(nSig,1000) + "  Distance " + distanceToEdge);
    	return computeGaussInt(nSig,1000);
    }
    
    private double computeGaussInt(double nSig, int nSteps){
    	double mean = 0;
    	double sigma = 1;
    	double dx = sigma*nSig/(double) nSteps;
    	double integral = 0;
    	for(int i = 0; i < nSteps; i++){
    		double x = dx*(i+0.5) + mean;
    		integral += dx * Gauss(x,mean,sigma);
    	}
    	return integral + 0.5;
    }
    
    private double Gauss(double x, double mean, double sigma){
    	return 1/(Math.sqrt(2*Math.PI*Math.pow(sigma,2)))*Math.exp(-Math.pow(x-mean,2)/(2*Math.pow(sigma,2)));
    }
    
    private boolean firstChanIsEdge(HpsSiSensor sensor){
    	int layer = (sensor.getLayerNumber() + 1)/2;
    	if(layer > 0 && layer < 4){
    		if(sensor.isAxial()) return false;
    		else return true;
    	}
    	else{
    		if(sensor.isAxial()){
    			if(sensor.getSide().matches("ELECTRON")) return false;
    			else return true;
    		}
    		else{
    			if(!sensor.getSide().matches("ELECTRON")) return false;
    			else return true;
    		}
    	}
    }
    
    private boolean isGBL(Track track){
    	return track.getType() > 25;
    }
    
    private Hep3Vector globalToSensor(Hep3Vector trkpos, HpsSiSensor sensor){
        SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
    	if(electrodes == null){
    		electrodes = sensor.getReadoutElectrodes(ChargeCarrier.ELECTRON);
    		System.out.println("Charge Carrier is NULL");
    	}
        return electrodes.getGlobalToLocal().transformed(trkpos);
    }
    
    private Hep3Vector getHitPos(Track track, int unusedLay){
    	if(unusedLay == 1){
        	return toHep3(track.getTrackStates().get(0).getReferencePoint());
        }
    	List<TrackerHit> hits = track.getTrackerHits();
    	for(TrackerHit hit:hits){
    		HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();  	            
            int layer = (sensor.getLayerNumber() + 1)/2;
            if(unusedLay - 1 == layer){
                return toHep3(hit.getPosition());
            }
    	}
    	return null;
    }
    
    private TrackState getTrackState(Track track, int unusedLay){
    	int layer = -1;
    	boolean isTop = track.getTrackStates().get(0).getTanLambda() > 0;
	    if(unusedLay == 1){
	    	return track.getTrackStates().get(0);
	    }
	    else{
	    	layer = unusedLay - 1;
	    }
	    HpsSiSensor sensorHole = getSensor(track,layer,isTop,true);
	    HpsSiSensor sensorSlot = getSensor(track,layer,isTop,false);
	    TrackState tState = TrackStateUtils.getTrackStateAtSensor(track,sensorHole.getMillepedeId());
	    if(tState == null){
	    	tState = TrackStateUtils.getTrackStateAtSensor(track,sensorSlot.getMillepedeId());
	    }
	    return tState;
    }
    
    private int getChan(Hep3Vector pos, HpsSiSensor sensor){
    	double readoutPitch = sensor.getReadoutStripPitch();
    	int nChan = sensor.getNumberOfChannels();
    	double height = readoutPitch * nChan;
    	return (int) ((height/2-pos.x())/readoutPitch);
    }
    
    private Hep3Vector toHep3(double[] arr) {
    	return new BasicHep3Vector(arr[0], arr[1], arr[2]);
    }
    
    private HpsSiSensor getSensor(Track track, int layer, boolean isAxial, boolean isHole) {
    	double tanLambda = track.getTrackStates().get(0).getTanLambda();
    	for(HpsSiSensor sensor: sensors){
    		int senselayer = (sensor.getLayerNumber() + 1)/2;
    		if(senselayer != layer) continue;
    		if((tanLambda > 0 && !sensor.isTopLayer()) || (tanLambda < 0 && sensor.isTopLayer())) continue;
    		if((isAxial && !sensor.isAxial()) || (!isAxial && sensor.isAxial())) continue;
    		if(layer < 4 && layer > 0){
    			return sensor;
    		}
    		else{
    		    if((!sensor.getSide().matches("ELECTRON") && isHole) || (sensor.getSide().matches("ELECTRON") && !isHole)) continue;
    		    return sensor;
    		}
    	}
    	return null;
    }
    
    /*public static Hep3Vector extrapolateTrackPosition(TrackState track, double startPositionX, double endPositionX, double stepSize, FieldMap fieldMap) {

        // Start by extrapolating the track to the approximate point where the
        // fringe field begins.
        Hep3Vector currentPosition = TrackUtils.extrapolateHelixToXPlane(track, startPositionX);
        // System.out.println("Track position at start of fringe: " +
        // currentPosition.toString());

        // Get the HelicalTrackFit object associated with the track. This will
        // be used to calculate the path length to the start of the fringe and
        // to find the initial momentum of the track.
        HelicalTrackFit helicalTrackFit = TrackUtils.getHTF(track);

        // Calculate the path length to the start of the fringe field.
        double pathToStart = HelixUtils.PathToXPlane(helicalTrackFit, startPositionX, 0., 0).get(0);

        // Get the momentum of the track and calculate the magnitude. The
        // momentum can be calculate using the track curvature and magnetic
        // field strength in the middle of the analyzing magnet.
        // FIXME: The position of the middle of the analyzing magnet should
        // be retrieved from the compact description.
        double bFieldY = fieldMap.getField(new BasicHep3Vector(0, 0, 500.0)).y();
        double p = Math.abs(helicalTrackFit.p(bFieldY));

        // Get a unit vector giving the track direction at the start of the of
        // the fringe field
        Hep3Vector helixDirection = HelixUtils.Direction(helicalTrackFit, pathToStart);
        // Calculate the momentum vector at the start of the fringe field
        Hep3Vector currentMomentum = VecOp.mult(p, helixDirection);
        // System.out.println("Track momentum vector: " +
        // currentMomentum.toString());

        // Get the charge of the track.
        double q = Math.signum(track.getOmega());
        // HACK: LCSim doesn't deal well with negative fields so they are
        // turned to positive for tracking purposes. As a result,
        // the charge calculated using the B-field, will be wrong
        // when the field is negative and needs to be flipped.
        if (bFieldY < 0)
            q = q * (-1);

        // Swim the track through the B-field until the end point is reached.
        // The position of the track will be incremented according to the step
        // size up to ~90% of the final position. At this point, a finer
        // track size will be used.
        boolean stepSizeChange = false;
        while (currentPosition.x() < endPositionX) {

            // The field map coordinates are in the detector frame so the
            // extrapolated track position needs to be transformed from the
            // track frame to detector.
            Hep3Vector currentPositionDet = CoordinateTransformations.transformVectorToDetector(currentPosition);

            // Get the field at the current position along the track.
            bFieldY = fieldMap.getField(currentPositionDet).y();
            // System.out.println("Field along y (z in detector): " + bField);

            // Get a tracjectory (Helix or Line objects) created with the
            // track parameters at the current position.
            Trajectory trajectory = TrackUtils.getTrajectory(currentMomentum, new org.lcsim.spacegeom.SpacePoint(currentPosition), q, bFieldY);

            // Using the new trajectory, extrapolated the track by a step and
            // update the extrapolated position.
            currentPosition = trajectory.getPointAtDistance(stepSize);
            // System.out.println("Current position: " + ((Hep3Vector)
            // currentPosition).toString());

            // Calculate the momentum vector at the new position. This will
            // be used when creating the trajectory that will be used to
            // extrapolate the track in the next iteration.
            currentMomentum = VecOp.mult(currentMomentum.magnitude(), trajectory.getUnitTangentAtLength(stepSize));

            // If the position of the track along X (or z in the detector frame)
            // is at 90% of the total distance, reduce the step size.
            if (currentPosition.x() / endPositionX > .80 && !stepSizeChange) {
                stepSize /= 10;
                // System.out.println("Changing step size: " + stepSize);
                stepSizeChange = true;
            }
        }
        
        return CoordinateTransformations.transformVectorToDetector(currentPosition);
    }*/
    
    private double[] computeExtrapErrorY(Track track, TrackState tState, HpsSiSensor sensor, Hep3Vector hitPos){
    	// Extract the corrections to the track parameters and the covariance matrix from the GBL trajectory
    	//System.out.println(sensor.getName());
    	Hep3Vector sensorPos = sensor.getGeometry().getPosition();
    	double bfac = Constants.fieldConversion * bfield;
    	double[] cov = tState.getCovMatrix();
    	HelicalTrackFit htf = TrackUtils.getHTF(tState);
    	//System.out.println("HTF " + htf);
    	//SymMatrix locCov = buildCovMatrix(cov);
    	SymmetricMatrix LocCov = new SymmetricMatrix(5,cov,true);
    	Matrix locCov = new Matrix(5,5);
    	for(int i = 0; i < 5; i++){
    		for(int j = 0; j < 5; j++){
    			locCov.set(i, j, LocCov.e(i, j));
    		}
    	}
    	
    	//System.out.println("Local Covariance " + locCov);
    	
    	int id = 1;
    	GBLStripClusterData strip = new GBLStripClusterData(id);
    	Hep3Vector a = strip.getU();
        Hep3Vector b = strip.getV();
        Hep3Vector c = strip.getW();
        
        //System.out.println("U " + a);
        //System.out.println("V " + b);
        //System.out.println("W " + c);
    	//GBLStripClusterData strip = hits.get(istrip);
    	//FittedGblTrajectory f = fit(hits,bfac,0);
    	
    	//BasicMatrix jacPointToPoint = new Matrix(5, 5);
        //jacPointToPoint.UnitMatrix();
        
        //int istrip = 0;
    	//GBLStripClusterData strip = new GBLStripClusterData(istrip);
    	
    	Hep3Vector point_on_plane = sensor.getGeometry().getPosition();
    	Hep3Vector pointInTrackingFrame = CoordinateTransformations.transformVectorToTracking(point_on_plane);
        Hep3Vector w = sensor.getGeometry().getLocalToGlobal().rotated(new BasicHep3Vector(0, 0, 1));
        Hep3Vector wInTrackingFrame = CoordinateTransformations.transformVectorToTracking(w);
        wInTrackingFrame = VecOp.unit(wInTrackingFrame);

        // get helix intercept: in global tracking coordinates
        //double step = HelixUtils.PathToXPlane(htf, pointInTrackingFrame.x(), 0., 0).get(0);
        //System.out.println("Step " + step);
        //Hep3Vector theInt = TrackUtils.getHelixPlaneIntercept(htf, wInTrackingFrame, pointInTrackingFrame, bfield, s_origin);


        // get helix intercept: in global tracking coordinates
        //double s_origin = HelixUtils.PathToXPlane(htf, pointInTrackingFrame.x(), 0., 0).get(0);
        //Hep3Vector theInt = TrackUtils.getHelixPlaneIntercept(htf, wInTrackingFrame, pointInTrackingFrame, bfield, s_origin);
        
    	// get measurement frame unit vectors
        //Hep3Vector u = new BasicHep3Vector(0,0,1);
        //Hep3Vector v = new BasicHep3Vector(0,-1,0);
        Hep3Vector u = new BasicHep3Vector(0,1,0);
        Hep3Vector v = new BasicHep3Vector(-1,0,0);
        

        // Measurement direction (perpendicular and parallel to strip direction)
        Matrix mDir = new Matrix(2, 3);
        mDir.set(0, 0, u.x());
        mDir.set(0, 1, u.y());
        mDir.set(0, 2, u.z());
        mDir.set(1, 0, v.x());
        mDir.set(1, 1, v.y());
        mDir.set(1, 2, v.z());
        
        Matrix mDirT = mDir.copy().transpose();


        // Track direction
        double sinLambda = sin(htf.slope());// ->GetLambda());
        double cosLambda = sqrt(1.0 - sinLambda * sinLambda);
        double phi = TrackUtils.getPhi(tState, w);
        double sinPhi = sin(phi);// ->GetPhi());
        double cosPhi = sqrt(1.0 - sinPhi * sinPhi);
        
     // Track direction in curvilinear frame (U,V,T)
        // U = Z x T / |Z x T|, V = T x U
        Matrix uvDir = new Matrix(2, 3);
        uvDir.set(0, 0, -sinPhi);
        uvDir.set(0, 1, cosPhi);
        uvDir.set(0, 2, 0.);
        uvDir.set(1, 0, -sinLambda * cosPhi);
        uvDir.set(1, 1, -sinLambda * sinPhi);
        uvDir.set(1, 2, cosLambda);
        
        /*Vector meas = new Vector(2);
        // double uRes = strip->GetUmeas() - strip->GetTrackPos().x(); // how can this be correct?
        double uRes = strip.getMeas() - strip.getTrackPos().x();
        meas.set(0, uRes);
        meas.set(1, 0.);
        Vector measErr = new Vector(2);
        measErr.set(0, strip.getMeasErr());
        measErr.set(1, 0.);
        Vector measPrec = new Vector(2);
        measPrec.set(0, 1.0 / (measErr.get(0) * measErr.get(0)));
        measPrec.set(1, 0.);*/
        

        // projection from measurement to local (curvilinear uv) directions (duv/dm)
        Matrix proM2l = uvDir.times(mDirT);

        // projection from local (uv) to measurement directions (dm/duv)
        Matrix proL2m = proM2l.copy();
        proL2m = proL2m.inverse();
        
        double step1 =  HelixUtils.PathToXPlane(htf,hitPos.x(),0,0).get(0);
        double step2 =  HelixUtils.PathToXPlane(htf,sensorPos.z(),0,0).get(0);
        double step = step2 - step1;
        
        BasicMatrix jacPointToPoint = GblUtils.gblSimpleJacobianLambdaPhi(step, cosLambda, abs(bfac));
    	
    	//System.out.println("Local Cov " + locCov);
        //System.out.println("Jacobian Simple " + jacPointToPoint);
       
        //double cosLambda = sqrt(1.0 - sin(tState.getTanLambda()) * sin(tState.getTanLambda()));  
    	//double cosLambda = Math.cos(Math.atan(tState.getTanLambda()));  
        //double step1 =  HelixUtils.PathToXPlane(htf,hitPos.x(),0,0).get(0);
        //double step2 =  HelixUtils.PathToXPlane(htf,sensorPos.z(),0,0).get(0);
        //double step = step2 - step1;
    	//BasicMatrix jacPointToPoint = GblUtils.gblSimpleJacobianLambdaPhi(step, cosLambda, Math.abs(bfac));
    	Matrix jacobian = new Matrix(5,5);
    	for(int i = 0; i < 5; i++){
    		for(int j = 0; j < 5; j++){
    			jacobian.set(i, j, jacPointToPoint.e(i,j));
    		}
    	}
    	
    	//Matrix ClToPerJac = GblUtils.getCLToPerigeeJacobian(htf,new HpsHelicalTrackFit(TrackUtils.getHTF(tState)),bfield);
    	//Matrix ClToPerJac = getCLToPerigeeJacobian(TrackUtils.getHTF(track.getTrackStates().get(0)),new HpsHelicalTrackFit(TrackUtils.getHTF(tState)),bfield);
    	Matrix ClToPerJac = GblUtils.getCLToPerigeeJacobian(htf,new HpsHelicalTrackFit(TrackUtils.getHTF(tState)),bfield);
    	Matrix PerToClJac = ClToPerJac.inverse();
    	Matrix MsCov = jacobian.times(PerToClJac.times(locCov.times(PerToClJac.transpose())).times(jacobian.transpose()));
    	Matrix helixCovariance = ClToPerJac.times(MsCov.times(ClToPerJac.transpose()));
    	//Matrix helixCovariance = locCov;
    	Matrix MsCov2 = new Matrix(3,3);
    	MsCov2.set(0,0,MsCov.get(3,3));
    	MsCov2.set(0,1,MsCov.get(3,4));
    	MsCov2.set(1,0,MsCov.get(4,3));
    	MsCov2.set(1,1,MsCov.get(4,4));
    	
    	//System.out.println("PerToClJac " + PerToClJac);
    	
    	//System.out.println("MsCov " + MsCov);
        
    	
    	//Matrix measMsCov = ClToPerJac.times(MsCov.times(ClToPerJac.transpose()));
    	//Matrix measMsCov = proL2m.times(MsCov2.times(proL2m.transpose()));
    	
    	SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
        Matrix rot = Hep3ToMatrix(electrodes.getGlobalToLocal().getRotation().getRotationMatrix());
        Matrix measMsCov = rot.times(MsCov2.times(rot.transpose()));
    	
    	//System.out.println("Global Y Error " + Math.sqrt(MsCov2.get(1,1)) + "  Y error = " + Math.sqrt(measMsCov.get(0, 0)));
    	//System.out.println("");
    	//System.out.println("");
    	
    	double d0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.dcaIndex,HelicalTrackFit.dcaIndex));
		double z0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.z0Index,HelicalTrackFit.z0Index));
		double tanlambda_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.slopeIndex,HelicalTrackFit.slopeIndex));
		double phi0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.phi0Index,HelicalTrackFit.phi0Index));
		double omega_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.curvatureIndex,HelicalTrackFit.curvatureIndex));
		
    	String sensorName = sensor.getName();
		D0_err.get(sensorName).fill(d0_err);
		Z0_err.get(sensorName).fill(z0_err);
		Tanlambda_err.get(sensorName).fill(tanlambda_err);
		Phi0_err.get(sensorName).fill(phi0_err);
		Omega_err.get(sensorName).fill(omega_err);
    	
    	//return new double[]{Math.sqrt(measMsCov.get(0, 0)),Math.sqrt(measMsCov.get(1,1))};
		return new double[]{Math.sqrt(measMsCov.get(0, 0)),Math.sqrt(measMsCov.get(1,1))};
		
    	/*System.out.println("Step " + step + "  Step 1 " + step1 + "  Step 2 " + step2);
    	
    	System.out.println("Jacobian " + jacobian);
    	
    	//Matrix perToCl = getPerToCl(htf, bfield);
    	//Matrix CltoPer = perToCl.inverse();
    	
    	//Matrix newjacobian = CltoPer.times(jacobian.times(CltoPer.transpose()));
    	
    	Matrix CltoPer = GblUtils.getCLToPerigeeJacobian(htf,new HpsHelicalTrackFit(TrackUtils.getHTF(track.getTrackStates().get(0))), bfield);
    	Matrix PertoCl = CltoPer.inverse();
    	
    	Matrix newloccov = PertoCl.times(locCov.times(PertoCl.transpose()));
    	
    	//Matrix newjacobian = CltoPer.times(jacobian.times(CltoPer.inverse()));
    	
        Matrix oldhelixCovariance = jacobian.times((newloccov).times(jacobian.transpose()));
        Matrix helixCovariance = CltoPer.times(oldhelixCovariance.times(CltoPer.transpose()));
    	//Matrix helixCovariance = newjacobian.times((locCov).times(newjacobian.transpose()));
        
        Matrix uvDir = new Matrix(2,3);
        uvDir.set(0,0,-sin(tState.getPhi()));
        uvDir.set(0,1,Math.cos(tState.getPhi()));
        uvDir.set(1,0,-sin(Math.atan(tState.getTanLambda()))*Math.cos(tState.getPhi()));
        uvDir.set(1,1,-sin(Math.atan(tState.getTanLambda()))*sin(tState.getPhi()));
        uvDir.set(1,2,Math.cos(Math.atan(tState.getTanLambda())));
        
        Matrix mDir = new Matrix(1,2);
        mDir.set(0,0,1);
        mDir.set(0,1,1);
        
        //Matrix proL2M = uvDir.times(mDir).inverse();
        
        //Matrix measMsCo = proL2M.times(helixCovariance.times(proL2M.transpose()));
        
        //System.out.println("Y error = " + Math.sqrt(measMsCo.get(0, 0)));
		System.out.println("");
		System.out.println("");
        
        //return new double[]{Math.sqrt(measMsCo.get(0, 0)),Math.sqrt(measMsCo.get(1,1))};
        return new double[]{1,1};*/
        
        /*double phi = htf.phi0() - step / htf.R();
        double x = htf.xc() -htf.R() * Math.sin(phi);
        double y = htf.yc() + htf.R() * Math.cos(phi);
        double d0 = Math.signum(tState.getD0()) * Math.sqrt(Math.pow(x,2)+Math.pow(y,2));
        //double z0 = 0;
        Matrix helixJacobian = new Matrix(3,3);
        helixJacobian.set(0,0,-Math.sin(phi));
        helixJacobian.set(1,0,-Math.cos(phi));
        helixJacobian.set(0,1,-Math.cos(phi)*d0);
        helixJacobian.set(1,1,-Math.sin(phi)*d0);
        helixJacobian.set(2,2,1);
        
        Matrix perCov = new Matrix(3,3);
        perCov.set(0,0,helixCovariance.get(0,0));
        perCov.set(0,1,helixCovariance.get(0,1));
        perCov.set(0,2,helixCovariance.get(0,3));
        perCov.set(1,0,helixCovariance.get(1,0));
        perCov.set(1,1,helixCovariance.get(1,1));
        perCov.set(1,2,helixCovariance.get(1,3));
        perCov.set(2,0,helixCovariance.get(3,0));
        perCov.set(2,1,helixCovariance.get(3,1));
        perCov.set(2,2,helixCovariance.get(3,3));
        
        Matrix trkToDet = Hep3ToMatrix(initializeInverse().getRotation().getRotationMatrix());
        
        Matrix trackCov = trkToDet.times((helixJacobian.times(perCov.times(helixJacobian.transpose()))));
              
        System.out.println("Helix Covariance " + helixCovariance);
        
        SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
        Matrix rot = Hep3ToMatrix(electrodes.getGlobalToLocal().getRotation().getRotationMatrix());
        
        System.out.println("Rotation Matrix " + rot);
        
        Matrix sensorCov = rot.times(trackCov.times(rot.transpose()));
        
        System.out.println("Sensor Cov " + sensorCov);
        
    	double de0 = d0;
		double tanlambda = tState.getTanLambda();
		double z0 = tState.getZ0() + tState.getReferencePoint()[2] + + step*tanlambda;
		double phi0 = phi;
		double omega = tState.getOmega();
    	
    	double d0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.dcaIndex,HelicalTrackFit.dcaIndex));
		double z0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.z0Index,HelicalTrackFit.z0Index));
		double tanlambda_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.slopeIndex,HelicalTrackFit.slopeIndex));
		double phi0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.phi0Index,HelicalTrackFit.phi0Index));
		double omega_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.curvatureIndex,HelicalTrackFit.curvatureIndex));
		
    	String sensorName = sensor.getName();
    	D0.get(sensorName).fill(de0);
		Z0.get(sensorName).fill(z0);
		Tanlambda.get(sensorName).fill(tanlambda);
		Phi0.get(sensorName).fill(phi0);
		Omega.get(sensorName).fill(omega);
		
		D0_err.get(sensorName).fill(d0_err);
		Z0_err.get(sensorName).fill(z0_err);
		Tanlambda_err.get(sensorName).fill(tanlambda_err);
		Phi0_err.get(sensorName).fill(phi0_err);
		Omega_err.get(sensorName).fill(omega_err);
		
		System.out.println("Y error = " + Math.sqrt(sensorCov.get(0, 0)));
		System.out.println("");
		System.out.println("");*/
        //return new double[]{Math.sqrt(sensorCov.get(0, 0)),Math.sqrt(sensorCov.get(1,1))};

    	//return 1.0;
    }
    
    /**
     * Calculate the Jacobian from Curvilinear to Perigee frame.
     * 
     * @param helicalTrackFit - original helix
     * @param helicalTrackFitAtIPCorrected - corrected helix at this point
     * @param bfield - magnitude of B-field
     * @return the Jacobian matrix from Curvilinear to Perigee frame
     */
    public static Matrix getCLToPerigeeJacobian(HelicalTrackFit helicalTrackFit,
            HpsHelicalTrackFit helicalTrackFitAtIPCorrected, double bfield) {

    	  double lambda_gbl = Math.atan(helicalTrackFitAtIPCorrected.slope());
          double qOverP_gbl = helicalTrackFitAtIPCorrected.curvature();
          //This part is taken from: // Strandlie, Wittek, NIMA 566, 2006 
    	  Matrix covariance_gbl = new Matrix(5, 5);
          //helpers
    	  double Bz = -Constants.fieldConversion * Math.abs(bfield); // TODO sign convention and should it be
          //it scaled from Telsa? 
          double p = Math.abs(1 / qOverP_gbl);
          double q = Math.signum(qOverP_gbl);
          double tanLambda = Math.tan(lambda_gbl);
          double cosLambda = Math.cos(lambda_gbl); // 
          Hep3Vector B = new BasicHep3Vector(0, 0, Bz); // TODO sign convention? 
          Hep3Vector H = new BasicHep3Vector(0, 0, 1); 
          Hep3Vector T = HelixUtils.Direction(helicalTrackFit, 0.); 
          Hep3Vector HcrossT = VecOp.cross(H, T); 
          double alpha = HcrossT.magnitude(); // this should be Bvec cross TrackDir/|B| 
          double Q = Bz * q / p; 
          Hep3Vector Z = new BasicHep3Vector(0, 0, 1); 
          Hep3Vector J = VecOp.mult(1. / VecOp.cross(T, Z).magnitude(), VecOp.cross(T, Z));
          Hep3Vector K = Z; 
          Hep3Vector U = VecOp.mult(-1, J); 
          Hep3Vector V = VecOp.cross(T, U); 
          Hep3Vector I = VecOp.cross(J, K); 
          Hep3Vector N = VecOp.mult(1 / alpha, VecOp.cross(H, T)); 
          double UdotI = VecOp.dot(U, I);
          double NdotV = VecOp.dot(N, V); 
          double NdotU = VecOp.dot(N, U); 
          double TdotI = VecOp.dot(T, I); 
          double VdotI = VecOp.dot(V, I);
          double VdotK = VecOp.dot(V, K); 
          covariance_gbl.set(HelicalTrackFit.dcaIndex,
          FittedGblTrajectory.GBLPARIDX.XT.getValue(), VdotK / TdotI); 
          covariance_gbl.set(HelicalTrackFit.phi0Index,
          FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 1); 
          covariance_gbl.set(HelicalTrackFit.phi0Index,
          FittedGblTrajectory.GBLPARIDX.XT.getValue(), -alpha * Q * UdotI * NdotU / (cosLambda * TdotI));
          covariance_gbl.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -alpha * Q * VdotI
          * NdotU / (cosLambda * TdotI)); 
          covariance_gbl.set(HelicalTrackFit.curvatureIndex,
          FittedGblTrajectory.GBLPARIDX.QOVERP.getValue(), -1 * Bz / cosLambda); //
          covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 0);
          covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), -1 * q *
          Bz * tanLambda / (p * cosLambda)); 
          covariance_gbl.set(HelicalTrackFit.curvatureIndex,
          FittedGblTrajectory.GBLPARIDX.XT.getValue(), q * Bz * alpha * Q * tanLambda * UdotI * NdotV / (p * cosLambda
          * TdotI)); 
          covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(), q
          * Bz * alpha * Q * tanLambda * VdotI * NdotV / (p * cosLambda * TdotI));
          covariance_gbl.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -1 / TdotI);
          covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), -1);
          covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), alpha * Q * UdotI
          * NdotV / TdotI); 
          covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(),
          alpha * Q * VdotI * NdotV / TdotI); 
          covariance_gbl.print(15, 13);
          
          return covariance_gbl;
         

        // Sho's magic below

        // Use projection matrix
        // TODO should this not be the corrected helix?
        /*Hep3Matrix perToClPrj = getPerToClPrj(helicalTrackFit);
        Hep3Matrix clToPerPrj = VecOp.inverse(perToClPrj);
        double C_gbl = helicalTrackFitAtIPCorrected.curvature();
        double lambda_gbl = Math.atan(helicalTrackFitAtIPCorrected.slope());
        double qOverP_gbl = helicalTrackFitAtIPCorrected.curvature()
                / (Constants.fieldConversion * Math.abs(bfield) * Math.sqrt(1 + Math.pow(
                        helicalTrackFitAtIPCorrected.slope(), 2)));

        Matrix jacobian = new Matrix(5, 5);
        jacobian.set(HelicalTrackFit.dcaIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), -clToPerPrj.e(1, 0));
        jacobian.set(HelicalTrackFit.dcaIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -clToPerPrj.e(1, 1));
        jacobian.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 1.0);
        jacobian.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), clToPerPrj.e(0, 1) * C_gbl);
        jacobian.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.QOVERP.getValue(),
                Constants.fieldConversion * Math.abs(bfield) / Math.cos(lambda_gbl));
        jacobian.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(),
                Constants.fieldConversion * Math.abs(bfield) * qOverP_gbl * Math.tan(lambda_gbl) / Math.cos(lambda_gbl));
        jacobian.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.XT.getValue(), clToPerPrj.e(2, 0));
        jacobian.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), clToPerPrj.e(2, 1));
        jacobian.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(),
                Math.pow(Math.cos(lambda_gbl), -2.0));

        return jacobian;*/
    }
    
    /*private double[] computeExtrapErrorY(Track track, TrackState tState, HpsSiSensor sensor, Hep3Vector hitPos){
    	// Extract the corrections to the track parameters and the covariance matrix from the GBL trajectory
    	Hep3Vector sensorPos = sensor.getGeometry().getPosition();
    	double bfac = Constants.fieldConversion * bfield;
    	double[] cov = tState.getCovMatrix();
    	HelicalTrackFit htf = TrackUtils.getHTF(tState);
    	System.out.println("HTF " + htf);
    	//SymMatrix locCov = buildCovMatrix(cov);
    	SymmetricMatrix LocCov = new SymmetricMatrix(5,cov,true);
    	Matrix locCov = new Matrix(5,5);
    	for(int i = 0; i < 5; i++){
    		for(int j = 0; j < 5; j++){
    			locCov.set(i, j, LocCov.e(i, j));
    		}
    	}
    	
    	System.out.println("Local Cov " + locCov);
       
        //double cosLambda = sqrt(1.0 - sin(tState.getTanLambda()) * sin(tState.getTanLambda()));  
    	double cosLambda = Math.cos(Math.atan(tState.getTanLambda()));  
        double step1 =  HelixUtils.PathToXPlane(htf,hitPos.x(),0,0).get(0);
        double step2 =  HelixUtils.PathToXPlane(htf,sensorPos.z(),0,0).get(0);
        double step = step2 - step1;
    	BasicMatrix jacPointToPoint = GblUtils.gblSimpleJacobianLambdaPhi(step, cosLambda, Math.abs(bfac));
    	Matrix jacobian = new Matrix(5,5);
    	for(int i = 0; i < 5; i++){
    		for(int j = 0; j < 5; j++){
    			jacobian.set(i, j, jacPointToPoint.e(i,j));
    		}
    	}
    	
    	System.out.println("Step " + step + "  Step 1 " + step1 + "  Step 2 " + step2);
    	
    	System.out.println("Jacobian " + jacobian);
    	
    	//Matrix perToCl = getPerToCl(htf, bfield);
    	//Matrix CltoPer = perToCl.inverse();
    	
    	//Matrix newjacobian = CltoPer.times(jacobian.times(CltoPer.transpose()));
    	
    	Matrix CltoPer = GblUtils.getCLToPerigeeJacobian(htf,new HpsHelicalTrackFit(TrackUtils.getHTF(track.getTrackStates().get(0))), bfield);
    	
    	Matrix newjacobian = CltoPer.times(jacobian.times(CltoPer.inverse()));
    	
        //Matrix helixCovariance = jacobian.times((locCov).times(jacobian.transpose()));
    	Matrix helixCovariance = newjacobian.times((locCov).times(newjacobian.transpose()));
        double phi = htf.phi0() - step / htf.R();
        double x = htf.xc() -htf.R() * Math.sin(phi);
        double y = htf.yc() + htf.R() * Math.cos(phi);
        double d0 = Math.signum(tState.getD0()) * Math.sqrt(Math.pow(x,2)+Math.pow(y,2));
        //double z0 = 0;
        Matrix helixJacobian = new Matrix(3,3);
        helixJacobian.set(0,0,-Math.sin(phi));
        helixJacobian.set(1,0,-Math.cos(phi));
        helixJacobian.set(0,1,-Math.cos(phi)*d0);
        helixJacobian.set(1,1,-Math.sin(phi)*d0);
        helixJacobian.set(2,2,1);
        
        Matrix perCov = new Matrix(3,3);
        perCov.set(0,0,helixCovariance.get(0,0));
        perCov.set(0,1,helixCovariance.get(0,1));
        perCov.set(0,2,helixCovariance.get(0,3));
        perCov.set(1,0,helixCovariance.get(1,0));
        perCov.set(1,1,helixCovariance.get(1,1));
        perCov.set(1,2,helixCovariance.get(1,3));
        perCov.set(2,0,helixCovariance.get(3,0));
        perCov.set(2,1,helixCovariance.get(3,1));
        perCov.set(2,2,helixCovariance.get(3,3));
        
        Matrix trkToDet = Hep3ToMatrix(initializeInverse().getRotation().getRotationMatrix());
        
        Matrix trackCov = trkToDet.times((helixJacobian.times(perCov.times(helixJacobian.transpose()))));
              
        System.out.println("Helix Covariance " + helixCovariance);
        
        SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
        Matrix rot = Hep3ToMatrix(electrodes.getGlobalToLocal().getRotation().getRotationMatrix());
        
        System.out.println("Rotation Matrix " + rot);
        
        Matrix sensorCov = rot.times(trackCov.times(rot.transpose()));
        
        System.out.println("Sensor Cov " + sensorCov);
        
    	double de0 = d0;
		double tanlambda = tState.getTanLambda();
		double z0 = tState.getZ0() + tState.getReferencePoint()[2] + + step*tanlambda;
		double phi0 = phi;
		double omega = tState.getOmega();
    	
    	double d0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.dcaIndex,HelicalTrackFit.dcaIndex));
		double z0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.z0Index,HelicalTrackFit.z0Index));
		double tanlambda_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.slopeIndex,HelicalTrackFit.slopeIndex));
		double phi0_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.phi0Index,HelicalTrackFit.phi0Index));
		double omega_err = Math.sqrt(helixCovariance.get(HelicalTrackFit.curvatureIndex,HelicalTrackFit.curvatureIndex));
		
    	String sensorName = sensor.getName();
    	D0.get(sensorName).fill(de0);
		Z0.get(sensorName).fill(z0);
		Tanlambda.get(sensorName).fill(tanlambda);
		Phi0.get(sensorName).fill(phi0);
		Omega.get(sensorName).fill(omega);
		
		D0_err.get(sensorName).fill(d0_err);
		Z0_err.get(sensorName).fill(z0_err);
		Tanlambda_err.get(sensorName).fill(tanlambda_err);
		Phi0_err.get(sensorName).fill(phi0_err);
		Omega_err.get(sensorName).fill(omega_err);
		
		System.out.println("Y error = " + Math.sqrt(sensorCov.get(0, 0)));
		System.out.println("");
		System.out.println("");
        return new double[]{Math.sqrt(sensorCov.get(0, 0)),Math.sqrt(sensorCov.get(1,1))};

    	//return 1.0;
    }*/
    
    /*private Matrix getPerToCl(HelicalTrackFit htf, double Bz) {
        Hep3Vector Z = new BasicHep3Vector(0, 0, 1);
        Hep3Vector T = HelixUtils.Direction(htf, 0.);
        Hep3Vector J = VecOp.mult(1. / VecOp.cross(T, Z).magnitude(), VecOp.cross(T, Z));
        Hep3Vector K = Z;
        Hep3Vector U = VecOp.mult(-1, J);
        Hep3Vector V = VecOp.cross(T, U);
        Hep3Vector I = VecOp.cross(J, K);
        double cosl = Math.cos(htf.slope());
        double phi = htf.phi0();
        double qOverp = htf.curvature()/Bz;
        double Q = -Bz*qOverp;
        Hep3Vector H = new BasicHep3Vector(0, 1, 0);
        double alpha = VecOp.cross(H,T).magnitude();
        Hep3Vector N = VecOp.mult(1/alpha,VecOp.cross(H,T));

        Matrix trans = new Matrix(5,5);     
        trans.set(0, 0, -Math.sin(phi)/Bz);
        trans.set(0, 1, qOverp/(Math.tan(phi)));
        trans.set(1, 1, -1);
        trans.set(1, 3, -alpha*Q*VecOp.dot(T,J)*VecOp.dot(V,N));
        trans.set(1, 4, -alpha*Q*VecOp.dot(T,K)*VecOp.dot(V,N));
        trans.set(2, 3, -alpha*Q*VecOp.dot(T,J)*VecOp.dot(U,N)/cosl);
        trans.set(2, 4, -alpha*Q*VecOp.dot(T,K)*VecOp.dot(U,N)/cosl);
        trans.set(3, 3, -1);
        trans.set(1, 3, VecOp.dot(V,K));
        
        /*trans.setElement(1, 1, VecOp.dot(T,I)*VecOp.dot(V,J));
        trans.setElement(1, 2, VecOp.dot(T,I)*VecOp.dot(V,K));
        trans.setElement(1, 3, -alpha*Q*VecOp.dot(T,J)*VecOp.dot(V,N));
        trans.setElement(1, 4, -alpha*Q*VecOp.dot(T,K)*VecOp.dot(V,N));
        trans.setElement(2, 1, VecOp.dot(T,I)*VecOp.dot(V,J)/cosl);
        trans.setElement(2, 2, VecOp.dot(T,I)*VecOp.dot(V,K)/cosl);
        trans.setElement(2, 3, -alpha*Q*VecOp.dot(T,J)*VecOp.dot(V,N)/cosl);
        trans.setElement(2, 4, -alpha*Q*VecOp.dot(T,K)*VecOp.dot(V,N)/cosl);
        trans.setElement(3, 3, VecOp.dot(U,J));
        trans.setElement(3, 4, VecOp.dot(U,K));
        trans.setElement(4, 3, VecOp.dot(V,J));
        trans.setElement(4, 4, VecOp.dot(V,K));
        return trans;
    }*/
    
    /*private Transform3D initialize() {
        BasicHep3Matrix tmp = new BasicHep3Matrix();
        tmp.setElement(0, 2, 1);
        tmp.setElement(1, 0, 1);
        tmp.setElement(2, 1, 1);
        return new Transform3D(new Rotation3D(tmp));
    }

    private Transform3D initializeInverse() {
        return initialize().inverse();
    }*/
    
    private Matrix Hep3ToMatrix(Hep3Matrix mat){
    	int Nrows = mat.getNRows();
    	int Ncolumns = mat.getNColumns();
    	Matrix matrix = new Matrix(Nrows,Ncolumns);
    	for(int i = 0; i < Nrows; i++){
    		for(int j = 0; j < Ncolumns; j++){	
    			matrix.set(i, j, mat.e(i,j));
    		}
    	}
    	return matrix;
    }
    
    public static Hep3Vector getLocationAtSensor(TrackState ts, HpsSiSensor sensor, double bfield) {
        if (ts == null || sensor == null)
            return null;
        if ((ts.getTanLambda() > 0 && sensor.isTopLayer()) || (ts.getTanLambda() < 0 && sensor.isBottomLayer()))
            return getLocationAtSensor(TrackUtils.getHTF(ts), sensor, bfield);
        return null;
    }

    public static Hep3Vector getLocationAtSensor(HelicalTrackFit htf, HpsSiSensor sensor, double bfield) {
        if (htf == null || sensor == null)
            return null;

        // get origin and normal of sensor, in global tracking coordinates
        Hep3Vector point_on_plane = sensor.getGeometry().getPosition();
        if (point_on_plane == null)
            return null;
        Hep3Vector pointInTrackingFrame = CoordinateTransformations.transformVectorToTracking(point_on_plane);
        Hep3Vector w = sensor.getGeometry().getLocalToGlobal().rotated(new BasicHep3Vector(0, 0, 1));
        Hep3Vector wInTrackingFrame = CoordinateTransformations.transformVectorToTracking(w);
        wInTrackingFrame = VecOp.unit(wInTrackingFrame);

        // get helix intercept: in global tracking coordinates
        double s_origin = HelixUtils.PathToXPlane(htf, pointInTrackingFrame.x(), 0., 0).get(0);
        Hep3Vector theInt = TrackUtils.getHelixPlaneIntercept(htf, wInTrackingFrame, pointInTrackingFrame, bfield, s_origin);

        // return in global detector coordinates
        if (theInt != null)
            return CoordinateTransformations.transformVectorToDetector(theInt);
        return null;
    }
    
    
    /*private SymMatrix buildCovMatrix(double[] cov){
    	SymMatrix mat = new SymMatrix(5);
    	SymmetricMatrix sym2 = new SymmetricMatrix(nBins, cov, cleanTridents);
    	mat.set(HelicalTrackFit.dcaIndex, HelicalTrackFit.dcaIndex, cov[0]);
    	mat.set(HelicalTrackFit.phi0Index, HelicalTrackFit.dcaIndex, cov[1]);
    	mat.set(HelicalTrackFit.phi0Index, HelicalTrackFit.phi0Index, cov[2]);
    	mat.set(HelicalTrackFit.curvatureIndex, HelicalTrackFit.dcaIndex, cov[3]);
    	mat.set(HelicalTrackFit.curvatureIndex, HelicalTrackFit.phi0Index, cov[4]);
    	mat.set(HelicalTrackFit.curvatureIndex, HelicalTrackFit.curvatureIndex, cov[5]);
    	mat.set(HelicalTrackFit.z0Index, HelicalTrackFit.dcaIndex, cov[6]);
    	mat.set(HelicalTrackFit.z0Index, HelicalTrackFit.phi0Index, cov[7]);
    	mat.set(HelicalTrackFit.z0Index, HelicalTrackFit.curvatureIndex, cov[8]);
    	mat.set(HelicalTrackFit.z0Index, HelicalTrackFit.z0Index, cov[9]);
    	mat.set(HelicalTrackFit.slopeIndex, HelicalTrackFit.dcaIndex, cov[10]);
    	mat.set(HelicalTrackFit.slopeIndex, HelicalTrackFit.phi0Index, cov[11]);
    	mat.set(HelicalTrackFit.slopeIndex, HelicalTrackFit.curvatureIndex, cov[12]);
    	mat.set(HelicalTrackFit.slopeIndex, HelicalTrackFit.z0Index, cov[13]);
    	mat.set(HelicalTrackFit.slopeIndex, HelicalTrackFit.slopeIndex, cov[14]);
    	return mat;
    }*/

    private int getUnusedSvtLayer(List<TrackerHit> stereoHits) {      
        int[] svtLayer = new int[6];
        
        // Loop over all of the stereo hits associated with the track
        for (TrackerHit stereoHit : stereoHits) {
            
            // Retrieve the sensor associated with one of the hits.  This will
            // be used to retrieve the layer number
            HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) stereoHit.getRawHits().get(0)).getDetectorElement();
            
            // Retrieve the layer number by using the sensor
            int layer = (sensor.getLayerNumber() + 1)/2;
           
            // If a hit is associated with that layer, increment its 
            // corresponding counter
            svtLayer[layer - 1]++;
        }
        
        // Loop through the layer counters and find which layer has not been
        // incremented i.e. is unused by the track
        for(int layer = 0; layer < svtLayer.length; layer++){
            if(svtLayer[layer] == 0) { 
                return (layer + 1);
            }
        }
        return -1;
    }
    
    private Pair<HpsSiSensor,Integer> isWithinSensorAcceptance(Track track, TrackState tState, int layer, boolean axial, Hep3Vector p,FieldMap fieldMap) {
  	   
        HpsSiSensor axialSensorHole = getSensor(track,layer,true,true);
        HpsSiSensor axialSensorSlot = getSensor(track,layer,true,false);
        HpsSiSensor stereoSensorHole = getSensor(track,layer,false,true);
        HpsSiSensor stereoSensorSlot = getSensor(track,layer,false,false);
        
    	Hep3Vector axialSensorHolePosition = axialSensorHole.getGeometry().getPosition();
    	Hep3Vector axialSensorSlotPosition = axialSensorSlot.getGeometry().getPosition();
        Hep3Vector stereoSensorHolePosition = stereoSensorHole.getGeometry().getPosition();
        Hep3Vector stereoSensorSlotPosition = stereoSensorSlot.getGeometry().getPosition();
        
        Hep3Vector endPosition = getHitPos(track, layer);
        
        //Hep3Vector axialTrackHolePos = extrapolateTrackPosition(tState,endPosition.x(),axialSensorHolePosition.z(),5,fieldMap);
        //Hep3Vector axialTrackSlotPos = extrapolateTrackPosition(tState,endPosition.x(),axialSensorSlotPosition.z(),5,fieldMap);
        //Hep3Vector stereoTrackHolePos = extrapolateTrackPosition(tState,endPosition.x(),stereoSensorHolePosition.z(),5,fieldMap);
        //Hep3Vector stereoTrackSlotPos = extrapolateTrackPosition(tState,endPosition.x(),stereoSensorSlotPosition.z(),5,fieldMap);
        
        HelicalTrackFit htf = TrackUtils.getHTF(tState);
        
        Hep3Vector axialTrackHolePos = getLocationAtSensor(htf,axialSensorHole,bfield);
        Hep3Vector axialTrackSlotPos = getLocationAtSensor(htf,axialSensorSlot,bfield);
        Hep3Vector stereoTrackHolePos = getLocationAtSensor(htf,stereoSensorHole,bfield);
        Hep3Vector stereoTrackSlotPos = getLocationAtSensor(htf,stereoSensorSlot,bfield);
        
        Pair<Boolean,Integer> axialHolePair = this.sensorContainsTrack(axialTrackHolePos, axialSensorHole);
        Pair<Boolean,Integer> axialSlotPair = this.sensorContainsTrack(axialTrackSlotPos, axialSensorSlot);
        Pair<Boolean,Integer> stereoHolePair = this.sensorContainsTrack(stereoTrackHolePos, stereoSensorHole);
        Pair<Boolean,Integer> stereoSlotPair = this.sensorContainsTrack(stereoTrackSlotPos, stereoSensorSlot);
        
        if(axialHolePair.getFirst() && axial){
        	return new Pair<> (axialSensorHole,axialHolePair.getSecond());
        }
    	
    	if(axialSlotPair.getFirst() && axial){
        	return new Pair<> (axialSensorSlot,axialSlotPair.getSecond());
        }
       		
        if(stereoHolePair.getFirst() && !axial){
            return new Pair<> (stereoSensorHole,stereoHolePair.getSecond());
        }
        
        if(stereoSlotPair.getFirst() && !axial){
            return new Pair<> (stereoSensorSlot,stereoSlotPair.getSecond());
        }
        
        return null;
    }
    
    public Pair<Boolean,Integer> sensorContainsTrack(Hep3Vector trackPosition, HpsSiSensor sensor){
    	Hep3Vector pos = globalToSensor(trackPosition, sensor);
    	int nChan = sensor.getNumberOfChannels();
    	int chan = getChan(pos,sensor);
    	//double width = 100;
    	
    	double width = getSensorLength(sensor);
    	
    	//System.out.println(sensor.getName() + " pos " + sensor.getGeometry().getPosition() + " sensor pos " + globalToSensor(sensor.getGeometry().getPosition(),sensor) + " width " + width);

    	if(chan < this.chanExtd || chan > (nChan + this.chanExtd)){
    		return new Pair<>(false,chan);
    	}
    	if(Math.abs(pos.y())>width/2){
    		return new Pair<>(false,chan);
    	}
    	return new Pair<>(true,chan);
    }

    protected double getSensorLength(HpsSiSensor sensor) {

        double length = 0;

        // Get the faces normal to the sensor
        final List<Polygon3D> faces = ((Box) sensor.getGeometry().getLogicalVolume().getSolid())
                .getFacesNormalTo(new BasicHep3Vector(0, 0, 1));
        for (final Polygon3D face : faces) {

            // Loop through the edges of the sensor face and find the longest
            // one
            final List<LineSegment3D> edges = face.getEdges();
            for (final LineSegment3D edge : edges) {
                if (edge.getLength() > length) {
                    length = edge.getLength();
                }
            }
        }
        return length;
    }
    
    public void endOfData(){
        System.out.println("End of Data. Computing Hit Efficiencies");
        PrintWriter out = null;
        try {
			out = new PrintWriter(outputFileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

        for(HpsSiSensor sensor:sensors){
        	String sensorName = sensor.getName();
        	int nChan = sensor.getNumberOfChannels() + chanExtd*2;
        	for(int i = 0; i < nChan; i++){ 
        		int chan = i - chanExtd;
        	    if(numberOfTracksChannel.get(sensorName).binHeight(i) != 0 && numberOfTracksChannel.get(sensorName).binHeight(i) != 0){
        	        hitEfficiencyChannel.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i)/(double)numberOfTracksChannel.get(sensorName).binHeight(i));
        	        hitEfficiencyChannelerr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannel.get(sensorName).binHeight(i)));
        	        hitEfficiencyChannelCorrected.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrected.get(sensorName).binHeight(i));
        	        hitEfficiencyChannelCorrectederr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannel.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksChannelEle.get(sensorName).binHeight(i) != 0 && numberOfTracksChannelEle.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyChannelEle.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i)/(double)numberOfTracksChannelEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelEleerr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelEle.get(sensorName).binHeight(i)));
        	    	hitEfficiencyChannelCorrectedEle.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrectedEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelCorrectedEleerr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelCorrectedEle.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksChannelPos.get(sensorName).binHeight(i) != 0 && numberOfTracksChannelPos.get(sensorName).binHeight(i)!= 0){
        	    	hitEfficiencyChannelPos.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i)/(double)numberOfTracksChannelPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelPoserr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelPos.get(sensorName).binHeight(i)));
        	    	hitEfficiencyChannelCorrectedPos.get(sensorName).fill(chan,numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrectedPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelCorrectedPoserr.get(sensorName).fill(chan,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelCorrectedPos.get(sensorName).binHeight(i)));
        	    }
        	}
        	for(int i = 0; i < nBins; i++){
        		double yMax = sensor.getNumberOfChannels() * sensor.getReadoutStripPitch() / 2;
        		double yMin = -yMax;
        		double y = (yMax-yMin)/(double) nBins * (i + 0.5) + yMin;
        		double pMax = 1.3 * ebeam;
        		double pMin = 0;
        		double p = (pMax-pMin)/(double) nBins * (i + 0.5) + pMin;
        	    if(numberOfTracksY.get(sensorName).binHeight(i) != 0 && numberOfTracksY.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyY.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerY.get(sensorName).binHeight(i)/(double)numberOfTracksY.get(sensorName).binHeight(i));
        	    	hitEfficiencyYerr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerY.get(sensorName).binHeight(i) + 1/(double)numberOfTracksY.get(sensorName).binHeight(i)));
        	    	hitEfficiencyYCorrected.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerY.get(sensorName).binHeight(i)/(double)numberOfTracksYCorrected.get(sensorName).binHeight(i));
        	    	hitEfficiencyYCorrectederr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerY.get(sensorName).binHeight(i) + 1/(double)numberOfTracksYCorrected.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksYEle.get(sensorName).binHeight(i) != 0 && numberOfTracksYEle.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyYEle.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerYEle.get(sensorName).binHeight(i)/(double)numberOfTracksYEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyYEleerr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerYEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksYEle.get(sensorName).binHeight(i)));
        	    	hitEfficiencyYCorrectedEle.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerYEle.get(sensorName).binHeight(i)/(double)numberOfTracksYCorrectedEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyYCorrectedEleerr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerYEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksYCorrectedEle.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksYPos.get(sensorName).binHeight(i) != 0 && numberOfTracksYPos.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyYPos.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerYPos.get(sensorName).binHeight(i)/(double)numberOfTracksYPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyYPoserr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerYPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksYPos.get(sensorName).binHeight(i)));
        	    	hitEfficiencyYCorrectedPos.get(sensorName).fill(y,numberOfTracksWithHitOnMissingLayerYPos.get(sensorName).binHeight(i)/(double)numberOfTracksYCorrectedPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyYCorrectedPoserr.get(sensorName).fill(y,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerYPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksYCorrectedPos.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksP.get(sensorName).binHeight(i) != 0 && numberOfTracksP.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyP.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerP.get(sensorName).binHeight(i)/(double)numberOfTracksP.get(sensorName).binHeight(i));
        	    	hitEfficiencyPerr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerP.get(sensorName).binHeight(i) + 1/(double)numberOfTracksP.get(sensorName).binHeight(i)));
        	    	hitEfficiencyPCorrected.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerP.get(sensorName).binHeight(i)/(double)numberOfTracksPCorrected.get(sensorName).binHeight(i));
        	    	hitEfficiencyPCorrectederr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerP.get(sensorName).binHeight(i) + 1/(double)numberOfTracksPCorrected.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksPEle.get(sensorName).binHeight(i) != 0 && numberOfTracksPEle.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyPEle.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerPEle.get(sensorName).binHeight(i)/(double)numberOfTracksPEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyPEleerr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerPEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksPEle.get(sensorName).binHeight(i)));
        	    	hitEfficiencyPCorrectedEle.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerPEle.get(sensorName).binHeight(i)/(double)numberOfTracksPCorrectedEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyPCorrectedEleerr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerPEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksPCorrectedEle.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksPPos.get(sensorName).binHeight(i) != 0 && numberOfTracksPPos.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyPPos.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerPPos.get(sensorName).binHeight(i)/(double)numberOfTracksPPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyPPoserr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerPPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksPPos.get(sensorName).binHeight(i)));
        	    	hitEfficiencyPCorrectedPos.get(sensorName).fill(p,numberOfTracksWithHitOnMissingLayerPPos.get(sensorName).binHeight(i)/(double)numberOfTracksPCorrectedPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyPCorrectedPoserr.get(sensorName).fill(p,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerPPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksPCorrectedPos.get(sensorName).binHeight(i)));
        	    }
        	}
        	double totalEff = 0;
        	double totalEffEle = 0;
        	double totalEffPos = 0;
        	double totalEfferr = 0;
        	double totalEffEleerr = 0;
        	double totalEffPoserr = 0;
        	
        	double totalEffCorrected = 0;
        	double totalEffEleCorrected = 0;
        	double totalEffPosCorrected = 0;
        	double totalEfferrCorrected = 0;
        	double totalEffEleerrCorrected = 0;
        	double totalEffPoserrCorrected = 0;
        	
        	if(numberOfTracksChannel.get(sensorName).sumAllBinHeights() != 0 && numberOfTracksChannel.get(sensorName).sumAllBinHeights() != 0){
        		totalEff = numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannel.get(sensorName).sumAllBinHeights();
        		totalEfferr = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannel.get(sensorName).sumAllBinHeights());
        		totalEffCorrected = numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannelCorrected.get(sensorName).sumAllBinHeights();
        		totalEfferrCorrected = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannelCorrected.get(sensorName).sumAllBinHeights());
        	}
        	if(numberOfTracksChannelEle.get(sensorName).sumAllBinHeights() != 0 && numberOfTracksChannelEle.get(sensorName).sumAllBinHeights() != 0){
        		totalEffEle = numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannelEle.get(sensorName).sumAllBinHeights();
        		totalEffEleerr = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannelEle.get(sensorName).sumAllBinHeights());
        		totalEffEleCorrected = numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannelCorrectedEle.get(sensorName).sumAllBinHeights();
        		totalEffEleerrCorrected = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannelCorrectedEle.get(sensorName).sumAllBinHeights());
        	}
        	if(numberOfTracksChannelPos.get(sensorName).sumAllBinHeights() != 0 && numberOfTracksChannelPos.get(sensorName).sumAllBinHeights() != 0){
        		totalEffPos = numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannelPos.get(sensorName).sumAllBinHeights();
        		totalEffPoserr = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannelPos.get(sensorName).sumAllBinHeights());
        		totalEffPosCorrected = numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).sumAllBinHeights()/(double)numberOfTracksChannelCorrectedPos.get(sensorName).sumAllBinHeights();
        		totalEffPoserrCorrected = Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).sumAllBinHeights() + 1/(double)numberOfTracksChannelCorrectedPos.get(sensorName).sumAllBinHeights());
        	}
        	
        	TotalEff.get(sensorName).fill(0,totalEff);
        	TotalEffEle.get(sensorName).fill(0,totalEffEle);
        	TotalEffPos.get(sensorName).fill(0,totalEffPos);
        	TotalEfferr.get(sensorName).fill(0,totalEfferr);
        	TotalEffEleerr.get(sensorName).fill(0,totalEffEleerr);
        	TotalEffPoserr.get(sensorName).fill(0,totalEffPoserr);
        	
        	TotalCorrectedEff.get(sensorName).fill(0,totalEffCorrected);
        	TotalCorrectedEffEle.get(sensorName).fill(0,totalEffEleCorrected);
        	TotalCorrectedEffPos.get(sensorName).fill(0,totalEffPosCorrected);
        	TotalCorrectedEfferr.get(sensorName).fill(0,totalEfferrCorrected);
        	TotalCorrectedEffEleerr.get(sensorName).fill(0,totalEffEleerrCorrected);
        	TotalCorrectedEffPoserr.get(sensorName).fill(0,totalEffPoserrCorrected);
        	
        	System.out.println(sensorName + " Total Efficiency = " + totalEff + " with error " + totalEfferr);
        	System.out.println(sensorName + " Total Corrected Efficiency = " + totalEffCorrected + " with error " + totalEfferrCorrected);
        	
        	org.hps.util.Pair<Integer, Integer> daqPair = getDaqPair(daqMap, sensor);
        	Collection<SvtChannel> channels = channelMap.find(daqPair);
        	for (SvtChannel channel : channels) {
                int chanID = channel.getChannelID();
                int chan = channel.getChannel();
        	    double eff = hitEfficiencyChannel.get(sensorName).binHeight(chan);
        	    out.println(chanID + ", " + eff);
        	}
        }
        out.close();
    }
    
    static org.hps.util.Pair<Integer, Integer> getDaqPair(SvtDaqMappingCollection daqMap, HpsSiSensor sensor) {

        final String svtHalf = sensor.isTopLayer() ? AbstractSvtDaqMapping.TOP_HALF : AbstractSvtDaqMapping.BOTTOM_HALF;
        for (final SvtDaqMapping object : daqMap) {

            if (svtHalf.equals(object.getSvtHalf()) && object.getLayerNumber() == sensor.getLayerNumber()
                    && object.getSide().equals(sensor.getSide())) {

                return new org.hps.util.Pair<Integer, Integer>(object.getFebID(), object.getFebHybridID());
            }
        }
        return null;
    }
}