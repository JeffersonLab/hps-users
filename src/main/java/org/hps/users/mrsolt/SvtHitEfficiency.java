/**
 * Driver used to compute SVT hit efficiencies at each sensor
 * as a function of strip and y
 */
/**
 * @author mrsolt
 *
 */
package org.hps.users.mrsolt;

import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotterFactory;
import hep.aida.IPlotter;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.ITree;
import hep.physics.matrix.BasicMatrix;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import org.lcsim.constants.Constants;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.IRotation3D;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.Rotation3D;
import org.lcsim.detector.Transform3D;
import org.lcsim.detector.converter.compact.subdetector.HpsTracker2;
import org.lcsim.detector.converter.compact.subdetector.SvtStereoLayer;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.solids.Point3D;
import org.lcsim.detector.solids.Polygon3D;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCIOParameters;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseTrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.geometry.compact.Subdetector;
import org.lcsim.lcio.LCIOConstants;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.swim.Trajectory;
import org.apache.commons.math3.util.Pair;
import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.svt.AbstractSvtDaqMapping;
import org.hps.conditions.svt.SvtChannel;
import org.hps.conditions.svt.SvtConditions;
import org.hps.conditions.svt.SvtDaqMapping;
import org.hps.conditions.svt.SvtDaqMapping.SvtDaqMappingCollection;
import org.hps.conditions.svt.TestRunSvtChannel;
import org.hps.conditions.svt.TestRunSvtDaqMapping;
import org.hps.conditions.svt.TestRunSvtChannel.TestRunSvtChannelCollection;
import org.hps.conditions.svt.TestRunSvtConditions;
import org.hps.conditions.svt.SvtChannel.SvtChannelCollection;
import org.hps.conditions.svt.TestRunSvtDaqMapping.TestRunSvtDaqMappingCollection;
import org.hps.detector.svt.SvtDetectorSetup;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.FittedRawTrackerHit;
import org.hps.recon.tracking.TrackStateUtils;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackerHitUtils;
import org.hps.recon.tracking.gbl.FittedGblTrajectory;
import org.hps.recon.tracking.gbl.GblPoint;
import org.hps.recon.tracking.gbl.GblTrajectory;
import org.hps.recon.tracking.gbl.GblUtils;
import org.hps.recon.tracking.gbl.HelicalTrackStripGbl;
import org.hps.recon.tracking.gbl.HpsGblRefitter;
import org.hps.recon.tracking.gbl.FittedGblTrajectory.GBLPOINT;
import org.hps.recon.tracking.gbl.matrix.Matrix;
import org.hps.recon.tracking.gbl.matrix.SymMatrix;
import org.hps.recon.tracking.gbl.matrix.Vector;

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
    Map<String, IHistogram1D> residualYerrorYDiff = new HashMap<String,IHistogram1D>();
    
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
    double maxTLambdaErr = 0.5;
    double minTLambdaErr = 0;
    double maxPhi0Err = 0.01;
    double minPhi0Err = 0;
    double maxOmegaErr = 0.01;
    double minOmegaErr = 0;
    
    String atIP = "IP";
    
    //Collection Strings
    private String stripHitOutputCollectionName = "StripClusterer_SiTrackerHitStrip1D";
    
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
        	int nChan = sensor.getNumberOfChannels();
        	double readoutPitch = sensor.getReadoutStripPitch();
        	double maxY = nChan * readoutPitch / 2;
        	numberOfTracksChannel.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel " + sensorName, nChan, 0, nChan));
        	numberOfTracksWithHitOnMissingLayerChannel.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannel.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel " + sensorName, nChan, 0, nChan));
        	numberOfTracksY.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksWithHitOnMissingLayerY.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyY.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksP.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerP.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyP.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Ele " + sensorName, nChan, 0, nChan));
        	numberOfTracksWithHitOnMissingLayerChannelEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel Ele " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannelEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Ele " + sensorName, nChan, 0, nChan));
        	numberOfTracksYEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Ele " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksWithHitOnMissingLayerYEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y Ele " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyYEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Ele " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksPEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Ele " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerPEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyPEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Pos " + sensorName, nChan, 0, nChan));
        	numberOfTracksWithHitOnMissingLayerChannelPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Channel Pos " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannelPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Pos " + sensorName, nChan, 0, nChan));
        	numberOfTracksYPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Pos " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksWithHitOnMissingLayerYPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit Y Pos " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyYPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Pos " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksPPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Pos " + sensorName, nBins, 0, 1.3*ebeam));
        	numberOfTracksWithHitOnMissingLayerPPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks With Hit P Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	hitEfficiencyPPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannelCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected " + sensorName, nChan, 0, nChan));
        	numberOfTracksYCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyYCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksPCorrected.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrected.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected Ele " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannelCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Ele " + sensorName, nChan, 0, nChan));
        	numberOfTracksYCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected Ele " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyYCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Ele " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksPCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected Ele " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrectedEle.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Ele " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	numberOfTracksChannelCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Channel Corrected Pos " + sensorName, nChan, 0, nChan));
        	hitEfficiencyChannelCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Pos " + sensorName, nChan, 0, nChan));
        	numberOfTracksYCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks Y Corrected Pos " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyYCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Pos " + sensorName, nBins, -maxY, maxY));
        	numberOfTracksPCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("Number of Tracks P Corrected Pos " + sensorName, nBins, 0, 1.3*ebeam));
        	hitEfficiencyPCorrectedPos.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Pos " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	TotalEff.put(sensorName,histogramFactory.createHistogram1D("Total Eff " + sensorName, 1, 0, 1));
        	TotalEffEle.put(sensorName,histogramFactory.createHistogram1D("Total Eff Ele " + sensorName, 1, 0, 1));
        	TotalEffPos.put(sensorName,histogramFactory.createHistogram1D("Total Eff Pos " + sensorName, 1, 0, 1));
        	
        	TotalCorrectedEff.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff " + sensorName, 1, 0, 1));
        	TotalCorrectedEffEle.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Ele " + sensorName, 1, 0, 1));
        	TotalCorrectedEffPos.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Pos " + sensorName, 1, 0, 1));
        	
        	hitEfficiencyChannelerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPCorrectederr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Ele Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Ele Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Ele Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Ele Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Ele Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPCorrectedEleerr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Ele Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Pos Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Pos Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Pos Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	hitEfficiencyChannelCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Channel Corrected Pos Error " + sensorName, nChan, 0, nChan));
        	hitEfficiencyYCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency Y Corrected Pos Error " + sensorName, nBins, -maxY, maxY));
        	hitEfficiencyPCorrectedPoserr.put(sensorName,histogramFactory.createHistogram1D("HitEfficiency P Corrected Pos Error " + sensorName, nBins, 0,  1.3*ebeam));
        	
        	TotalEfferr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Error " + sensorName, 1, 0, 1));
        	TotalEffEleerr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Ele Error " + sensorName, 1, 0, 1));
        	TotalEffPoserr.put(sensorName,histogramFactory.createHistogram1D("Total Eff Pos Error " + sensorName, 1, 0, 1));
        	
        	TotalCorrectedEfferr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Error " + sensorName, 1, 0, 1));
        	TotalCorrectedEffEleerr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Ele Error " + sensorName, 1, 0, 1));
        	TotalCorrectedEffPoserr.put(sensorName,histogramFactory.createHistogram1D("Total Corrected Eff Pos Error " + sensorName, 1, 0, 1));
        	
        	residualY.put(sensorName,histogramFactory.createHistogram1D("Residual Y " + sensorName, nBins, -10, 10));
        	errorY.put(sensorName,histogramFactory.createHistogram1D("Error Y " + sensorName, nBins, 0, 0.5));
        	residualYerrorYDiff.put(sensorName,histogramFactory.createHistogram1D("Residual Y - Error Y " + sensorName, nBins, -10, 10));
        	
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
        
        //Grab all the clusters in the event
        List<SiTrackerHitStrip1D> stripHits = event.get(SiTrackerHitStrip1D.class, stripHitOutputCollectionName);
        
        for(List<Track> tracks:trackCollections){
        	
    		if(cleanTridents){
    			// Require an event to have exactly two tracks
    			if (tracks.size() != 2) continue;

    			// Require the two tracks to be in opposite volumes
    			if (tracks.get(0).getTrackStates().get(0).getTanLambda()*tracks.get(1).getTrackStates().get(0).getTanLambda() >= 0) continue;

    			// Require the two tracks to be oppositely charged
    			if (tracks.get(0).getTrackStates().get(0).getOmega()*tracks.get(1).getTrackStates().get(0).getOmega() >= 0) continue;
    		}
    		
            for(Track track:tracks){
            	boolean gbl = isGBL(track);//fix this function
            	if(!gbl) continue;
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
    	    	//HelicalTrackFit htf = TrackUtils.getHTF(track.getTrackStates().get(0));
    	    	//Hep3Matrix matrix = getPerToClPrj(htf);
    	    	
    	    	Hep3Vector axialExtrapPos = extrapolateTrackPosition(tState,endHitPos.x(),axialPos.z(),5,bFieldMap);
    	    	Hep3Vector stereoExtrapPos = extrapolateTrackPosition(tState,endHitPos.x(),stereoPos.z(),5,bFieldMap);
    	    	
    	    	Hep3Vector axialExtrapPosSensor = globalToSensor(axialExtrapPos,axialSensor);
    	    	Hep3Vector stereoExtrapPosSensor = globalToSensor(stereoExtrapPos,stereoSensor);
    	    	
    	    	//Compute the extrapolation errors
    	    	//compute error correctly
    	    	System.out.println("Layer " + unusedLay);
    	    	double yErrorAxial = computeExtrapErrorY(tState,axialSensor,endHitPos)[0];
    	    	double yErrorStereo = computeExtrapErrorY(tState,stereoSensor,endHitPos)[0];
    	    	
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
	    		
    	    	//Loop over all tracker hits in the event
    	    	boolean fillAxial = false;
    	    	boolean fillStereo = false;
    	    	for(SiTrackerHitStrip1D hit:stripHits){
    	    		if(fillAxial && fillStereo) break;
    	    		//Get the sensor and position of the hit
    	    		HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();
            	    double[] hitPos = hit.getPosition();
            	    //Check to see if the sensor of this hit is the same sensor you expect to see an axial hit
            	    if(sensorAxialName == sensor.getName() && !fillAxial){
            	    	//Compute the residual between extrapolated track and hit position and fill histo
            	    	double residual = axialExtrapPos.y() - hitPos[1];
            	    	residualY.get(sensorAxialName).fill(residual);
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
            	    	double residual = stereoExtrapPos.y() - hitPos[1];
            	    	residualY.get(sensorStereoName).fill(residual);
            	    	//Check to see if the residual is within desired error values
            	    	//If it is, then fill numerator of efficiency histos
            	    	residualYerrorYDiff.get(sensorStereoName).fill(Math.abs(residual) - yErrorStereo);
            	    	if((Math.abs(residual) > this.nSig * yErrorAxial)) continue;
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
    	    	}
            }
        }
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
	    if(unusedLay == 1){
	    	return track.getTrackStates().get(0);
	    }
	    else{
	    	layer = unusedLay - 1;
	    }
	    HpsSiSensor sensorHole = getSensor(track,layer,true,true);
	    HpsSiSensor sensorSlot = getSensor(track,layer,true,false);
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
    
    public static Hep3Vector extrapolateTrackPosition(TrackState track, double startPositionX, double endPositionX, double stepSize, FieldMap fieldMap) {

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
    }
    
    private double[] computeExtrapErrorY(TrackState tState, HpsSiSensor sensor, Hep3Vector hitPos){
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
       
        double cosLambda = sqrt(1.0 - sin(tState.getTanLambda()) * sin(tState.getTanLambda()));       
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
    	
    	
        Matrix helixCovariance = jacobian.times(( locCov).times(jacobian.transpose()));
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
    }
    
    private Transform3D initialize() {
        BasicHep3Matrix tmp = new BasicHep3Matrix();
        tmp.setElement(0, 2, 1);
        tmp.setElement(1, 0, 1);
        tmp.setElement(2, 1, 1);
        return new Transform3D(new Rotation3D(tmp));
    }

    private Transform3D initializeInverse() {
        return initialize().inverse();
    }
    
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
    
    private SymMatrix buildCovMatrix(double[] cov){
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
    }

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
        
        Hep3Vector axialTrackHolePos = extrapolateTrackPosition(tState,endPosition.x(),axialSensorHolePosition.z(),5,fieldMap);
        Hep3Vector axialTrackSlotPos = extrapolateTrackPosition(tState,endPosition.x(),axialSensorSlotPosition.z(),5,fieldMap);
        Hep3Vector stereoTrackHolePos = extrapolateTrackPosition(tState,endPosition.x(),stereoSensorHolePosition.z(),5,fieldMap);
        Hep3Vector stereoTrackSlotPos = extrapolateTrackPosition(tState,endPosition.x(),stereoSensorSlotPosition.z(),5,fieldMap);
        
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
    	double width = 100;

    	if(chan < 0 || chan > nChan){
    		return new Pair<>(false,chan);
    	}
    	if(Math.abs(pos.y())>width/2){
    		return new Pair<>(false,chan);
    	}
    	return new Pair<>(true,chan);
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
        	int nChan = sensor.getNumberOfChannels();
        	for(int i = 0; i < nChan; i++){ 
        	    if(numberOfTracksChannel.get(sensorName).binHeight(i) != 0 && numberOfTracksChannel.get(sensorName).binHeight(i) != 0){
        	        hitEfficiencyChannel.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i)/(double)numberOfTracksChannel.get(sensorName).binHeight(i));
        	        hitEfficiencyChannelerr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannel.get(sensorName).binHeight(i)));
        	        hitEfficiencyChannelCorrected.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrected.get(sensorName).binHeight(i));
        	        hitEfficiencyChannelCorrectederr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannel.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannel.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksChannelEle.get(sensorName).binHeight(i) != 0 && numberOfTracksChannelEle.get(sensorName).binHeight(i) != 0){
        	    	hitEfficiencyChannelEle.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i)/(double)numberOfTracksChannelEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelEleerr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelEle.get(sensorName).binHeight(i)));
        	    	hitEfficiencyChannelCorrectedEle.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrectedEle.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelCorrectedEleerr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelEle.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelCorrectedEle.get(sensorName).binHeight(i)));
        	    }
        	    if(numberOfTracksChannelPos.get(sensorName).binHeight(i) != 0 && numberOfTracksChannelPos.get(sensorName).binHeight(i)!= 0){
        	    	hitEfficiencyChannelPos.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i)/(double)numberOfTracksChannelPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelPoserr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelPos.get(sensorName).binHeight(i)));
        	    	hitEfficiencyChannelCorrectedPos.get(sensorName).fill(i,numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i)/(double)numberOfTracksChannelCorrectedPos.get(sensorName).binHeight(i));
        	    	hitEfficiencyChannelCorrectedPoserr.get(sensorName).fill(i,Math.sqrt(1/(double)numberOfTracksWithHitOnMissingLayerChannelPos.get(sensorName).binHeight(i) + 1/(double)numberOfTracksChannelCorrectedPos.get(sensorName).binHeight(i)));
        	    }
        	}
        	for(int i = 0; i < nBins; i++){
        		double yMax = sensor.getNumberOfChannels() * sensor.getReadoutStripPitch() / 2;
        		double yMin = -yMax;
        		double y = (yMax-yMin)/(double) nBins * i + yMin;
        		double pMax = 1.3 * ebeam;
        		double pMin = 0;
        		double p = (pMax-pMin)/(double) nBins * i + pMin;
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