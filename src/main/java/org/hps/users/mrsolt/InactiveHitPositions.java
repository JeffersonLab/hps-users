/**
 * InactiveHitPositions Driver
 */
/**
 * @author mrsolt
 *
 */
package org.hps.users.mrsolt;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.IHistogram2D;
import hep.aida.ITree;

import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;


public class InactiveHitPositions extends Driver {

    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree; 
    IHistogramFactory histogramFactory; 

    //List of Sensors
    private List<HpsSiSensor> sensors = null;

    Map<String, IHistogram2D> ActiveHitMap = new HashMap<String,IHistogram2D>();
    Map<String, IHistogram2D> InActiveHitMap = new HashMap<String,IHistogram2D>();
    
    IHistogram1D Chi2;
    IHistogram1D D0;
    IHistogram1D Pz;
    IHistogram1D Omega;
    IHistogram1D Phi;
    IHistogram1D TanLambda;
    IHistogram1D Z0;
    
    //Histogram Settings
    double minX = -85;
    double maxX = 85;
    double minY = 0;
    double maxY = 60;
    int nBins = 50;
    
    //Collection Strings
    private String activeSimHitsCollectionName = "TrackerHits";
    private String inactiveSimHitsCollectionName = "TrackerHits_Inactive";
   
    //Constants
    private static final String SUBDETECTOR_NAME = "Tracker";
    
    public void detectorChanged(Detector detector){
    	
    	aida.tree().cd("/");
    	tree = aida.tree();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);
        // Get the HpsSiSensor objects from the tracker detector element
        sensors = detector.getSubdetector(SUBDETECTOR_NAME)
                          .getDetectorElement().findDescendants(HpsSiSensor.class);
   
        // If the detector element had no sensors associated with it, throw
        // an exception
        if (sensors.size() == 0) {
            throw new RuntimeException("No sensors were found in this detector.");
        }
        
        Chi2 = aida.histogram1D("chi2", 50, 0, 50);
        D0 = aida.histogram1D("D0", 50, -5, 5);
        Pz = aida.histogram1D("Pz", 50, 0, 2);
        Omega = aida.histogram1D("Omega", 50, 0, 0.02);
        Phi = aida.histogram1D("Phi", 50, 0, 0.2);
        TanLambda = aida.histogram1D("TanLambda", 50, 0, 0.1);
        Z0 = aida.histogram1D("Z0", 50, -5, 5);
        
        
        for(HpsSiSensor sensor:sensors){
        	if(sensor.isTopLayer()){
        		ActiveHitMap.put(sensor.getName(),histogramFactory.createHistogram2D("Active Hits " + sensor.getName(), nBins, minX, maxX, nBins, 0, maxY));
        		InActiveHitMap.put(sensor.getName(),histogramFactory.createHistogram2D("Inactive Hits " + sensor.getName(), nBins, minX, maxX, nBins, 0, maxY));
        	}
        	else{
        		ActiveHitMap.put(sensor.getName(),histogramFactory.createHistogram2D("Active Hits " + sensor.getName(), nBins, minX, maxX, nBins,-maxY, 0));
        		InActiveHitMap.put(sensor.getName(),histogramFactory.createHistogram2D("Inactive Hits " + sensor.getName(), nBins, minX, maxX, nBins, -maxY, 0));
        	}
        }

    }

    public void process(EventHeader event){
		aida.tree().cd("/");	
        
        List<SimTrackerHit> activetrackerHits = event.get(SimTrackerHit.class, activeSimHitsCollectionName);
        
        List<List<Track>> trackCollections = event.get(Track.class);
        
        for (List<Track> tracks : trackCollections) {
        	for(Track track:tracks){
        		double chi2 = track.getChi2();
        		TrackState state = track.getTrackStates().get(0);
        		Chi2.fill(chi2);
        		D0.fill(state.getD0());
        		Pz.fill(state.getMomentum()[2]);
        		Omega.fill(state.getOmega());
        		Phi.fill(state.getPhi());
        		TanLambda.fill(state.getTanLambda());
        		Z0.fill(state.getZ0());
        	}
        }
        
        for(SimTrackerHit hit:activetrackerHits){
        	double[] pos = hit.getPosition();
        	for(HpsSiSensor sensor:sensors){
        		boolean sameLay = sensor.getLayerNumber() == hit.getLayer();
        		boolean sameVolume = (sensor.isTopLayer() && pos[1] > 0) || (!sensor.isTopLayer() && pos[1] < 0);
        		if(sameLay && sameVolume){
        		    ActiveHitMap.get(sensor.getName()).fill(pos[0],pos[1]);	
        		}
        	}
        }
        if (event.hasCollection(SimTrackerHit.class, inactiveSimHitsCollectionName)){
            List<SimTrackerHit> inactivetrackerHits = event.get(SimTrackerHit.class, inactiveSimHitsCollectionName);
            for(SimTrackerHit hit:inactivetrackerHits){
        	    double[] pos = hit.getPosition();
        	    for(HpsSiSensor sensor:sensors){
        		    boolean sameLay = sensor.getLayerNumber() == hit.getLayer();
        		    boolean sameVolume = (sensor.isTopLayer() && pos[1] > 0) || (!sensor.isTopLayer() && pos[1] < 0);
        		    if(sameLay && sameVolume){
        		        InActiveHitMap.get(sensor.getName()).fill(pos[0],pos[1]);	
        		    }
        	    }
            }
        }
        	

    }
    
    public void endOfData(){
        System.out.println("End of Data");
    }
}