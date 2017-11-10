/**
 * Analysis driver to calculate hit efficiencies in the SVT
 */
/**
 * @author mrsolt
 *
 */

package org.hps.users.mrsolt;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import hep.aida.IHistogram1D;
import hep.aida.ITree;
import hep.physics.vec.Hep3Vector;

import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.apache.commons.math3.util.Pair;
import org.hps.analysis.MC.MCFullDetectorTruth;

public class SimpAcceptance extends Driver {
	
    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree; 

    
    int FindableDecaysWithDaughtersSameVol = 0;
    int FindableDecaysWithDaughtersOppoVol = 0;
    int FindableDecaysWithDaughtersTot = 0;
    int FindableDecaysWithDaughtersSameVolWithEcal = 0;
    int FindableDecaysWithDaughtersOppoVolWithEcal = 0;
    int FindableDecaysWithDaughtersTotWithEcal = 0;
    
    int FindableDecaysWithDaughtersWithRecoilSameVol = 0;
    int FindableDecaysWithDaughtersWithRecoilOppoVol = 0;
    int FindableDecaysWithDaughtersWithRecoilTot = 0;
    int FindableDecaysWithDaughtersWithRecoilSameVolWithEcal = 0;
    int FindableDecaysWithDaughtersWithRecoilOppoVolWithEcal = 0;
    int FindableDecaysWithDaughtersWithRecoilTotWithEcal = 0;
    
    int events = 0;
    int nLay = 6;
    int nHitsReq = 5;
    

    IHistogram1D findableDecaysWithDaughters;
    IHistogram1D findableDecaysWithDaughtersSameVolWithEcal;
    IHistogram1D findableDecaysWithDaughtersOppoVolWithEcal;
    IHistogram1D findableDecaysWithDaughtersTotWithEcal;
   

    IHistogram1D findableDecaysWithDaughtersWithRecoil;
    IHistogram1D findableDecaysWithDaughtersWithRecoilSameVolWithEcal;
    IHistogram1D findableDecaysWithDaughtersWithRecoilOppoVolWithEcal;
    IHistogram1D findableDecaysWithDaughtersWithRecoilTotWithEcal;
    
    IHistogram1D NumberofEvents;

    
    
    public void detectorChanged(Detector detector){
    	
    	aida.tree().cd("/");
    	tree = aida.tree();
        
        findableDecaysWithDaughters = aida.histogram1D("findableDecaysWithDaughters", 1, 0, 1);
        findableDecaysWithDaughtersSameVolWithEcal = aida.histogram1D("findableDecaysWithDaughtersSameVolWithEcal", 1, 0, 1);
        findableDecaysWithDaughtersOppoVolWithEcal = aida.histogram1D("findableDecaysWithDaughtersOppoVolWithEcal", 1, 0, 1);
        findableDecaysWithDaughtersTotWithEcal = aida.histogram1D("findableDecaysWithDaughtersTotWithEcal", 1, 0, 1);
       
        findableDecaysWithDaughtersWithRecoil = aida.histogram1D("findableDecaysWithDaughtersWithRecoil", 1, 0, 1);
        findableDecaysWithDaughtersWithRecoilSameVolWithEcal = aida.histogram1D("findableDecaysWithDaughtersWithRecoilSameVolWithEcal", 1, 0, 1);
        findableDecaysWithDaughtersWithRecoilOppoVolWithEcal = aida.histogram1D("findableDecaysWithDaughtersWithRecoilOppoVolWithEcal", 1, 0, 1);
        findableDecaysWithDaughtersWithRecoilTotWithEcal = aida.histogram1D("findableDecaysWithDaughtersWithRecoilTotWithEcal", 1, 0, 1);
        
        NumberofEvents = aida.histogram1D("Number of Events", 1, 0, 1);
        
    }
    
    public void setNLay(int nLay) { 
        this.nLay = nLay;
    }
    
    public void setNHitsReq(int nHitsReq) { 
        this.nHitsReq = nHitsReq;
    }
    
    private Map<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> constructHitMap(Map<MCParticle, List<SimTrackerHit>> trackerHitMap,Map<MCParticle, List<SimCalorimeterHit>> calHitMap){
    	Map<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> hitMap = new HashMap<MCParticle,Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>>();
    	for (Entry<MCParticle, List<SimTrackerHit>> entry : trackerHitMap.entrySet()){
            MCParticle p = entry.getKey();
            List<SimTrackerHit> hits = entry.getValue();
            if (hitMap.get(p) == null) {
                hitMap.put(p, new Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>(new ArrayList<SimTrackerHit>(), new ArrayList<SimCalorimeterHit>()));
            }
            hitMap.get(p).getFirst().addAll(hits);
    	}
        for (Entry<MCParticle, List<SimCalorimeterHit>> entry : calHitMap.entrySet()){
            MCParticle p = entry.getKey();
            List<SimCalorimeterHit> hits = entry.getValue();
            if (hitMap.get(p) == null) {
                hitMap.put(p, new Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>(new ArrayList<SimTrackerHit>(), new ArrayList<SimCalorimeterHit>()));
            }
            hitMap.get(p).getSecond().addAll(hits);
        }
    	return hitMap;
    }
    
    private boolean IsVDaughterEle(MCParticle p){
    	return IsVDaughter(p) && p.getPDGID() == 11;
    }
    
    private boolean IsVDaughterPos(MCParticle p){
    	return IsVDaughter(p) && p.getPDGID() == -11;
    }
    
    private boolean IsVDaughter(MCParticle p){
    	return true;
    	//return p.getParents() == null;
    }
    
    private boolean IsRecoil(MCParticle p){
    	return false;
    	//return p.getParents() == null && p.getOriginZ() < 0.0001;
    }
    
    private boolean IsTrackFindable(List<SimTrackerHit> trackerhits){
    	boolean[] hitsOnSensor = new boolean[this.nLay*2];
    	int nHits = 0;
    	for(int i = 0; i < 2*this.nLay; i++){
    		hitsOnSensor[i] = false;
    	}
    	for(SimTrackerHit hit:trackerhits){
    		hitsOnSensor[hit.getLayer()-1] = true;
    	}
    	for(int i = 0; i < this.nLay; i++){
    		if(hitsOnSensor[2*i] == true && hitsOnSensor[2*i+1] == true){
    			nHits++;
    		}
    	}
    	if(nHits >= this.nHitsReq){
    		return true;	
    	}
    	else {
    		return false;
    	}
    }
    
    
    @Override
    public void process(EventHeader event){
		//aida.tree().cd("/");
		
		events++;
        
        List<SimTrackerHit> trackerHits = event.get(SimTrackerHit.class, "TrackerHits");
        List<SimCalorimeterHit> calHits = event.get(SimCalorimeterHit.class, "EcalHits");
        
        Map<MCParticle, List<SimTrackerHit>> trackerHitMap = MCFullDetectorTruth.BuildTrackerHitMap(trackerHits);
        Map<MCParticle, List<SimCalorimeterHit>> calHitMap = MCFullDetectorTruth.BuildCalHitMap(calHits);
        Map<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> hitMap = constructHitMap(trackerHitMap, calHitMap);
        Map<String,boolean[]> findable = new HashMap<>();
        boolean[] dummy = {false,false,false};
        
        findable.put("ele", dummy);
        findable.put("pos", dummy);
        findable.put("eleRec", dummy);       
        
        for (Entry<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> entry : hitMap.entrySet()) {
        	//System.out.println("Pass 1");
            
            MCParticle p = entry.getKey();
            boolean isVDaughterPos = IsVDaughterPos(p);
            boolean isVDaughterEle = IsVDaughterEle(p);
            boolean isRecoil = IsRecoil(p);
            if(!isVDaughterEle && !isVDaughterPos && !isRecoil) continue;
            List<SimTrackerHit> trackerhits = entry.getValue().getFirst();
            List<SimCalorimeterHit> ecalhits = entry.getValue().getSecond();
            //System.out.println("Pass 2");

            boolean findableTrack = IsTrackFindable(trackerhits);          
            if(!findableTrack) continue;
            //System.out.println("Pass 3");
            
            if(isVDaughterEle) findable.get("ele")[0] = true;
            if(isVDaughterPos) findable.get("pos")[0] = true;
            if(isRecoil) findable.get("eleRec")[0] = true;
            
            if(ecalhits.size() != 0){
            	Hep3Vector ecalPos = ecalhits.get(0).getPositionVec();
            	if(ecalPos.y() > 0){
            		if(isVDaughterEle) findable.get("ele")[1] = true;
                	if(isVDaughterPos) findable.get("pos")[1] = true;
                	if(isRecoil) findable.get("eleRec")[1] = true;
            	}
            	else{
            		if(isVDaughterEle) findable.get("ele")[2] = true;
                	if(isVDaughterPos) findable.get("pos")[2] = true;
                	if(isRecoil) findable.get("eleRec")[2] = true;
            	}
            }
        }
        
        if(findable.get("ele")[0] && findable.get("pos")[0]){
        	findableDecaysWithDaughters.fill(0);
        	if(findable.get("eleRec")[0]){
        		findableDecaysWithDaughtersWithRecoil.fill(0);
        	}
        	if((findable.get("ele")[1] && findable.get("pos")[1]) || (findable.get("ele")[2] && findable.get("pos")[2])){
        		findableDecaysWithDaughtersSameVolWithEcal.fill(0);
        		if(findable.get("eleRec")[0] && (findable.get("eleRec")[1] || findable.get("eleRec")[2])){
        			findableDecaysWithDaughtersWithRecoilSameVolWithEcal.fill(0);
        		}
        	}
        	if((findable.get("ele")[1] && findable.get("pos")[2]) || (findable.get("ele")[2] && findable.get("pos")[1])){
        		findableDecaysWithDaughtersOppoVolWithEcal.fill(0);
        		if(findable.get("eleRec")[0] && (findable.get("eleRec")[1] || findable.get("eleRec")[2])){
        			findableDecaysWithDaughtersWithRecoilOppoVolWithEcal.fill(0);
        		}
        	}
        }
    }
    
    @Override
    public void endOfData(){
    	findableDecaysWithDaughtersTotWithEcal.fill(0,findableDecaysWithDaughtersSameVolWithEcal.binHeight(0)+findableDecaysWithDaughtersOppoVolWithEcal.binHeight(0));
    	findableDecaysWithDaughtersWithRecoilTotWithEcal.fill(0,findableDecaysWithDaughtersWithRecoilSameVolWithEcal.binHeight(0)+findableDecaysWithDaughtersWithRecoilOppoVolWithEcal.binHeight(0));
        
        System.out.println("findableDecaysWithDaughters = " + findableDecaysWithDaughters.binHeight(0));
        System.out.println("findableDecaysWithDaughtersSameVolWithEcal = " + findableDecaysWithDaughtersSameVolWithEcal.binHeight(0));
        System.out.println("findableDecaysWithDaughtersOppoVolWithEcal = " + findableDecaysWithDaughtersOppoVolWithEcal.binHeight(0));
        System.out.println("findableDecaysWithDaughtersTotWithEcal = " + findableDecaysWithDaughtersTotWithEcal.binHeight(0));
       
        System.out.println("findableDecaysWithDaughtersWithRecoil = " + findableDecaysWithDaughtersWithRecoil.binHeight(0));
        System.out.println("findableDecaysWithDaughtersWithRecoilSameVolWithEcal = " + findableDecaysWithDaughtersWithRecoilSameVolWithEcal.binHeight(0));
        System.out.println("findableDecaysWithDaughtersWithRecoilOppoVolWithEcal = " + findableDecaysWithDaughtersWithRecoilOppoVolWithEcal.binHeight(0));
        System.out.println("findableDecaysWithDaughtersWithRecoilTotWithEcal = " + findableDecaysWithDaughtersWithRecoilTotWithEcal.binHeight(0));
    	
    	NumberofEvents.fill(0,events);
    	System.out.println("Total Number of Events = " + events);
    	
    	System.out.println("%===================================================================%");
    	System.out.println("%======================  Simp Acceptance Complete==========================%");
    	System.out.println("%===================================================================% \n%");
    	System.out.println("% \n%===================================================================%");
    }
}
