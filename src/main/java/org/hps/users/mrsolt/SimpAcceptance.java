/**
 * Analysis driver to calculate acceptances of Simps
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
import hep.aida.IHistogram2D;
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
    
    IHistogram2D ecalHitPosEle;
    IHistogram2D ecalHitPosPos;
    IHistogram2D ecalHitPosEleRec;
    
    IHistogram2D ecalHitPosEle_PosTop;
    IHistogram2D ecalHitPosEle_PosBot;
    IHistogram2D ecalHitPosPos_EleTop;
    IHistogram2D ecalHitPosPos_EleBot;
    
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
        
        ecalHitPosEle = aida.histogram2D("Ecal Hit Position Electron", 50, -400, 400, 50, -100, 100);
        ecalHitPosPos = aida.histogram2D("Ecal Hit Position Positron", 50, -400, 400, 50, -100, 100);
        ecalHitPosEleRec = aida.histogram2D("Ecal Hit Position Recoil Electron", 50, -400, 400, 50, -100, 100);
        
        ecalHitPosEle_PosTop = aida.histogram2D("Ecal Hit Position Electron (Positron Top)", 50, -400, 400, 50, -100, 100);
        ecalHitPosEle_PosBot = aida.histogram2D("Ecal Hit Position Electron (Positron Bot)", 50, -400, 400, 50, -100, 100);
        ecalHitPosPos_EleTop = aida.histogram2D("Ecal Hit Position Positron (Electron Top)", 50, -400, 400, 50, -100, 100);
        ecalHitPosPos_EleBot = aida.histogram2D("Ecal Hit Position Positron (Electron Bot)", 50, -400, 400, 50, -100, 100);

        
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
            if(!hitMap.get(p).getSecond().contains(hits)){
                hitMap.get(p).getSecond().addAll(hits);
            }
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
    	List<MCParticle> parents = p.getParents();
    	if(parents == null) return false;
    	if(parents.size() != 2) return false;
    	boolean V0Parent = false;
    	boolean apParent = false;
    	for(MCParticle parent:parents){
    		if (parent.getPDGID() == 622) apParent = true;
    		if (parent.getPDGID() == 625) V0Parent = true;
    	}
    	return apParent && V0Parent;
    }
    
    private boolean IsRecoil(MCParticle p){
    	List<MCParticle> parents = p.getParents();
    	if(parents == null) return false;
    	if(parents.size() != 1) return false;
    	boolean apParent = parents.get(0).getPDGID() == 622;
    	return apParent && p.getPDGID() == 11;
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
		events++;
        
        List<SimTrackerHit> trackerHits = event.get(SimTrackerHit.class, "TrackerHits");
        List<SimCalorimeterHit> calHits = event.get(SimCalorimeterHit.class, "EcalHits");
        
        //Build the tracker hit and ecal hit maps to MCParticles
        Map<MCParticle, List<SimTrackerHit>> trackerHitMap = MCFullDetectorTruth.BuildTrackerHitMap(trackerHits);
        Map<MCParticle, List<SimCalorimeterHit>> calHitMap = MCFullDetectorTruth.BuildCalHitMap(calHits);
        Map<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> hitMap = constructHitMap(trackerHitMap, calHitMap);
        
        boolean[] findableEle = {false,false,false};
        boolean[] findablePos = {false,false,false};
        boolean[] findableEleRec = {false,false,false};
        Hep3Vector posEle = null;
        Hep3Vector posPos = null;

        for (Entry<MCParticle, Pair<List<SimTrackerHit>, List<SimCalorimeterHit>>> entry : hitMap.entrySet()) {
            MCParticle p = entry.getKey();
            boolean isVDaughterPos = IsVDaughterPos(p);
            boolean isVDaughterEle = IsVDaughterEle(p);
            boolean isRecoil = IsRecoil(p);
            
            //Is particle one of the particles of interest?
            if(!isVDaughterEle && !isVDaughterPos && !isRecoil) continue;
            List<SimTrackerHit> trackerhits = entry.getValue().getFirst();
            List<SimCalorimeterHit> ecalhits = entry.getValue().getSecond();

            //Is particle within tracker acceptance?
            boolean findableTrack = IsTrackFindable(trackerhits);   
            if(!findableTrack) continue;
            
            if(isVDaughterEle) findableEle[0] = true;
            if(isVDaughterPos) findablePos[0] = true;
            if(isRecoil) findableEleRec[0] = true;
            
            if(ecalhits.size() != 0){
            	double energy = 0;
            	//Find position in ecal of particle hit
            	//This needs to be checked
            	Hep3Vector ecalPos = ecalhits.get(0).getPositionVec();
            	for(SimCalorimeterHit hit:ecalhits){
            		for(int i = 0; i < hit.getMCParticleCount(); i++){
            			if(p != hit.getMCParticle(i)) continue;
            			if(hit.getCorrectedEnergy() > energy){
            				energy = hit.getCorrectedEnergy();
            				ecalPos = hit.getPositionVec();
            			}
            		}
            	}
            	if(p.getMomentum().y() > 0){
                	if(isVDaughterEle){
                		findableEle[1] = true;
                		ecalHitPosEle.fill(ecalPos.x(),ecalPos.y());
                		posEle = ecalPos;
                	}
                	if(isVDaughterPos){
                		findablePos[1] = true;
                		ecalHitPosPos.fill(ecalPos.x(),ecalPos.y());
                		posPos = ecalPos;
                	}
                	if(isRecoil){
                		findableEleRec[1] = true;
                		ecalHitPosEleRec.fill(ecalPos.x(),ecalPos.y());
                	}
            	}
            	else{
            		if(isVDaughterEle){
            			findableEle[2] = true;
            			ecalHitPosEle.fill(ecalPos.x(),ecalPos.y());
            			posEle = ecalPos;
            		}
                	if(isVDaughterPos){
                		findablePos[2] = true;
                		ecalHitPosPos.fill(ecalPos.x(),ecalPos.y());
                		posPos = ecalPos;
                	}
                	if(isRecoil){
                		findableEleRec[2] = true;
                		ecalHitPosEleRec.fill(ecalPos.x(),ecalPos.y());
                	}
            	}
            }
        }
        
        if(findableEle[1] && posPos != null){
        	ecalHitPosPos_EleTop.fill(posPos.x(),posPos.y());
        }
        
        if(findableEle[2] && posPos != null){
        	ecalHitPosPos_EleBot.fill(posPos.x(),posPos.y());
        }
        
        if(findablePos[1] && posEle != null){
        	ecalHitPosEle_PosTop.fill(posEle.x(),posEle.y());
        }
        
        if(findablePos[2] && posEle != null){
        	ecalHitPosEle_PosBot.fill(posEle.x(),posEle.y());
        }

        //Fill appropriate histograms with the number of particles within acceptance
        if(findableEle[0] && findablePos[0]){
        	findableDecaysWithDaughters.fill(0);
        	if(findableEleRec[0]){
        		findableDecaysWithDaughtersWithRecoil.fill(0);
        	}
        	if((findableEle[1] && findablePos[1]) || (findableEle[2] && findablePos[2])){
        		findableDecaysWithDaughtersSameVolWithEcal.fill(0);
        		if(findableEleRec[0] && (findableEleRec[1] || findableEleRec[2])){
        			findableDecaysWithDaughtersWithRecoilSameVolWithEcal.fill(0);
        		}
        	}
        	if((findableEle[1] && findablePos[2]) || (findableEle[2] && findablePos[1])){
        		findableDecaysWithDaughtersOppoVolWithEcal.fill(0);
        		if(findableEleRec[0] && (findableEleRec[1] || findableEleRec[2])){
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
