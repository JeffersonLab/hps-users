package org.hps.users.mrsolt;

import java.util.Arrays;
import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;

/**
 * This is the tuple template driver
 * Use this to add your code and variables to make a tuple
 * Run the GeneralTupleDriver to output info into a text file
 * Change the steering file to include this driver
 * Run "makeTree.py" on text file to create a root tuple
 *
 * @author mrsolt on Aug 31, 2017
 */

public class PhotonPositionDriver extends GeneralTupleDriverPhotonPosition {
	
    public static void fillVariables(EventHeader event, MCParticle p) {
        
        if(p.getPDGID() == 22){
        	tupleMap.put("run/I", (double) event.getRunNumber());
            tupleMap.put("event/I", (double) event.getEventNumber());
            tupleMap.put("energy/D",p.getEnergy());
        	tupleMap.put("startX/D",p.getOrigin().x());
        	tupleMap.put("startY/D",p.getOrigin().y());
        	tupleMap.put("startZ/D",p.getOrigin().z());
        	tupleMap.put("Px/D",p.getPX());
        	tupleMap.put("Py/D",p.getPY());
        	tupleMap.put("Pz/D",p.getPZ());
        	tupleMap.put("genStatus/I", (double) p.getGeneratorStatus());
        }
    }
    
    public static void addVariables() {
    	
        String[] newVars = new String[]{"run/I", "event/I",
        		"energy/D","startX/D","startY/D","startZ/D","endX/D","endY/D","endZ/D",
        		"Px/D","Py/D","Pz/D","parentID/I","genStatus/I"};
        tupleVariables.addAll(Arrays.asList(newVars));
    }
}