package org.hps.users.mrsolt;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;

/**
 * This is the general tuple driver
 * Use TupleTemplateDriver to add you code and variables
 * Run this driver to output variables into a text file
 * Specify your driver name in the steering file
 * Run "makeTree.py" to create a root tuple
 *
 * @author mrsolt on Aug 31, 2017
 */

public class GeneralTupleDriverPhotonPosition extends Driver {

    protected String tupleFile = null;
    protected PrintWriter tupleWriter = null;
    protected final static List<String> tupleVariables = new ArrayList<String>();
    protected final static Map<String, Double> tupleMap = new HashMap<String, Double>();
    protected String driverName = "TupleTemplateDriver";

    //abstract protected void setupVariables();

    @Override
    protected void detectorChanged(Detector detector) {
        setupVariables();
        if (tupleFile != null) {
            try {
                tupleWriter = new PrintWriter(tupleFile);
            } catch (FileNotFoundException e) {
                tupleWriter = null;
            }
            tupleWriter.println(StringUtils.join(tupleVariables, ":"));
        }
    }

    public void setDriverName(String driverName) {
        this.driverName = driverName;
    }
    
    @Override
    public void endOfData() {
        if (tupleWriter != null) {
            tupleWriter.close();
        }
    }

    protected void writeTuple() {
        for (String variable : tupleVariables) {
            Double value = tupleMap.get(variable);
            if (value == null || Double.isNaN(value)) {
                value = -9999.0;
            }
            if (variable.endsWith("/I") || variable.endsWith("/B")) {
                tupleWriter.format("%d\t", Math.round(value));
            } else {
                tupleWriter.format("%g\t", value);
            }
        }
        tupleWriter.println();
    }

    public void setTupleFile(String tupleFile) {
        this.tupleFile = tupleFile;
    }
    
    protected void addGeneralVariables() {
    	PhotonPositionDriver.addVariables();
    }
    
    protected void fillGeneralVariables(EventHeader event, MCParticle p) {
    	PhotonPositionDriver.fillVariables(event,p);
    }
    
    protected void setupVariables() {
        tupleVariables.clear();
        addGeneralVariables();
    }
    
    @Override
    public void process(EventHeader event) {         
        
        List<MCParticle> particles = event.get(MCParticle.class, "MCParticle"); 
        
        for(MCParticle p:particles){
        	tupleMap.clear();
        	fillGeneralVariables(event,p);
        	if (tupleWriter != null) writeTuple(); 
        }
    }
}
