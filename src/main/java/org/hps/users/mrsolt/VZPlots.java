/**
 * Plotting VZ variables
 */
/**
 * @author mrsolt
 *
 */
package org.hps.users.mrsolt;

import java.util.List;
import hep.aida.IAnalysisFactory;
import hep.aida.IHistogramFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.ITree;
import hep.physics.vec.Hep3Vector;

import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.recon.particle.ReconParticleDriver;

public class VZPlots extends Driver {


    // Use JFreeChart as the default plotting backend
    static { 
        hep.aida.jfree.AnalysisFactory.register();
    }

    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree; 
    IHistogramFactory histogramFactory; 
    
    //Sample histograms
    IHistogram1D uncVZ;
    IHistogram1D uncV0prob;
    IHistogram1D eleTrackChisq;
    IHistogram1D posTrackChisq;
    IHistogram2D uncVZ_uncV0prob; 
    IHistogram2D uncVZ_eleTrackChisq; 
    IHistogram2D uncVZ_posTrackChisq; 
    
    //Histogram Settings
    double minVZ = 0;
    double maxVZ = 100;
    double maxV0prob = 0.01;
    int nBins = 100;
    
    //Collection Strings
    private String unconstrainedV0CandidatesColName = "UnconstrainedV0Candidates";
    
    //Configurable Variables
    boolean isMC = false;

    public void setIsMC(boolean isMC) { 
        this.isMC = isMC;
    }
    
    public void setMaxV0prob(double maxV0prob) { 
        this.maxV0prob = maxV0prob;
    }
    
    //Beam Energy
    double ebeam;
    
    public void detectorChanged(Detector detector){
        aida.tree().cd("/");
        tree = aida.tree();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);
    
        //Set Beam Energy
        BeamEnergyCollection beamEnergyCollection = 
                this.getConditionsManager().getCachedConditions(BeamEnergyCollection.class, "beam_energies").getCachedData();        
        ebeam = beamEnergyCollection.get(0).getBeamEnergy();      
   
        
        //Setup the Histograms
        uncVZ = aida.histogram1D("Reconstructed Z", nBins, minVZ, maxVZ);
        uncV0prob = aida.histogram1D("Vertex Probability", nBins, 0, maxV0prob);
        eleTrackChisq = aida.histogram1D("Electron Track Chisq", nBins, 0, 100);
        posTrackChisq = aida.histogram1D("Positron Track Chisq", nBins, 0, 100);
        uncVZ_uncV0prob = aida.histogram2D("Reconstructed Z vs Vertex Proability", nBins, 0, maxV0prob, nBins, minVZ, maxVZ);
        uncVZ_eleTrackChisq = aida.histogram2D("Reconstructed Z vs Electron Track Chisq", nBins, 0, 100, nBins, minVZ, maxVZ);
        uncVZ_posTrackChisq = aida.histogram2D("Reconstructed Z vs Positron Track Chisq", nBins, 0, 100, nBins, minVZ, maxVZ);
    }

    public void process(EventHeader event){
        aida.tree().cd("/");	
        
        List<ReconstructedParticle> unConstrainedV0List = event.get(ReconstructedParticle.class, unconstrainedV0CandidatesColName);
        
        for(ReconstructedParticle v0 : unConstrainedV0List){
            Hep3Vector v0Pos = v0.getStartVertex().getPosition();
            double v0prob = v0.getStartVertex().getProbability();
            ReconstructedParticle electron = v0.getParticles().get(ReconParticleDriver.ELECTRON);
            ReconstructedParticle positron = v0.getParticles().get(ReconParticleDriver.POSITRON);
            double eleTrkChisq = electron.getTracks().get(0).getChi2();
            double posTrkChisq = positron.getTracks().get(0).getChi2();
            
            uncVZ.fill(v0Pos.z());
            uncV0prob.fill(v0prob);
            eleTrackChisq.fill(eleTrkChisq);
            posTrackChisq.fill(posTrkChisq);
            uncVZ_uncV0prob.fill(v0prob,v0Pos.z());
            uncVZ_eleTrackChisq.fill(eleTrkChisq,v0Pos.z());
            uncVZ_posTrackChisq.fill(posTrkChisq,v0Pos.z());
        }
    }

    protected void startOfData() { 
    	System.out.println("Print Start of Data Variables");
    }
    
    public void endOfData(){
        System.out.println("Print End of Data Variables");
    }
}