package org.hps.users.spaul.skims;



import java.util.List;

import org.hps.recon.filtering.EventReconFilter;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;

public class V0SkimTest extends EventReconFilter{
    private String v0CollectionName = "TargetConstrainedV0Candidates";

    public void setV0CollectionName(String val){
        this.v0CollectionName = val;
    }

    private boolean requireOppositeClusters = false;

    private double maximumTrackChisquared = Double.MAX_VALUE;
    public void setRequireOppositeClusters(boolean val){
        this.requireOppositeClusters = val;
    }

    public void setMaximumTrackChisquared(double val){
        this.maximumTrackChisquared = val;
    }

    public void process(EventHeader event){
        incrementEventProcessed();
        List<ReconstructedParticle> v0s = event.get(ReconstructedParticle.class, v0CollectionName);
        if(v0s.size() == 0)
            skipEvent();



        boolean foundGoodV0 = false;
        for(ReconstructedParticle v0 : v0s){
            if(v0.getType()<32)
                continue; //it's 2017, and we still have to do this sort of thing to check if a member of a collection is GBL or not! why???
            if(requireOppositeClusters){
                //require that both tracks be matched to a cluster,
                //and that they be on opposite sides of the Ecal.  
                if(v0.getParticles().get(0).getClusters().isEmpty())
                    continue;
                if(v0.getParticles().get(1).getClusters().isEmpty())
                    continue;

                if(v0.getParticles().get(0).getClusters().get(0).getPosition()[1]
                        *v0.getParticles().get(1).getClusters().get(0).getPosition()[1] > 0)
                    continue;
            }
            if(v0.getParticles().get(0).getTracks().get(0).getChi2() > maximumTrackChisquared 
                    || v0.getParticles().get(1).getTracks().get(0).getChi2() > maximumTrackChisquared )
                continue;
            foundGoodV0 = true;
        }
        if(!foundGoodV0)
            skipEvent();



        incrementEventPassed();
    }
}
