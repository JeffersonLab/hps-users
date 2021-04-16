package org.hps.users.ngraf.skim;

import java.util.List;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimFeeClusterMissingTheTrack extends Driver {

    int _numberOfEventsSelected;
    private AIDA aida = AIDA.defaultInstance();

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
        for (ReconstructedParticle rp : rpList) {
            //photons
            if (rp.getParticleIDUsed().getPDG() == 22) {
                Cluster c = rp.getClusters().get(0);
                double e = rp.getEnergy();
                double[] cPos = c.getPosition();
                String topOrBottom = cPos[1] > 0 ? " top " : " bottom ";
                CalorimeterHit seed = c.getCalorimeterHits().get(0);
                int ix = seed.getIdentifierFieldValue("ix");
                int iy = seed.getIdentifierFieldValue("iy");
                aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                //           aida.histogram1D(ix + " " + iy + " cluster energy", 100, 0., 5.0).fill(c.getEnergy());
                aida.histogram1D(topOrBottom + " cluster energy", 100, 0., 5.0).fill(c.getEnergy());
            }
        }
    }

    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " events");
    }

}
