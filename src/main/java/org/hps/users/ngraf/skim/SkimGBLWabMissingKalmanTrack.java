package org.hps.users.ngraf.skim;

import java.util.List;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Select clean WAB candidates missing a Kalman Track
 *
 * @author Norman A. Graf
 */
public class SkimGBLWabMissingKalmanTrack extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        // require two and only two ECal clusters
        if (clusters.size() == 2) {
            // check that we have one electron and one photon using GBL tracks
            List<ReconstructedParticle> rpListGBL = event.get(ReconstructedParticle.class, "FinalStateParticles");
            if (rpListGBL.size() == 2) {
                if (rpListGBL.get(0).getParticleIDUsed().getPDG() == 11 && rpListGBL.get(1).getParticleIDUsed().getPDG() == 22) {
                    Cluster electronCluster = rpListGBL.get(0).getClusters().get(0);
                    double electronmomentum = rpListGBL.get(0).getMomentum().magnitude();
                    Cluster photonCluster = rpListGBL.get(1).getClusters().get(0);
                    // check on energy sum and energy momentum sum
                    double esum = electronCluster.getEnergy() + photonCluster.getEnergy();
                    double epsum = electronmomentum + photonCluster.getEnergy();
                    aida.histogram1D("wab e1 + e2", 100, 0., 5.).fill(esum);
                    aida.histogram1D("wab p1 + e2", 100, 0., 5.).fill(epsum);
                    String half = electronCluster.getPosition()[1] < 0 ? " bottom " : " top ";
                    aida.histogram1D("wab e1 + e2 electron " + half, 100, 0., 5.).fill(esum);
                    aida.histogram1D("wab p1 + e2 electron " + half, 100, 0., 5.).fill(epsum);

                    List<ReconstructedParticle> rpListKF = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
                    boolean electronClusterHasKfTrack = false;
                    //check whether any electron built from a KF track is associated to the GBL track electron cluster
                    for (ReconstructedParticle rp : rpListKF) {
                        if (rp.getParticleIDUsed().getPDG() == 11) {
                            for (Cluster c : rp.getClusters()) {
                                if (c == electronCluster) {
                                    electronClusterHasKfTrack = true;
                                }
                            }
                        }
                    }
                    // write out the event if the GBL electron cluster does not have an associated KF track
                    if (!electronClusterHasKfTrack) {
                        skipEvent = false;
                    }
                }
            }
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }
}
