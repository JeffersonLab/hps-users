package org.hps.users.ngraf.skim.golden;

import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimGoldenFee2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed;

    private int _maxNClusters = 1;
    private double _minClusterEnergy = 3.0;
    private double _minSeedHitEnergy = 1.9;
    private int _minHitsOnTrack = 13;
    private boolean _requireFiducialFee = true;

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        boolean isGoldenFee = false;
        _numberOfEventsProcessed++;
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        if (ecalClusters.size() <= _maxNClusters) {
            List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
            if (rps.size() <= _maxNClusters) {
                for (ReconstructedParticle rp : rps) {
                    if (rp.getParticleIDUsed().getPDG() == 11) {
                        if (!rp.getClusters().isEmpty()) {
                            Cluster c = rp.getClusters().get(0);
                            if (c.getEnergy() >= _minClusterEnergy) {
                                double seedEnergy = ClusterUtilities.findSeedHit(c).getCorrectedEnergy();
                                if (seedEnergy > _minSeedHitEnergy) {
                                    Track t = rp.getTracks().get(0);
                                    int nHits = t.getTrackerHits().size();
                                    if (nHits >= _minHitsOnTrack) {
                                        isGoldenFee = true;
                                        boolean clusterIsFiducial = TriggerModule.inFiducialRegion(c);
                                        if (_requireFiducialFee && !clusterIsFiducial) {
                                            isGoldenFee = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (isGoldenFee) {
            skipEvent = false;
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " golden FEEs of " + _numberOfEventsProcessed + "  events processed");
    }
}
