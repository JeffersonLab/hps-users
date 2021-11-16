package org.hps.users.ngraf.skim.golden;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
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
public class SkimGoldenTrident2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed;

    private int _minHitsOnTrack = 13;
    private double _minTridentClusterEnergy = 0.5;
    private double _maxTridentClusterEnergy = 2.5;

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        boolean isGoldenTrident = false;
        _numberOfEventsProcessed++;

        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        if (ecalClusters.size() == 3) {
            boolean allFiducial = true;
            for (Cluster cluster : ecalClusters) {
                if (!TriggerModule.inFiducialRegion(cluster)) {
                    allFiducial = false;
                }
            }
            if (allFiducial) {
                List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
                if (rps.size() == 3) {
                    ReconstructedParticle electron1 = null;
                    ReconstructedParticle electron2 = null;
                    ReconstructedParticle positron = null;
                    for (ReconstructedParticle rp : rps) {
                        if (rp.getParticleIDUsed().getPDG() == 11) {
                            if (electron1 == null) {
                                electron1 = rp;
                            } else {
                                electron2 = rp;
                            }
                        }
                        if (rp.getParticleIDUsed().getPDG() == -11) {
                            positron = rp;
                        }
                    }
                    boolean allGoodTracks = true;
                    boolean allInTime = true;
                    if (electron1 != null && electron2 != null && positron != null) {
                        double posT = ClusterUtilities.getSeedHitTime(positron.getClusters().get(0));
                        double ele1T = ClusterUtilities.getSeedHitTime(electron1.getClusters().get(0));
                        double ele2T = ClusterUtilities.getSeedHitTime(electron2.getClusters().get(0));
                        if (abs(posT - ele1T) > 2.) {
                            allInTime = false;
                        }
                        if (abs(posT - ele2T) > 2.) {
                            allInTime = false;
                        }
                        if (abs(ele1T - ele2T) > 2.) {
                            allInTime = false;
                        }
                        if (allInTime) {
                            if (electron1.getTracks().get(0).getTrackerHits().size() < 11) {
                                allGoodTracks = false;
                            }
                            if (electron2.getTracks().get(0).getTrackerHits().size() < 11) {
                                allGoodTracks = false;
                            }
                            if (positron.getTracks().get(0).getTrackerHits().size() < 11) {
                                allGoodTracks = false;
                            }
                        }
                        if (allGoodTracks && allInTime) {
                            double pysum = 0.;
                            double esum = 0.;
                            pysum += electron1.getMomentum().y();
                            pysum += electron2.getMomentum().y();
                            pysum += positron.getMomentum().y();
                            esum += electron1.getEnergy();
                            esum += electron2.getEnergy();
                            esum += positron.getEnergy();
                            if (esum > 2. && esum < 5. && abs(pysum) < 0.2) {
                                isGoldenTrident = true;
                            }
                        }
                    }
                } // end of check on 3 rps
            }// end of check on all fiducial clusters
        } // end of check on 3 clusters
        if (isGoldenTrident) {
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
        System.out.println("Selected " + _numberOfEventsSelected + " golden Tridents of " + _numberOfEventsProcessed + "  events processed");
    }
}
