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
public class SkimGoldenWab2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed;

    private int _minHitsOnTrack = 13;
    private double _cluster1MinEnergy = 1.8;
    private double _cluster1MaxEnergy = 2.8;
    private double _cluster2MinEnergy = 0.5;
    double _esumCut = 3.0;
    private boolean _requireFiducialWab = true;

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        boolean isGoldenWab = false;
        _numberOfEventsProcessed++;
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        if (ecalClusters.size() == 2) {
            Cluster c1 = ecalClusters.get(0);
            double e1 = c1.getEnergy();
            Cluster c2 = ecalClusters.get(1);
            double e2 = c2.getEnergy();
            // check on cluster energies
            if (e1 < _cluster1MaxEnergy && e1 > _cluster1MinEnergy) {
                if (e2 > _cluster2MinEnergy) {
                    Hep3Vector pos1 = new BasicHep3Vector(c1.getPosition());
                    double t1 = ClusterUtilities.getSeedHitTime(c1);
                    double t2 = ClusterUtilities.getSeedHitTime(c2);
                    //
                    double deltaT = t1 - t2;
                    if (abs(deltaT) < 2.) {
                        double esum = e1 + e2;
                        Hep3Vector pos2 = new BasicHep3Vector(c2.getPosition());
                        // opposite hemispheres
                        if (pos1.x() * pos2.x() < 0. && pos1.y() * pos2.y() < 0.) {

                            if (esum > _esumCut) {
                                // check on one electron and one photon
                                List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
                                ReconstructedParticle electron = null;
                                ReconstructedParticle photon = null;
                                if (rps.size() == 2) {
                                    for (ReconstructedParticle rp : rps) {
                                        if (rp.getParticleIDUsed().getPDG() == 11) {
                                            electron = rp;
                                        }
                                        if (rp.getParticleIDUsed().getPDG() == 22) {
                                            photon = rp;
                                        }
                                    }
                                    if (electron != null && photon != null) {
                                        Track t = electron.getTracks().get(0);
                                        int nHits = t.getTrackerHits().size();
                                        if (nHits >= _minHitsOnTrack) {
                                            isGoldenWab = true;
                                            if (_requireFiducialWab) {
                                                boolean e1IsFiducial = TriggerModule.inFiducialRegion(c1);
                                                boolean e2IsFiducial = TriggerModule.inFiducialRegion(c2);
                                                boolean wabIsFiducial = e1IsFiducial && e2IsFiducial;
                                                if (!wabIsFiducial) {
                                                    isGoldenWab = false;
                                                }
                                            }
                                        }
                                    }
                                }
                            } // end of check on esum
                        } // end of check on opposite hemispheres
                    } // end of check on deltaT
                } // end of check on cluster2 energy
            } // end of check on cluster1 energy
        } // end of check on 2 clusters
        if (isGoldenWab) {
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
        System.out.println("Selected " + _numberOfEventsSelected + " golden WABs of " + _numberOfEventsProcessed + "  events processed");
    }
}
