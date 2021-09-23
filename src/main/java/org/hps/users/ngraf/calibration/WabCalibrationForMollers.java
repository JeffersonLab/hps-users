package org.hps.users.ngraf.calibration;

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
 * Analyze WAB events with the electron in the Moller fiducial region to derive
 * calorimeter cluster energy corrections and track momentum corrections.
 *
 * @author Norman A. Graf
 */
public class WabCalibrationForMollers extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;

    protected void process(EventHeader event) {
        _numberOfEventsProcessed++;
        boolean skipEvent = true;
        // select events with two and only two clusters, one of which is an
        // electron in the Moller fiducial region and the other is a photon in
        // the calorimeter fiducial region.
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        if (ecalClusters.size() == 2) {
            List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
            if (rpList.size() == 2) {
                ReconstructedParticle electron = null;
                ReconstructedParticle photon = null;
                for (ReconstructedParticle rp : rpList) {
                    int pdgId = rp.getParticleIDUsed().getPDG();
                    if (abs(pdgId) == 11 && !rp.getClusters().isEmpty()) {
                        CalorimeterHit seed = rp.getClusters().get(0).getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        // is electron cluster in the Moller fiducial region?
                        if (abs(iy) <= 2 && ix > -16 && ix < -9) {
                            electron = rp;
                        }
                    }
                    if (abs(pdgId) == 22) {
                        // is photon cluster in the calorimeter fiducial region?
                        if (TriggerModule.inFiducialRegion(rp.getClusters().get(0))) {
                            photon = rp;
                        }
                    }
                }
                if (electron != null && photon != null) {
                    Cluster eclus = electron.getClusters().get(0);
                    Cluster pclus = photon.getClusters().get(0);
                    // are the clusters in-time?
                    double t1 = ClusterUtilities.getSeedHitTime(eclus);
                    double t2 = ClusterUtilities.getSeedHitTime(pclus);
                    if (abs(t1 - t2) < 2.) {
                        Hep3Vector pos1 = new BasicHep3Vector(eclus.getPosition());
                        Hep3Vector pos2 = new BasicHep3Vector(pclus.getPosition());
                        // opposite hemispheres
                        if (pos1.x() * pos2.x() < 0. && pos1.y() * pos2.y() < 0.) {
                            skipEvent = false;
                            double electronEnergy = electron.getEnergy();
                            double electronMomentum = electron.getMomentum().magnitude();
                            boolean electronIsTop = eclus.getPosition()[1] > 0 ? true : false;
                            String torb = electronIsTop ? "top" : "bottom";
                            Track t = electron.getTracks().get(0);
                            aida.histogram1D("electron energy", 100, 0., 5.).fill(electron.getEnergy());

                            aida.histogram1D("electron energy " + torb, 100, 0., 5.).fill(electron.getEnergy());
                            aida.histogram1D("electron momentum " + torb, 100, 0., 5.).fill(electron.getMomentum().magnitude());

                            double EoverP = electronEnergy / electronMomentum;
                            aida.profile1D("EoverP vs p profile " + torb, 50, 1., 4.).fill(electronMomentum, EoverP);

                            double photonEnergy = photon.getEnergy();
                            int ix = eclus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                            int iy = eclus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
                            aida.histogram2D("cluster ix vs iy electron", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            ix = pclus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                            iy = pclus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
                            aida.histogram2D("cluster ix vs iy photon", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            aida.histogram1D("photon energy", 100, 0., 5.).fill(photon.getEnergy());

                            double eesum = electron.getEnergy() + photon.getEnergy();
                            aida.histogram1D("esum " + torb, 100, 2., 5.).fill(eesum);
                            double pesum = electron.getMomentum().magnitude() + photon.getEnergy();
                            aida.histogram1D("pesum " + torb, 100, 2., 5.).fill(pesum);

                            // let's see if we can correct a few things
                            // functions are fit to the profile plots.
                            // this gives us E/p as a function of the measured p
                            // multiplying by the measured p brings us to the correct E
                            double[] parTop = {1.26, -0.15};
                            double[] parBottom = {1.16, -0.138};

                            double correctedTrackMomentum = 0.;
                            if (electronIsTop) {
                                correctedTrackMomentum = (parTop[0] + parTop[1] * electronMomentum) * electronMomentum;
                            } else {
                                correctedTrackMomentum = (parBottom[0] + parBottom[1] * electronMomentum) * electronMomentum;
                            }
                            // correction for the ECal energy scale 3.74/3.5
                            double eCorr = 1.07;
                            aida.histogram1D("electron momentum " + torb + " corrected", 100, 0., 5.).fill(correctedTrackMomentum * eCorr);
                            aida.histogram1D("pesum " + torb + " corrected", 100, 2., 5.).fill((correctedTrackMomentum + photon.getEnergy()) * eCorr);
                        }    // end of check on opposite hemispheres
                    } // end of check on cluster times
                }
            }
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }
}
