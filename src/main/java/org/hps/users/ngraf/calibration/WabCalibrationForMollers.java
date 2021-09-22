package org.hps.users.ngraf.calibration;

import static java.lang.Math.abs;
import java.util.List;
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
                    skipEvent = false;
                    _numberOfEventsSelected++;
                    Cluster eclus = electron.getClusters().get(0);
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
                    Cluster pclus = photon.getClusters().get(0);
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
                }
            }
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }
}
