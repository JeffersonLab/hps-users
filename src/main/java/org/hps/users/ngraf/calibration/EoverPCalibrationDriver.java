package org.hps.users.ngraf.calibration;

import static java.lang.Math.abs;
import java.util.List;
import org.hps.record.triggerbank.TriggerModule;
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
public class EoverPCalibrationDriver extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    protected void process(EventHeader event) {
        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
        for (ReconstructedParticle rp : rpList) {
            // is this RP an electron or positron?
            if (abs(rp.getParticleIDUsed().getPDG()) == 11) {
                boolean isElectron = rp.getParticleIDUsed().getPDG() == 11;
                boolean isPositron = rp.getParticleIDUsed().getPDG() == -11;
                String type = "";
                if (isElectron) {
                    type = "electron ";
                }
                if (isPositron) {
                    type = "positron ";
                }
                //does this RP have an associated cluster?
                if (!rp.getClusters().isEmpty()) {
                    Cluster cluster = rp.getClusters().get(0);
                    boolean isFiducial = TriggerModule.inFiducialRegion(cluster);
                    if (isFiducial) {
                        String fid = isFiducial ? "fiducial" : "";
                        CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double e = rp.getEnergy();
                        double p = rp.getMomentum().magnitude();
                        double eOverP = e / p;
                        int nHits = rp.getTracks().get(0).getTrackerHits().size();
                        String topOrBottom = rp.getMomentum().y() > 0 ? " top " : " bottom ";
                        aida.histogram1D(type + " " + topOrBottom + " nHits", 15, -0.5, 14.5).fill(nHits);
                        if (nHits > 13) {
                            aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                            aida.histogram2D(type + "cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            aida.histogram2D(type + "Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                            aida.histogram2D(type + "" + ix + " " + iy + " E vs EoverP " + fid, 50, 0., 5., 100, 0.5, 1.5).fill(e, eOverP);
                            aida.histogram1D(type + "" + ix + " " + iy + " EoverP " + fid, 100, 0.5, 1.5).fill(eOverP);
                            aida.histogram1D(type + "" + ix + " " + iy + " E " + fid, 50, 0., 5.).fill(e);
                            aida.histogram1D(type + "" + ix + " " + iy + " P " + fid, 50, 0., 5.).fill(p);
                        }
                    }
                }

            }
        }

    }
}
