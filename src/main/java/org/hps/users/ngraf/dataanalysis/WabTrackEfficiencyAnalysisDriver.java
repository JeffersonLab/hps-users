package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import static java.lang.Math.atan2;
import java.util.List;
import org.hps.recon.tracking.TrackType;
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
public class WabTrackEfficiencyAnalysisDriver extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    String[] finalStateParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF"};

    protected void process(EventHeader event) {
        // quick loop over the two ReconstructedParticle collections.
        // kill events with a reconstructed positron associated with the "photon" cluster
        boolean skipit = false;
        for (String s : finalStateParticleCollectionNames) {
            if (event.hasCollection(ReconstructedParticle.class, s)) {
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, s);
                for (ReconstructedParticle rp : rpList) {
                    if (!rp.getParticleIDs().isEmpty()) { // why do I have to check here?
                        int pdgId = rp.getParticleIDUsed().getPDG();
                        if (pdgId == -11 && !rp.getClusters().isEmpty()) {
                            skipit = true;
                        }
                    } else {
                        skipit = true; // bail if I can't use pdgId
                    }
                }
            }
        }
        if (!skipit) {
            List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
            int nClusters = clusters.size();
            aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
            if (nClusters == 2) {
                Cluster electronCluster = null;
                double esum = 0.;
                for (Cluster cluster : clusters) {
                    esum += cluster.getEnergy();
                    aida.histogram2D("Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                    if (cluster.getPosition()[1] > 0.) {
                        aida.histogram1D("Top cluster energy ", 100, 0., 5.5).fill(cluster.getEnergy());
                    } else {
                        aida.histogram1D("Bottom cluster energy ", 100, 0., 5.5).fill(cluster.getEnergy());
                    }
                    // identify electron cluster as the one on the electron side of the calorimeter
                    if (cluster.getPosition()[0] < 0.) {
                        electronCluster = cluster;
                    }
                }
                aida.histogram1D("Cluster esum ", 100, 0., 5.5).fill(esum);

                for (String s : finalStateParticleCollectionNames) {
                    String trackType = s.contains("_KF") ? " KF " : " GBL ";
                    if (event.hasCollection(ReconstructedParticle.class, s)) {
                        String dir = s + " WAB Tracking Analysis";
                        aida.tree().mkdirs(dir);
                        aida.tree().cd(dir);
                        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, s);
                        boolean foundElectronTrack = false;
                        Track electronTrack = null;
                        ReconstructedParticle foundElectron = null;
                        for (ReconstructedParticle rp : rpList) {
                            // is this rp associated with the "electron" cluster
                            List<Cluster> rpClusters = rp.getClusters();
                            if (!rpClusters.isEmpty()) {
                                if (rpClusters.get(0) == electronCluster) {
                                    // did we find an associated track?
                                    if (!rp.getTracks().isEmpty()) {
                                        foundElectronTrack = true;
                                        electronTrack = rp.getTracks().get(0);
                                        foundElectron = rp;
                                    }
                                }
                            }
                        }
                        String type = foundElectronTrack ? " track found " : " track not found ";
                        aida.histogram2D("Cluster x vs y" + trackType + type, 320, -270.0, 370.0, 90, -90.0, 90.0).fill(electronCluster.getPosition()[0], electronCluster.getPosition()[1]);
                        aida.histogram1D("Cluster esum" + trackType + type, 100, 0., 5.5).fill(esum);

                        aida.tree().cd("..");
                    }
                }// end of loop over final state collections
            } // end of check on two clusters
        } // end of check on skipit
    }
}
