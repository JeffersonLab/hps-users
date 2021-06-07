package org.hps.users.ngraf.dataanalysis;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import java.util.List;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class WabTrackEfficiencyAnalysisDriver extends Driver {

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    private AIDA aida = AIDA.defaultInstance();
    String[] finalStateParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF"};

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        // quick loop over the two ReconstructedParticle collections.
        // kill events with a reconstructed positron associated with the "photon" cluster
        boolean skipit = false;
        String gblEventType = "";
        String kfEventType = "";

        ReconstructedParticle gblElectron = null;

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
                }
                aida.histogram1D("Cluster esum ", 100, 0., 5.5).fill(esum);

                Cluster c1 = clusters.get(0);
                Cluster c2 = clusters.get(1);
                ReconstructedParticle p1 = null;
                ReconstructedParticle p2 = null;

                String eventType = "null";
                // loop over the final state collections in the event
                for (String s : finalStateParticleCollectionNames) {
                    String trackType = s.contains("_KF") ? " KF " : " GBL ";
                    int nRpsWithCluster = 0;
                    if (event.hasCollection(ReconstructedParticle.class, s)) {
                        String dir = s + " WAB Tracking Analysis";
                        aida.tree().mkdirs(dir);
                        aida.tree().cd(dir);
                        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, s);
                        aida.histogram1D("Number of ReconstructedParticles in event", 10, -0.5, 9.5).fill(rpList.size());
                        // loop over the RPs and associate them with the two clusters
                        for (ReconstructedParticle rp : rpList) {
                            // is this rp associated with the "electron" cluster
                            List<Cluster> rpClusters = rp.getClusters();
                            if (!rpClusters.isEmpty()) {
                                if (rpClusters.get(0) == c1) {
                                    p1 = rp;
                                    nRpsWithCluster++;
                                }
                                if (rpClusters.get(0) == c2) {
                                    p2 = rp;
                                    nRpsWithCluster++;
                                }
                            }
                        }
                        // should now have two RPs that are associated with the two clusters.
                        if (p1 == null || p2 == null) {
                            System.out.println("didn't find RP associated to clusters");
                        } else {
                            // can either have e-e-, e-gamma or gammagamma
                            int pdgId1 = 0;
                            if (!p1.getParticleIDs().isEmpty()) { // why do I have to check here?
                                pdgId1 = p1.getParticleIDUsed().getPDG();
                            }
                            int pdgId2 = 0;
                            if (!p2.getParticleIDs().isEmpty()) { // why do I have to check here?
                                pdgId2 = p2.getParticleIDUsed().getPDG();
                            }
                            if (pdgId1 == 0 || pdgId2 == 0) {
                                System.out.println("didn't find pdgId for RPs");
                            }
                            if (pdgId1 == 22 && pdgId2 == 22) {
                                eventType = "gg";
                            }
                            if (pdgId1 == 22 && pdgId2 == 11) {
                                eventType = "ge";
                            }
                            if (pdgId1 == 11 && pdgId2 == 22) {
                                eventType = "eg";
                            }
                            if (pdgId1 == 11 && pdgId2 == 11) {
                                eventType = "ee";
                            }

                            aida.histogram2D("Cluster1 x vs y" + trackType + eventType, 320, -270.0, 370.0, 90, -90.0, 90.0).fill(c1.getPosition()[0], c1.getPosition()[1]);
                            aida.histogram2D("Cluster2 x vs y" + trackType + eventType, 320, -270.0, 370.0, 90, -90.0, 90.0).fill(c2.getPosition()[0], c2.getPosition()[1]);

                            aida.histogram1D("Cluster esum" + trackType + eventType, 100, 0., 5.5).fill(esum);
                            aida.histogram1D("Number of RPs with Cluster" + trackType, 10, -0.5, 9.5).fill(nRpsWithCluster);
                        }
                        aida.tree().cd("..");
                    }
                    if (trackType.equals(" GBL ")) {
                        gblEventType = eventType;
                        if (eventType.equals("ge")) {
                            gblElectron = p2;
                        }
                        if (eventType.equals("eg")) {
                            gblElectron = p1;
                        }
                    }
                    if (trackType.equals(" KF ")) {
                        kfEventType = eventType;
                    }
                }// end of loop over final state collections
                //
                // analyze events where GBL found track, KF did not.
                String dir = "GBL track without KF track";
                if (gblEventType.equals("ge") || gblEventType.equals("eg")) {
                    if (kfEventType == "gg") {
                        skipEvent = false;
                        aida.tree().mkdirs(dir);
                        aida.tree().cd(dir);
                        aida.histogram1D("Cluster esum", 100, 0., 5.5).fill(esum);
                        Cluster gblElectronCluster = gblElectron.getClusters().get(0);
                        String topOrBottom = gblElectronCluster.getPosition()[1] > 0. ? " top " : " bottom ";
                        double e = gblElectron.getEnergy();
                        double p = gblElectron.getMomentum().magnitude();
                        aida.histogram2D("Cluster with GBL Track missing KF Track x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(gblElectronCluster.getPosition()[0], gblElectronCluster.getPosition()[1]);
                        aida.histogram1D("Cluster with GBL Track missing KF Track energy", 100, 0., 5.5).fill(gblElectronCluster.getEnergy());
                        aida.histogram1D("Cluster with GBL Track missing KF Track energy" + topOrBottom, 100, 0., 5.5).fill(gblElectronCluster.getEnergy());
                        aida.histogram1D("Cluster with GBL Track missing KF Track EoverP" + topOrBottom, 100, 0., 2.).fill(e / p);
                        aida.tree().cd("..");
                    }
                }
            } // end of check on two clusters
        } // end of check on skipit
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        IAnalysisFactory af = IAnalysisFactory.create();
        IHistogramFactory hf = af.createHistogramFactory(af.createTreeFactory().create());
        ITree tree = aida.tree();
        IHistogram1D gblGG = (IHistogram1D) aida.tree().find("/FinalStateParticles WAB Tracking Analysis/Cluster esum GBL gg");
        IHistogram1D kfGG = (IHistogram1D) aida.tree().find("/FinalStateParticles_KF WAB Tracking Analysis/Cluster esum KF gg");
        IHistogram1D GGesum = (IHistogram1D) aida.tree().find("Cluster esum ");
        System.out.println("GBL gg all entries " + gblGG.allEntries());
        IHistogram1D gblInefficiency = hf.divide("GBL tracking inefficiency", gblGG, GGesum);
        IHistogram1D kfInefficiency = hf.divide("KF tracking inefficiency", kfGG, GGesum);

        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");

    }
}
