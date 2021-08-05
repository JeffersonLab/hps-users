package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.atan2;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.vertexing.BilliorTrack;
import org.hps.recon.vertexing.BilliorVertex;
import org.hps.recon.vertexing.BilliorVertexer;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TSData2019;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class MollerAnalysis2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    BilliorVertexer vtxFitter;

    String[] finalStateParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF"};

    protected void detectorChanged(Detector detector) {
        double bfield = detector.getFieldMap().getField(new BasicHep3Vector(0., 0., -7.5)).y();
        vtxFitter = new BilliorVertexer(bfield);
        System.out.println("using bfield " + bfield);
    }

    @Override
    protected void process(EventHeader event) {

        boolean triggered = false;
        boolean inMollerRegion = false;
        if (event.hasCollection(GenericObject.class, "TSBank")) {
            List<GenericObject> triggerList = event.get(GenericObject.class, "TSBank");
            for (GenericObject data : triggerList) {
                if (AbstractIntData.getTag(data) == TSData2019.BANK_TAG) {
                    TSData2019 triggerData = new TSData2019(data);
                    triggered = triggerData.isSingle0BotTrigger() || triggerData.isSingle0TopTrigger();
                }
            }
        }
        aida.histogram1D("Triggered singles0", 2, -0.5, 1.5).fill(triggered ? 1 : 0);
        if (triggered) {
            if (event.hasCollection(Cluster.class, "EcalClustersCorr")) {
                List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
                aida.histogram1D("Number of Clusters in triggered Event", 10, 0., 10.).fill(clusters.size());
                if (clusters.size() == 1) {
                    aida.tree().mkdirs("Singles0 single cluster");
                    aida.tree().cd("Singles0 single cluster");
                    for (Cluster cluster : clusters) {
                        CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double[] cPos = cluster.getPosition();
                        double energy = cluster.getEnergy();
                        aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram2D("cluster ix vs energy", 47, -23.5, 23.5, 100, 0., 5.0).fill(ix, energy);
                        aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cPos[0], cPos[1]);
                        String half = cluster.getPosition()[1] > 0. ? "Top " : "Bottom ";
                        aida.histogram1D(half + "Cluster Energy", 100, 0., 5.).fill(cluster.getEnergy());
                        aida.histogram1D("Cluster Energy", 100, 0., 5.0).fill(cluster.getEnergy());
                        //Mollers?
                        if (energy > 0.6 && ix < -10) {
                            inMollerRegion = true;
                            aida.histogram1D("Cluster Energy signal region", 100, 0., 5.0).fill(cluster.getEnergy());
                        }
                        if (energy > 1. && energy < 2.) {
                            aida.histogram2D("cluster ix vs iy 1. GeV < e <2. GeV ", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            aida.histogram2D("cluster ix vs energy 1. GeV < e <2. GeV ", 47, -23.5, 23.5, 100, 0., 5.0).fill(ix, energy);
                            aida.histogram2D("cluster x vs y 1. GeV < e <2. GeV ", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cPos[0], cPos[1]);
                            aida.histogram1D("Cluster Energy 1. GeV < e <2. GeV ", 100, 0., 5.0).fill(cluster.getEnergy());
                        }
                    }
                    aida.tree().cd("..");
                }
            }

            if (inMollerRegion) {
                // let's look for events with at least two electrons, one in the top and the other in the bottom
                for (String s : finalStateParticleCollectionNames) {
                    int nTopElectrons = 0;
                    List<ReconstructedParticle> topElectrons = new ArrayList<>();
                    List<ReconstructedParticle> bottomElectrons = new ArrayList<>();
                    int nBottomElectrons = 0;
                    String trackType = "No Track ";
                    if (event.hasCollection(ReconstructedParticle.class, s)) {
                        String dir = s + " ReconstructedParticle Analysis";
                        aida.tree().mkdirs(dir);
                        aida.tree().cd(dir);
                        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, s);
                        for (ReconstructedParticle rp : rpList) {
                            int pdgId = rp.getParticleIDUsed().getPDG();
                            if (pdgId == 11) {

                                if (TrackType.isGBL(rp.getType())) {
                                    trackType = "GBL ";
                                }
                                if (rp.getType() == 1) {
                                    trackType = "Kalman ";
                                }
                                Track t = rp.getTracks().get(0);
                                int trackFinderType = t.getType();
                                int nHits = t.getTrackerHits().size();
                                String id = "electron";
                                Hep3Vector pmom = rp.getMomentum();
                                double thetaY = atan2(pmom.y(), pmom.z());//asin(pmom.y() / pmom.magnitude());
                                double z0 = t.getTrackStates().get(0).getZ0();
                                String torb = pmom.y() > 0. ? " top " : " bottom ";
                                aida.histogram1D(trackType + torb + id + " track momentum", 100, 0., 7.).fill(rp.getMomentum().magnitude());
                                aida.histogram1D(trackType + torb + id + " track nHits", 20, -0.5, 19.5).fill(nHits);
                                double theta = Math.acos(pmom.z() / pmom.magnitude());
                                if (pmom.y() > 0) {
                                    nTopElectrons++;
                                    topElectrons.add(rp);
                                    aida.histogram2D("p vs theta top all", 100, 0., 2.5, 100, 0.01, 0.08).fill(pmom.magnitude(), theta);
                                    if (!rp.getClusters().isEmpty()) {
                                        aida.histogram2D("p vs theta top with cluster", 100, 0., 2.5, 100, 0.01, 0.08).fill(pmom.magnitude(), theta);
                                        aida.histogram2D("e vs theta top with cluster", 100, 0., 2.5, 100, 0.01, 0.08).fill(rp.getClusters().get(0).getEnergy(), theta);
                                    }
                                } else {
                                    nBottomElectrons++;
                                    bottomElectrons.add(rp);
                                    aida.histogram2D("p vs theta bottom all", 100, 0., 2.5, 100, 0.01, 0.08).fill(pmom.magnitude(), theta);
                                    if (!rp.getClusters().isEmpty()) {
                                        aida.histogram2D("p vs theta bottom with cluster", 100, 0., 2.5, 100, 0.01, 0.08).fill(pmom.magnitude(), theta);
                                        aida.histogram2D("e vs theta bottom with cluster", 100, 0., 2.5, 100, 0.01, 0.08).fill(rp.getClusters().get(0).getEnergy(), theta);
                                    }
                                }
                            }
                        }

                        aida.histogram2D(trackType + "nTopElectron vs nBottomElectron", 5, -0.5, 4.5, 5, -0.5, 4.5).fill(nTopElectrons, nBottomElectrons);
                        // let's start clean 
                        ReconstructedParticle electronWithCluster = null;
                        if (nTopElectrons == 1 && nBottomElectrons == 1) {
                            ReconstructedParticle topElectron = topElectrons.get(0);
                            ReconstructedParticle bottomElectron = bottomElectrons.get(0);
                            if (!topElectron.getClusters().isEmpty()) {
                                electronWithCluster = topElectron;
                            }
                            if (!bottomElectron.getClusters().isEmpty()) {
                                electronWithCluster = bottomElectron;
                            }
                            if (electronWithCluster != null) {
                                Hep3Vector momTop = topElectrons.get(0).getMomentum();
                                Hep3Vector momBottom = bottomElectrons.get(0).getMomentum();
                                double pTop = momTop.magnitude();
                                double pBottom = momBottom.magnitude();
                                // cut on track momenta to remove FEEs and poorly-reconstructed low momentum tracks
                                if (pTop > 0.5 && pTop < 2.5 && pBottom > 0.5 && pBottom < 2.5) {
                                    double psum = pTop + pBottom;
                                    double thetaTop = Math.acos(momTop.z() / momTop.magnitude());
                                    double thetaBottom = Math.acos(momBottom.z() / momBottom.magnitude());

                                    aida.histogram2D("top momentum vs bottom momentum", 100, 0., 2.5, 100, 0., 2.5).fill(pTop, pBottom);
                                    aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                                    aida.histogram1D("p Top", 100, 0., 2.5).fill(pTop);
                                    aida.histogram1D("p Bottom", 100, 0., 2.5).fill(pBottom);
                                    aida.histogram2D("thetaTop vs thetaBottom", 100, 0.01, 0.08, 100, 0.01, 0.08).fill(thetaTop, thetaBottom);
                                    aida.histogram2D("p vs theta top", 100, 0., 2.5, 100, 0.01, 0.08).fill(pTop, thetaTop);
                                    aida.histogram2D("p vs theta bottom", 100, 0., 2.5, 100, 0.01, 0.08).fill(pBottom, thetaBottom);

                                    // let's try to fit a Moller vertex here...
                                    // TODO add a few more track quality cuts
                                    // TODO when the top track is associated to the ECal cluster, use the cluster energy instead of the track momentum
                                    List<BilliorTrack> tracksToVertex = new ArrayList<>();
                                    tracksToVertex.add(new BilliorTrack(topElectron.getTracks().get(0)));
                                    tracksToVertex.add(new BilliorTrack(bottomElectron.getTracks().get(0)));
                                    BilliorVertex vtx = vtxFitter.fitVertex(tracksToVertex);
                                    double mass = vtx.getInvMass();
                                    Hep3Vector vertexPos = vtx.getPosition();
                                    aida.histogram2D("psum vs Moller mass", 50, 0., 7., 50, 0., 0.2).fill(psum, mass);
                                    aida.histogram1D("Moller mass", 100, 0., 0.2).fill(mass);
                                    aida.histogram1D("Moller vertex z", 100, -25., 25.).fill(vertexPos.z());
                                }
                            }
                        }
                        aida.tree().cd("..");
                    }
                } // loop over ReconstructedParticle collections
            } // if Ecal cluster in Moller region
        } // if triggered
    }
}
