package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.HepLorentzVector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackData;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.Vertex;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.fourvec.Lorentz4Vector;
import org.lcsim.util.fourvec.Momentum4Vector;

/**
 *
 * @author Norman A. Graf
 */
public class MollerAnalysis2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private boolean _analyzeMollerCollections = false;
    private String[] mollerCollectionNames = {"UnconstrainedMollerVertices", "UnconstrainedMollerVertices_KF"};

    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles_KF"};//, "FinalStateParticles_KF"};

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean _skipEvent = true;

    final double pScale = 1.761 / 2.09;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    @Override
    protected void process(EventHeader event) {
        _skipEvent = true;
        _numberOfEventsProcessed++;
        // simple cluster analysis
        aida.tree().mkdirs("Cluster analysis ");
        aida.tree().cd("Cluster analysis ");
        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        int nClusters = clusters.size();
        Cluster clus = null;
        aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
        int nFiducialClusters = 0;
        for (Cluster cluster : clusters) {
            aida.histogram2D("Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            if (cluster.getPosition()[1] > 0.) {
                aida.histogram1D("Top cluster energy ", 100, 0., 5.).fill(cluster.getEnergy());
            } else {
                aida.histogram1D("Bottom cluster energy ", 100, 0., 5.).fill(cluster.getEnergy());
            }
            if (nClusters < 4) {
                aida.histogram1D("cluster energy " + nClusters + " clusters", 100, 0., 5.).fill(cluster.getEnergy());
            }
            CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
            int ix = seed.getIdentifierFieldValue("ix");
            int iy = seed.getIdentifierFieldValue("iy");
            double energy = cluster.getEnergy();
            aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            if (abs(iy) == 1 && ix > -15 && ix < -9) {
                clus = cluster;
                nFiducialClusters++;
                aida.histogram1D("cluster energy fiducial", 100, 0., 3.0).fill(energy);
                aida.histogram2D("cluster ix vs iy fiducial", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            }
        }
        aida.histogram1D("number of Fiducial Clusters", 5, -0.5, 4.5).fill(nFiducialClusters);
        aida.tree().cd("..");
        // continue only with events with one fiducial cluster with energy > 0.8GeV
        if (nFiducialClusters == 1 && clus.getEnergy() > 0.8) {

            // my own event analysis
            for (String rpCollectionName : ReconstructedParticleCollectionNames) {
//                System.out.println("analyzing " + rpCollectionName);
                aida.tree().mkdirs("Moller analysis " + rpCollectionName);
                aida.tree().cd("Moller analysis " + rpCollectionName);
                aida.histogram1D("number of clusters", 5, -0.5, 4.5).fill(nClusters);
                aida.histogram1D("cluster energy fiducial", 100, 0., 3.0).fill(clus.getEnergy());

                // electron with cluster
                ReconstructedParticle electron1 = null;
                // electron with only track, will have to match time and py
                ReconstructedParticle electron2 = null;
                List<ReconstructedParticle> electrons = new ArrayList<>();
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, rpCollectionName);
                for (ReconstructedParticle rp : rpList) {
                    int pdgId = rp.getParticleIDUsed().getPDG();
                    if (pdgId == 11) {
                        // let's find the electron of interest that triggered...
                        if (!rp.getClusters().isEmpty()) {
                            if (rp.getClusters().get(0) == clus) {
                                // require the track to have 12 or more hits...
                                if (rp.getTracks().get(0).getTrackerHits().size() >= 11) {
                                    electron1 = rp;
                                }
                            }
                        }
                        if (rp != electron1) {
                            if (rp.getTracks().get(0).getTrackerHits().size() >= 12) {
                                electrons.add(rp);
                            }
                        }
                    }
                }
                // if we found a track associated to the cluster let's see if we can find another electron...
                if (electron1 != null) {

                    double clusterEnergy = clus.getEnergy();
                    boolean isTop = clus.getPosition()[1] > 0 ? true : false;
                    double clusterTime = ClusterUtilities.findSeedHit(clus).getTime();
//                    System.out.println("found electron matched to cluster");
//                    System.out.println("found " + electrons.size() + " other electrons");
                    RelationalTable trackToData = getKFTrackDataRelations(event);
                    Track track1 = electron1.getTracks().get(0);
                    double track1_time = ((GenericObject) trackToData.from(track1)).getFloatVal(0);
                    double track1_momentum = electron1.getMomentum().magnitude();
                    // require track and cluster to be in time...
                    aida.histogram1D("electron cluster time", 100, 30., 50.).fill(clusterTime);
                    aida.histogram1D("track1 time", 100, -40., 40.).fill(track1_time);
                    aida.histogram1D("electron1 track - cluster time", 100, -70., -20.).fill(track1_time - clusterTime);

//                    if (track1_time - clusterTime < -45.) { // track-cluster timing can differ by 24ns depending on trigger synch.
                    aida.histogram1D("number of additional electrons found", 5, -0.5, 4.5).fill(electrons.size());
                    aida.histogram1D("electron1 track time - cluster time", 100, -100., 0.).fill(track1_time - clusterTime);
                    aida.histogram1D("electron1 cluster energy", 100, 0., 3.0).fill(clus.getEnergy());
                    aida.histogram1D("electron1 track momentum", 100, 0., 3.0).fill(track1_momentum);
                    aida.histogram2D("electron1 cluster energy vs track momentum", 100, 0., 3.0, 100, 0., 3.0).fill(clus.getEnergy(), track1_momentum);
                    aida.histogram1D("electron1 number of hits on track", 20, -0.5, 19.5).fill(track1.getTrackerHits().size());
                    if (isTop) {
                        aida.histogram1D("electron1 track time - cluster time top ", 100, -100., 0.).fill(track1_time - clusterTime);
                    } else {
                        aida.histogram1D("electron1 track time - cluster time bottom", 100, -100., 0.).fill(track1_time - clusterTime);
                    }

                    // let's look for the other electron
                    if (electrons.size() == 1) {
                        for (ReconstructedParticle rp2 : electrons) {
                            double psum = electron1.getMomentum().magnitude() + rp2.getMomentum().magnitude();
                            double psumCorr = 0.;
                            if (isTop) {
                                psumCorr = electron1.getMomentum().magnitude() + rp2.getMomentum().magnitude() * pScale;
                            } else {
                                psumCorr = electron1.getMomentum().magnitude() * pScale + rp2.getMomentum().magnitude();
                            }
                            aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                            aida.histogram1D("psumCorr", 100, 0., 7.).fill(psumCorr);
                            Track track2 = rp2.getTracks().get(0);
                            double track2_time = ((GenericObject) trackToData.from(track2)).getFloatVal(0);
                            double track2_momentum = rp2.getMomentum().magnitude();
                            aida.histogram1D("electron2 track momentum", 100, 0., 5.0).fill(track2_momentum);
                            aida.histogram1D("track2 time", 100, -40., 40.).fill(track2_time);
                            aida.histogram2D("track1 time vs track2 time", 50, -40., 40., 50, -40., 40.).fill(track1_time, track2_time);
                            aida.histogram1D("track delta time", 100, -100., 100.).fill(track1_time - track2_time);
                            aida.histogram1D("track sum pY", 100, -1., 1.).fill(electron1.getMomentum().y() + rp2.getMomentum().y());
                            aida.histogram1D("electron2 number of hits on track", 20, -0.5, 19.5).fill(track2.getTrackerHits().size());
                            aida.histogram1D("electron1 track time - cluster time after match", 100, -100., 0.).fill(track1_time - clusterTime);
                            aida.histogram2D("track delta time vs track1-cluster delta time", 100, -100., 0., 100, -100., 100.).fill(track1_time - clusterTime, track1_time - track2_time);
                            aida.histogram2D("track delta time vs track sum pY", 50, -10., 10, 50, -0.1, 0.1).fill(track1_time - track2_time, electron1.getMomentum().y() + rp2.getMomentum().y());
                            // cut on sumPY and track delta time...
                            if (abs(electron1.getMomentum().y() + rp2.getMomentum().y()) < 0.04) {
                                if (track2_momentum > 1.6 && track2_momentum < 3.2) {
                                    aida.histogram1D("track delta time after sumPy and track2 momentum cuts", 100, -100., 100.).fill(track1_time - track2_time);
                                    // signal should be in-time
                                    if (abs(track1_time - track2_time) < 4.0) {
                                        aida.histogram1D("electron1 track momentum final", 100, 0., 5.0).fill(track1_momentum);
                                        aida.histogram1D("electron2 track momentum final", 100, 0., 5.0).fill(track2_momentum);
                                        aida.histogram2D("electron1 track momentum vs electron2 track momentum final", 100, 0., 5.0, 100, 0., 5.0).fill(track1_momentum, track2_momentum);
                                        aida.histogram1D("psum final", 100, 0., 7.).fill(psum);
                                        aida.histogram1D("psumCorr final", 100, 0., 7.).fill(psumCorr);
                                        aida.histogram1D("track delta time final", 100, -100., 100.).fill(track1_time - track2_time);
                                        aida.histogram1D("track delta time final finescale", 50, -10., 10.).fill(track1_time - track2_time);
                                        aida.histogram1D("track sum pY final", 50, -0.1, 0.1).fill(electron1.getMomentum().y() + rp2.getMomentum().y());
                                        aida.histogram1D("track sum pX final", 50, 0.0, 0.2).fill(electron1.getMomentum().x() + rp2.getMomentum().x());
                                        if (isTop) {
                                            analyzeTwoElectrons(electron1, rp2);
                                        } else {
                                            analyzeTwoElectrons(rp2, electron1);
                                        }
                                        _skipEvent = false;
                                    }
                                }
                            }
                            //let's now look at out-of-time tracks keeping all the other cuts the same
                            if (abs(electron1.getMomentum().y() + rp2.getMomentum().y()) < 0.04 && (track1_time - track2_time) < -10.0) {
                                if (track2_momentum > 1.6 && track2_momentum < 3.2) {
                                    aida.tree().mkdirs("out of time track2");
                                    aida.tree().cd("out of time track2");
                                    aida.histogram1D("electron1 track momentum final", 100, 0., 5.0).fill(track1_momentum);
                                    aida.histogram1D("electron2 track momentum final", 100, 0., 5.0).fill(track2_momentum);
                                    aida.histogram2D("electron1 track momentum vs electron2 track momentum final", 100, 0., 5.0, 100, 0., 5.0).fill(track1_momentum, track2_momentum);
                                    aida.histogram1D("psum final", 100, 0., 7.).fill(psum);
                                    aida.histogram1D("psumCorr final", 100, 0., 7.).fill(psumCorr);
                                    aida.histogram1D("track delta time final", 100, -100., 100.).fill(track1_time - track2_time);
                                    aida.histogram1D("track sum pY final", 50, -0.1, 0.1).fill(electron1.getMomentum().y() + rp2.getMomentum().y());
                                    aida.histogram1D("track sum pX final", 50, 0.0, 0.2).fill(electron1.getMomentum().x() + rp2.getMomentum().x());
                                    if (isTop) {
                                        analyzeTwoElectrons(electron1, rp2);
                                    } else {
                                        analyzeTwoElectrons(rp2, electron1);
                                    }
                                    aida.tree().cd("..");
                                }// end of cut on track2 momentum
                            }// end of out-of-time analysis
                        }//loop over other electrons
                    }//end of check on one and only one other electron
                }
//                }
                aida.tree().cd("..");
            }
        }// end of cut on fiducial cluster
        // analysis of Moller collections
        if (_analyzeMollerCollections) {
            for (String mollerCollectionName : mollerCollectionNames) {
                aida.tree().mkdirs("Moller analysis " + mollerCollectionName);
                aida.tree().cd("Moller analysis " + mollerCollectionName);
                List<Vertex> vertices = event.get(Vertex.class, mollerCollectionName);
                for (Vertex v : vertices) {
                    ReconstructedParticle rp = v.getAssociatedParticle();

                    List<ReconstructedParticle> parts = rp.getParticles();
                    ReconstructedParticle rp1 = parts.get(0);
                    ReconstructedParticle rp2 = parts.get(1);
                    double p1 = rp1.getMomentum().magnitude();
                    double p2 = rp2.getMomentum().magnitude();
                    int nclusters = 0;
                    Cluster c1 = null;
                    Cluster c2 = null;
                    if (!rp1.getClusters().isEmpty()) {
                        nclusters++;
                        c1 = rp1.getClusters().get(0);
                        CalorimeterHit seed = c1.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double energy = c1.getEnergy();
                        aida.histogram2D("cluster 1 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram1D("cluster 1 energy", 100, 0., 3.0).fill(energy);
                    }
                    if (!rp2.getClusters().isEmpty()) {
                        nclusters++;
                        c2 = rp2.getClusters().get(0);
                        CalorimeterHit seed = c2.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double energy = c2.getEnergy();
                        aida.histogram2D("cluster 2 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram1D("cluster 2 energy", 100, 0., 3.0).fill(energy);
                    }

                    //rotate into physics frame of reference
                    Hep3Vector rprot = VecOp.mult(beamAxisRotation, rp.getMomentum());

                    Hep3Vector p1rot = VecOp.mult(beamAxisRotation, rp1.getMomentum());
                    Hep3Vector p2rot = VecOp.mult(beamAxisRotation, rp2.getMomentum());
                    double theta1 = Math.acos(p1rot.z() / p1rot.magnitude());
                    double theta2 = Math.acos(p2rot.z() / p2rot.magnitude());
                    double thetasum = theta1 + theta2;

                    aida.histogram1D("invariant mass", 100, 0.02, 0.1).fill(rp.getMass());
                    aida.histogram1D("px", 100, -0.1, 0.1).fill(rprot.x());
                    aida.histogram1D("py", 100, -0.05, 0.05).fill(rprot.y());
                    aida.histogram1D("pz", 100, 2.5, 5.0).fill(rprot.z());
                    aida.histogram1D("p", 100, 2.5, 5.0).fill(rprot.magnitude());
                    aida.histogram1D("number of clusters", 3, -0.5, 2.5).fill(nclusters);

                    aida.histogram1D("theta sum", 100, 0.02, 0.06).fill(thetasum);
                    aida.histogram2D("p1 vs theta1", 100, 1.2, 3., 100, 0.012, 0.024).fill(p1, theta1);
                    aida.histogram2D("p2 vs theta2", 100, 1.2, 3., 100, 0.012, 0.024).fill(p2, theta2);
                    aida.histogram2D("p1 vs p2", 100, 1.2, 3., 100, 1.2, 3.).fill(p1, p2);
                    aida.histogram2D("theta1 vs theta2", 100, 0.012, 0.024, 100, 0.012, 0.024).fill(theta1, theta2);
                }
                aida.tree().cd("..");
            }
        }
        if (_skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }

    private void analyzeTwoElectrons(ReconstructedParticle ele1, ReconstructedParticle ele2) {
        // note that e1 is the top electron
        double emass = 0.000511;
        double emass2 = emass * emass;
        double[] p1 = ele1.getMomentum().v();
        double[] p2 = ele2.getMomentum().v();
        // apply ad hoc correction to bottom momentum
        for (int i = 0; i < 2; ++i) {
            p2[i] *= p2[i] * pScale;
        }
        double e1 = 0;
        double e2 = 0.;
        for (int i = 0; i < 3; ++i) {
            e1 += p1[i] * p1[i];
            e2 += p2[i] * p2[i];
        }
        e1 = sqrt(e1 + emass2);
        e2 = sqrt(e2 + emass2);
        Momentum4Vector evec1 = new Momentum4Vector(p1[0], p1[1], p1[2], e1);
        Momentum4Vector evec2 = new Momentum4Vector(p2[0], p2[1], p2[2], e2);
        Lorentz4Vector eesum = evec1.plus(evec2);
        double eemass = eesum.mass();
        aida.histogram1D("invariant mass", 50, 0.0, 0.15).fill(eemass);
    }

    public RelationalTable getKFTrackDataRelations(EventHeader event) {

        List<TrackData> TrackData;
        RelationalTable trackToData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackRelations;
        TrackData trackdata;
        TrackData = event.get(TrackData.class, "KFTrackData");
        trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
        for (LCRelation relation : trackRelations) {
            if (relation != null && relation.getTo() != null) {
                trackToData.add(relation.getFrom(), relation.getTo());
            }
        }
        return trackToData;
    }
}
