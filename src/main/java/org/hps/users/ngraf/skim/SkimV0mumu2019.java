package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.List;
import java.util.Map;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackType;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.fourvec.Lorentz4Vector;
import org.lcsim.util.fourvec.Momentum4Vector;

/**
 * Skim off events consistent with a Vertex object composed of mu+ mu- muon ID
 * relies on a MIP cluster in the ECal, so by definition both
 * ReconstructedParticles must have an Ecal cluster
 *
 * @author Norman A. Graf
 */
public class SkimV0mumu2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;
    String[] vertexCollectionNames = {"UnconstrainedV0Vertices", "UnconstrainedV0Vertices_KF"};
    private static final double muMass = 0.10566;
    private static final double mumass2 = muMass * muMass;
    private double _maxMuClusterEnergy = 0.45;
    private double _clusterDeltaTimeCut = 5.0;
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private boolean _selectOnMass = true;
    private double _minMassCut = 0.5;
    private double _maxMassCut = 0.85;

    double[] p1 = new double[4];
    double[] p2 = new double[4];
    double[] pV = new double[3];
    double[] particleMasses = {0.000511, 0.10566, 0.13957, 0.493667};
    String[] particleNames = {"electron", "muon", "pion", "kaon"};

    protected void detectorChanged(Detector detector) {
        if (_minMassCut >= _maxMassCut) {
            throw new RuntimeException("Incorrect mass selection range");
        }
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    public void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;
        for (String vertexCollectionName : vertexCollectionNames) {
//            System.out.println(vertexCollectionName);
            if (event.hasCollection(Vertex.class, vertexCollectionName)) {

//                System.out.println("found " + vertexCollectionName + " collection");
                List<Vertex> vertices = event.get(Vertex.class, vertexCollectionName);
//                System.out.println("found " + vertices.size() + " vertices");
                aida.tree().mkdirs(vertexCollectionName + "_eemumu");
                aida.tree().cd(vertexCollectionName + "_eemumu");
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                for (Vertex v : vertices) {
                    ReconstructedParticle v0 = v.getAssociatedParticle();
                    int minNhits = 5;
//                    String trackType = "SeedTrack ";
//                    if (TrackType.isGBL(v0.getType())) {
//                        trackType = "GBL ";
//                    }
                    if (v0.getType() == 1) {
//                        trackType = "Kalman ";
                        minNhits = 10;
                    }
                    Vertex uncVert = v0.getStartVertex();
//                    Hep3Vector pVtxRot = VecOp.mult(beamAxisRotation, v0.getMomentum());
                    Hep3Vector vtxPosRot = VecOp.mult(beamAxisRotation, uncVert.getPosition());
//                    double theta = Math.acos(pVtxRot.z() / pVtxRot.magnitude());
//                    double phi = Math.atan2(pVtxRot.y(), pVtxRot.x());

                    // this always has 2 tracks.
                    List<ReconstructedParticle> trks = v0.getParticles();
                    ReconstructedParticle neg = trks.get(0);
                    ReconstructedParticle pos = trks.get(1);
                    // let's also require both to have clusters, since this will distinguish e from mu
                    if (!neg.getClusters().isEmpty() && !pos.getClusters().isEmpty()) {
                        Cluster negClus = neg.getClusters().get(0);
                        int negIX = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                        int negIY = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

                        Cluster posClus = pos.getClusters().get(0);
                        int posIX = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                        int posIY = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

                        // in time
                        double negTime = ClusterUtilities.getSeedHitTime(negClus);
                        double posTime = ClusterUtilities.getSeedHitTime(posClus);
                        double deltaTime = negTime - posTime;
                        double topMinusBottomTime = deltaTime;
                        if (negClus.getPosition()[1] < 0) {
                            topMinusBottomTime = -deltaTime;
                        }
                        aida.histogram1D("top - bottom cluster delta time", 100, -5., 5.).fill(deltaTime);
                        aida.histogram1D("cluster pair delta time", 100, -5., 5.).fill(deltaTime);
                        double negE = negClus.getEnergy();
                        double posE = posClus.getEnergy();
                        int negNhits = neg.getTracks().get(0).getTrackerHits().size();
                        int posNhits = pos.getTracks().get(0).getTrackerHits().size();
                        double negMom = neg.getMomentum().magnitude();
                        double posMom = pos.getMomentum().magnitude();
                        double negEoverP = negE / negMom;
                        double posEoverP = posE / posMom;

                        if (negNhits >= minNhits && posNhits >= minNhits && abs(deltaTime) < _clusterDeltaTimeCut) {
                            aida.histogram1D("negative track nHits", 20, 0., 20.).fill(negNhits);
                            aida.histogram1D("positive track nHits", 20, 0., 20.).fill(posNhits);
                            aida.histogram1D("negative momentum", 100, 0., 6.0).fill(negMom);
                            aida.histogram1D("positive momentum", 100, 0., 6.0).fill(posMom);
                            aida.histogram1D("v0 x", 50, -5., 5.).fill(vtxPosRot.x());
                            aida.histogram1D("v0 y", 50, -2., 2.).fill(vtxPosRot.y());
                            aida.histogram1D("v0 z", 50, -25., 0.).fill(vtxPosRot.z());
                            aida.histogram1D("v0 x neg " + negNhits + " pos " + posNhits + " hits on track", 50, -5., 5.).fill(vtxPosRot.x());
                            aida.histogram1D("v0 y neg " + negNhits + " pos " + posNhits + " hits on track", 50, -2., 2.).fill(vtxPosRot.y());
                            aida.histogram1D("v0 z neg " + negNhits + " pos " + posNhits + " hits on track", 50, -25., 0.).fill(vtxPosRot.z());
                            aida.histogram1D("v0 energy", 100, 0., 10.).fill(v0.getEnergy());
                            aida.histogram1D("v0 mass", 50, 0., 0.5).fill(v0.getMass());
                            aida.histogram2D("v0 mass vs Z vertex", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
                            aida.histogram2D("v0 mass vs Z vertex neg " + negNhits + " pos " + posNhits + " hits on track", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
                            aida.profile1D("v0 mass vs Z vertex profile", 50, 0.05, 0.25).fill(v0.getMass(), vtxPosRot.z());
                            aida.histogram1D("psum", 100, 0., 6.0).fill(negMom + posMom);
                            aida.histogram1D("psum neg " + negNhits + " pos " + posNhits + " hits on track", 100, 0., 6.0).fill(negMom + posMom);
                            aida.histogram2D("negative vs positive momentum", 100, 0., 6.0, 100, 0., 6.).fill(negMom, posMom);
                            aida.histogram1D("psum both ECal Clusters", 100, 0., 6.0).fill(negMom + posMom);
                            aida.histogram1D("esum both ECal Clusters", 100, 0., 6.0).fill(neg.getClusters().get(0).getEnergy() + pos.getClusters().get(0).getEnergy());

                            aida.histogram2D("negative E vs positive E all", 100, 0., 5., 100, 0., 5.).fill(negE, posE);
                            aida.histogram2D("negative E vs positive E lowE", 100, 0., 0.5, 100, 0., 0.5).fill(negE, posE);
                            aida.histogram1D("negative E over P", 100, 0., 2.).fill(negEoverP);
                            aida.histogram1D("positive E over P", 100, 0., 2.).fill(posEoverP);
                            aida.histogram2D("negative vs positive E over P", 100, 0., 2., 100, 0., 2.).fill(negEoverP, posEoverP);

                            aida.histogram2D("negative cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(negIX, negIY);
                            aida.histogram2D("positive cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(posIX, posIY);

                            // define mu+mu- sample...
                            String v0type = "bkgnd";
                            if (negE < _maxMuClusterEnergy && posE < _maxMuClusterEnergy) {
                                v0type = "mu+mu-";
                                Map<String, Double> vals = v.getParameters();
                                // System.out.println(vals);
                                double[] p1 = new double[4];
                                double[] p2 = new double[4];
                                double[] pV = new double[3];
                                p1[0] = vals.get("p1X");
                                p1[1] = vals.get("p1Y");
                                p1[2] = vals.get("p1Z");
                                p2[0] = vals.get("p2X");
                                p2[1] = vals.get("p2Y");
                                p2[2] = vals.get("p2Z");
                                double v0p = vals.get("V0P");
                                pV[0] = vals.get("V0Px");
                                pV[1] = vals.get("V0Py");
                                pV[2] = vals.get("V0Pz");
                                double mu1 = 0;
                                double mu2 = 0.;
                                for (int i = 0; i < 3; ++i) {
                                    mu1 += p1[i] * p1[i];
                                    mu2 += p2[i] * p2[i];
                                }
                                mu1 = sqrt(mu1 + mumass2);
                                mu2 = sqrt(mu2 + mumass2);
                                Momentum4Vector fourVec1 = new Momentum4Vector(p1[0], p1[1], p1[2], mu1);
                                Momentum4Vector fourVec2 = new Momentum4Vector(p2[0], p2[1], p2[2], mu2);
                                Lorentz4Vector mumusum = fourVec1.plus(fourVec2);
                                double mumumass = mumusum.mass();
                                aida.histogram1D("v0 mu+mu- mass " + v0type, 100, 0., 1.0).fill(mumumass);
                                aida.histogram1D("v0 mu+mu- high mass " + v0type, 100, 0.275, 1.0).fill(mumumass);

                                // analyze with a few different particle hypotheses
                                for (int j = 0; j < particleMasses.length; ++j) {
                                    double mass2 = particleMasses[j] * particleMasses[j];
                                    double k1 = 0;
                                    double k2 = 0.;
                                    for (int i = 0; i < 3; ++i) {
                                        k1 += p1[i] * p1[i];
                                        k2 += p2[i] * p2[i];
                                    }
                                    k1 = sqrt(k1 + mass2);
                                    k2 = sqrt(k2 + mass2);
                                    fourVec1 = new Momentum4Vector(p1[0], p1[1], p1[2], k1);
                                    fourVec2 = new Momentum4Vector(p2[0], p2[1], p2[2], k2);
                                    Lorentz4Vector vertex4vecSum = fourVec1.plus(fourVec2);
                                    double vertexMass = vertex4vecSum.mass();
                                    aida.histogram1D("vertex invariant mass phi search two clusters below " + _maxMuClusterEnergy + " " + particleNames[j] + " hypothesis", 200, 0., 2.0).fill(vertexMass);
                                } //end of loop over particle hypotheses
                                skipEvent = false;
                                if (_selectOnMass) {
                                    // reject events outside of mass region of interest
                                    if (mumumass < _minMassCut || mumumass > _maxMassCut) {
                                        skipEvent = true;
                                    }
                                    // analyze "background" region
                                    if (mumumass < _minMassCut) {
                                        analyzeVertex(v, "background", p1, fourVec1, p2, fourVec2, pV, mumusum);
                                    }
                                    // analyze "signal" region
                                    if (mumumass > _minMassCut && mumumass < _maxMassCut) {
                                        analyzeVertex(v, "signal", p1, fourVec1, p2, fourVec2, pV, mumusum);

                                    }
                                }
                            }
                            // define e+e- sample...
                            if (negE > _maxMuClusterEnergy && posE > _maxMuClusterEnergy) {
                                v0type = "e+e-";
                            }
                            aida.histogram1D("v0 mass " + v0type, 50, 0., 0.5).fill(v0.getMass());
                            aida.histogram1D("cluster pair delta time " + v0type, 100, -5., 5.).fill(deltaTime);
                            aida.histogram1D("negative momentum " + v0type, 100, 0., 6.0).fill(negMom);
                            aida.histogram1D("positive momentum " + v0type, 100, 0., 6.0).fill(posMom);
                            aida.histogram2D("negative vs positive momentum " + v0type, 100, 0., 6.0, 100, 0., 6.).fill(negMom, posMom);
                            aida.histogram2D("negative E vs positive E all " + v0type, 100, 0., 5., 100, 0., 5.).fill(negE, posE);
                            aida.histogram1D("negative E over P " + v0type, 100, 0., 2.).fill(negEoverP);
                            aida.histogram1D("positive E over P " + v0type, 100, 0., 2.).fill(posEoverP);
                            aida.histogram2D("negative vs positive E over P " + v0type, 100, 0., 2., 100, 0., 2.).fill(negEoverP, posEoverP);

                            aida.histogram2D("negative cluster ix vs iy " + v0type, 47, -23.5, 23.5, 11, -5.5, 5.5).fill(negIX, negIY);
                            aida.histogram2D("positive cluster ix vs iy " + v0type, 47, -23.5, 23.5, 11, -5.5, 5.5).fill(posIX, posIY);

                        }//end of check on track having greater than minHits
                    }// end of check on both clusters
                }// end of loop over vertices
                aida.tree().cd("..");
            } // end of loop over check on collection
        } // end of loop over vertex collections 
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

    private void analyzeVertex(Vertex v, String type, double[] p1, Momentum4Vector fourVec1, double[] p2, Momentum4Vector fourVec2, double[] pv, Lorentz4Vector mumusum) {
        aida.tree().mkdirs(type);
        aida.tree().cd(type);

        ReconstructedParticle v0 = v.getAssociatedParticle();
        // this always has 2 tracks.
        List<ReconstructedParticle> trks = v0.getParticles();
        ReconstructedParticle neg = trks.get(0);
        ReconstructedParticle pos = trks.get(1);
        Vertex uncVert = v0.getStartVertex();
        Hep3Vector vtxPosRot = VecOp.mult(beamAxisRotation, uncVert.getPosition());
//                  
        // let's also require both to have clusters, since this will distinguish e from mu
        if (!neg.getClusters().isEmpty() && !pos.getClusters().isEmpty()) {
            Cluster negClus = neg.getClusters().get(0);
            int negIX = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
            int negIY = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

            Cluster posClus = pos.getClusters().get(0);
            int posIX = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
            int posIY = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

            // in time
            double negTime = ClusterUtilities.getSeedHitTime(negClus);
            double posTime = ClusterUtilities.getSeedHitTime(posClus);
            double deltaTime = negTime - posTime;
            double topMinusBottomTime = deltaTime;
            if (negClus.getPosition()[1] < 0) {
                topMinusBottomTime = -deltaTime;
            }
            aida.histogram1D("top - bottom cluster delta time", 100, -5., 5.).fill(deltaTime);
            aida.histogram1D("cluster pair delta time", 100, -5., 5.).fill(deltaTime);
            double negE = negClus.getEnergy();
            double posE = posClus.getEnergy();
            int negNhits = neg.getTracks().get(0).getTrackerHits().size();
            int posNhits = pos.getTracks().get(0).getTrackerHits().size();
            double negMom = neg.getMomentum().magnitude();
            double posMom = pos.getMomentum().magnitude();
            double negEoverP = negE / negMom;
            double posEoverP = posE / posMom;

            aida.histogram1D("negative track nHits", 20, 0., 20.).fill(negNhits);
            aida.histogram1D("positive track nHits", 20, 0., 20.).fill(posNhits);
            aida.histogram1D("negative momentum", 100, 0., 6.0).fill(negMom);
            aida.histogram1D("positive momentum", 100, 0., 6.0).fill(posMom);
            aida.histogram1D("v0 x", 50, -5., 5.).fill(vtxPosRot.x());
            aida.histogram1D("v0 y", 50, -2., 2.).fill(vtxPosRot.y());
            aida.histogram2D("v0 x vs v0 y", 50, -5., 5., 50, -1., 1.).fill(vtxPosRot.x(), vtxPosRot.y());
            aida.histogram1D("v0 z", 50, -25., 0.).fill(vtxPosRot.z());
            aida.histogram1D("v0 x neg " + negNhits + " pos " + posNhits + " hits on track", 50, -5., 5.).fill(vtxPosRot.x());
            aida.histogram1D("v0 y neg " + negNhits + " pos " + posNhits + " hits on track", 50, -2., 2.).fill(vtxPosRot.y());
            aida.histogram1D("v0 z neg " + negNhits + " pos " + posNhits + " hits on track", 50, -25., 0.).fill(vtxPosRot.z());
            aida.histogram1D("v0 energy", 100, 0., 10.).fill(v0.getEnergy());
            aida.histogram1D("v0 mass", 50, 0., 0.5).fill(v0.getMass());
            aida.histogram2D("v0 mass vs Z vertex", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
            aida.histogram2D("v0 mass vs Z vertex neg " + negNhits + " pos " + posNhits + " hits on track", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
            aida.profile1D("v0 mass vs Z vertex profile", 50, 0.05, 0.25).fill(v0.getMass(), vtxPosRot.z());
            aida.histogram1D("psum", 100, 0., 6.0).fill(negMom + posMom);
            aida.histogram1D("psum neg " + negNhits + " pos " + posNhits + " hits on track", 100, 0., 6.0).fill(negMom + posMom);
            aida.histogram2D("negative vs positive momentum", 100, 0., 6.0, 100, 0., 6.).fill(negMom, posMom);
            aida.histogram1D("psum both ECal Clusters", 100, 0., 6.0).fill(negMom + posMom);
            aida.histogram1D("esum both ECal Clusters", 100, 0., 6.0).fill(neg.getClusters().get(0).getEnergy() + pos.getClusters().get(0).getEnergy());

            aida.histogram2D("negative E vs positive E all", 100, 0., 5., 100, 0., 5.).fill(negE, posE);
            aida.histogram2D("negative E vs positive E lowE", 100, 0., 0.5, 100, 0., 0.5).fill(negE, posE);
            aida.histogram1D("negative E over P", 100, 0., 2.).fill(negEoverP);
            aida.histogram1D("positive E over P", 100, 0., 2.).fill(posEoverP);
            aida.histogram2D("negative vs positive E over P", 100, 0., 2., 100, 0., 2.).fill(negEoverP, posEoverP);

            aida.histogram2D("negative cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(negIX, negIY);
            aida.histogram2D("positive cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(posIX, posIY);
        } // end of check on clusters

        aida.histogram2D("p1x vs p1y", 100, -0.5, 1.0, 100, -0.5, 0.5).fill(p1[0], p1[1]);
        aida.histogram2D("p2x vs p2y", 100, -0.5, 1.0, 100, -0.5, 0.5).fill(p2[0], p2[1]);
        aida.histogram2D("pvx vs pvy", 100, -0.5, 1.0, 100, -0.5, 0.5).fill(pv[0], pv[1]);
        aida.histogram1D("vertex pz", 100, 0., 7.).fill(pv[2]);

        aida.tree().cd("..");
    }

    public void setMaxMuClusterEnergy(double d) {
        _maxMuClusterEnergy = d;
    }

    public void setClusterDeltaTimeCut(double d) {
        _clusterDeltaTimeCut = d;
    }

    public void setSelectOnMass(boolean b) {
        _selectOnMass = b;
    }

    public void setMinMassCut(double d) {
        _minMassCut = d;
    }

    public void set_maxMassCut(double d) {
        _maxMassCut = d;
    }
}
