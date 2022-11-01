package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.BasicHepLorentzVector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.HepLorentzVector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimEcalFeeWabTrident2016 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    double _esumCut = 1.75;
    private int _numberOfEventsSelected;
    private int _numberOfFeesSelected;
    private int _numberOfWabsSelected;
    private int _numberOfTridentsSelected;
    private int _numberOfEventsProcessed = 0;

    private int _maxNClusters = 3;
    private double _minClusterEnergy = 1.75;

    private double _minSeedHitEnergy = 0.;
    private boolean _requireFiducialFee = false;
    private boolean _requireFiducialWab = false;

    private double _minTridentClusterEnergy = 0.4;
    private double _maxTridentClusterEnergy = 1.5;

    private boolean _skimFee = true;
    private boolean _skimWab = true;
    private boolean _skimTrident = true;

    double emass = 0.000511;
    double emass2 = emass * emass;

    //20220609 empirically determine beam rotation to be 28.5mr from 2021 Mollers
    private double _beamAxisRotationY = -0.025;//-0.0305;
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, _beamAxisRotationY, -Math.PI / 2);
    }

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        boolean isWabCandidate = false;
        if (_skimWab) {
            isWabCandidate = isWabCandidate(ecalClusters);
        }
        boolean isFeeCandidate = false;
        if (_skimFee) {
            isFeeCandidate = isFeeCandidate(ecalClusters);
        }

        if (isWabCandidate || isFeeCandidate) {
            skipEvent = false;
        }

        boolean isTridentCandidate = false;
        if (_skimTrident) {
            isTridentCandidate = isTridentCandidate(ecalClusters, event);
        }

        if (isFeeCandidate) {
            _numberOfFeesSelected++;
            skipEvent = false;
        }
        if (isWabCandidate) {
            _numberOfWabsSelected++;
            skipEvent = false;
        }
        if (isTridentCandidate) {
            _numberOfTridentsSelected++;
            skipEvent = false;
        }

        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    private boolean isFeeCandidate(List<Cluster> ecalClusters) {
        boolean isFeeCandidate = false;
        aida.tree().mkdirs("fee candidate analysis");
        aida.tree().cd("fee candidate analysis");
        aida.histogram1D("number of clusters", 10, -0.5, 9.5).fill(ecalClusters.size());
        int nClusters = ecalClusters.size();
        if (ecalClusters.size() > 0 && nClusters <= _maxNClusters) {
            for (Cluster cluster : ecalClusters) {
                //System.out.println("cluster energy "+cluster.getEnergy());
                if (cluster.getEnergy() > _minClusterEnergy) {
                    double seedEnergy = ClusterUtilities.findSeedHit(cluster).getCorrectedEnergy();
                    if (seedEnergy > _minSeedHitEnergy) {
                        isFeeCandidate = true;
                        boolean clusterIsFiducial = TriggerModule.inFiducialRegion(cluster);
                        String fid = clusterIsFiducial ? " fiducial " : "";
                        if (_requireFiducialFee && !clusterIsFiducial) {
                            isFeeCandidate = false;
                        }
                        if (isFeeCandidate) {
                            CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
                            int ix = seed.getIdentifierFieldValue("ix");
                            int iy = seed.getIdentifierFieldValue("iy");
                            aida.histogram2D("cluster ix vs iy" + fid, 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                            aida.histogram2D("Cluster x vs y" + fid, 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                            aida.histogram1D("Cluster energy" + fid, 100, 0.0, 3.0).fill(cluster.getEnergy());
                            aida.histogram2D("Cluster energy vs seed hit energy" + fid, 100, 0., 3., 100, 0., 3.).fill(cluster.getEnergy(), seedEnergy);

                            if (cluster.getPosition()[1] > 0.) {
                                aida.histogram1D("Top cluster energy" + fid, 100, 0.0, 3.0).fill(cluster.getEnergy());
                                aida.histogram1D("Top cluster energy " + nClusters + " clusters" + fid, 100, 0.0, 3.0).fill(cluster.getEnergy());

                            } else {
                                aida.histogram1D("Bottom cluster energy" + fid, 100, 0.0, 3.0).fill(cluster.getEnergy());
                                aida.histogram1D("Bottom cluster energy " + nClusters + " clusters" + fid, 100, 0.0, 3.0).fill(cluster.getEnergy());
                            }
                        }
                    }
                }
            }
        }
        aida.tree().cd("..");
        return isFeeCandidate;
    }

    private boolean isWabCandidate(List<Cluster> ecalClusters) {
        boolean isWabCandidate = false;
        // let's start by requiring two and only two clusters, in opposite hemispheres,
        // whose energies sum to the beam energy
        aida.tree().mkdirs("wab candidate analysis");
        aida.tree().cd("wab candidate analysis");
        aida.histogram1D("number of clusters", 10, -0.5, 9.5).fill(ecalClusters.size());
        if (ecalClusters.size() == 2) {
            Cluster c1 = ecalClusters.get(0);
            double e1 = c1.getEnergy();
            // eliminate FEEs...
            if (e1 < 3.5) {
                Hep3Vector pos1 = new BasicHep3Vector(c1.getPosition());
                double t1 = ClusterUtilities.getSeedHitTime(c1);

                Cluster c2 = ecalClusters.get(1);
                double e2 = c2.getEnergy();
                double t2 = ClusterUtilities.getSeedHitTime(c2);
                //
                double deltaT = t1 - t2;
                if (abs(deltaT) < 2.) {
                    double esum = e1 + e2;
                    Hep3Vector pos2 = new BasicHep3Vector(c2.getPosition());
//            aida.histogram2D("two cluster e1 vs e2", 100, 0., 5., 100, 0., 5.).fill(e1, e2);
//            aida.histogram1D("two cluster e1 + e2", 100, 0., 5.).fill(esum);
                    // opposite hemispheres
                    if (pos1.x() * pos2.x() < 0. && pos1.y() * pos2.y() < 0.) {
//                aida.histogram2D("two opposite cluster e1 vs e2", 100, 0., 5., 100, 0., 5.).fill(e1, e2);
//                aida.histogram1D("two opposite cluster e1 + e2", 100, 0., 5.).fill(esum);
                        if (esum > _esumCut) {
//                    aida.histogram2D("two opposite esum > " + _esumCut + " cluster e1 vs e2", 100, 0., 5., 100, 0.,
//                            5.).fill(e1, e2);
//                    aida.histogram1D("two opposite esum > " + _esumCut + " cluster e1 + e2", 100, _esumCut, 5.)
//                            .fill(esum);
//                    aida.histogram2D("two opposite esum > " + _esumCut + " cluster1 x vs y", 320, -270.0, 370.0, 90,
//                            -90.0, 90.0).fill(pos1.x(), pos1.y());
//                    aida.histogram2D("two opposite esum > " + _esumCut + " cluster2 x vs y", 320, -270.0, 370.0, 90,
//                            -90.0, 90.0).fill(pos2.x(), pos2.y());
                            boolean e1IsFiducial = TriggerModule.inFiducialRegion(c1);
                            boolean e2IsFiducial = TriggerModule.inFiducialRegion(c2);
                            isWabCandidate = true;
                            boolean wabIsFiducial = e1IsFiducial && e2IsFiducial;
                            if (_requireFiducialWab && !wabIsFiducial) {
                                isWabCandidate = false;
                            }
                            if (isWabCandidate) {
                                if (wabIsFiducial) {
                                    aida.tree().mkdirs("fiducial");
                                    aida.tree().cd("fiducial");
                                }
                                aida.histogram1D("opposite esum > " + _esumCut + " cluster e1", 100, 0.0, 3.0)
                                        .fill(e1);
                                aida.histogram1D("opposite esum > " + _esumCut + " cluster e2", 100, 0.0, 3.0)
                                        .fill(e2);
                                aida.histogram2D("opposite esum > " + _esumCut + " cluster e1 vs e2", 100, 0.,
                                        3., 100, 0.0, 3.0).fill(e1, e2);
                                aida.histogram1D("opposite esum > " + _esumCut + " cluster e1 + e2", 100,
                                        _esumCut, 3.).fill(esum);
                                aida.histogram2D("opposite esum > " + _esumCut + " cluster1 x vs y", 320,
                                        -270.0, 370.0, 90, -90.0, 90.0).fill(pos1.x(), pos1.y());
                                aida.histogram2D("opposite esum > " + _esumCut + " cluster2 x vs y", 320,
                                        -270.0, 370.0, 90, -90.0, 90.0).fill(pos2.x(), pos2.y());
                                aida.histogram1D("cluster delta time", 100, -5., 5.).fill(deltaT);
                                CalorimeterHit seed = c1.getCalorimeterHits().get(0);
                                int ix = seed.getIdentifierFieldValue("ix");
                                int iy = seed.getIdentifierFieldValue("iy");
                                aida.histogram2D("cluster1 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                                seed = c2.getCalorimeterHits().get(0);
                                ix = seed.getIdentifierFieldValue("ix");
                                iy = seed.getIdentifierFieldValue("iy");
                                aida.histogram2D("cluster2 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                                if (wabIsFiducial) {
                                    aida.tree().cd("..");
                                }
                            }
                        } // end of check on esum
                    } // end of check on opposite hemispheres
                } // end of check on deltaT
            } // end of check on e1 < 3.% to get rid of FEEs
        } // end of check on 2 clusters
        aida.tree().cd("..");
        return isWabCandidate;
    }

    private boolean isTridentCandidate(List<Cluster> ecalClusters, EventHeader event) {
        boolean isTridentCandidate = false;
        // let's start by requiring two and only two clusters, in opposite hemispheres,
        // whose energies sum to the beam energy
        aida.tree().mkdirs("trident candidate analysis");
        aida.tree().cd("trident candidate analysis");
        aida.histogram1D("number of clusters", 10, -0.5, 9.5).fill(ecalClusters.size());
        if (ecalClusters.size() > 2) {
            // positrons defined as clusters with x>100
            List<Cluster> positrons = new ArrayList<>();
            //electrons defined as clusters with x<0
            List<Cluster> electrons = new ArrayList<>();
            for (Cluster cluster : ecalClusters) {
                aida.histogram1D("all cluster energy", 100, 0., 3.).fill(cluster.getEnergy());
                // remove FEEs and low energy dross
                if (cluster.getEnergy() > _minTridentClusterEnergy && cluster.getEnergy() < _maxTridentClusterEnergy) {
                    if (cluster.getPosition()[0] > 100.) {
                        aida.histogram1D("all cluster positron energy", 100, 0., 3.).fill(cluster.getEnergy());
                        positrons.add(cluster);
                    }
                    if (cluster.getPosition()[0] < 0.) {
                        aida.histogram1D("all cluster electron energy", 100, 0., 3.).fill(cluster.getEnergy());
                        electrons.add(cluster);
                    }
                }
            }
            // do we have at least one positron and two electron clusters?

            aida.histogram1D("number of electrons", 10, -0.5, 9.5).fill(electrons.size());
            aida.histogram1D("number of positrons", 10, -0.5, 9.5).fill(positrons.size());
            aida.histogram2D("number of positrons vs electrons", 10, -0.5, 9.5, 10, -0.5, 9.5).fill(positrons.size(), electrons.size());

            if (positrons.size() > 0 && electrons.size() > 1) {
                List<Cluster> tridentClusters = new ArrayList<>();
                Cluster positron = null;
                Cluster electron1 = null;
                Cluster electron2 = null;
                // check that they are in time
                for (Cluster pos : positrons) {
                    double t1 = ClusterUtilities.getSeedHitTime(pos);
                    List<Cluster> inTimeElectronClusters = new ArrayList<>();
                    for (Cluster ele : electrons) {
                        double t2 = ClusterUtilities.getSeedHitTime(ele);
                        if (abs(t1 - t2) < 2.) {
                            inTimeElectronClusters.add(ele);
                        }
                    }
                    // do we have two electron clusters in time with the positron?
                    aida.histogram1D("number of in-time electrons", 10, -0.5, 9.5).fill(inTimeElectronClusters.size());
                    if (inTimeElectronClusters.size() == 2) {
                        tridentClusters.add(pos);

                        tridentClusters.addAll(inTimeElectronClusters);
                        positron = pos;
                        electron1 = inTimeElectronClusters.get(0);
                        electron2 = inTimeElectronClusters.get(1);

                    }
                }
                // do we have three in-time clusters?
                if (tridentClusters.size() == 3) {
                    double esum = 0.;
                    double pysum = 0.;
                    //let's look at the kinematics.
                    for (Cluster c : tridentClusters) {
                        esum += c.getEnergy();
                        pysum += c.getEnergy() * sinTheta(c.getPosition());
                    }
                    aida.histogram1D("trident positron energy", 100, 0.0, 3.0).fill(positron.getEnergy());
                    aida.histogram1D("trident electron energy 1", 100, 0.0, 3.0).fill(electron1.getEnergy());
                    aida.histogram1D("trident electron energy 2", 100, 0.0, 3.0).fill(electron2.getEnergy());
                    aida.histogram1D("trident esum", 100, 0.0, 3.0).fill(esum);
                    aida.histogram2D("trident esum vs electron energy 1", 50, 0.0, 3.0, 50, 0.0, 1.5).fill(esum, electron1.getEnergy());
                    aida.histogram2D("trident esum vs electron energy 2", 50, 0.0, 3.0, 50, 0.0, 1.5).fill(esum, electron2.getEnergy());

                    aida.histogram1D("trident pysum", 100, -1.0, 1.0).fill(pysum);
                    //let's see if we can clean things up...

                    if (TriggerModule.inFiducialRegion(positron) && TriggerModule.inFiducialRegion(electron1) && TriggerModule.inFiducialRegion(electron2)) {
                        aida.histogram1D("trident positron energy all fiducial", 100, 0.0, 3.0).fill(positron.getEnergy());
                        aida.histogram1D("trident electron energy 1 all fiducial", 100, 0.0, 3.0).fill(electron1.getEnergy());
                        aida.histogram1D("trident electron energy 2 all fiducial", 100, 0.0, 3.0).fill(electron2.getEnergy());
                        aida.histogram1D("trident esum all fiducial", 100, 0.0, 3.0).fill(esum);
                        aida.histogram1D("trident pysum all fiducial", 100, -1.0, 1.0).fill(pysum);
                    }
                    if (esum > .5 && esum < 5.5 && abs(pysum) < 0.2) {
                        isTridentCandidate = true;
                        analyzeTridentEvent(event, tridentClusters);
                    }
                }
            }
        } // end of check on 2 clusters
        aida.tree().cd("..");
        return isTridentCandidate;
    }

    // look at tracking for events that pass ECal-only trident selection
    private void analyzeTridentEvent(EventHeader event, List<Cluster> ecalClusters) {
        aida.tree().mkdirs("trident event analysis");
        aida.tree().cd("trident event analysis");
        if (event.hasCollection(ReconstructedParticle.class, "FinalStateParticles_KF")) {
            List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
            int nRPs = rps.size();
            int nClusters = event.get(Cluster.class, "EcalClustersCorr").size();
            ReconstructedParticle positron = null;
            boolean positronIsPositron = false;
            ReconstructedParticle electron1 = null;
            boolean electron1IsElectron = false;
            ReconstructedParticle electron2 = null;
            boolean electron2IsElectron = false;
            aida.histogram2D("number of ReconstructedParticles vs number of Clusters all events", 10, -0.5, 9.5, 10, -0.5, 9.5).fill(nRPs, nClusters);
            // loop over RPs to find which ones are associated with the Ecal clusters
            // note that cluster list is: positron, electron1, electron2
            for (ReconstructedParticle rp : rps) {
                if (!rp.getClusters().isEmpty()) {
                    Cluster c = rp.getClusters().get(0);
                    if (c == ecalClusters.get(0)) {
                        positron = rp;
                        if (rp.getParticleIDUsed().getPDG() == -11) {
                            positronIsPositron = true;
                        }
                    }
                    if (c == ecalClusters.get(1)) {
                        electron1 = rp;
                        if (rp.getParticleIDUsed().getPDG() == 11) {
                            electron1IsElectron = true;
                        }
                    }
                    if (c == ecalClusters.get(2)) {
                        electron2 = rp;
                        if (rp.getParticleIDUsed().getPDG() == 11) {
                            electron2IsElectron = true;
                        }
                    }
                }
            }
            double esum = 0.;
            for (Cluster c : ecalClusters) {
                esum += c.getEnergy();
            }
            aida.histogram1D("cluster esum all events", 100, 0.0, 3.0).fill(esum);
            for (int i = 0; i < 3; ++i) {
                aida.histogram1D("cluster " + (i + 1) + " energy all events", 100, 0., 3.).fill(ecalClusters.get(i).getEnergy());
                aida.histogram2D("cluster esum vs cluster " + (i + 1) + " energy all events", 100, 0., 10., 100, 0., 3.).fill(esum, ecalClusters.get(i).getEnergy());
            }
            // sanity check...
            if (positron != null && electron1 != null && electron2 != null) {
                aida.histogram2D("number of ReconstructedParticles vs number of Clusters all matched to RP", 10, -0.5, 9.5, 10, -0.5, 9.5).fill(nRPs, nClusters);
                aida.histogram1D("cluster esum all matched to RP", 100, 0.0, 3.0).fill(esum);
                aida.histogram1D("trident positron cluster energy", 100, 0.0, 3.0).fill(positron.getEnergy());
                aida.histogram1D("trident electron1 cluster energy", 100, 0.0, 3.0).fill(electron1.getEnergy());
                aida.histogram1D("trident electron2 cluster energy", 100, 0.0, 3.0).fill(electron2.getEnergy());
                if (positronIsPositron) {
                    aida.histogram1D("trident positron RP energy", 100, 0.0, 3.0).fill(positron.getEnergy());
                }
                if (electron1IsElectron) {
                    aida.histogram1D("trident electron1 RP energy", 100, 0.0, 3.0).fill(electron1.getEnergy());
                }
                if (electron2IsElectron) {
                    aida.histogram1D("trident electron2 RP energy", 100, 0.0, 3.0).fill(electron2.getEnergy());
                }
                // all the clusters are associated with a correctly-signed track
                if (positronIsPositron && electron1IsElectron && electron2IsElectron) {
                    aida.histogram2D("number of ReconstructedParticles vs number of Clusters all matched to correct sign track", 10, -0.5, 9.5, 10, -0.5, 9.5).fill(nRPs, nClusters);
                    aida.histogram1D("cluster esum all matched to correct sign track", 100, 0.0, 3.0).fill(esum);
                    // let's make a trident...
                    HepLorentzVector fourVec1 = makeFourVector(positron);
                    HepLorentzVector fourVec2 = makeFourVector(electron1);
                    HepLorentzVector fourVec3 = makeFourVector(electron2);
                    HepLorentzVector trident = VecOp.add(fourVec1, fourVec2);
                    trident = VecOp.add(trident, fourVec3);
                    Hep3Vector tridentVecRot = VecOp.mult(beamAxisRotation, trident.v3());
                    aida.histogram1D("trident energy", 100, 0.0, 3.).fill(trident.t());
                    aida.histogram1D("trident pX", 100, -0.2, 0.2).fill(trident.v3().x());
                    aida.histogram1D("trident pY", 100, -0.2, 0.2).fill(trident.v3().y());
                    aida.histogram1D("trident pZ", 100, 0., 10.).fill(trident.v3().z());
                    aida.histogram1D("trident energy rot " + _beamAxisRotationY, 100, 0.0, 3.).fill(trident.t());
                    aida.histogram1D("trident pX rot " + _beamAxisRotationY, 100, -0.2, 0.2).fill(tridentVecRot.x());
                    aida.histogram1D("trident pY rot " + _beamAxisRotationY, 100, -0.2, 0.2).fill(tridentVecRot.y());
                    aida.histogram1D("trident pZ rot " + _beamAxisRotationY, 100, 0., 10.).fill(tridentVecRot.z());

                    for (int i = 0; i < 3; ++i) {
                        aida.histogram1D("cluster " + (i + 1) + " energy all matched to correct sign track", 100, 0., 3.).fill(ecalClusters.get(i).getEnergy());
                        aida.histogram2D("cluster esum vs cluster " + (i + 1) + " energy all matched to correct sign track", 100, 0., 10., 100, 0., 3.).fill(esum, ecalClusters.get(i).getEnergy());
                        aida.histogram2D("trident energy vs cluster " + (i + 1) + " energy all matched to correct sign track", 100, 0., 10., 100, 0., 3.).fill(trident.t(), ecalClusters.get(i).getEnergy());
                    }
                }
            }
        }
        aida.tree().cd("..");
    }

    private BasicHepLorentzVector makeFourVector(ReconstructedParticle rp) {
        double eSquared = emass2;
        double[] mom = rp.getMomentum().v();
        for (int i = 0; i < 3; ++i) {
            eSquared += mom[i] * mom[i];
        }
        return new BasicHepLorentzVector(sqrt(eSquared), mom);
    }

    private double sinTheta(double[] p) {
        return p[1] / sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
        System.out.println("Selected " + _numberOfFeesSelected + " FEEs");
        System.out.println("Selected " + _numberOfWabsSelected + " WABs");
        System.out.println("Selected " + _numberOfTridentsSelected + " Tridents");
    }

    public void setMinClusterEnergy(double d) {
        _minClusterEnergy = d;
    }

    public void setMinSeedHitEnergy(double d) {
        _minSeedHitEnergy = d;
    }

    public void setMaxNClusters(int i) {
        _maxNClusters = i;
    }

    public void setRequireFiducialFee(boolean b) {
        _requireFiducialFee = b;
    }

    public void setRequireFiducialWab(boolean b) {
        _requireFiducialWab = b;
    }

    public void setSkimFee(boolean b) {
        _skimFee = b;
    }

    public void setSkimWab(boolean b) {
        _skimWab = b;
    }

    public void setBeamAxisRotationY(double d) {
        _beamAxisRotationY = d;
    }
}
