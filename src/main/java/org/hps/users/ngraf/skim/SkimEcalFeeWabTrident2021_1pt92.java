package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimEcalFeeWabTrident2021_1pt92 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    double _esumCut = 1.5;//2.0;
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;

    private int _maxNClusters = 3;
    private double _minClusterEnergy = 1.0;//2.0;
    private double _minSeedHitEnergy = 0.95;//1.9;
    private boolean _requireFiducialFee = true;
    private boolean _requireFiducialWab = true;
    private double _cluster1MinEnergy = 0.0;//0.9;//1.8;
    private double _cluster1MaxEnergy = 1.7;//1.4;//2.8;
    private double _cluster2MinEnergy = 0.0;//0.5;
    private double _minTridentClusterEnergy = 0.0;//0.2;//0.5;
    private double _maxTridentClusterEnergy = 1.25;//2.5;

    private boolean _skimFee = true;
    private boolean _skimWab = true;
    private boolean _skimTrident = true;

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
        boolean isTridentCandidate = false;
        if (_skimTrident) {
            isTridentCandidate = isTridentCandidate(ecalClusters);
        }
        //TODO implement calorimeter-only trident selection
        if (isWabCandidate || isFeeCandidate || isTridentCandidate) {
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
                            aida.histogram1D("Cluster energy" + fid, 100, 0.5, 5.).fill(cluster.getEnergy());
                            aida.histogram2D("Cluster energy vs seed hit energy" + fid, 100, 0.5, 5., 100, 0.5, 3.5).fill(cluster.getEnergy(), seedEnergy);

                            if (cluster.getPosition()[1] > 0.) {
                                aida.histogram1D("Top cluster energy" + fid, 100, 0.5, 5.).fill(cluster.getEnergy());
                                aida.histogram1D("Top cluster energy " + nClusters + " clusters" + fid, 100, 0.5, 5.).fill(cluster.getEnergy());

                            } else {
                                aida.histogram1D("Bottom cluster energy" + fid, 100, 0.5, 5.).fill(cluster.getEnergy());
                                aida.histogram1D("Bottom cluster energy " + nClusters + " clusters" + fid, 100, 0.5, 5.).fill(cluster.getEnergy());
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
//            aida.histogram2D("two cluster e1 vs e2", 100, 0., 5., 100, 0., 5.).fill(e1, e2);
//            aida.histogram1D("two cluster e1 + e2", 100, 0., 5.).fill(esum);
                        // opposite hemispheres
                        if (pos1.x() * pos2.x() < 0. && pos1.y() * pos2.y() < 0.) {
//                aida.histogram2D("two opposite cluster e1 vs e2", 100, 0., 5., 100, 0., 5.).fill(e1, e2);
//                aida.histogram1D("two opposite cluster e1 + e2", 100, 0., 5.).fill(esum);
//let's look at the kinematics.
                            double pysum = c1.getEnergy() * sinTheta(c1.getPosition());
                            pysum += c2.getEnergy() * sinTheta(c2.getPosition());
                            aida.histogram1D("wab pysum", 100, -0.1, 0.1).fill(pysum);
                            double coplanarity = TriggerModule.getValueCoplanarity(new Cluster[]{c1, c2});
                            aida.histogram1D("wab pair coplanarity", 90, 0.0, 180.0).fill(coplanarity);
                            aida.histogram2D("esum vs wab pair coplanarity", 100, 0., 5., 90, 0.0, 180.0).fill(esum, coplanarity);

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
                                    aida.histogram1D("opposite esum > " + _esumCut + " cluster e1", 100, 0., 5.)
                                            .fill(e1);
                                    aida.histogram1D("opposite esum > " + _esumCut + " cluster e2", 100, 0., 5.)
                                            .fill(e2);
                                    aida.histogram2D("opposite esum > " + _esumCut + " cluster e1 vs e2", 100, 0.,
                                            5., 100, 0., 5.).fill(e1, e2);
                                    aida.histogram1D("opposite esum > " + _esumCut + " cluster e1 + e2", 100,
                                            0., 5.).fill(esum);
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
                } // end of check on cluster2 energy
            } // end of check on cluster1 energy
        } // end of check on 2 clusters
        aida.tree().cd("..");
        return isWabCandidate;
    }

    private boolean isTridentCandidate(List<Cluster> ecalClusters) {
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
                // remove FEEs and low energy dross
                if (cluster.getEnergy() > _minTridentClusterEnergy && cluster.getEnergy() < _maxTridentClusterEnergy) {
                    if (cluster.getPosition()[0] > 100.) {
                        positrons.add(cluster);
                    }
                    if (cluster.getPosition()[0] < 0.) {
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
                    }
                }
                // do we have three in-time clusters?
                if (tridentClusters.size() == 3) {
                    double esum = 0.;
                    double pysum = 0.;
                    //let's look at the kinematics.
                    for (Cluster c : tridentClusters) {
                        aida.histogram1D("cluster energy", 100, 0., 5.).fill(c.getEnergy());
                        esum += c.getEnergy();
                        pysum += c.getEnergy() * sinTheta(c.getPosition());
                    }
                    aida.histogram1D("trident esum", 100, 0., 5.).fill(esum);
                    aida.histogram1D("trident pysum", 100, -1.0, 1.0).fill(pysum);
                    if (esum > 0.5 && esum < 5. && abs(pysum) < 0.2) {
                        isTridentCandidate = true;
                    }
                }
            }
        } // end of check on 2 clusters
        aida.tree().cd("..");
        return isTridentCandidate;
    }

    private double sinTheta(double[] p) {
        return p[1] / sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
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

    public void setSkimTrident(boolean b) {
        _skimTrident = b;
    }

    public void setCluster1MinEnergy(double d) {
        _cluster1MinEnergy = d;
    }

    public void setCluster1MaxEnergy(double d) {
        _cluster1MaxEnergy = d;
    }

    public void setCluster2MinEnergy(double d) {
        _cluster2MinEnergy = d;
    }
}
