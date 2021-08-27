package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
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
public class SkimEcalFeeWabTrident2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    double _esumCut = 3.5;
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;

    private int _maxNClusters = 3;
    private double _minClusterEnergy = 4.0;
    private double _minSeedHitEnergy = 2.6;
    private boolean _requireFiducialFee = false;
    private boolean _requireFiducialWab = true;

    private boolean _skimFee = true;
    private boolean _skimWab = true;

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
        //TODO implement calorimeter-only trident selection
        if (isWabCandidate || isFeeCandidate) {
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
                            aida.histogram1D("Cluster energy" + fid, 100, 3.5, 5.5).fill(cluster.getEnergy());
                            aida.histogram2D("Cluster energy vs seed hit energy" + fid, 100, 2.5, 5., 100, 0.5, 3.5).fill(cluster.getEnergy(), seedEnergy);

                            if (cluster.getPosition()[1] > 0.) {
                                aida.histogram1D("Top cluster energy" + fid, 100, 3.5, 5.5).fill(cluster.getEnergy());
                                aida.histogram1D("Top cluster energy " + nClusters + " clusters" + fid, 100, 3.5, 5.5).fill(cluster.getEnergy());

                            } else {
                                aida.histogram1D("Bottom cluster energy" + fid, 100, 3.5, 5.5).fill(cluster.getEnergy());
                                aida.histogram1D("Bottom cluster energy " + nClusters + " clusters" + fid, 100, 3.5, 5.5).fill(cluster.getEnergy());
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
                            aida.histogram1D("opposite esum > " + _esumCut + " cluster e1", 100, 0., 5.)
                                    .fill(e1);
                            aida.histogram1D("opposite esum > " + _esumCut + " cluster e2", 100, 0., 5.)
                                    .fill(e2);
                            aida.histogram2D("opposite esum > " + _esumCut + " cluster e1 vs e2", 100, 0.,
                                    5., 100, 0., 5.).fill(e1, e2);
                            aida.histogram1D("opposite esum > " + _esumCut + " cluster e1 + e2", 100,
                                    _esumCut, 5.).fill(esum);
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
                    }
                }
            }
        }
        aida.tree().cd("..");
        return isWabCandidate;
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
}
