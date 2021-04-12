package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimTwoFiducialEcalClusters2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsWritten = 0;
    double _esumCut = 3.5;
    private int _numberOfEventsSelected = 0;

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        // let's start by requiring two and only two clusters, in opposite hemispheres,
        // whose energies sum to the beam energy
        aida.tree().mkdirs("two cluster analysis");
        aida.tree().cd("two cluster analysis");
        aida.histogram1D("number of clusters", 10, 0., 10.).fill(ecalClusters.size());
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
                        boolean e1IsFiducial = isFiducial(ClusterUtilities.findSeedHit(c1));
                        boolean e2IsFiducial = isFiducial(ClusterUtilities.findSeedHit(c2));
                        if (e1IsFiducial && e2IsFiducial) {
                            aida.histogram1D("two fiducial opposite esum > " + _esumCut + " cluster e1", 100, 0., 5.)
                                    .fill(e1);
                            aida.histogram1D("two fiducial opposite esum > " + _esumCut + " cluster e2", 100, 0., 5.)
                                    .fill(e2);
                            aida.histogram2D("two fiducial opposite esum > " + _esumCut + " cluster e1 vs e2", 100, 0.,
                                    5., 100, 0., 5.).fill(e1, e2);
                            aida.histogram1D("two fiducial opposite esum > " + _esumCut + " cluster e1 + e2", 100,
                                    _esumCut, 5.).fill(esum);
                            aida.histogram2D("two fiducial opposite esum > " + _esumCut + " cluster1 x vs y", 320,
                                    -270.0, 370.0, 90, -90.0, 90.0).fill(pos1.x(), pos1.y());
                            aida.histogram2D("two fiducial opposite esum > " + _esumCut + " cluster2 x vs y", 320,
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
                        }
                        skipEvent = false;
                        if (!e1IsFiducial) {
                            skipEvent = true;
                        }
                        if (!e2IsFiducial) {
                            skipEvent = true;
                        }
                    }
                }
            }
        }
        aida.tree().cd("..");
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " events");
    }

    public boolean isFiducial(CalorimeterHit hit) {
        int ix = hit.getIdentifierFieldValue("ix");
        int iy = hit.getIdentifierFieldValue("iy");
        // Get the x and y indices for the cluster.
        int absx = Math.abs(ix);
        int absy = Math.abs(iy);

        // Check if the cluster is on the top or the bottom of the
        // calorimeter, as defined by |y| == 5. This is an edge cluster
        // and is not in the fiducial region.
        if (absy == 5) {
            return false;
        }

        // Check if the cluster is on the extreme left or right side
        // of the calorimeter, as defined by |x| == 23. This is also
        // an edge cluster and is not in the fiducial region.
        if (absx == 23) {
            return false;
        }

        // Check if the cluster is along the beam gap, as defined by
        // |y| == 1. This is an internal edge cluster and is not in the
        // fiducial region.
        if (absy == 1) {
            return false;
        }

        // Lastly, check if the cluster falls along the beam hole, as
        // defined by clusters with -11 <= x <= -1 and |y| == 2. This
        // is not the fiducial region.
        if (absy == 2 && ix <= -1 && ix >= -11) {
            return false;
        }

        // If all checks fail, the cluster is in the fiducial region.
        return true;
    }

    public void setEsumCut(double d) {
        _esumCut = d;
    }
}
