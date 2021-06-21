package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.Hep3Vector;
import hep.physics.vec.HepLorentzVector;
import hep.physics.vec.VecOp;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class RhoPiPiAnalysisDriver extends Driver {

    AIDA aida = AIDA.defaultInstance();

    @Override
    protected void process(EventHeader event) {
        List<MCParticle> mcparticles = event.get(MCParticle.class).get(0);
        double mass = 0.;
        if (mcparticles.size() > 1) {
            MCParticle piplus = mcparticles.get(0);
            MCParticle piminus = mcparticles.get(1);
            HepLorentzVector rho = VecOp.add(piplus.asFourVector(), piminus.asFourVector());
            mass = sqrt(rho.t() * rho.t() - rho.v3().magnitudeSquared());
            aida.histogram1D("rho mass all events", 150, 0., 1.5).fill(mass);
            aida.histogram1D("pi+ momentum all events", 100, 0., 5.0).fill(piplus.getMomentum().magnitude());
            aida.histogram1D("pi- momentum all events", 100, 0., 5.0).fill(piminus.getMomentum().magnitude());
            if (piplus.getSimulatorStatus().isDecayedInCalorimeter() && piminus.getSimulatorStatus().isDecayedInCalorimeter()) {
                aida.histogram1D("rho mass in calorimeter acceptance", 150, 0., 1.5).fill(mass);
                Hep3Vector momplus = piplus.getMomentum();
                Hep3Vector momminus = piminus.getMomentum();
                if (momplus.x() * momminus.x() < 0 && momplus.y() * momminus.y() < 0) {
                    aida.histogram1D("rho mass in calorimeter acceptance opposite quadrants", 150, 0., 1.5).fill(mass);
                }
            }

            if (event.hasCollection(Cluster.class, "EcalClusters")) {
                List<Cluster> clusters = event.get(Cluster.class, "EcalClusters");
                aida.histogram1D("Number of Clusters in Event", 10, 0., 10.).fill(clusters.size());
                if (clusters.size() > 1) {
                    List<Cluster> mipClusters = new ArrayList<>();
                    for (Cluster cluster : clusters) {
                        double[] cPos = cluster.getPosition();
                        aida.histogram2D("Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cPos[0], cPos[1]);
                        aida.histogram1D("Cluster Energy", 100, 0., 0.5).fill(cluster.getEnergy());
                        aida.histogram1D("Cluster number of crystals", 20, 0., 20.).fill(cluster.getCalorimeterHits().size());
                        aida.histogram2D("Cluster energy vs number of crystals", 100, 0., 3., 20, 0., 20.).fill(cluster.getEnergy(), cluster.getCalorimeterHits().size());
                        if (cluster.getEnergy() > 0.125 && cluster.getEnergy() < 0.200) {
                            mipClusters.add(cluster);
                        }
                    }
                    if (clusters.size() == 2) {
                        Cluster c1 = clusters.get(0);
                        Cluster c2 = clusters.get(1);
                        aida.histogram1D("Two Clusters Cluster Energy", 100, 0., 0.5).fill(c1.getEnergy());
                        aida.histogram1D("Two Clusters Cluster Energy", 100, 0., 0.5).fill(c2.getEnergy());
                        aida.histogram2D("Two Clusters e1 vs e2", 100, 0., 0.5, 100, 0., 0.5).fill(c1.getEnergy(), c2.getEnergy());
                    }
                    aida.histogram1D("Number of MIP Clusters in Event", 10, 0., 10.).fill(mipClusters.size());
                    if (mipClusters.size() == 2) {
                        Cluster c1 = mipClusters.get(0);
                        Cluster c2 = mipClusters.get(1);
                        aida.histogram2D("Two MIP Clusters e1 vs e2", 100, 0., 0.25, 100, 0., 0.25).fill(c1.getEnergy(), c2.getEnergy());
                        aida.histogram1D("Two MIP Clusters rho mass", 150, 0., 1.5).fill(mass);

                        // see if we can identify the MC parentage of the clusters
                        // check whether each cluster has contributions only from pi+ or pi- 
                        //cluster 1
                        boolean c1HasPiPlus = false;
                        boolean c1HasPiMinus = false;
                        List<CalorimeterHit> c1Hits = c1.getCalorimeterHits();
                        Set<MCParticle> c1mcContributors = new HashSet<>();
                        for (CalorimeterHit chit : c1Hits) {
                            if (chit instanceof SimCalorimeterHit) {
                                SimCalorimeterHit simHit = (SimCalorimeterHit) chit;
                                for (int i = 0; i < simHit.getMCParticleCount(); ++i) {
                                    c1mcContributors.add(simHit.getMCParticle(i));
                                }
                            }
                        }
                        aida.histogram1D("number of MCParticles contributing to cluster 1", 10, 0., 10.).fill(c1mcContributors.size());
                        for (MCParticle mcp : c1mcContributors) {
//                            System.out.println("cluster 1 mcpdgid: " + mcp.getPDGID());
                            if (mcp.getPDGID() == 211) {
                                c1HasPiPlus = true;
                            }
                            if (mcp.getPDGID() == -211) {
                                c1HasPiMinus = true;
                            }
                        }
                        //cluster 2
                        boolean c2HasPiPlus = false;
                        boolean c2HasPiMinus = false;
                        List<CalorimeterHit> c2Hits = c2.getCalorimeterHits();
                        Set<MCParticle> c2mcContributors = new HashSet<>();
                        for (CalorimeterHit chit : c2Hits) {
                            if (chit instanceof SimCalorimeterHit) {
                                SimCalorimeterHit simHit = (SimCalorimeterHit) chit;
                                for (int i = 0; i < simHit.getMCParticleCount(); ++i) {
                                    c2mcContributors.add(simHit.getMCParticle(i));
                                }
                            }
                        }
                        aida.histogram1D("number of MCParticles contributing to cluster 2", 10, 0., 10.).fill(c2mcContributors.size());
                        for (MCParticle mcp : c2mcContributors) {
//                            System.out.println("cluster 2mcpdgid: " + mcp.getPDGID());
                            if (mcp.getPDGID() == 211) {
                                c2HasPiPlus = true;
                            }
                            if (mcp.getPDGID() == -211) {
                                c2HasPiMinus = true;
                            }
                        }

                        // do we have any "golden" two cluster events where one cluster has only pi+
                        // and the other cluster only has pi-?
                        boolean c1IsPiPlus = c1HasPiPlus & !c1HasPiMinus;
                        boolean c1IsPiMinus = !c1HasPiPlus & c1HasPiMinus;

                        boolean c2IsPiPlus = c2HasPiPlus & !c2HasPiMinus;
                        boolean c2IsPiMinus = !c2HasPiPlus & c2HasPiMinus;

                        if ((c1IsPiPlus && c2IsPiMinus) || (c2IsPiPlus && c1IsPiMinus)) {

//                            for (MCParticle mcp : c1mcContributors) {
//                                System.out.println("cluster 1 mcpdgid: " + mcp.getPDGID());
//                            }
//                            for (MCParticle mcp : c2mcContributors) {
//                                System.out.println("cluster 2 mcpdgid: " + mcp.getPDGID());
//                            }
                            aida.histogram1D("Cluster Energy two MIP clusters pi+ pi-", 100, 0., 0.5).fill(c1.getEnergy());
                            aida.histogram1D("Cluster Energy two MIP clusters pi+ pi-", 100, 0., 0.5).fill(c2.getEnergy());

                            aida.histogram1D("rho mass two MIP clusters pi+ pi-", 150, 0., 1.5).fill(mass);
                            aida.histogram1D("pi+ momentum two MIP clusters pi+ pi-", 100, 0., 5.0).fill(piplus.getMomentum().magnitude());
                            aida.histogram1D("pi- momentum two MIP clusters pi+ pi-", 100, 0., 5.0).fill(piminus.getMomentum().magnitude());
                            aida.histogram2D("pi+ vs pi- momentum two MIP clusters pi+ pi-", 100, 0., 5.0, 100, 0., 5.0).fill(piplus.getMomentum().magnitude(), piminus.getMomentum().magnitude());

                            double[] pipluspos = null;
                            double[] piminuspos = null;
                            if (c1IsPiPlus) {
                                pipluspos = c1.getPosition();
                                piminuspos = c2.getPosition();
                            }
                            if (c2IsPiPlus) {
                                pipluspos = c2.getPosition();
                                piminuspos = c1.getPosition();
                            }
                            aida.histogram2D("pi+ cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(pipluspos[0], pipluspos[1]);
                            aida.histogram2D("pi- cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(piminuspos[0], piminuspos[1]);
                        }

                    }
                }
            }
        }
    }
}
