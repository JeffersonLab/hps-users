package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.BasicHepLorentzVector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.HepLorentzVector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Skim off events containing Vertex object(s) that are consistent with
 * fully-reconstructed tridents (e+e-e-) or muon final states (mu+ mu-) This
 * Driver uses simple, robust cuts, to be refined later in the analysis stage.
 *
 * For tridents I require the positron to have an associated cluster, since it
 * must have triggered.
 *
 * For muons, I require both legs to have an associated cluster and that the
 * cluster energy be consistent with a MIP deposit.
 *
 * Loop over both GBL and KF collections
 *
 * @author Norman A. Graf
 */
public class SkimV0TridentAndMuMu2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _numberOfEventsSelected = 0;
    private int _numberOfTridentsSelected = 0;
    private int _numberOfV0mumuSelected = 0;
    private int _numberOfEventsProcessed = 0;
    boolean _skipEvent = true;
    String[] _vertexCollectionNames = {"UnconstrainedV0Vertices", "UnconstrainedV0Vertices_KF"};

    // mumu analysis
    private double _clusterDeltaTimeCut = 5.0;
    private double _muonClusterEnergyCut = 0.45;

    // trident analysis
    double _emass = 0.000511;
    double _emass2 = _emass * _emass;
    private final BasicHep3Matrix _beamAxisRotation = new BasicHep3Matrix();
    private int _minNumberOfElectronsWithClusters = 2;

    protected void detectorChanged(Detector detector) {
        _beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    public void process(EventHeader event) {
        _skipEvent = true;
        _numberOfEventsProcessed++;
        for (String vertexCollectionName : _vertexCollectionNames) {
            aida.tree().mkdirs(vertexCollectionName + " mu+ mu- analysis");
            aida.tree().cd(vertexCollectionName + " mu+ mu- analysis");
            if (event.hasCollection(Vertex.class, vertexCollectionName)) {
                List<Vertex> vertices = event.get(Vertex.class, vertexCollectionName);
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                if (vertices.size() > 0) {
                    if (skimV0mumu(vertices)) {
                        _numberOfV0mumuSelected++;
                    }
                }
            } // end of loop over check on collection
            aida.tree().cd("..");
        } // end of loop over vertex collections 
        if (skimTrident(event)) {
            _numberOfTridentsSelected++;
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
        System.out.println("Found " + _numberOfTridentsSelected + " trident candidates");
        System.out.println("Found " + _numberOfV0mumuSelected + " V0 mu+mu- candidates");
    }

    private boolean skimV0mumu(List<Vertex> vertices) {
        boolean foundOne = false;
        for (Vertex v : vertices) {
            ReconstructedParticle v0 = v.getAssociatedParticle();
            int minNhits = 5; // GBL tracks have 2D hits
            if (v0.getType() == 1) {
                minNhits = 10; // Kalman tracks have 1D hits
            }
            // These vertices always have 2 tracks.
            List<ReconstructedParticle> trks = v0.getParticles();
            ReconstructedParticle neg = trks.get(0);
            ReconstructedParticle pos = trks.get(1);
            // let's require both particles to have clusters, since cluster energy will distinguish e from mu
            if (!neg.getClusters().isEmpty() && !pos.getClusters().isEmpty()) {
                Cluster negClus = neg.getClusters().get(0);
                Cluster posClus = pos.getClusters().get(0);
                // in time
                double negTime = ClusterUtilities.getSeedHitTime(negClus);
                double posTime = ClusterUtilities.getSeedHitTime(posClus);
                double topMinusBottomTime = negTime - posTime;
                if (negClus.getPosition()[1] < 0) {
                    topMinusBottomTime = -topMinusBottomTime;
                }
                double negE = negClus.getEnergy();
                double posE = posClus.getEnergy();
                aida.histogram2D("all negative cluster energy vs positive cluster energy", 100, 0., 5., 100, 0., 5.).fill(negE, posE);
                aida.histogram1D("all top - bottom cluster delta time", 100, -5., 5.).fill(topMinusBottomTime);

                // are Ecal clusters consistent with MIP deposition?
                if (negE <= _muonClusterEnergyCut && posE <= _muonClusterEnergyCut) {
                    // are ECal clusters in time?
                    if (abs(topMinusBottomTime) <= _clusterDeltaTimeCut) {
                        int negNhits = neg.getTracks().get(0).getTrackerHits().size();
                        int posNhits = pos.getTracks().get(0).getTrackerHits().size();
                        if (negNhits >= minNhits && posNhits >= minNhits && abs(topMinusBottomTime) < _clusterDeltaTimeCut) {
                            _skipEvent = false;
                            foundOne = true;
                            aida.histogram2D("negative cluster energy vs positive cluster energy", 100, 0., 1., 100, 0., 1.).fill(negE, posE);
                            aida.histogram1D("top - bottom cluster delta time", 100, -5., 5.).fill(topMinusBottomTime);
                        }//end of check on track having greater than minHits
                    }// end of check on cluster delta time
                }// end of check on cluster energies
            }// end of check on both clusters
        }// end of loop over vertices
        return foundOne;
    }

    private boolean skimTrident(EventHeader event) {
        boolean foundOne = false;
        // note that this algorithm is quick, but only picks up events with the recoil in the same half of the detector as the positron.
        String[] collectionType = {"", "_KF"};

        for (String s : collectionType) {
            aida.tree().mkdirs("trident analysis" + s);
            aida.tree().cd("trident analysis" + s);
            List<ReconstructedParticle> V0List = event.get(ReconstructedParticle.class, "UnconstrainedV0Candidates" + s);
            int nV0 = V0List.size();
            List<ReconstructedParticle> VCList = event.get(ReconstructedParticle.class, "UnconstrainedVcCandidates" + s);
            int nVc = VCList.size();
            aida.histogram2D("number of V0 vs Vc ", 5, 0., 5., 5, 0., 5.).fill(nV0, nVc);
            if (nV0 == 1 && nVc == 1) {
                ReconstructedParticle v0 = V0List.get(0);
                // electron comes first
                List<ReconstructedParticle> v0parts = v0.getParticles();
                ReconstructedParticle vc = VCList.get(0);
                List<ReconstructedParticle> vcparts = vc.getParticles();
//                System.out.println("V0 particles" + s);
//                for (ReconstructedParticle p : v0parts) {
//                    System.out.println(p);
//                }
//                System.out.println("Vc particles" + s);
//                for (ReconstructedParticle p : vcparts) {
//                    System.out.println(p);
//                }
                List<ReconstructedParticle> electrons = new ArrayList<>();
                List<ReconstructedParticle> positrons = new ArrayList<>();

                positrons.add(v0parts.get(1));
                electrons.add(v0parts.get(0));
                electrons.add(vcparts.get(0));

                // positron must have associated Ecal cluster has to have triggered...
                if (!positrons.get(0).getClusters().isEmpty()) {
                    Cluster pc1 = positrons.get(0).getClusters().get(0);
                    double posClusTime = ClusterUtilities.findSeedHit(pc1).getTime();
                    String positorb = pc1.getPosition()[1] > 0 ? "top " : "bot ";
                    BasicHepLorentzVector p1 = makeFourVector(electrons.get(0));
                    BasicHepLorentzVector p2 = makeFourVector(electrons.get(1));
                    BasicHepLorentzVector p3 = makeFourVector(positrons.get(0));

                    HepLorentzVector trident = VecOp.add(p1, p2);
                    trident = VecOp.add(trident, p3);
                    Hep3Vector tridentVec = VecOp.mult(_beamAxisRotation, trident.v3());

                    aida.histogram1D("trident energy", 100, 0.0, 10.).fill(trident.t());

                    double dc1 = 0.;  //positron electron1
                    double dc2 = 0.;  //positron electron2
                    double dc3 = 0.;  //electron1 electron2
                    double ele1ClusTime = 0.;
                    double ele2ClusTime = 0.;
                    List<Cluster> electronClusters = new ArrayList<>();
                    if (!electrons.get(0).getClusters().isEmpty()) {
                        Cluster c = electrons.get(0).getClusters().get(0);
                        ele1ClusTime = ClusterUtilities.findSeedHit(c).getTime();
                        electronClusters.add(c);
                        dc1 = posClusTime - ele1ClusTime;
                        aida.histogram1D("delta cluster time ele1 pos", 50, -5., 5.).fill(ele1ClusTime - posClusTime);
                    }
                    if (!electrons.get(1).getClusters().isEmpty()) {
                        Cluster c = electrons.get(1).getClusters().get(0);
                        ele2ClusTime = ClusterUtilities.findSeedHit(c).getTime();
                        electronClusters.add(c);
                        dc2 = posClusTime - ele2ClusTime;
                        aida.histogram1D("delta cluster time ele2 pos", 50, -5., 5.).fill(ele2ClusTime - posClusTime);
                    }
                    // require at least one of the electrons to have a cluster?
                    int nClusters = electronClusters.size();
                    if (nClusters == 2) {
                        dc3 = ele1ClusTime - ele2ClusTime;
                        aida.histogram1D("delta cluster time ele1 ele2", 50, -5., 5.).fill(ele1ClusTime - ele2ClusTime);
                        double ele1e = electrons.get(0).getClusters().get(0).getEnergy();
                        double ele2e = electrons.get(1).getClusters().get(0).getEnergy();
                        double pose = positrons.get(0).getClusters().get(0).getEnergy();
                        aida.histogram1D("ele1 cluster energy", 100, 0., 5.).fill(ele1e);
                        aida.histogram1D("ele2 cluster energy", 100, 0., 5.).fill(ele2e);
                        aida.histogram1D("pos cluster energy", 100, 0., 5.).fill(pose);
                        aida.histogram2D("ele1 vs ele2 cluster energy", 100, 0., 3., 100, 0., 3.).fill(ele1e, ele2e);
                        aida.histogram2D("ele1 vs pose cluster energy", 100, 0., 3., 100, 0., 3.).fill(ele1e, pose);
                        aida.histogram2D("ele2 vs pose cluster energy", 100, 0., 3., 100, 0., 3.).fill(ele2e, pose);

                        aida.histogram2D("ele1 cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(electrons.get(0).getClusters().get(0).getPosition()[0], electrons.get(0).getClusters().get(0).getPosition()[1]);
                        aida.histogram2D("ele2 cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(electrons.get(1).getClusters().get(0).getPosition()[0], electrons.get(1).getClusters().get(0).getPosition()[1]);
                        aida.histogram2D("pos cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(positrons.get(0).getClusters().get(0).getPosition()[0], positrons.get(0).getClusters().get(0).getPosition()[1]);
                    }
                    // require clusters to be in time...
                    if (abs(dc1) < _clusterDeltaTimeCut && abs(dc2) < _clusterDeltaTimeCut && abs(dc3) < _clusterDeltaTimeCut) {
                        aida.histogram1D("trident energy after cluster timing cuts " + nClusters + " electron clusters", 100, 0.0, 10.).fill(trident.t());
                        aida.histogram1D("Number of V0s in the event", 10, -0.5, 9.5).fill(nV0);
                        aida.histogram1D("Number of VCs in the event", 10, -0.5, 9.5).fill(nVc);
                        aida.histogram1D("trident energy after cluster timing cuts " + nClusters + " electron clusters " + nV0 + " V0 " + nVc + " Vc", 100, 0.0, 10.).fill(trident.t());
                        // keep this event
                        if (nClusters >= _minNumberOfElectronsWithClusters) {
                            _skipEvent = false;
                            foundOne = true;
                            aida.histogram1D("good trident run number", 750, 10000, 10750).fill(event.getRunNumber());
                        }
                    }
                }
            }
            aida.tree().cd("..");
        }
        return foundOne;
    }

    private BasicHepLorentzVector makeFourVector(ReconstructedParticle rp) {
        double emass2 = _emass * _emass;
        double eSquared = emass2;
        double[] mom = rp.getMomentum().v();
        for (int i = 0; i < 3; ++i) {
            eSquared += mom[i] * mom[i];
        }
        return new BasicHepLorentzVector(sqrt(eSquared), mom);
    }

    public void setClusterDeltaTimeCut(double d) {
        _clusterDeltaTimeCut = d;
    }

    public void setMinNumberOfElectronsWithClusters(int i) {
        _minNumberOfElectronsWithClusters = i;
    }

    public void setMuonClusterEnergyCut(double d) {
        _muonClusterEnergyCut = d;
    }
}
