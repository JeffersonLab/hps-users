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
import org.hps.recon.tracking.TrackData;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimTridents2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _numberOfEventsSelected = 0;
    boolean skipEvent = true;
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private int minNumberOfElectronsWithClusters = 2;
    String[] reconstructedParticleCollectionNames = {"FinalStateParticles"};//, "FinalStateParticles_KF"};
    double trackDeltaTimeCut = 4.;
    double clusterDeltaTimeCut = 2.;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    double emass = 0.000511;
    double emass2 = emass * emass;

    @Override
    protected void process(EventHeader event) {
        skipEvent = true;
        List<ReconstructedParticle> V0List = event.get(ReconstructedParticle.class, "UnconstrainedV0Candidates");
        int nV0 = V0List.size();
        List<ReconstructedParticle> VCList = event.get(ReconstructedParticle.class, "UnconstrainedVcCandidates");
        int nVc = VCList.size();

        if (nV0 == 1 && nVc == 1) {
            for (String reconstructedParticleCollectionName : reconstructedParticleCollectionNames) {
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, reconstructedParticleCollectionName);
                if (rpList.size() == 3) {
                    aida.tree().mkdirs(reconstructedParticleCollectionName + " trident track analysis");
                    aida.tree().cd(reconstructedParticleCollectionName + " trident track analysis");

                    List<ReconstructedParticle> electrons = new ArrayList<>();
                    List<ReconstructedParticle> positrons = new ArrayList<>();
                    for (ReconstructedParticle rp : rpList) {
                        if (!rp.getParticleIDs().isEmpty()) {
                            int pdgId = rp.getParticleIDUsed().getPDG();
                            if (pdgId == 11) {
                                electrons.add(rp);
                            }
                            if (pdgId == -11) {
                                positrons.add(rp);
                            }
                        }
                    }
                    // for now, let's be selective and require three and only three tracks in the event
                    // cleans things up and doesn't require complicated looping over combinatorics...
                    // should I further require no photons in the event?
                    if (electrons.size() == 2 && positrons.size() == 1) {
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
                            Hep3Vector tridentVec = VecOp.mult(beamAxisRotation, trident.v3());

                            aida.histogram1D("trident energy", 100, 0.0, 10.).fill(trident.t());

                            // let's look at track times...
                            Track posTrack = positrons.get(0).getTracks().get(0);
                            double posTrackTime = TrackData.getTrackTime(TrackData.getTrackData(event, posTrack));

                            Track ele1Track = electrons.get(0).getTracks().get(0);
                            double ele1TrackTime = TrackData.getTrackTime(TrackData.getTrackData(event, ele1Track));

                            Track ele2Track = electrons.get(1).getTracks().get(0);
                            double ele2TrackTime = TrackData.getTrackTime(TrackData.getTrackData(event, ele2Track));
                            // before cuts
                            aida.histogram1D("delta track time ele1 ele2 before time cuts", 50, -5., 5.).fill(ele1TrackTime - ele2TrackTime);
                            aida.histogram1D("delta track time ele1 pos1 before time cuts", 50, -5., 5.).fill(ele1TrackTime - posTrackTime);
                            aida.histogram1D("delta track time ele2 pos1 before time cuts", 50, -5., 5.).fill(ele2TrackTime - posTrackTime);

                            double dt1 = ele1TrackTime - ele2TrackTime;
                            double dt2 = ele1TrackTime - posTrackTime;
                            double dt3 = ele2TrackTime - posTrackTime;
                            // require all three tracks to be in time
                            if (abs(dt1) < trackDeltaTimeCut && abs(dt2) < trackDeltaTimeCut && abs(dt3) < trackDeltaTimeCut) {
                                aida.histogram1D("trident energy after track timing cuts", 100, 0.0, 10.).fill(trident.t());
                                aida.histogram1D("delta track time ele1 ele2 after time cuts", 50, -5., 5.).fill(ele1TrackTime - ele2TrackTime);
                                aida.histogram1D("delta track time ele1 pos after time cuts", 50, -5., 5.).fill(ele1TrackTime - posTrackTime);
                                aida.histogram1D("delta track time ele2 pos after time cuts", 50, -5., 5.).fill(ele2TrackTime - posTrackTime);
                                aida.histogram1D("positron track time - cluster time", 100, -50., -30.).fill(posTrackTime - posClusTime);
                                aida.histogram1D(positorb + "positron track time - cluster time", 100, -50., -30.).fill(posTrackTime - posClusTime);
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
                                    aida.histogram1D("electron1 track time - cluster time", 100, -50., -30.).fill(ele1TrackTime - ele1ClusTime);
                                    aida.histogram1D("delta cluster time ele1 pos", 50, -5., 5.).fill(ele1ClusTime - posClusTime);
                                }
                                if (!electrons.get(1).getClusters().isEmpty()) {
                                    Cluster c = electrons.get(1).getClusters().get(0);
                                    ele2ClusTime = ClusterUtilities.findSeedHit(c).getTime();
                                    electronClusters.add(c);
                                    dc2 = posClusTime - ele2ClusTime;
                                    aida.histogram1D("electron2 track time - cluster time", 100, -50., -30.).fill(ele2TrackTime - ele2ClusTime);
                                    aida.histogram1D("delta cluster time ele2 pos", 50, -5., 5.).fill(ele2ClusTime - posClusTime);
                                }
                                // require at least one of the electrons to have a cluster?
                                int nClusters = electronClusters.size();
                                if (nClusters == 2) {
                                    dc3 = ele1ClusTime - ele2ClusTime;
                                    aida.histogram1D("delta cluster time ele1 ele2", 50, -5., 5.).fill(ele1ClusTime - ele2ClusTime);
                                    aida.histogram1D("ele1 cluster energy", 100, 0., 5.).fill(electrons.get(0).getClusters().get(0).getEnergy());
                                    aida.histogram1D("ele2 cluster energy", 100, 0., 5.).fill(electrons.get(1).getClusters().get(0).getEnergy());
                                    aida.histogram1D("pos cluster energy", 100, 0., 5.).fill(positrons.get(0).getClusters().get(0).getEnergy());
                                    aida.histogram2D("ele1 cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(electrons.get(0).getClusters().get(0).getPosition()[0], electrons.get(0).getClusters().get(0).getPosition()[1]);
                                    aida.histogram2D("ele2 cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(electrons.get(1).getClusters().get(0).getPosition()[0], electrons.get(1).getClusters().get(0).getPosition()[1]);
                                    aida.histogram2D("pos cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(positrons.get(0).getClusters().get(0).getPosition()[0], positrons.get(0).getClusters().get(0).getPosition()[1]);
                                }
                                aida.histogram1D("trident energy after track timing cuts " + nClusters + " electron clusters", 100, 0.0, 10.).fill(trident.t());
                                // require clusters to be in time...
                                if (abs(dc1) < clusterDeltaTimeCut && abs(dc2) < clusterDeltaTimeCut && abs(dc3) < clusterDeltaTimeCut) {
                                    aida.histogram1D("trident energy after track and cluster timing cuts " + nClusters + " electron clusters", 100, 0.0, 10.).fill(trident.t());
                                    aida.histogram1D("Number of V0s in the event", 10, -0.5, 9.5).fill(nV0);
                                    aida.histogram1D("Number of VCs in the event", 10, -0.5, 9.5).fill(nVc);
                                    aida.histogram1D("trident energy after track and cluster timing cuts " + nClusters + " electron clusters " + nV0 + " V0 " + nVc + " Vc", 100, 0.0, 10.).fill(trident.t());
                                    // keep this event
                                    if (nClusters >= minNumberOfElectronsWithClusters) {
                                        skipEvent = false;
                                        aida.histogram1D("good trident run number", 750, 10000, 10750).fill(event.getRunNumber());
                                    }
                                }
                            }
                        }
                    }
                    aida.tree().cd("..");
                }// end of cut on only 3 RPs
            } // end of loop over RP collections
        }
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

    private BasicHepLorentzVector makeFourVector(ReconstructedParticle rp) {
        double emass2 = emass * emass;
        double eSquared = emass2;
        double[] mom = rp.getMomentum().v();
        for (int i = 0; i < 3; ++i) {
            eSquared += mom[i] * mom[i];
        }
        return new BasicHepLorentzVector(sqrt(eSquared), mom);
    }

    public void setMinNumberOfElectronsWithClusters(int i) {
        minNumberOfElectronsWithClusters = i;
    }

    public void setTrackDeltaTimeCut(double d) {
        trackDeltaTimeCut = d;
    }

    public void setClusterDeltaTimeCut(double d) {
        clusterDeltaTimeCut = d;
    }

}
