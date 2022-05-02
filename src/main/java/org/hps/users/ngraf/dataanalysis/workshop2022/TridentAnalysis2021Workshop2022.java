package org.hps.users.ngraf.dataanalysis.workshop2022;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackData;
import org.hps.recon.vertexing.BilliorTrack;
import org.hps.recon.vertexing.BilliorVertex;
import org.hps.recon.vertexing.BilliorVertexer;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class TridentAnalysis2021Workshop2022 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    double _esumCut = 2.0;
    private int _numberOfEventsSelected;
    private int _numberOfTridentsSelected;
    private int _numberOfEventsProcessed = 0;

    private double _minTridentClusterEnergy = 0.5;
    private double _maxTridentClusterEnergy = 2.5;
    String rpCollectionName = "FinalStateParticles_KF";
    String[] particles = {"positron", "electron1", "electron2"};
    Map<String, ReconstructedParticle> particleMap = new HashMap<>();

    protected double bField;
    private FieldMap fieldMap;
    BilliorVertexer vtxFitter;

    List<Cluster> tridentClusters = new ArrayList<>();
    private boolean _skimTrident = true;

    Cluster FeeCluster = null;

    protected void detectorChanged(Detector detector) {

        // Set the magnetic field parameters to the appropriate values.
        // we want the field near the target...
        Hep3Vector ip = new BasicHep3Vector(0., 0., 0.0);
        fieldMap = detector.getFieldMap();
        bField = fieldMap.getField(ip).y();
        for (int i = -50; i < 500; ++i) {
            double by = fieldMap.getField(new BasicHep3Vector(0., 0., i)).y();
//            System.out.println(i + " : " + by);
            aida.histogram1D("Fieldmap By at z", 550, -50., 500.).fill(i, -by);
        }

        vtxFitter = new BilliorVertexer(bField);
    }

    protected void process(EventHeader event) {
        if (event.getRunNumber() > 14623 && event.getRunNumber() < 14674) {
            _esumCut = 1.25;
            _minTridentClusterEnergy = 0.0;//0.2;//0.5;
            _maxTridentClusterEnergy = 1.25;//2.5;
        }
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");

        boolean isTridentCandidate = isTridentCandidate(ecalClusters);

        if (isTridentCandidate) {
            if (event.hasCollection(ReconstructedParticle.class, rpCollectionName)) {
                aida.tree().mkdirs(rpCollectionName);
                aida.tree().cd(rpCollectionName);
                aida.histogram1D("number of trident clusters", 10, -0.5, 9.5).fill(tridentClusters.size());
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, rpCollectionName);
                ReconstructedParticle positron = null;
                ReconstructedParticle electron1 = null;
                ReconstructedParticle electron2 = null;
                for (Cluster c : tridentClusters) {
                    for (ReconstructedParticle rp : rpList) {
                        //let's find the positron...
                        if (rp.getParticleIDUsed().getPDG() == -11) {
                            if (!rp.getClusters().isEmpty()) {
                                if (rp.getClusters().get(0) == c) {
                                    positron = rp;
                                }
                            }
                        }
                        //let's find the electrons...
                        if (rp.getParticleIDUsed().getPDG() == 11) {
                            if (!rp.getClusters().isEmpty()) {
                                if (rp.getClusters().get(0) == c) {
                                    if (electron1 == null) {
                                        electron1 = rp;
                                    } else {
                                        electron2 = rp;
                                    }
                                }
                            }
                        }
                    }
                }
                if (positron != null && electron1 != null && electron2 != null) {
                    particleMap.clear();
                    particleMap.put("positron", positron);
                    particleMap.put("electron1", electron1);
                    particleMap.put("electron2", electron2);
                    for (int i = 0; i < 3; ++i) {
                        aida.histogram1D(particles[i] + " energy", 100, 0., 5.).fill(particleMap.get(particles[i]).getEnergy());
                        aida.histogram1D(particles[i] + " momentum", 100, 0., 5.).fill(particleMap.get(particles[i]).getMomentum().magnitude());
                        aida.histogram1D(particles[i] + " eOverP", 100, 0., 2.).fill(particleMap.get(particles[i]).getEnergy() / particleMap.get(particles[i]).getMomentum().magnitude());
                        aida.histogram2D(particles[i] + " x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(particleMap.get(particles[i]).getClusters().get(0).getPosition()[0], particleMap.get(particles[i]).getClusters().get(0).getPosition()[1]);
                    }
                    double esum = 0.;
                    double pysum = 0.;
                    //let's look at the kinematics.
                    for (Cluster c : tridentClusters) {
                        esum += c.getEnergy();
                        pysum += c.getEnergy() * sinTheta(c.getPosition());
                    }
                    aida.histogram1D("trident esum", 100, 0., 10.).fill(esum);
                    aida.histogram1D("trident esum_y", 100, -1.0, 1.0).fill(pysum);
                    aida.histogram1D("trident pysum", 100, -1.0, 1.0).fill(positron.getMomentum().y() + electron1.getMomentum().y() + electron2.getMomentum().y());

                    //let's try vertexing same-side vs opposite side electrons with the positron...
                    String topOrBottom = positron.getMomentum().y() > 0 ? " top " : " bottom ";

                    //the tracks to vertex...
                    List<BilliorTrack> tracksToVertex = new ArrayList<BilliorTrack>();
                    // add the positron
                    tracksToVertex.add(new BilliorTrack(positron.getTracks().get(0)));
                    //is electron1 in the same hemisphere as the positron?
                    if (positron.getMomentum().y() * electron1.getMomentum().y() > 0) {
                        tracksToVertex.add(new BilliorTrack(electron1.getTracks().get(0)));
                    }
                    //is electron2 in the same hemisphere as the positron?
                    if (positron.getMomentum().y() * electron2.getMomentum().y() > 0) {
                        tracksToVertex.add(new BilliorTrack(electron2.getTracks().get(0)));
                    }
                    if (tracksToVertex.size() == 2) {
                        BilliorVertex vtx = vtxFitter.fitVertex(tracksToVertex);
                        aida.tree().mkdirs("positron" + topOrBottom + " same side");
                        aida.tree().cd("positron" + topOrBottom + " same side");
                        aida.histogram1D("vtx_x", 100, -5, 5).fill(vtx.getPosition().x());
                        aida.histogram1D("vtx_y", 100, -5, 5).fill(vtx.getPosition().y());
                        aida.histogram1D("vtx_z", 100, -20., 20.).fill(vtx.getPosition().z());
                        aida.histogram2D("vtx_x_y", 200, -5, 5, 200, -5, 5).fill(vtx.getPosition().x(), vtx.getPosition().y());
                        aida.histogram2D("vtx_x_z", 100, -5, 5, 100, -20., 20.).fill(vtx.getPosition().x(), vtx.getPosition().z());
                        aida.histogram2D("vtx_y_z", 100, -5, 5, 100, -20., 20.).fill(vtx.getPosition().y(), vtx.getPosition().z());
                        aida.histogram1D("vtx_chi2", 100, 0, 500).fill(vtx.getChi2());
                        aida.tree().cd("..");
                    }

                    //the tracks to vertex...
                    tracksToVertex.clear();
                    // add the positron
                    tracksToVertex.add(new BilliorTrack(positron.getTracks().get(0)));
                    //is electron1 in the same hemisphere as the positron?
                    if (positron.getMomentum().y() * electron1.getMomentum().y() < 0) {
                        tracksToVertex.add(new BilliorTrack(electron1.getTracks().get(0)));
                    }
                    //is electron2 in the same hemisphere as the positron?
                    if (positron.getMomentum().y() * electron2.getMomentum().y() < 0) {
                        tracksToVertex.add(new BilliorTrack(electron2.getTracks().get(0)));
                    }
                    if (tracksToVertex.size() == 2) {
                        BilliorVertex vtx = vtxFitter.fitVertex(tracksToVertex);
                        aida.tree().mkdirs("positron" + topOrBottom + " opposite side");
                        aida.tree().cd("positron" + topOrBottom + " opposite side");
                        aida.histogram1D("vtx_x", 100, -5, 5).fill(vtx.getPosition().x());
                        aida.histogram1D("vtx_y", 100, -5, 5).fill(vtx.getPosition().y());
                        aida.histogram1D("vtx_z", 100, -20., 20.).fill(vtx.getPosition().z());
                        aida.histogram2D("vtx_x_y", 200, -5, 5, 200, -5, 5).fill(vtx.getPosition().x(), vtx.getPosition().y());
                        aida.histogram2D("vtx_x_z", 100, -5, 5, 100, -20., 20.).fill(vtx.getPosition().x(), vtx.getPosition().z());
                        aida.histogram2D("vtx_y_z", 100, -5, 5, 100, -20., 20.).fill(vtx.getPosition().y(), vtx.getPosition().z());
                        aida.histogram1D("vtx_chi2", 100, 0, 500).fill(vtx.getChi2());
                        aida.tree().cd("..");
                    }
                }

                aida.tree().cd("..");
            }
            _numberOfTridentsSelected++;
            skipEvent = false;
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    private boolean isTridentCandidate(List<Cluster> ecalClusters) {
        boolean isTridentCandidate = false;
        tridentClusters.clear();
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
                        esum += c.getEnergy();
                        pysum += c.getEnergy() * sinTheta(c.getPosition());
                    }
                    aida.histogram1D("trident esum", 100, 0., 10.).fill(esum);
                    aida.histogram1D("trident esum_y", 100, -1.0, 1.0).fill(pysum);
                    if (esum > _esumCut && esum < 5. && abs(pysum) < 0.2) {
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

    private void analyzeReconstructedParticles(EventHeader event, Cluster FeeCluster) {
        List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class,
                "FinalStateParticles_KF");
        Track FeeTrack = null;
        // find the electron associated with this cluster (if any)
        for (ReconstructedParticle rp : rps) {
            if (rp.getParticleIDUsed().getPDG() == 11) {
                if (!rp.getClusters().isEmpty()) {
                    if (rp.getClusters().get(0) == FeeCluster) {
                        FeeTrack = rp.getTracks().get(0);
                    }
                }
            }
        }
        // Found a track associated with the FEE cluster
//        if (FeeTrack != null) {
//            RelationalTable trackToData = getKFTrackDataRelations(event);
//            double feeTrackTime = ((GenericObject) trackToData.from(FeeTrack)).getFloatVal(0);
//            aida.histogram1D("FEE Track time", 100, -30., 30.).fill(feeTrackTime);
//            // get all the tracks in the event, and plot delta time
//            if (FeeTrack.getTrackerHits().size() > 10) {
//                List<Track> tracks = event.get(Track.class, "KalmanFullTracks");
//                for (Track t : tracks) {
//                    if (t != FeeTrack) {
//                        double trackTime = ((GenericObject) trackToData.from(t)).getFloatVal(0);
//                        int nHits = t.getTrackerHits().size();
//                        aida.histogram1D("FEE Track time - track time", 100, -100., 100.).fill(feeTrackTime - trackTime);
//                        aida.histogram1D("FEE Track time - track time " + nHits, 100, -100., 100.).fill(feeTrackTime - trackTime);
//                    }
//                }
//            }
//        }
    }

    public RelationalTable getKFTrackDataRelations(EventHeader event) {

        List<TrackData> TrackData;
        RelationalTable trackToData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackRelations;
        TrackData trackdata;
        TrackData = event.get(TrackData.class,
                "KFTrackData");
        trackRelations = event.get(LCRelation.class,
                "KFTrackDataRelations");
        for (LCRelation relation : trackRelations) {
            if (relation != null && relation.getTo() != null) {
                trackToData.add(relation.getFrom(), relation.getTo());
            }
        }
        return trackToData;
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");

        System.out.println("Selected " + _numberOfTridentsSelected + " Tridents");

    }

}
