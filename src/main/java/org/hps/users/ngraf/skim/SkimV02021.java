package org.hps.users.ngraf.skim;

import static java.lang.Math.abs;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackData;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Vertex;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Skim off events containing Vertex object(s) Loop over both GBL and KF
 * collections
 *
 * @author Norman A. Graf
 */
public class SkimV02021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;
    String[] vertexCollectionNames = {"UnconstrainedV0Vertices_KF"};//, "UnconstrainedV0Vertices"};
    private double _maxNumberOfVertices = 1;
    private double _deltaTrackTimeCut = 20.;
    private int _minNhitsOnTrack = 11;
    private double _minMomentumCut = 0.5;
    private double _maxMomentumCut = 10.;
    private boolean _requireClusterMatch = true;
    private boolean _requireFiducialClusters = true;
    private double _deltaClusterTimeCut = 2.0;

    public void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;
        for (String vertexCollectionName : vertexCollectionNames) {
            aida.tree().mkdirs(vertexCollectionName);
            aida.tree().cd(vertexCollectionName);
            if (event.hasCollection(Vertex.class, vertexCollectionName)) {
                List<Vertex> vertices = event.get(Vertex.class, vertexCollectionName);
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                if (vertices.size() > 0 && vertices.size() <= _maxNumberOfVertices) {
                    RelationalTable trackToData = getKFTrackDataRelations(event);
                    for (Vertex v : vertices) {
                        ReconstructedParticle v0 = v.getAssociatedParticle();
                        // this always has 2 tracks.
                        List<ReconstructedParticle> trks = v0.getParticles();
                        ReconstructedParticle neg = trks.get(0);
                        ReconstructedParticle pos = trks.get(1);

                        int negNhits = neg.getTracks().get(0).getTrackerHits().size();
                        int posNhits = pos.getTracks().get(0).getTrackerHits().size();
                        double negMom = neg.getMomentum().magnitude();
                        double posMom = pos.getMomentum().magnitude();
                        double negTrackTime = ((GenericObject) trackToData.from(neg.getTracks().get(0))).getFloatVal(0);
                        double posTrackTime = ((GenericObject) trackToData.from(pos.getTracks().get(0))).getFloatVal(0);
                        double deltaTrackTime = negTrackTime - posTrackTime;

                        boolean posTrackIsTop = pos.getMomentum().y() > 0.;
                        boolean negTrackIsTop = neg.getMomentum().y() > 0.;
                        String posTorB = posTrackIsTop ? " top " : " bottom ";
                        String negTorB = negTrackIsTop ? " top " : " bottom ";

                        if (negNhits >= _minNhitsOnTrack && posNhits >= _minNhitsOnTrack && abs(deltaTrackTime) < _deltaTrackTimeCut) {
                            if (negMom >= _minMomentumCut && negMom <= _maxMomentumCut) {
                                if (posMom >= _minMomentumCut && posMom <= _maxMomentumCut) {
                                    aida.histogram1D("delta track time", 100, -20., 20.).fill(deltaTrackTime);
                                    aida.histogram1D("negative momentum", 100, 0., 6.0).fill(negMom);
                                    aida.histogram1D("positive momentum", 100, 0., 6.0).fill(posMom);

                                    aida.histogram1D("negative momentum" + negTorB, 100, 0., 6.0).fill(negMom);
                                    aida.histogram1D("positive momentum" + posTorB, 100, 0., 6.0).fill(posMom);
                                    aida.histogram1D("negative nHits" + negTorB, 20, -0.5, 19.5).fill(negNhits);
                                    aida.histogram1D("positive nHits" + posTorB, 20, -0.5, 19.5).fill(posNhits);

                                    skipEvent = false;
                                    if (_requireClusterMatch) {
                                        skipEvent = true;
                                        // let's also require both to have clusters
                                        if (!neg.getClusters().isEmpty() && !pos.getClusters().isEmpty()) {
                                            Cluster negClus = neg.getClusters().get(0);
                                            Cluster posClus = pos.getClusters().get(0);
                                            // in time
                                            double negTime = ClusterUtilities.getSeedHitTime(negClus);
                                            double posTime = ClusterUtilities.getSeedHitTime(posClus);
                                            double deltaClusterTime = negTime - posTime;
                                            aida.histogram1D("delta cluster time", 100, -20., 20.).fill(deltaClusterTime);
                                            if (abs(deltaClusterTime) < _deltaClusterTimeCut) {
                                                skipEvent = false;
                                                // fiducial
                                                boolean negIsFiducial = TriggerModule.inFiducialRegion(negClus);
                                                boolean posIsFiducial = TriggerModule.inFiducialRegion(posClus);
                                                if (_requireFiducialClusters) {
                                                    if (!negIsFiducial) {
                                                        skipEvent = true;
                                                    }
                                                    if (!posIsFiducial) {
                                                        skipEvent = true;
                                                    }
                                                    if (negIsFiducial && posIsFiducial) {
                                                        aida.histogram1D("delta track time two fiducial clusters", 100, -20., 20.).fill(deltaTrackTime);
                                                        aida.histogram1D("negative momentum two fiducial clusters", 100, 0., 6.0).fill(negMom);
                                                        aida.histogram1D("positive momentum two fiducial clusters", 100, 0., 6.0).fill(posMom);
                                                        aida.histogram1D("delta cluster time two fiducial clusters", 100, -20., 20.).fill(deltaClusterTime);

                                                        aida.histogram1D("negative momentum" + negTorB + " two fiducial clusters", 100, 0., 6.0).fill(negMom);
                                                        aida.histogram1D("positive momentum" + posTorB + " two fiducial clusters", 100, 0., 6.0).fill(posMom);
                                                        aida.histogram1D("negative nHits" + negTorB + " two fiducial clusters", 20, -0.5, 19.5).fill(negNhits);
                                                        aida.histogram1D("positive nHits" + posTorB + " two fiducial clusters", 20, -0.5, 19.5).fill(posNhits);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } // end of check on number of track hits, delta time and momentum
                    }//end of loop over vertices
                }// end of check on existence of vertices
            } // end of loop over check on collection
            aida.tree().cd("..");
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

    public void setMinNhitsOnTrack(int i) {
        _minNhitsOnTrack = i;
    }

    public void setDeltaTrackTimeCut(double d) {
        _deltaTrackTimeCut = d;
    }

    public void setMinMomentumCut(double d) {
        _minMomentumCut = d;
    }

    public void setMaxMomentumCut(double d) {
        _maxMomentumCut = d;
    }

    public void setRequireClusterMatch(boolean b) {
        _requireClusterMatch = b;
    }

    public void setRequireFiducialClusters(boolean b) {
        _requireFiducialClusters = b;
    }

    public void setDeltaClusterTimeCut(double d) {
        _deltaClusterTimeCut = d;
    }

    public void setMaxNumberOfVertices(int i) {
        _maxNumberOfVertices = i;
    }
}
