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

    private double _trackDeltaTimeCut = 20.;
    private int _minNhitsOnTrack = 11;
    private double _minMomentumCut = 0.5;
    private double _maxMomentumCut = 10.;

    public void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;
        for (String vertexCollectionName : vertexCollectionNames) {
            aida.tree().mkdirs(vertexCollectionName);
            aida.tree().cd(vertexCollectionName);
            if (event.hasCollection(Vertex.class, vertexCollectionName)) {
                List<Vertex> vertices = event.get(Vertex.class, vertexCollectionName);
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                if (vertices.size() > 0) {
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
                        double deltaTime = negTrackTime - posTrackTime;

                        if (negNhits >= _minNhitsOnTrack && posNhits >= _minNhitsOnTrack && abs(deltaTime) < _trackDeltaTimeCut) {
                            if (negMom >= _minMomentumCut && negMom <= _maxMomentumCut) {
                                if (posMom >= _minMomentumCut && posMom <= _maxMomentumCut) {
                                    aida.histogram1D("delta track time", 100, -20., 20.).fill(deltaTime);
                                    aida.histogram1D("negative momentum", 100, 0., 6.0).fill(negMom);
                                    aida.histogram1D("positive momentum", 100, 0., 6.0).fill(posMom);
                                    skipEvent = false;
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

    public void setMinNhitsOnTrack(int i)
    {
        _minNhitsOnTrack = i;
    }
    public void setTrackDeltaTimeCut(double d) {
        _trackDeltaTimeCut = d;
    }

    public void setMinMomentumCut(double d) {
        _minMomentumCut = d;
    }

    public void setMaxMomentumCut(double d) {
        _maxMomentumCut = d;
    }

}
