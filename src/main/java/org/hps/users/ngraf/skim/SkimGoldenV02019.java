package org.hps.users.ngraf.skim;

import static java.lang.Math.abs;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Skim off events containing Vertex object(s) Loop over both GBL and KF
 * collections
 *
 * @author Norman A. Graf
 */
public class SkimGoldenV02019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;
    String[] vertexCollectionNames = {"UnconstrainedV0Vertices", "UnconstrainedV0Vertices_KF"};

    private double _clusterDeltaTimeCut = 5.0;
    private double _trackDeltaTimeCut = 20.;
    private double _minMomentumCut = 0.5;
    private double _maxMomentumCut = 10.;
    private double _minEnergyCut = 0.5;
    private double _maxEnergyCut = 10.;

    @Override
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
                    for (Vertex v : vertices) {
                        ReconstructedParticle v0 = v.getAssociatedParticle();
                        int minNhits = 7;
                        if (v0.getType() == 1) {
                            minNhits = 14;
                        }
                        // this always has 2 tracks.
                        List<ReconstructedParticle> trks = v0.getParticles();
                        ReconstructedParticle neg = trks.get(0);
                        ReconstructedParticle pos = trks.get(1);
                        // let's also require both to have clusters
                        if (!neg.getClusters().isEmpty() && !pos.getClusters().isEmpty()) {
                            Cluster negClus = neg.getClusters().get(0);
                            int negIX = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                            int negIY = negClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
                            boolean negIsFiducial = TriggerModule.inFiducialRegion(negClus);
                            Cluster posClus = pos.getClusters().get(0);
                            int posIX = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
                            int posIY = posClus.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
                            boolean posIsFiducial = TriggerModule.inFiducialRegion(posClus);

                            // in time
                            double negTime = ClusterUtilities.getSeedHitTime(negClus);
                            double posTime = ClusterUtilities.getSeedHitTime(posClus);
                            double deltaTime = negTime - posTime;
                            double topMinusBottomTime = deltaTime;
                            if (negClus.getPosition()[1] < 0) {
                                topMinusBottomTime = -deltaTime;
                            }
                            aida.histogram1D("top - bottom cluster delta time", 100, -5., 5.).fill(deltaTime);
                            aida.histogram1D("cluster pair delta time", 100, -5., 5.).fill(deltaTime);
                            double negE = negClus.getEnergy();
                            double posE = posClus.getEnergy();
                            int negNhits = neg.getTracks().get(0).getTrackerHits().size();
                            int posNhits = pos.getTracks().get(0).getTrackerHits().size();
                            double negMom = neg.getMomentum().magnitude();
                            double posMom = pos.getMomentum().magnitude();
                            double negEoverP = negE / negMom;
                            double posEoverP = posE / posMom;

                            if (negNhits >= minNhits && posNhits >= minNhits && abs(deltaTime) < _clusterDeltaTimeCut) {
                                if (negMom >= _minMomentumCut && negMom <= _maxMomentumCut) {
                                    if (posMom >= _minMomentumCut && posMom <= _maxMomentumCut) {
                                        if (negE >= _minEnergyCut && negE <= _maxEnergyCut) {
                                            if (posE >= _minEnergyCut && posE <= _maxEnergyCut) {
                                                aida.histogram1D("negative momentum", 100, 0., 6.0).fill(negMom);
                                                aida.histogram1D("positive momentum", 100, 0., 6.0).fill(posMom);
                                                aida.histogram2D("negative cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(negIX, negIY);
                                                aida.histogram2D("positive cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(posIX, posIY);
                                                aida.histogram1D("negative E over P", 100, 0., 2.).fill(negEoverP);
                                                aida.histogram1D("positive E over P", 100, 0., 2.).fill(posEoverP);
                                                if (negIsFiducial) {
                                                    aida.histogram1D("negative momentum fiducial", 100, 0., 6.0).fill(negMom);
                                                    aida.histogram1D("negative E over P fiducial", 100, 0., 2.).fill(negEoverP);
                                                }
                                                if (posIsFiducial) {
                                                    aida.histogram1D("positive momentum fiducial", 100, 0., 6.0).fill(posMom);
                                                    aida.histogram1D("positive E over P fiducial", 100, 0., 2.).fill(posEoverP);
                                                }
                                                skipEvent = false;
                                            }
                                        }
                                    }
                                }
                            } // end of check on number of track hits
                        } // end of check on clusters
                    }//end of loop over vertices
                } // end of loop over check on collection
            } // end of loop over vertex collections 
            aida.tree().cd("..");
        }
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

    public void setTrackDeltaTimeCut(double d) {
        _trackDeltaTimeCut = d;
    }
}
