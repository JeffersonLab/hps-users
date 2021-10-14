package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import static java.lang.Math.abs;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.tracking.TrackData;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class MollerTrackOnlyAnalysis2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    private boolean _skipEvent = true;
    private boolean _skimEvents = false;
    final double pScale = 1.761 / 2.09;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    @Override
    protected void process(EventHeader event) {
        _skipEvent = true;
        _numberOfEventsProcessed++;
        // simple cluster analysis
        aida.tree().mkdirs("Moller Track-only");
        aida.tree().cd("Moller Track-only");
        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
        List<ReconstructedParticle> topElectrons = new ArrayList<>();
        List<ReconstructedParticle> bottomElectrons = new ArrayList<>();

        for (ReconstructedParticle rp : rpList) {
            int pdgId = rp.getParticleIDUsed().getPDG();
            if (pdgId == 11) {
                //let's find Moller electron candidates with minimal cuts
//                //get rid of FEEs
//                if (rp.getMomentum().magnitude() < 2.8) {
                // require the track to have 12 or more hits...
                if (rp.getTracks().get(0).getTrackerHits().size() >= 11) {
                    if (rp.getMomentum().y() > 0) {
                        topElectrons.add(rp);
                    } else {
                        bottomElectrons.add(rp);
                    }
                }
//                }
            }
        }
        aida.histogram2D("number of top electrons vs number of bottom electrons", 5, -0.5, 4.5, 5, -0.5, 4.5).fill(topElectrons.size(), bottomElectrons.size());
        //let's start with one top and one bottom electron
        if (topElectrons.size() == 1 && bottomElectrons.size() == 1) {
            ReconstructedParticle topElectron = topElectrons.get(0);
            ReconstructedParticle bottomElectron = bottomElectrons.get(0);

            RelationalTable trackToData = getKFTrackDataRelations(event);
            Track topTrack = topElectron.getTracks().get(0);
            Track bottomTrack = bottomElectron.getTracks().get(0);
            double deltaTrackTime = ((GenericObject) trackToData.from(topTrack)).getFloatVal(0) - ((GenericObject) trackToData.from(bottomTrack)).getFloatVal(0);

            double pTop = topElectron.getMomentum().magnitude();
            double pBottom = bottomElectron.getMomentum().magnitude();
            double pBottomCorr = pBottom * pScale;
            double psum = pTop + pBottom;
            double pdiff = pTop - pBottom;
            double psumCorr = pTop + pBottomCorr;
            double pdiffCorr = pTop - pBottomCorr;

            aida.histogram1D("track delta time", 100, -100., 100.).fill(deltaTrackTime);
            aida.histogram1D("top track number of hits on track", 20, -0.5, 19.5).fill(topTrack.getTrackerHits().size());
            aida.histogram1D("bottom track number of hits on track", 20, -0.5, 19.5).fill(bottomTrack.getTrackerHits().size());
            aida.histogram1D("top electron track momentum", 100, 0., 6.0).fill(pTop);
            aida.histogram1D("bottom electron track momentum", 100, 0., 6.0).fill(pBottom);
            aida.histogram2D("top momentum vs bottom momentum", 100, 0., 3.0, 100, 0., 3.0).fill(pTop, pBottom);
            aida.histogram1D("bottom electron track momentum corr", 100, 0., 6.0).fill(pBottomCorr);
            aida.histogram2D("top momentum vs bottom momentum corr", 100, 0., 3.0, 100, 0., 3.0).fill(pTop, pBottomCorr);
            aida.histogram1D("psum", 100, 0., 7.).fill(psum);
            aida.histogram1D("pdiff", 100, -3.0, 3.0).fill(pdiff);
            aida.histogram1D("psumCorr", 100, 0., 7.).fill(psumCorr);
            aida.histogram1D("pdiffCorr", 100, -3.0, 3.0).fill(pdiffCorr);

            if (abs(deltaTrackTime) < 4) {
                aida.tree().mkdirs("delta time cut");
                aida.tree().cd("delta time cut");

                aida.histogram1D("track delta time", 100, -100., 100.).fill(deltaTrackTime);
                aida.histogram1D("top track number of hits on track", 20, -0.5, 19.5).fill(topTrack.getTrackerHits().size());
                aida.histogram1D("bottom track number of hits on track", 20, -0.5, 19.5).fill(bottomTrack.getTrackerHits().size());
                aida.histogram1D("top electron track momentum", 100, 0., 6.0).fill(pTop);
                aida.histogram1D("bottom electron track momentum", 100, 0., 6.0).fill(pBottom);
                aida.histogram2D("top momentum vs bottom momentum", 100, 0., 3.0, 100, 0., 3.0).fill(pTop, pBottom);
                aida.histogram1D("bottom electron track momentum corr", 100, 0., 6.0).fill(pBottomCorr);
                aida.histogram2D("top momentum vs bottom momentum corr", 100, 0., 3.0, 100, 0., 3.0).fill(pTop, pBottomCorr);
                aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                aida.histogram1D("pdiff", 100, -3.0, 3.0).fill(pdiff);
                aida.histogram1D("psumCorr", 100, 0., 7.).fill(psumCorr);
                aida.histogram1D("pdiffCorr", 100, -3.0, 3.0).fill(pdiffCorr);
                _skipEvent = false;
                aida.tree().cd("..");
            } else {
                aida.histogram1D("psumCorr out of time", 100, 0., 7.).fill(psumCorr);
            }
        }
        aida.tree().cd("..");
        if (_skipEvent) {
            if(_skimEvents) throw new Driver.NextEventException();
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

    public void setSkimEvents(boolean b) {
        _skimEvents = b;
    }
}
