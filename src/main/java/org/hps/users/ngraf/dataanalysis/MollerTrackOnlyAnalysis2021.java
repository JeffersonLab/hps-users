package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
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
import org.lcsim.util.fourvec.Lorentz4Vector;
import org.lcsim.util.fourvec.Momentum4Vector;

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
    private double pScale = 1.761 / 2.09;
    private double fee_cut = 2.8;
    private double _maxDeltaTrackTime = 4.;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    @Override
    protected void process(EventHeader event) {

        int run = event.getRunNumber();
        if (run > 14624 && run < 14685) {
            fee_cut = 1.5;
            pScale = 0.853 / 1.13; //1.92 / 2.87;
        }
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
                        if (rp.getMomentum().magnitude() < fee_cut) {
                            topElectrons.add(rp);
                        }
                    } else {
                        if (rp.getMomentum().magnitude() * pScale < fee_cut) {
                            bottomElectrons.add(rp);
                        }
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
            aida.histogram1D("top electron track momentum", 100, 0., 5.0).fill(pTop);
            aida.histogram1D("bottom electron track momentum", 100, 0., 5.0).fill(pBottom);
            aida.histogram2D("top momentum vs bottom momentum", 100, 0., 5.0, 100, 0., 5.0).fill(pTop, pBottom);
            aida.histogram1D("bottom electron track momentum corr", 100, 0., 5.0).fill(pBottomCorr);
            aida.histogram2D("top momentum vs bottom momentum corr", 100, 0., 5.0, 100, 0., 5.0).fill(pTop, pBottomCorr);
            aida.histogram1D("psum", 100, 0., 7.).fill(psum);
            aida.histogram1D("pdiff", 100, -3.0, 3.0).fill(pdiff);
            aida.histogram1D("psumCorr", 100, 0., 7.).fill(psumCorr);
            aida.histogram1D("pdiffCorr", 100, -3.0, 3.0).fill(pdiffCorr);

            if (abs(deltaTrackTime) < _maxDeltaTrackTime) {
                aida.tree().mkdirs("delta time cut");
                aida.tree().cd("delta time cut");

                aida.histogram1D("track delta time", 100, -100., 100.).fill(deltaTrackTime);
                aida.histogram1D("top track number of hits on track", 20, -0.5, 19.5).fill(topTrack.getTrackerHits().size());
                aida.histogram1D("bottom track number of hits on track", 20, -0.5, 19.5).fill(bottomTrack.getTrackerHits().size());
                aida.histogram1D("top electron track momentum", 100, 0., 5.0).fill(pTop);
                aida.histogram1D("bottom electron track momentum", 100, 0., 5.0).fill(pBottom);
                aida.histogram2D("top momentum vs bottom momentum", 100, 0., 5.0, 100, 0., 5.0).fill(pTop, pBottom);
                aida.histogram1D("bottom electron track momentum corr", 100, 0., 5.0).fill(pBottomCorr);
                aida.histogram2D("top momentum vs bottom momentum corr", 100, 0., 5.0, 100, 0., 5.0).fill(pTop, pBottomCorr);
                aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                aida.histogram1D("pdiff", 100, -3.0, 3.0).fill(pdiff);
                aida.histogram1D("psumCorr", 100, 0., 7.).fill(psumCorr);
                aida.histogram1D("pdiffCorr", 100, -3.0, 3.0).fill(pdiffCorr);
                _skipEvent = false;
                boolean topHasCluster = !topElectron.getClusters().isEmpty();
                boolean bottomHasCluster = !bottomElectron.getClusters().isEmpty();
                if (topHasCluster) {
                    double trackTime = ((GenericObject) trackToData.from(topTrack)).getFloatVal(0);
                    double clusterTime = ClusterUtilities.findSeedHit(topElectron.getClusters().get(0)).getTime();
                    aida.histogram1D("top track-cluster deltaT", 100, -70., -20.).fill(trackTime - clusterTime);
                    aida.histogram1D("psumCorr top track has cluster", 100, 0., 7.).fill(psumCorr);
                }
                if (bottomHasCluster) {
                    double trackTime = ((GenericObject) trackToData.from(bottomTrack)).getFloatVal(0);
                    double clusterTime = ClusterUtilities.findSeedHit(bottomElectron.getClusters().get(0)).getTime();
                    aida.histogram1D("bottom track-cluster deltaT", 100, -70., -20.).fill(trackTime - clusterTime);
                    aida.histogram1D("psumCorr bottom track has cluster", 100, 0., 7.).fill(psumCorr);
                }
                if (topHasCluster && bottomHasCluster) {
                    aida.histogram1D("psumCorr top and bottom tracks have cluster", 100, 0., 7.).fill(psumCorr);
                }
                if (topHasCluster || bottomHasCluster) {
                    aida.histogram1D("psumCorr top or bottom track has cluster", 100, 0., 7.).fill(psumCorr);
                }
                if (!topHasCluster && !bottomHasCluster) {
                    aida.histogram1D("psumCorr neither top nor bottom track has cluster", 100, 0., 7.).fill(psumCorr);
                }
                double invMass = analyzeTwoElectrons(topElectron, bottomElectron);
                aida.histogram2D("psumCorr vs invariant mass", 100, 0., 3., 100, 0.0, 0.1).fill(psumCorr, invMass);
                aida.tree().cd("..");
            } else {
                aida.histogram1D("psumCorr out of time", 100, 0., 7.).fill(psumCorr);
            }
        }
        aida.tree().cd("..");
        if (_skipEvent) {
            if (_skimEvents) {
                throw new Driver.NextEventException();
            }
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

    private double analyzeTwoElectrons(ReconstructedParticle ele1, ReconstructedParticle ele2) {
        // note that e1 is the top electron
        double emass = 0.000511;
        double emass2 = emass * emass;
        double[] p1 = ele1.getMomentum().v();
        double[] p2 = ele2.getMomentum().v();
        // apply ad hoc correction to bottom momentum
        for (int i = 0; i < 2; ++i) {
            p2[i] *= p2[i] * pScale;
        }
        double e1 = 0;
        double e2 = 0.;
        for (int i = 0; i < 3; ++i) {
            e1 += p1[i] * p1[i];
            e2 += p2[i] * p2[i];
        }
        e1 = sqrt(e1 + emass2);
        e2 = sqrt(e2 + emass2);
        Momentum4Vector evec1 = new Momentum4Vector(p1[0], p1[1], p1[2], e1);
        Momentum4Vector evec2 = new Momentum4Vector(p2[0], p2[1], p2[2], e2);
        Lorentz4Vector eesum = evec1.plus(evec2);
        double eemass = eesum.mass();
        aida.histogram1D("invariant mass", 100, 0.0, 0.1).fill(eemass);
        return eemass;
    }

    public void setSkimEvents(boolean b) {
        _skimEvents = b;
    }

    public void setMaxDeltaTrackTime(double d) {
        _maxDeltaTrackTime = d;
    }
}
