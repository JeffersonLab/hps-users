package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackData;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;
import org.lcsim.math.chisq.ChisqProb;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class FeeAnalysis2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private boolean _analyzeTracksByNhits = true;
    private boolean _analyzeHitsInFit = false;
    private boolean _analyzeMomentumByCalorimeterCell = true;
    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF", "OtherElectrons"};//, "OtherElectrons", "OtherElectrons_KF"};

    private int _maxNumberOfClusters = 3;
    private double _minClusterEnergy = 2.5;
    RelationalTable _hitToStrips;
    RelationalTable _hitToRotated;
    RelationalTable _trackToData;

    Map<String, DescriptiveStatistics> _channelMaps = new HashMap<>();

    private SvtSensorTimeCorrector timeCorrector;

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed;
    private boolean _requireFiducial = true;
    private boolean _skipEvent;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
        timeCorrector = new SvtSensorTimeCorrector();
    }

    protected void process(EventHeader event) {
        _hitToStrips = TrackUtils.getHitToStripsTable(event);
        _hitToRotated = TrackUtils.getHitToRotatedTable(event);

        _skipEvent = true;
        _numberOfEventsProcessed++;

        // 2021 has data at 3.74 and 1.92GeV
        // adjust cuts accordingly
        if (event.getRunNumber() > 14623 && event.getRunNumber() < 14674) {
            _minClusterEnergy = 1.2;
        }
        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        int nClusters = clusters.size();
        aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
        Cluster maxEnergyCluster = null;
        double maxEnergy = 0.;
        // select the maximum energy cluster in the event
        for (Cluster cluster : clusters) {
            if (cluster.getEnergy() > maxEnergy) {
                maxEnergy = cluster.getEnergy();
                maxEnergyCluster = cluster;
            }
            aida.histogram2D("All Clusters x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
            int ix = seed.getIdentifierFieldValue("ix");
            int iy = seed.getIdentifierFieldValue("iy");
            aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            if (cluster.getPosition()[1] > 0.) {
                aida.histogram1D("Top cluster energy ", 100, 0.0, 5.0).fill(cluster.getEnergy());
                if (nClusters == 1) {
                    aida.histogram1D("Top single cluster energy ", 100, 0.0, 5.0).fill(cluster.getEnergy());
                }
            } else {
                aida.histogram1D("Bottom cluster energy ", 100, 0.0, 5.0).fill(cluster.getEnergy());
                if (nClusters == 1) {
                    aida.histogram1D("Bottom single cluster energy ", 100, 0.0, 5.0).fill(cluster.getEnergy());
                }
            }
        }
        aida.histogram1D("Maximum cluster energy ", 100, 0.0, 5.0).fill(maxEnergy);

        setupSensors(event);
        _trackToData = getKFTrackDataRelations(event);

        // let's look at ReconstructedParticles and find the matching RP...
        for (String rpCollectionName : ReconstructedParticleCollectionNames) {
            if (event.hasCollection(ReconstructedParticle.class, rpCollectionName)) {
                aida.tree().mkdirs(rpCollectionName);
                aida.tree().cd(rpCollectionName);
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, rpCollectionName);
                ReconstructedParticle matchedRP = null;
                for (ReconstructedParticle rp : rpList) {
                    if (abs(rp.getParticleIDUsed().getPDG()) == 11) {
                        if (!rp.getClusters().isEmpty()) {
                            if (rp.getClusters().get(0) == maxEnergyCluster) {
                                matchedRP = rp;
                            }
                        }
                    }
                }
                if (matchedRP != null) {
                    if (nClusters == 1 && matchedRP.getClusters().get(0).getEnergy() > _minClusterEnergy) {
                        analyzeReconstructedParticle(matchedRP);
                    }
                }
                aida.tree().cd("..");
            }
        }
        if (_skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    void analyzeReconstructedParticle(ReconstructedParticle rp) {
        boolean isElectron = rp.getParticleIDUsed().getPDG() == 11;
        boolean isPositron = rp.getParticleIDUsed().getPDG() == -11;
        boolean isPhoton = rp.getParticleIDUsed().getPDG() == 22;
        String type = "";
        if (isElectron) {
            type = "electron ";
        }
        if (isPositron) {
            type = "positron ";
        }
        if (isPhoton) {
            type = "photon ";
        }
        if (isElectron || isPositron) {
            if (rp.getClusters().size() == 1) {
                analyzeTrack(rp);
            }
        }
    }

    void analyzeTrack(ReconstructedParticle rp) {
        boolean isGBL = TrackType.isGBL(rp.getType());
        int minNhits = isGBL ? 6 : 13;
        String trackDir = isGBL ? "gbl" : "htf";
        if (rp.getType() == 1) {
            trackDir = "kf";
        }
        aida.tree().mkdirs(trackDir);
        aida.tree().cd(trackDir);

        Track t = rp.getTracks().get(0);
        int nHits = t.getTrackerHits().size();
        if (nHits >= minNhits) {
            // ensure that there are no crossovers
            boolean[] holeOrSlot = isHoleorSlotTrack(t);
            boolean isHole = holeOrSlot[0] && !holeOrSlot[1];
            boolean isSlot = !holeOrSlot[0] && holeOrSlot[1];
            String side = " overlap ";
            if (isHole) {
                side = " hole ";
            }
            if (isSlot) {
                side = " slot ";
            }
//        aida.cloud1D("ReconstructedParticle Type").fill(rp.getType());
//        aida.cloud1D("Track Type").fill(t.getType());
            //rotate into physiscs frame of reference
            Hep3Vector rprot = VecOp.mult(beamAxisRotation, rp.getMomentum());
            double theta = Math.acos(rprot.z() / rprot.magnitude());
            double chiSquared = t.getChi2();
            int ndf = t.getNDF();
            double chi2Ndf = t.getChi2() / t.getNDF();
            double chisqProb = 1.;
            if (ndf != 0) {
                chisqProb = ChisqProb.gammp(ndf, chiSquared);
            }
            double dEdx = t.getdEdx();
            double e = rp.getEnergy();
            double p = rp.getMomentum().magnitude();
//            String topOrBottom = isTopTrack(t) ? " top " : " bottom ";
            String topOrBottom = rp.getMomentum().y() > 0 ? " top " : " bottom ";
            aida.histogram1D("Track chisq per df" + topOrBottom, 100, 0., 50.).fill(chiSquared / ndf);
            aida.histogram1D("Track chisq prob" + topOrBottom, 100, 0., 1.).fill(chisqProb);
            aida.histogram1D("Track nHits" + topOrBottom, 14, 0.5, 14.5).fill(t.getTrackerHits().size());
            aida.histogram1D("Track momentum" + topOrBottom, 200, 0., 10.0).fill(p);
            aida.histogram1D("Track momentum" + topOrBottom + side, 200, 0., 10.0).fill(p);
            aida.histogram1D("Track deDx" + topOrBottom, 100, 0.00004, 0.00013).fill(t.getdEdx());
            aida.histogram1D("Track theta" + topOrBottom, 100, 0.010, 0.160).fill(theta);
            aida.histogram2D("Track theta vs p" + topOrBottom, 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
            aida.histogram1D("rp x0" + topOrBottom, 100, -0.50, 0.50).fill(TrackUtils.getX0(t));
            aida.histogram1D("rp y0" + topOrBottom, 100, -5.0, 5.0).fill(TrackUtils.getY0(t));
            aida.histogram1D("rp z0" + topOrBottom, 100, -1.0, 1.0).fill(TrackUtils.getZ0(t));

            double tanLambda = t.getTrackStates().get(0).getTanLambda();
            aida.histogram1D("track tanlambda" + topOrBottom + " " + nHits + " hits", 50, 0.01, 0.07).fill(abs(tanLambda));
            aida.histogram2D("Track tanlambda vs p" + topOrBottom + " " + nHits + " hits", 50, 0.01, 0.07, 100, 0., 7.0).fill(abs(tanLambda), p);
            aida.profile1D("Track tanlambda vs p profile" + topOrBottom + " " + nHits + " hits", 50, 0.01, 0.07).fill(abs(tanLambda), p);

            //
            if (_analyzeTracksByNhits) {
                aida.histogram1D("Track chisq per df" + topOrBottom + " " + nHits + " hits", 100, 0., 50.).fill(chiSquared / ndf);
                aida.histogram1D("Track chisq prob" + topOrBottom + " " + nHits + " hits", 100, 0., 1.).fill(chisqProb);
                aida.histogram1D("Track momentum" + topOrBottom + " " + nHits + " hits", 200, 0., 10.0).fill(p);
                aida.histogram1D("Track deDx" + topOrBottom + " " + nHits + " hits", 100, 0.00004, 0.00013).fill(t.getdEdx());
                aida.histogram1D("Track theta" + topOrBottom + " " + nHits + " hits", 100, 0.010, 0.160).fill(theta);
                aida.histogram2D("Track theta vs p" + topOrBottom + " " + nHits + " hits", 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
                aida.histogram1D("rp x0" + topOrBottom + " " + nHits + " hits", 100, -0.50, 0.50).fill(TrackUtils.getX0(t));
                aida.histogram1D("rp y0" + topOrBottom + " " + nHits + " hits", 100, -5.0, 5.0).fill(TrackUtils.getY0(t));
                aida.histogram1D("rp z0" + topOrBottom + " " + nHits + " hits", 100, -1.0, 1.0).fill(TrackUtils.getZ0(t));
            }
            boolean hasCluster = rp.getClusters().size() == 1;
            if (hasCluster) {
                //let's look at E/p and timing
                Cluster c = rp.getClusters().get(0);
                if (_analyzeHitsInFit && rp.getType() == 1) {
                    analyzeHitsInFit(rp);
                }
                boolean isFiducial = TriggerModule.inFiducialRegion(c);
                CalorimeterHit seed = c.getCalorimeterHits().get(0);
                int ix = seed.getIdentifierFieldValue("ix");
                int iy = seed.getIdentifierFieldValue("iy");
                aida.histogram2D("track momentum vs cluster x iy: " + iy, 320, -270.0, 370.0, 100, 0., 10.0).fill(c.getPosition()[0], p);
                aida.histogram1D("EoverP" + topOrBottom, 100, 0., 2.).fill(e / p);
                aida.histogram2D("E vs p" + topOrBottom, 100, 2., 7., 100, 2., 7.).fill(e, p);

                //Analyze track-cluster matching...
                //note that the reference point coordinates are listed as (z,x,y).
                double[] tpos = t.getTrackStates().get(t.getTrackStates().size() - 1).getReferencePoint();
                double tx = tpos[1];
                double ty = tpos[2];
                double tz = tpos[0];
                double[] cpos = c.getPosition();
                double cx = cpos[0];
                double cy = cpos[1];
                double cz = cpos[2];

                double dx = cx - tx;
                double dy = cy - ty;
                double thetaY = abs(Math.asin(rprot.y() / rprot.magnitude()));

                aida.cloud1D("cluster z").fill(cz);
                aida.cloud1D("track z").fill(tz);
                aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cx, cy);
                aida.histogram2D("track x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(tx, ty);

                aida.histogram1D("cluster x - track x", 100, -15., 15.).fill(dx);
                aida.histogram2D("cluster x - track x vs cluster x", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                if (cy > 0) {
                    aida.histogram1D("cluster x - track x top", 100, -15., 15.).fill(dx);
                    aida.histogram2D("cluster x - track x vs cluster x top", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                } else {
                    aida.histogram1D("cluster x - track x bottom", 100, -15., 15.).fill(dx);
                    aida.histogram2D("cluster x - track x vs cluster x bottom", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                }
                aida.histogram2D("cluster x - track x vs cluster y", 90, -90.0, 90.0, 100, -15., 15.).fill(cy, dx);
                aida.histogram1D("cluster y - track y", 100, -10., 10.).fill(dy);
                aida.histogram2D("cluster y - track y vs cluster x", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                aida.histogram2D("cluster y - track y vs cluster y", 90, -90.0, 90.0, 100, -10., 10.).fill(cy, dy);
                aida.histogram2D("cluster y - track y vs track y", 90, -90.0, 90.0, 100, -10., 10.).fill(ty, dy);
                if (cy > 0) {
                    aida.histogram1D("ThetaY top", 100, 0.01, 0.07).fill(thetaY);
                    aida.histogram1D("cluster y - track y top", 100, -10., 10.).fill(dy);
                    aida.histogram2D("cluster y - track y vs cluster x top", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                    aida.histogram2D("cluster y - track y vs thetaY top", 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                } else {
                    aida.histogram1D("ThetaY bottom", 100, 0.01, 0.07).fill(thetaY);
                    aida.histogram1D("cluster y - track y bottom", 100, -10., 10.).fill(dy);
                    aida.histogram2D("cluster y - track y vs cluster x bottom", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                    aida.histogram2D("cluster y - track y vs thetaY bottom", 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                }
                if (isFiducial) {
                    aida.tree().mkdirs("fiducial");
                    aida.tree().cd("fiducial");

                    aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cx, cy);
                    aida.histogram2D("track x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(tx, ty);
                    aida.histogram1D("Track momentum" + topOrBottom, 100, 0., 5.0).fill(p);

                    aida.histogram1D("cluster x - track x", 100, -15., 15.).fill(dx);
                    aida.histogram2D("cluster x - track x vs cluster x", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    if (cy > 0) {
                        aida.histogram1D("cluster energy top", 100, 0.0, 5.0).fill(c.getEnergy());
                        aida.histogram1D("cluster x - track x top", 100, -15., 15.).fill(dx);
                        aida.histogram2D("cluster x - track x vs cluster x top", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    } else {
                        aida.histogram1D("cluster energy bottom", 100, 0.0, 5.0).fill(c.getEnergy());
                        aida.histogram1D("cluster x - track x bottom", 100, -15., 15.).fill(dx);
                        aida.histogram2D("cluster x - track x vs cluster x bottom", 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    }
                    aida.histogram2D("cluster x - track x vs cluster y", 90, -90.0, 90.0, 100, -15., 15.).fill(cy, dx);
                    aida.histogram1D("cluster y - track y", 100, -10., 10.).fill(dy);
                    aida.histogram2D("cluster y - track y vs cluster x", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                    aida.histogram2D("cluster y - track y vs cluster y", 90, -90.0, 90.0, 100, -10., 10.).fill(cy, dy);
                    aida.histogram2D("cluster y - track y vs track y", 90, -90.0, 90.0, 100, -10., 10.).fill(ty, dy);
                    if (cy > 0) {
                        aida.histogram1D("ThetaY top", 100, 0.01, 0.07).fill(thetaY);
                        aida.histogram1D("cluster y - track y top", 100, -10., 10.).fill(dy);
                        aida.histogram2D("cluster y - track y vs cluster x top", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                        aida.histogram2D("cluster y - track y vs thetaY top", 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                    } else {
                        aida.histogram1D("ThetaY bottom", 100, 0.01, 0.07).fill(thetaY);
                        aida.histogram1D("cluster y - track y bottom", 100, -10., 10.).fill(dy);
                        aida.histogram2D("cluster y - track y vs cluster x bottom", 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                        aida.histogram2D("cluster y - track y vs thetaY bottom", 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                    }

                    // by side
                    aida.histogram1D("cluster x - track x" + side, 100, -15., 15.).fill(dx);
                    aida.histogram2D("cluster x - track x vs cluster x" + side, 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    if (cy > 0) {
                        aida.histogram1D("cluster x - track x top" + side, 100, -15., 15.).fill(dx);
                        aida.histogram2D("cluster x - track x vs cluster x top" + side, 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    } else {
                        aida.histogram1D("cluster x - track x bottom" + side, 100, -15., 15.).fill(dx);
                        aida.histogram2D("cluster x - track x vs cluster x bottom" + side, 320, -270.0, 370.0, 100, -15., 15.).fill(cx, dx);
                    }
                    aida.histogram2D("cluster x - track x vs cluster y" + side, 90, -90.0, 90.0, 100, -15., 15.).fill(cy, dx);
                    aida.histogram1D("cluster y - track y" + side, 100, -10., 10.).fill(dy);
                    aida.histogram2D("cluster y - track y vs cluster x" + side, 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                    aida.histogram2D("cluster y - track y vs cluster y" + side, 90, -90.0, 90.0, 100, -10., 10.).fill(cy, dy);
                    if (cy > 0) {
                        aida.histogram1D("ThetaY top" + side, 100, 0.01, 0.07).fill(thetaY);
                        aida.histogram2D("cluster y - track y vs cluster x top" + side, 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                        aida.histogram2D("cluster y - track y vs thetaY top" + side, 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                    } else {
                        aida.histogram1D("ThetaY bottom" + side, 100, 0.01, 0.07).fill(thetaY);
                        aida.histogram2D("cluster y - track y vs cluster x bottom" + side, 320, -270.0, 370.0, 100, -10., 10.).fill(cx, dy);
                        aida.histogram2D("cluster y - track y vs thetaY bottom" + side, 100, 0.01, 0.07, 100, -10., 10.).fill(thetaY, dy);
                    }
                    aida.tree().cd("..");
                }

                if (_analyzeMomentumByCalorimeterCell) {
                    aida.tree().mkdirs("calorimeter cell analysis");
                    aida.tree().cd("calorimeter cell analysis");
                    aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                    aida.histogram1D("Track momentum ix: " + ix + " iy: " + iy, 100, 0., 10.0).fill(p);
                    aida.histogram1D("Track theta ix: " + ix + " iy: " + iy, 100, 0.010, 0.160).fill(theta);
                    aida.histogram2D("Track theta vs p ix: " + ix + " iy: " + iy, 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
//cng                analyzeTrackHits(ix, iy, t);
                    aida.tree().cd("..");
                }
                // have a good FEE candidate here...
                _skipEvent = false;
                if (_requireFiducial && !isFiducial) {
                    _skipEvent = true;
                }
                if (!_skipEvent) {
                    _numberOfEventsSelected++;
                }
            } // end of check on whether RP has ECal cluster
        }// end of check on minimum number of hits on track
        aida.tree().cd("..");
    }

    private void analyzeHitsInFit(ReconstructedParticle rp) {

        Cluster clus = rp.getClusters().get(0);
        if (clus.getEnergy() > _minClusterEnergy) {
            aida.tree().mkdirs("hits in fit analysis");
            aida.tree().cd("hits in fit analysis");
            double clusterTime = ClusterUtilities.findSeedHit(clus).getTime();
            String torb = clus.getPosition()[1] > 0 ? "top" : "bottom";
            Track track = rp.getTracks().get(0);
            double trackTime = ((GenericObject) _trackToData.from(track)).getFloatVal(0);
            aida.histogram1D("Track time - cluster time", 100, -70., -30.).fill(trackTime - clusterTime);
            aida.histogram1D("Track time", 100, -20., 0.).fill(trackTime);
            aida.histogram1D("cluster time " + torb, 100, 30., 50.).fill(clusterTime);
            double correctedTrackTime = timeCorrector.correctTrackTime(track);
            aida.histogram1D("Track time corrected - cluster time", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
            aida.histogram1D("Track time corrected", 100, -20., 0.).fill(correctedTrackTime);
            aida.histogram1D("Track time - track time corrected", 50, -10., -40.).fill(trackTime - correctedTrackTime);

            if (clus.getPosition()[1] > 0) {
                aida.histogram1D("Top cluster energy ", 100, 0.0, 5.0).fill(clus.getEnergy());
                aida.histogram1D("Track time - cluster time top", 100, -70., -30.).fill(trackTime - clusterTime);
                aida.histogram1D("Track time top", 100, -20., 0.).fill(trackTime);
                aida.histogram1D("Track time corrected - cluster time top", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
                aida.histogram1D("Track time corrected top", 100, 0., 100.).fill(correctedTrackTime);
            } else {
                aida.histogram1D("Bottom cluster energy ", 100, 0.0, 5.0).fill(clus.getEnergy());
                aida.histogram1D("Track time - cluster time bottom", 100, -70., -30.).fill(trackTime - clusterTime);
                aida.histogram1D("Track time bottom", 100, -20., 0.).fill(trackTime);
                aida.histogram1D("Track time corrected - cluster time bottom", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
                aida.histogram1D("Track time corrected bottom", 100, 0., 100.).fill(correctedTrackTime);
            }

            for (TrackerHit hit : track.getTrackerHits()) {
                List<RawTrackerHit> rthits = hit.getRawHits();
                RawTrackerHit rawHit = rthits.get(0);
                HpsSiSensor sensor = ((HpsSiSensor) rawHit.getDetectorElement());
                String sensorName = sensor.getName();
                double time = hit.getTime();
                aida.histogram1D(sensorName + " hit time", 100, -30., 20.).fill(time);
                aida.histogram1D(sensorName + " hit time - cluster time", 100, -70., -30.).fill(time - clusterTime);
            }
            aida.tree().cd("..");
        }
    }

    private boolean[] isHoleorSlotTrack(Track t) {
        boolean isHole = false;
        boolean isSlot = false;
        int nHole = 0;
        int nSlot = 0;
        for (TrackerHit hit : t.getTrackerHits()) {
            List<RawTrackerHit> rawHits = hit.getRawHits();
            for (RawTrackerHit rawHit : rawHits) {
                String sensorName = ((HpsSiSensor) rawHit.getDetectorElement()).getName();
                if (sensorName.contains("_hole")) {
                    nHole++;
                }
                if (sensorName.contains("_slot")) {
                    nSlot++;
                }
            }
        }
        // require either all hole or all slot hits
        if (nHole > 0 && nHole * nSlot == 0) {
            isHole = true;
        }
        if (nSlot > 0 && nHole * nSlot == 0) {
            isSlot = true;
        }
        return new boolean[]{isHole, isSlot};
    }

    private void analyzeTrackHits(int ix, int iy, Track t) {
        List<TrackerHit> hitsOnTrack = TrackUtils.getStripHits(t, _hitToStrips, _hitToRotated);
        for (TrackerHit hit : hitsOnTrack) {
            double[] pnt = hit.getPosition();
//            System.out.format("    Hit global position: %10.6f %10.6f %10.6f\n", pnt[0], pnt[1], pnt[2]);
            List<RawTrackerHit> rawHits = hit.getRawHits();
            for (RawTrackerHit rawHit : rawHits) {
                int channel = rawHit.getIdentifierFieldValue("strip");
                HpsSiSensor sensor = (HpsSiSensor) rawHit.getDetectorElement();
                String sensorName = sensor.getName();
                int Layer = sensor.getLayerNumber();
                //accumulate some statistics here...
                if (_channelMaps.containsKey(ix + " " + iy + " " + sensorName)) {
                    _channelMaps.get(ix + " " + iy + " " + sensorName).addValue(channel);
                } else {
                    _channelMaps.put(ix + " " + iy + " " + sensorName, new DescriptiveStatistics());
                    _channelMaps.get(ix + " " + iy + " " + sensorName).addValue(channel);
                }
                aida.histogram1D(ix + " " + iy + " " + sensorName + " channel number", 150, 0., 750.).fill(channel);
//                System.out.format("      Raw hit in layer %d, channel %d\n", Layer, channel);
            }
        }
    }

    @Override
    protected void endOfData() {
//        List<String> sortedNames = new ArrayList<>(_channelMaps.keySet());
//        Collections.sort(sortedNames);
//        System.out.println(" Sensor    Min    Max  Mean  RMS");
//        for (String s : sortedNames) {
//            DescriptiveStatistics dat = _channelMaps.get(s);
//            System.out.format("%52s : %5.1f %5.1f %6.4f %6.4f %n", StringUtils.rightPad(s, 52), dat.getMin(), dat.getMax(), dat.getMean(), dat.getStandardDeviation());
////            System.out.println(s + " : " + dat.getMin()+" "+ dat.getMax() + " " +dat.getMean() + " " + dat.getStandardDeviation());
//        }
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }

    private boolean isTopTrack(Track t) {
        List<TrackerHit> hits = t.getTrackerHits();
        int n[] = {0, 0};
        int nHits = hits.size();
        for (TrackerHit h : hits) {
            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) h.getRawHits().get(0)).getDetectorElement());
            if (sensor.isTopLayer()) {
                n[0] += 1;
            } else {
                n[1] += 1;
            }
        }
        if (n[0] == nHits && n[1] == 0) {
            return true;
        }
        if (n[1] == nHits && n[0] == 0) {
            return false;
        }
        throw new RuntimeException("mixed top and bottom hits on same track");

    }

    private void setupSensors(EventHeader event) {
        List<RawTrackerHit> rawTrackerHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
        EventHeader.LCMetaData meta = event.getMetaData(rawTrackerHits);
        // Get the ID dictionary and field information.
        IIdentifierDictionary dict = meta.getIDDecoder().getSubdetector().getDetectorElement().getIdentifierHelper().getIdentifierDictionary();
        int fieldIdx = dict.getFieldIndex("side");
        int sideIdx = dict.getFieldIndex("strip");
        for (RawTrackerHit hit : rawTrackerHits) {
            // The "side" and "strip" fields needs to be stripped from the ID for sensor lookup.
            IExpandedIdentifier expId = dict.unpack(hit.getIdentifier());
            expId.setValue(fieldIdx, 0);
            expId.setValue(sideIdx, 0);
            IIdentifier strippedId = dict.pack(expId);
            // Find the sensor DetectorElement.
            List<IDetectorElement> des = DetectorElementStore.getInstance().find(strippedId);
            if (des == null || des.size() == 0) {
                throw new RuntimeException("Failed to find any DetectorElements with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            } else if (des.size() == 1) {
                hit.setDetectorElement((SiSensor) des.get(0));
            } else {
                // Use first sensor found, which should work unless there are sensors with duplicate IDs.
                for (IDetectorElement de : des) {
                    if (de instanceof SiSensor) {
                        hit.setDetectorElement((SiSensor) de);
                        break;
                    }
                }
            }
            // No sensor was found.
            if (hit.getDetectorElement() == null) {
                throw new RuntimeException("No sensor was found for hit with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            }
        }
    }

    private RelationalTable getKFTrackDataRelations(EventHeader event) {

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

    public void setAnalyzeHitsInFit(boolean b) {
        _analyzeHitsInFit = b;
    }

    public void setAnalyzeTracksByNhits(boolean b) {
        _analyzeTracksByNhits = b;
    }

    public void setAnalyzeMomentumByCalorimeterCell(boolean b) {
        _analyzeMomentumByCalorimeterCell = b;
    }

    public void setRequireFiducial(boolean b) {
        _requireFiducial = b;
    }

}
