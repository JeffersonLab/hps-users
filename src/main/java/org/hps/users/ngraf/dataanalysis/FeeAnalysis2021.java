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

    private boolean _analyzeTracksByNhits = false;
    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF"};//, "OtherElectrons", "OtherElectrons_KF"};

    RelationalTable _hitToStrips;
    RelationalTable _hitToRotated;
    RelationalTable _trackToData;

    Map<String, DescriptiveStatistics> _channelMaps = new HashMap<>();

    private SvtSensorTimeCorrector timeCorrector;

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
        timeCorrector = new SvtSensorTimeCorrector();
    }

    protected void process(EventHeader event) {
        _hitToStrips = TrackUtils.getHitToStripsTable(event);
        _hitToRotated = TrackUtils.getHitToRotatedTable(event);

        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        int nClusters = clusters.size();
        aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
        for (Cluster cluster : clusters) {
            aida.histogram2D("Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            if (cluster.getPosition()[1] > 0.) {
                aida.histogram1D("Top cluster energy ", 200, 0.0, 10.0).fill(cluster.getEnergy());
            } else {
                aida.histogram1D("Bottom cluster energy ", 200, 0.0, 10.0).fill(cluster.getEnergy());
            }
        }
        setupSensors(event);
        _trackToData = getKFTrackDataRelations(event);

        // let's look at ReconstructedParticles...
        for (String rpCollectionName : ReconstructedParticleCollectionNames) {
            if (event.hasCollection(ReconstructedParticle.class, rpCollectionName)) {
                aida.tree().mkdirs(rpCollectionName);
                aida.tree().cd(rpCollectionName);
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, rpCollectionName);
                int nTracksWithClusters = 0;
                int nTracks = 0;
                for (ReconstructedParticle rp : rpList) {
                    if (abs(rp.getParticleIDUsed().getPDG()) == 11) {
                        nTracks++;
                        if (rp.getClusters().size() == 1) {
                            nTracksWithClusters++;
                        }
                    }
                    analyzeReconstructedParticle(rp);
                }
                aida.histogram1D("number of tracks", 5, -0.5, 4.5).fill(nTracks);
                aida.histogram1D("number of tracks with clusters", 5, -0.5, 4.5).fill(nTracksWithClusters);
                aida.histogram2D("number of clusters vs number of tracks with clusters", 5, -0.5, 4.5, 5, -0.5, 4.5).fill(nClusters, nTracksWithClusters);
                aida.tree().cd("..");
            }
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
        if (rp.getClusters().size() == 1) {
            Cluster c = rp.getClusters().get(0);
            double e = rp.getEnergy();
            aida.histogram2D(type + "single Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(c.getPosition()[0], c.getPosition()[1]);
            if (c.getPosition()[1] > 0.) {
                aida.histogram1D(type + "single Top cluster energy ", 200, 0.0, 10.0).fill(c.getEnergy());
            } else {
                aida.histogram1D(type + "single Bottom cluster energy ", 200, 0.0, 10.0).fill(c.getEnergy());
            }
        }
    }

    void analyzeTrack(ReconstructedParticle rp) {
        boolean isGBL = TrackType.isGBL(rp.getType());
        int minNhits = isGBL ? 6 : 12;
        String trackDir = isGBL ? "gbl" : "htf";
        if (rp.getType() == 1) {
            trackDir = "kf";
        }
        aida.tree().mkdirs(trackDir);
        aida.tree().cd(trackDir);

        Track t = rp.getTracks().get(0);
        int nHits = t.getTrackerHits().size();
        if (nHits >= minNhits) {
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
            String topOrBottom = isTopTrack(t) ? " top " : " bottom ";
            aida.histogram1D("Track chisq per df" + topOrBottom, 100, 0., 50.).fill(chiSquared / ndf);
            aida.histogram1D("Track chisq prob" + topOrBottom, 100, 0., 1.).fill(chisqProb);
            aida.histogram1D("Track nHits" + topOrBottom, 14, 0.5, 14.5).fill(t.getTrackerHits().size());
            aida.histogram1D("Track momentum" + topOrBottom, 200, 0., 10.0).fill(p);
            aida.histogram1D("Track deDx" + topOrBottom, 100, 0.00004, 0.00013).fill(t.getdEdx());
            aida.histogram1D("Track theta" + topOrBottom, 100, 0.010, 0.160).fill(theta);
            aida.histogram2D("Track theta vs p" + topOrBottom, 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
            aida.histogram1D("rp x0" + topOrBottom, 100, -0.50, 0.50).fill(TrackUtils.getX0(t));
            aida.histogram1D("rp y0" + topOrBottom, 100, -5.0, 5.0).fill(TrackUtils.getY0(t));
            aida.histogram1D("rp z0" + topOrBottom, 100, -1.0, 1.0).fill(TrackUtils.getZ0(t));

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
                if (rp.getType() == 1) {
                    analyzeHitsInFit(rp);
                }
                aida.histogram1D("EoverP" + topOrBottom, 100, 0., 2.).fill(e / p);
                aida.histogram2D("E vs p" + topOrBottom, 100, 2., 7., 100, 2., 7.).fill(e, p);
                aida.tree().mkdirs("calorimeter cell analysis");
                aida.tree().cd("calorimeter cell analysis");
                CalorimeterHit seed = c.getCalorimeterHits().get(0);
                int ix = seed.getIdentifierFieldValue("ix");
                int iy = seed.getIdentifierFieldValue("iy");
                aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                aida.histogram1D("Track momentum ix: " + ix + " iy: " + iy, 100, 0., 10.0).fill(p);
//                aida.histogram1D("Track theta ix: " + ix + " iy: " + iy, 100, 0.010, 0.160).fill(theta);
//                aida.histogram2D("Track theta vs p ix: " + ix + " iy: " + iy, 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
//                aida.histogram2D("track momentum vs cluster x iy: " + iy, 320, -270.0, 370.0, 100, 0., 10.0).fill(c.getPosition()[0], p);
//cng                analyzeTrackHits(ix, iy, t);
                aida.tree().cd("..");
            } // end of check on whether RP has ECal cluster
        }// end of check on minimum number of hits on track
        aida.tree().cd("..");
    }

    private void analyzeHitsInFit(ReconstructedParticle rp) {

        Cluster clus = rp.getClusters().get(0);
        if (clus.getEnergy() > 2.5) {
            aida.tree().mkdirs("hits in fit analysis");
            aida.tree().cd("hits in fit analysis");
            double clusterTime = ClusterUtilities.findSeedHit(clus).getTime();
            String torb = clus.getPosition()[1] > 0 ? "top" : "bottom";
            Track track = rp.getTracks().get(0);
            double trackTime = ((GenericObject) _trackToData.from(track)).getFloatVal(0);
            aida.histogram1D("Track time - cluster time", 100, -70., -30.).fill(trackTime - clusterTime);
            aida.histogram1D("Track time",100, -20., 0.).fill(trackTime);
            aida.histogram1D("cluster time " + torb, 100, 30., 50.).fill(clusterTime);
            double correctedTrackTime = timeCorrector.correctTrackTime(track);
            aida.histogram1D("Track time corrected - cluster time", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
            aida.histogram1D("Track time corrected",100, -20., 0.).fill(correctedTrackTime);
            aida.cloud1D("Track time - track time corrected").fill(trackTime - correctedTrackTime);

            if (clus.getPosition()[1] > 0) {
                aida.histogram1D("Top cluster energy ", 200, 0.0, 10.0).fill(clus.getEnergy());
                aida.histogram1D("Track time - cluster time top", 100, -70., -30.).fill(trackTime - clusterTime);
                aida.histogram1D("Track time top",100, -20., 0.).fill(trackTime);
                aida.histogram1D("Track time corrected - cluster time top", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
                aida.cloud1D("Track time corrected top").fill(correctedTrackTime);
            } else {
                aida.histogram1D("Bottom cluster energy ", 200, 0.0, 10.0).fill(clus.getEnergy());
                aida.histogram1D("Track time - cluster time bottom", 100, -70., -30.).fill(trackTime - clusterTime);
                aida.histogram1D("Track time bottom", 100, -20., 0.).fill(trackTime);
                aida.histogram1D("Track time corrected - cluster time bottom", 100, -70., -30.).fill(correctedTrackTime - clusterTime);
                aida.cloud1D("Track time corrected bottom").fill(correctedTrackTime);
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
        List<String> sortedNames = new ArrayList<>(_channelMaps.keySet());
        Collections.sort(sortedNames);
        System.out.println(" Sensor    Min    Max  Mean  RMS");
        for (String s : sortedNames) {
            DescriptiveStatistics dat = _channelMaps.get(s);
            System.out.format("%52s : %5.1f %5.1f %6.4f %6.4f %n", StringUtils.rightPad(s, 52), dat.getMin(), dat.getMax(), dat.getMean(), dat.getStandardDeviation());
//            System.out.println(s + " : " + dat.getMin()+" "+ dat.getMax() + " " +dat.getMean() + " " + dat.getStandardDeviation());
        }
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

}
