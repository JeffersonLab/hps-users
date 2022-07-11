package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.math.chisq.ChisqProb;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class FeeAnalysis2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private boolean _analyzeTracksByNhits = false;
    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF", "OtherElectrons", "OtherElectrons_KF"};

    private boolean _skimFee = true;
    private boolean _isFeeCandidate;
    private double _minEnergy = 4.0;
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    // use this to select an even number of FEEs over the calorimeter face
    private int _maxClustersPerCrystal = 500;
    private Map<Long, Integer> crystalOccupancyMap = new HashMap<>();

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    protected void process(EventHeader event) {
        boolean skipEvent = false;
        _numberOfEventsProcessed++;
        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        int nClusters = clusters.size();
        aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
        for (Cluster cluster : clusters) {
            aida.histogram2D("All clusters x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
            int ix = seed.getIdentifierFieldValue("ix");
            int iy = seed.getIdentifierFieldValue("iy");
            double e = cluster.getEnergy();
            aida.histogram2D("All clusters ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            if (cluster.getPosition()[1] > 0.) {
                aida.histogram1D("Top cluster energy ", 200, 0.0, 10.0).fill(cluster.getEnergy());
            } else {
                aida.histogram1D("Bottom cluster energy ", 200, 0.0, 10.0).fill(cluster.getEnergy());
            }
            if (e > _minEnergy) {
                _isFeeCandidate = true;
                long cellId = seed.getCellID();
                if (crystalOccupancyMap.containsKey(cellId)) {
//                                System.out.println("found cell "+cellId+" with "+crystalOccupancyMap.get(cellId)+" hits ");
                    crystalOccupancyMap.put(cellId, crystalOccupancyMap.get(cellId) + 1);
                } else {
                    crystalOccupancyMap.put(cellId, 1);
                }
                if (crystalOccupancyMap.get(cellId) > _maxClustersPerCrystal) {
                    _isFeeCandidate = false;
                    skipEvent = true;
                } else {
                    aida.histogram2D("Selected clusters ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                }
            } else {
                _isFeeCandidate = false;
                skipEvent = true;
            }
        }
        if (_isFeeCandidate) {
            setupSensors(event);
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
        if (_skimFee) {
            if (skipEvent) {
                throw new Driver.NextEventException();
            } else {
                _numberOfEventsSelected++;
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
        String trackDir = isGBL ? "gbl" : "htf";
        if (rp.getType() == 1) {
            trackDir = "kf";
        }
        aida.tree().mkdirs(trackDir);
        aida.tree().cd(trackDir);

        Track t = rp.getTracks().get(0);

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
        int nHits = t.getTrackerHits().size();
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
            aida.histogram1D("EoverP" + topOrBottom, 100, 0., 2.).fill(e / p);
            aida.histogram2D("E vs p" + topOrBottom, 100, 2., 7., 100, 2., 7.).fill(e, p);
            aida.tree().mkdirs("calorimeter cell analysis");
            aida.tree().cd("calorimeter cell analysis");
            CalorimeterHit seed = c.getCalorimeterHits().get(0);
            int ix = seed.getIdentifierFieldValue("ix");
            int iy = seed.getIdentifierFieldValue("iy");
            aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            aida.histogram1D("Track momentum ix: " + ix + " iy: " + iy, 100, 0., 10.0).fill(p);
            aida.histogram1D("Track theta ix: " + ix + " iy: " + iy, 100, 0.010, 0.160).fill(theta);
            aida.histogram2D("Track theta vs p ix: " + ix + " iy: " + iy, 100, 0.010, 0.160, 100, 0., 10.0).fill(theta, p);
            aida.histogram2D("track momentum vs cluster x iy: " + iy, 320, -270.0, 370.0, 100, 0., 10.0).fill(c.getPosition()[0], p);
            aida.tree().cd("..");
        }

        aida.tree().cd("..");
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

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }

    public void setSkimFee(boolean b) {
        _skimFee = b;
    }

    public void setMaxClustersPerCrystal(int i) {
        _maxClustersPerCrystal = i;
    }
}
