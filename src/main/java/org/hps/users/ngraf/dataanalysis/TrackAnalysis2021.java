package org.hps.users.ngraf.dataanalysis;

import static java.lang.Math.abs;
import java.util.List;
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
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class TrackAnalysis2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    public void process(EventHeader event) {
        if (event.hasCollection(ReconstructedParticle.class, "FinalStateParticles_KF")) {
            List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
            for (ReconstructedParticle rp : rps) {
                int pdgId = rp.getParticleIDUsed().getPDG();
                // only consider charged particles
                if (abs(pdgId) == 11) {
                    // only consider tracks with all sensors hit
                    if (rp.getTracks().get(0).getTrackerHits().size() >= 12) {
                        //only consider particles with associated cluster
                        if (!rp.getClusters().isEmpty()) {
                            String type = pdgId == 11 ? " electron" : " positron";

                            Cluster c = rp.getClusters().get(0);
                            Track t = rp.getTracks().get(0);
                            double p = rp.getMomentum().magnitude();
                            double e = rp.getEnergy();
                            if (e > 0.5) {
                                aida.tree().mkdirs(type + " Analysis");
                                aida.tree().cd(type + " Analysis");
                                CalorimeterHit seed = c.getCalorimeterHits().get(0);
                                int ix = seed.getIdentifierFieldValue("ix");
                                int iy = seed.getIdentifierFieldValue("iy");
                                String fid = TriggerModule.inFiducialRegion(c) ? " fiducial " : " non fid ";
                                String topOrBottom = iy > 0 ? " top " : " bottom ";
                                aida.histogram2D("cluster ix vs iy" + type, 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
//                                aida.histogram1D("EoverP ix: " + ix + " iy: " + iy + type, 100, 0., 2.0).fill(e / p);
                                aida.histogram1D("Track momentum" + topOrBottom + type, 200, 0., 7.0).fill(p);
                                aida.histogram2D("E vs p" + topOrBottom + type, 100, 0., 7., 100, 0., 7.).fill(e, p);
                                aida.histogram2D("E vs p" + topOrBottom + fid + type, 100, 0., 7., 100, 0., 7.).fill(e, p);
                                aida.histogram1D("EoverP" + topOrBottom + fid + type, 100, 0., 2.0).fill(e / p);

                                setupSensors(event);
                                List<TrackerHit> hits = t.getTrackerHits();
                                for (TrackerHit hit : hits) {
                                    HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();
                                    String sensorName = sensor.getName();
                                    int layer = sensor.getLayerNumber();
                                    String sensorType = layer < 3 ? "thin" : "thick";
                                    int clusSize = hit.getRawHits().size();
                                    double dEdx = 1000000 * hit.getdEdx();
                                    double time = hit.getTime();
                                    aida.histogram1D("hit dEdx " + sensorType, 100, 0., 20.).fill(dEdx);
                                    aida.histogram1D("hit time " + sensorType, 100, -20., 20.).fill(time);
//                                    aida.histogram1D("hit cluster size " + sensorName, 10, 0.5, 10.5).fill(clusSize);
                                    if (type.equals(" electron") && !sensorName.contains("slot")) {
                                        aida.histogram1D(sensorName + " hit dEdx", 100, 0., 20.).fill(dEdx);
                                        if (clusSize < 3) {
                                            aida.histogram1D(sensorName + " hit dEdx " + clusSize + " strip cluster", 100, 0., 20.).fill(dEdx);
                                        }
                                        aida.histogram1D(sensorName + " cluster size", 10, -0.5, 9.5).fill(clusSize);
                                    }
                                    if (type.equals(" positron") && !sensorName.contains("hole")) {
                                        aida.histogram1D(sensorName + " hit dEdx", 100, 0., 20.).fill(dEdx);
                                        if (clusSize < 3) {
                                            aida.histogram1D(sensorName + " hit dEdx " + clusSize + " strip cluster", 100, 0., 20.).fill(dEdx);
                                        }
                                        aida.histogram1D(sensorName + " cluster size", 10, -0.5, 9.5).fill(clusSize);
                                    }
                                    // let's look at some individual strips, so select single strip clusters.
                                    if (clusSize == 1) {
                                        List<RawTrackerHit> clusterRawHits = hit.getRawHits();
                                        int channel = clusterRawHits.get(0).getIdentifierFieldValue("strip");
                                        aida.histogram1D(sensor.getName() + " channel", 640, 0, 640.).fill(channel);
                                        //let's start with the split-strip sensors...
                                        if (sensorName.contains("L1") || sensorName.contains("L2")) {
                                            if (channel < 20 || channel > 490) {
                                                aida.tree().mkdirs(sensor.getName());
                                                aida.tree().cd(sensor.getName());
                                                aida.histogram1D(sensorName + " hit dEdx channel "+channel, 100, 0., 20.).fill(dEdx);
                                                aida.tree().cd("..");
                                            }
                                        }
                                    }
//                                    aida.histogram1D(sensorName + " hit time", 100, -20., 20.).fill(time);
//                                    aida.histogram1D(sensorName + " hit cluster size", 10, -0.5, 9.5).fill(clusSize);
//                                    aida.histogram2D(sensorName + " dEdx vs time", 100, 0., 20., 100, -20., 20.).fill(dEdx, time);
                                }
                                aida.tree().cd("..");
                            }
                        }
                    }
                }
            }
        }
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
}
