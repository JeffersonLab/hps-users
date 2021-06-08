package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.hps.recon.tracking.FittedRawTrackerHit;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
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
public class KFTrackAnalysis extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    String[] collectionNames = {"FinalStateParticles_KF", "OtherElectrons_KF"};

    boolean _analyzeRawHits = true;
    boolean _analyzeBySensor = false;

    protected void process(EventHeader event) {
        setupSensors(event);
        List<LCRelation> fittedHits = event.get(LCRelation.class, "SVTFittedRawTrackerHits");
        List<RawTrackerHit> rawHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
//        if (_analyzeRawHits) {
        // Map the fitted hits to their corresponding raw hits
        Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();
        for (LCRelation fittedHit : fittedHits) {
            fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
        }
//        }
        for (String RPCollection : collectionNames) {
            if (event.hasCollection(ReconstructedParticle.class, RPCollection)) {
                aida.tree().mkdirs(RPCollection);
                aida.tree().cd(RPCollection);
                List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, RPCollection);
                int nRPs = 0;
                for (ReconstructedParticle rp : rpList) {
                    if (abs(rp.getParticleIDUsed().getPDG()) == 11) {
                        nRPs++;
                        Cluster ecalClus = null;
                        if (!rp.getClusters().isEmpty()) {
                            ecalClus = rp.getClusters().get(0);
                        }
                        Track t = rp.getTracks().get(0);
                        int nHitsOnTrack = t.getTrackerHits().size();

                        int pdgId = rp.getParticleIDUsed().getPDG();
                        Hep3Vector pmom = rp.getMomentum();
//                        double thetaY = asin(pmom.y() / pmom.magnitude());
//                        double z0 = t.getTrackStates().get(0).getZ0();
                        String torb = isTopTrack(t) ? "top " : "bottom ";
                        aida.histogram1D(torb + pdgId + " track nHits", 15, -0.5, 14.5).fill(nHitsOnTrack);
                        aida.histogram1D(torb + pdgId + " track momentum", 100, 0., 7.).fill(rp.getMomentum().magnitude());
//                        aida.cloud1D(torb + pdgId + "|thetaY|").fill(abs(thetaY));
//                        aida.histogram1D(torb + pdgId + "z0", 100, -2., 2.).fill(z0);
//                        aida.cloud2D(torb + pdgId + "|thetaY| vs z0").fill(abs(thetaY), z0);
//                        aida.profile1D(torb + pdgId + "|thetaY| vs z0 profile", 10, 0.01, 0.1).fill(abs(thetaY), z0);
                        aida.tree().mkdirs(nHitsOnTrack + " hit tracks");
                        aida.tree().cd(nHitsOnTrack + " hit tracks");
                        aida.histogram1D(torb + pdgId + " " + nHitsOnTrack + " hit track momentum", 100, 0., 7.).fill(rp.getMomentum().magnitude());
                        if (ecalClus != null) {
                            aida.histogram1D(torb + pdgId + " " + nHitsOnTrack + " hit track cluster energy", 100, 0., 7.).fill(ecalClus.getEnergy());
                            aida.histogram1D(torb + pdgId + " " + nHitsOnTrack + " hit track e over p", 100, 0., 7.).fill(ecalClus.getEnergy() / rp.getMomentum().magnitude());
                        }
//                        aida.cloud1D(torb + pdgId + " " + nHitsOnTrack + " hit |thetaY|").fill(abs(thetaY));
//                        aida.histogram1D(torb + pdgId + " " + nHitsOnTrack + " hit z0", 100, -2., 2.).fill(z0);
//                        aida.cloud2D(torb + pdgId + " " + nHitsOnTrack + " hit |thetaY| vs z0").fill(abs(thetaY), z0);
//                        aida.profile1D(torb + pdgId + " " + nHitsOnTrack + " hit |thetaY| vs z0 profile", 10, 0.01, 0.1).fill(abs(thetaY), z0);
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
                            aida.histogram1D("hit cluster size " + sensorName, 10, 0.5, 10.5).fill(clusSize);
                            if (_analyzeBySensor) {
                                aida.histogram1D(sensorName + " hit dEdx", 100, 0., 20.).fill(dEdx);
                                aida.histogram1D(sensorName + " hit time", 100, -20., 20.).fill(time);
                                aida.histogram1D(sensorName + " hit cluster size", 10, -0.5, 9.5).fill(clusSize);
                                aida.histogram2D(sensorName + " dEdx vs time", 100, 0., 20., 100, -20., 20.).fill(dEdx, time);
                                if (sensorName.contains("L1") || sensor.getName().contains("L1")) {
                                    if (clusSize < 3) {
                                        aida.histogram1D(sensorName + " " + clusSize + " hit dEdx", 100, 0., 20.).fill(dEdx);
                                        aida.histogram1D(sensorName + " " + clusSize + " hit time", 100, -20., 20.).fill(time);
                                        aida.histogram2D(sensorName + " " + clusSize + " dEdx vs time", 100, 0., 20., 100, -20., 20.).fill(dEdx, time);
                                    }
                                }
                            }
                            // look at the raw hits
                            List<RawTrackerHit> clusterRawHits = hit.getRawHits();
                            if (clusSize == 1) {
//                                if (1000000 * hit.getdEdx() < 2 || 1000000 * hit.getdEdx() > 3) {
//                                String lowHigh = 1000000 * hit.getdEdx() < 2 ? " low" : " high";
                                for (RawTrackerHit rth : clusterRawHits) {
                                    int channel = rth.getIdentifierFieldValue("strip");
                                    short[] adcs = rth.getADCValues();
                                    short lo = Short.MAX_VALUE;
                                    int loIndex = -1;
                                    short hi = Short.MIN_VALUE;
                                    int hiIndex = -1;
                                    for (int i = 0; i < 6; ++i) {
                                        if (adcs[i] < lo) {
                                            lo = adcs[i];
                                            loIndex = i;
                                        }
                                        if (adcs[i] > hi) {
                                            hi = adcs[i];
                                            hiIndex = i;
                                        }
//                                            aida.histogram1D(rth.getCellID() + lowHigh, 6, -0.5, 5.5).fill(i, adcs[i]);
                                        //TODO  figure out how to load the pedestals
//                                            System.out.println(sensorName + " " + i + " pedestal " + sensor.getPedestal(channel, i));
                                    }
                                    short deltaADC = (short) (hi - lo);
                                    int sign = loIndex < hiIndex ? 1 : -1;
                                    // get raw amplitude
                                    double amplitude = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rth));
                                    aida.histogram2D("amplitude vs dEdx " + sensorType, 250, 0., 5000., 100, 0., 20.).fill(amplitude, dEdx);
                                    aida.histogram1D(sensor.getName() + " channel", 640, 0, 640.).fill(channel);
                                    aida.histogram1D("deltaADC" + sensorType, 150, 200., 3000.).fill(sign * deltaADC);
                                    aida.histogram2D("deltaADC vs dEdx " + sensorType, 150, 200., 3000., 100, 0., 20.).fill(sign * deltaADC, dEdx);
                                    aida.histogram2D("deltaADC vs time " + sensorType, 150, 200., 3000., 100, -20., 20.).fill(sign * deltaADC, time);
                                    aida.histogram2D("time vs dEdx " + sensorType, 100, -20., 20., 100, 0., 20.).fill(time, dEdx);
                                    aida.histogram2D("deltaADC vs time wide " + sensorType, 150, -5000., 7000., 100, -250., 100.).fill(sign * deltaADC, time);
                                    aida.histogram2D("time wide vs dEdx " + sensorType, 100, -250., 100., 100, 0., 20.).fill(time, dEdx);
                                    int deltaADCCut = 700;
                                    if (sign * deltaADC > deltaADCCut) {
                                        for (int i = 0; i < 6; ++i) {
                                            aida.histogram2D("APV25 Samples " + sensorName, 6, 0., 6., 100, 4000., 8000.).fill(i, adcs[i]);
                                        }
                                        // quick check on split-strip sensors...(and one other just to compare)
                                        if (sensor.getLayerNumber() <= 6) {
                                            for (int i = 0; i < 6; ++i) {
                                                if (channel < 300) {
                                                    aida.histogram2D("APV25 Samples " + sensorName + " channel < 300", 6, 0., 6., 100, 4000., 8000.).fill(i, adcs[i]);
                                                } else {
                                                    aida.histogram2D("APV25 Samples " + sensorName + " channel > 300", 6, 0., 6., 100, 4000., 8000.).fill(i, adcs[i]);
                                                }
                                            }
                                        }
                                        aida.histogram1D("1 hit dEdx deltaADC> " + deltaADCCut + " " + sensorType, 100, 0., 20.).fill(dEdx);
                                        aida.histogram1D("1 hit time deltaADC> " + deltaADCCut + " " + sensorType, 100, -20., 20.).fill(time);
                                    } else {
                                        aida.histogram1D("1 hit dEdx deltaADC< " + deltaADCCut + " " + sensorType, 100, 0., 20.).fill(dEdx);
                                        aida.histogram1D("1 hit time deltaADC< " + deltaADCCut + " " + sensorType, 100, -20., 20.).fill(time);
                                    }
                                    aida.histogram1D("1 hit dEdx " + sensorType, 100, 0., 20.).fill(dEdx);
                                    aida.histogram1D("1 hit time " + sensorType, 100, -20., 20.).fill(time);
                                }
//                                }
                            } else if (clusSize == 2) {
                                RawTrackerHit s1 = clusterRawHits.get(0);
                                RawTrackerHit s2 = clusterRawHits.get(1);

                                double a1 = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(s1));
                                double a2 = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(s2));
                                double t1 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(s1));
                                double t2 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(s2));

                                aida.cloud1D(sensorName + " 2 strip cluster amplitude 1").fill(a1);
                                aida.cloud1D(sensorName + " 2 strip cluster amplitude 2").fill(a2);
                                aida.cloud2D(sensorName + " 2 strip cluster amplitude 1 vs amplitude 2").fill(a1, a2);
                                aida.cloud1D(sensorName + " 2 strip cluster time 1").fill(t1);
                                aida.cloud1D(sensorName + " 2 strip cluster time 2").fill(t2);
                                aida.cloud1D(sensorName + " 2 strip cluster delta time").fill(t1 - t2);
                                aida.cloud2D(sensorName + " 2 strip cluster time 1 vs time 2").fill(t1, t2);
                                if (abs(t1 - t2) > 8) {
                                    System.out.println(event.getRunNumber() + " " + event.getEventNumber() + " " + sensorName + " out of time");
                                    System.out.println("hit time: " + hit.getTime());
                                    System.out.println("strip 1: " + s1.getIdentifierFieldValue("strip") + " " + a1 + " " + t1);
                                    System.out.println(" FittedRawTrackerHit.getT0 " + t1 + " RawTrackerHit.getTime() " + s1.getTime());
                                    System.out.println("strip 2: " + s2.getIdentifierFieldValue("strip") + " " + a2 + " " + t2);
                                    System.out.println(" FittedRawTrackerHit.getT0 " + t2 + " RawTrackerHit.getTime() " + s2.getTime());

                                }
                            }
                        }
                        aida.tree().cd("..");
                    } // end of check on electron/positron
                }//end of loop over ReconstructedParticle collections
                aida.histogram1D("number of charged ReconstructedParticles", 10, 0., 10.).fill(nRPs);
                aida.tree().cd("..");
            }  //end of check on existence of ReconstructedParticle collections
        } //end of loop over collections
        if (_analyzeRawHits) {
            HashMap<String, Integer> sensorCountMap = new HashMap<>();
            aida.tree().mkdirs("SVTRawTrackerHits");
            aida.tree().cd("SVTRawTrackerHits");
//            List<LCRelation> fittedHits = event.get(LCRelation.class, "SVTFittedRawTrackerHits");
//            List<RawTrackerHit> rawHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
            // Map the fitted hits to their corresponding raw hits
//            Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();
//            for (LCRelation fittedHit : fittedHits) {
//                fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
//            }
            for (RawTrackerHit rawHit : rawHits) {
                // Access the sensor associated with the raw hit
                HpsSiSensor sensor = (HpsSiSensor) rawHit.getDetectorElement();
                String sensorName = sensor.getName();
                if (sensorCountMap.containsKey(sensorName)) {
                    sensorCountMap.put(sensorName, sensorCountMap.get(sensorName) + 1);
                } else {
                    sensorCountMap.put(sensorName, 1);
                }
                int layer = sensor.getLayerNumber();
                String sensorType = layer < 3 ? "thin" : "thick";
                // Get the hit amplitude
                double amplitude = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rawHit));

                aida.histogram1D(sensor.getName() + " amplitude", 250, 0., 5000.).fill(amplitude);
                // Get the t0 of the hit
                double t0 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(rawHit));
                aida.histogram1D(sensor.getName() + " t0", 250, -400., 100.).fill(t0);

                // Retrieve the channel ID of the raw hit
                int channel = rawHit.getIdentifierFieldValue("strip");
                aida.histogram1D(sensor.getName() + " channel", 640, 0, 640.).fill(channel);
                short[] adcs = rawHit.getADCValues();
                short lo = Short.MAX_VALUE;
                int loIndex = -1;
                short hi = Short.MIN_VALUE;
                int hiIndex = -1;
                for (int i = 0; i < 6; ++i) {
                    if (adcs[i] < lo) {
                        lo = adcs[i];
                        loIndex = i;
                    }
                    if (adcs[i] > hi) {
                        hi = adcs[i];
                        hiIndex = i;
                    }
//                                            aida.histogram1D(rth.getCellID() + lowHigh, 6, -0.5, 5.5).fill(i, adcs[i]);
                    //TODO  figure out how to load the pedestals
//                                            System.out.println(sensorName + " " + i + " pedestal " + sensor.getPedestal(channel, i));
                }
                short deltaADC = (short) (hi - lo);
                int sign = loIndex < hiIndex ? 1 : -1;
                aida.histogram1D("deltaADC " + sensorType, 150, -5000., 7000.).fill(sign * deltaADC);
                aida.histogram2D("deltaADC vs amplitude " + sensorType, 150, -5000., 7000., 250, 0., 5000.).fill(sign * deltaADC, amplitude);
                aida.histogram2D("deltaADC vs time " + sensorType, 150, -5000., 7000., 250, -400., 100.).fill(sign * deltaADC, t0);
                aida.histogram2D("time vs amplitude " + sensorType, 250, -400., 100., 250, 0., 5000.).fill(t0, amplitude);

                aida.histogram1D("amplitude " + sensorType, 250, 0., 5000.).fill(amplitude);
                aida.histogram1D("t0 " + sensorType, 250, -400., 100.).fill(t0);
                int deltaADCCut = 700;
                if (sign * deltaADC > deltaADCCut) {
                    aida.histogram1D("amplitude deltaADC> " + deltaADCCut + " " + sensorType, 250, 0., 5000.).fill(amplitude);
                    aida.histogram1D("t0 deltaADC> " + deltaADCCut + " " + sensorType, 250, -400., 100.).fill(t0);
                } else {
                    aida.histogram1D("amplitude deltaADC< " + deltaADCCut + " " + sensorType, 250, 0., 5000.).fill(amplitude);
                    aida.histogram1D("t0 deltaADC< " + deltaADCCut + " " + sensorType, 250, -400., 100.).fill(t0);
                }
            }
            for (Map.Entry entry : sensorCountMap.entrySet()) {
                aida.histogram1D(entry.getKey() + " hits per event", 500, 0., 500.).fill((int) entry.getValue());
            }
            aida.tree().cd("..");
        }
    }

    public void setAnalyzeRawHits(boolean b) {
        _analyzeRawHits = b;
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
}
