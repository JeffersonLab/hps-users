package org.hps.users.ngraf.dataanalysis;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class RawTrackerHitAnalysisDriver extends Driver {

    private boolean _cutOnSign = true;
    private int _deltaADCCut = 700;

    protected AIDA aida = AIDA.defaultInstance();

    @Override
    public void process(EventHeader event) {
        if (!event.hasCollection(RawTrackerHit.class, "SVTRawTrackerHits")) {
            // System.out.println(rawHitCollectionName + " does not exist; skipping event");
            return;
        }

        List<RawTrackerHit> rawHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
        if (rawHits == null) {
            throw new RuntimeException("Event is missing SVT hits collection!");
        }

        List<RawTrackerHit> goodRawHits = new ArrayList<>();
        List<RawTrackerHit> badRawHits = new ArrayList<>();

        setupSensors(event);
        // organize lists of hits keyed by strip number into a map keyed by sensor name
        Map<String, Map<Integer, RawTrackerHit>> rawHitMap = new TreeMap<>();
//        System.out.println("found " + rawHits.size() + " raw hits");
        for (RawTrackerHit hit : rawHits) {
            int strip = hit.getIdentifierFieldValue("strip");
            HpsSiSensor sensor = (HpsSiSensor) hit.getDetectorElement();
            String sensorName = sensor.getName();
//            System.out.println(sensor.getName() + " " + sensor.getLayerNumber() + " " + strip);
            //TODO create map with hard-code sensornames
            if (rawHitMap.containsKey(sensorName)) {
                rawHitMap.get(sensorName).put(strip, hit);
            } else {
                rawHitMap.put(sensorName, new TreeMap<Integer, RawTrackerHit>());
                rawHitMap.get(sensorName).put(strip, hit);
            }
        }
        // check
        int nMappedHits = 0;
        for (String s : rawHitMap.keySet()) {
//            System.out.println(s);
            Map<Integer, RawTrackerHit> hits = rawHitMap.get(s);
            // loop through strips
            for (Integer i : hits.keySet()) {
//                System.out.println(i);
                nMappedHits++;
            }
            // check whether isolated strips are "bad"
            boolean isBad = false;
            Iterator<Integer> it = hits.keySet().iterator();
//            System.out.println("input hits size " + hits.size());
            if (hits.size() == 1) {
                Integer firstChannel = it.next();
                isBad = isBadHit(hits.get(firstChannel));
                if (!isBad) {
                    goodRawHits.add(hits.get(firstChannel));
                } else {
                    badRawHits.add(hits.get(firstChannel));
                }
            } else {
                //get the first channel
                Integer firstChannel = it.next();
                int deltaBefore = 0;
                int deltaAfter = 0;
                int deltaChannel = 0;
                // now continue looping over the rest of the strips
                while (it.hasNext()) {
//                    System.out.println("first " + firstChannel);
                    Integer thisChannel = it.next();
//                    System.out.println("next " + thisChannel);
                    // are these channels nearest neighbors?
                    deltaChannel = thisChannel - firstChannel;
//                    System.out.println("delta " + deltaChannel + " deltaBefore " + deltaBefore);
                    if (deltaChannel != 1 && deltaBefore != 1) {
                        // check whether this isolated channel is bad
                        isBad = isBadHit(hits.get(firstChannel));
//                        System.out.println(isBad ? "is bad" : "is good");
                        // only add isolated hits if they are not bad
                        if (!isBad) {
                            goodRawHits.add(hits.get(firstChannel));
                        } else {
                            badRawHits.add(hits.get(firstChannel));
                        }
                    } else // add all non-isolated hits
                    {
                        goodRawHits.add(hits.get(firstChannel));
                    }
                    deltaBefore = deltaChannel;
                    firstChannel = thisChannel;
                }
                // check last channel
                if (deltaChannel != 1) {
                    // check whether this isolated channel is bad
                    isBad = isBadHit(hits.get(firstChannel));
                    // only add isolated hits if they are not bad
                    if (!isBad) {
                        goodRawHits.add(hits.get(firstChannel));
                    } else {
                        badRawHits.add(hits.get(firstChannel));
                    }
                } else // add all non-isolated hits
                {
                    goodRawHits.add(hits.get(firstChannel));
                }
            }
        }
        event.put("SVTRawTrackerHitsGood", goodRawHits, RawTrackerHit.class, 0);
        event.put("SVTRawTrackerHitsBad", badRawHits, RawTrackerHit.class, 0);
//        System.out.println("found " + rawHits.size() + " raw hits");
//        System.out.println("mapped " + nMappedHits + " hits");
//        System.out.println("ended up with " + goodRawHits.size() + " good raw hits");
        for (RawTrackerHit hit : goodRawHits) {
            int strip = hit.getIdentifierFieldValue("strip");
            HpsSiSensor sensor = (HpsSiSensor) hit.getDetectorElement();
            String sensorName = sensor.getName();
//            System.out.println(sensor.getName() + " " + sensor.getLayerNumber() + " " + strip);
        }
        aida.histogram1D("Number of SVTRawTrackerHits", 100, 0., 1000.).fill(rawHits.size());
        aida.histogram1D("Number of good SVTRawTrackerHits", 100, 0., 1000.).fill(goodRawHits.size());
        aida.histogram2D("Number of SVTRawTrackerHits vs number of good SVTRawTrackerHits", 100, 0., 1000., 100, 0., 1000.).fill(rawHits.size(), goodRawHits.size());

        //hits are now collected by sensor and sorted in ascending order by strip
        // loop over hits, identify isolated single-strip hits
        // check if max possible amplitude is below our hit threshold
        // remove so we don't waste time fitting the waveforms for hits we'll never use.
    }

    private boolean isBadHit(RawTrackerHit hit) {
//        System.out.println("checking hit " + hit.getIdentifierFieldValue("strip"));
        short[] adcs = hit.getADCValues();
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
        if (deltaADC < _deltaADCCut) {
            return true;
        }
        if (_cutOnSign) {
            return loIndex < hiIndex ? false : true;
        }
        return false;
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

    public void setDeltaADCCut(int i) {
        _deltaADCCut = i;
    }

    public void setCutOnSign(boolean b) {
        _cutOnSign = b;
    }
}
