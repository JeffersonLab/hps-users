package org.hps.users.ngraf.dataanalysis;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.StatUtils;
import org.hps.recon.ecal.cluster.ClusterUtilities;
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
import org.lcsim.geometry.Detector;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author ngraf
 */
public class SvtHitAnalysisDriver extends Driver {

    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree;
    IHistogramFactory histogramFactory;

    //List of Sensors
    private List<HpsSiSensor> sensors = null;

    Map<String, IHistogram1D> nRawHits = new HashMap<String, IHistogram1D>();
//    Map<String, IHistogram2D> adc = new HashMap<String, IHistogram2D>();
//    Map<String, IHistogram2D> adcMonster = new HashMap<String, IHistogram2D>();
//    Map<String, IHistogram1D> channelMonster = new HashMap<String, IHistogram1D>();
//    Map<String, IHistogram1D> channelAboveThresh = new HashMap<String, IHistogram1D>();
//    Map<String, IHistogram1D> channelBelowThresh = new HashMap<String, IHistogram1D>();

    //Histogram Settings
    double minX = 0;
    double maxX = 10;
    double minADC = -1000;
    double maxADC = 5000;
    int nBins = 120;
    int nHitMax = 100;

    //Collection Strings
    private String fittedHitsCollectionName = "SVTFittedRawTrackerHits";
    String rawTrackerHitCollectionName = "SVTRawTrackerHits";

    private static final String SUBDETECTOR_NAME = "Tracker";

    public void detectorChanged(Detector detector) {

        aida.tree().cd("/");
        tree = aida.tree();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

        //Set Beam Energy
//        BeamEnergy.BeamEnergyCollection beamEnergyCollection
//                = this.getConditionsManager().getCachedConditions(BeamEnergy.BeamEnergyCollection.class, "beam_energies").getCachedData();
        // Get the HpsSiSensor objects from the tracker detector element
        sensors = detector.getSubdetector(SUBDETECTOR_NAME)
                .getDetectorElement().findDescendants(HpsSiSensor.class);

        // If the detector element had no sensors associated with it, throw
        // an exception
        if (sensors.size() == 0) {
            throw new RuntimeException("No sensors were found in this detector.");
        }

        for (HpsSiSensor sensor : sensors) {
            nRawHits.put(sensor.getName(), histogramFactory.createHistogram1D("Number of SVT Raw Hits " + sensor.getName(), 640, 0, 640));
//            adc.put(sensor.getName(), histogramFactory.createHistogram2D("ADC Counts " + sensor.getName(), 6, 0, 6, nBins, minADC, maxADC));
//            adcMonster.put(sensor.getName(), histogramFactory.createHistogram2D("ADC Counts Monster Events " + sensor.getName(), 6, 0, 6, nBins, minADC, maxADC));
//            channelMonster.put(sensor.getName(), histogramFactory.createHistogram1D("Channels Hit Monster Events " + sensor.getName(), 640, 0, 640));
//            channelAboveThresh.put(sensor.getName(), histogramFactory.createHistogram1D("Channels Hit Above Threshold " + sensor.getName(), 640, 0, 640));
//            channelBelowThresh.put(sensor.getName(), histogramFactory.createHistogram1D("Channels Hit Below Threshold " + sensor.getName(), 640, 0, 640));
        }

    }

    public void process(EventHeader event) {
        aida.tree().cd("/");
        setupSensors(event);
        // Get the list of fitted hits from the event
        List<LCRelation> fittedHits = event.get(LCRelation.class, fittedHitsCollectionName);

        // Map the fitted hits to their corresponding raw hits
        Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();

        List<RawTrackerHit> rawHits = event.get(RawTrackerHit.class, rawTrackerHitCollectionName);

        aida.histogram1D("Number of RawTrackerHits in Event", 100, 0., 2000.).fill(rawHits.size());
        aida.histogram1D("Number of FittedRawTrackerHits in Event", 100, 0., 2000.).fill(fittedHits.size());
        
        /*for (RawTrackerHit rawHit : rawHits) {
             
            // Access the sensor associated with the raw hit
            HpsSiSensor sensor = (HpsSiSensor) rawHit.getDetectorElement();
             
            // Retrieve the channel ID of the raw hit
            int channel = rawHit.getIdentifierFieldValue("strip");
        } */
        for (LCRelation fittedHit : fittedHits) {
            fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
        }
        analyzeFittedRawTrackerHits(fittedHits);

        List<SiTrackerHitStrip1D> stripHits = event.get(SiTrackerHitStrip1D.class, "StripClusterer_SiTrackerHitStrip1D");
        analyzeSiTrackerHitStrip1D(stripHits);

        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
        analyzeHitsOnTracks(rpList);

        Map<String, Integer> hits = new HashMap<String, Integer>();
        for (HpsSiSensor sensor : sensors) {
            hits.put(sensor.getName(), 0);
        }

        int hitNumber = 0;
        for (RawTrackerHit rawHit : rawHits) {
            // Access the sensor associated with the raw hit
            HpsSiSensor sensor = (HpsSiSensor) rawHit.getDetectorElement();
            Integer nHits = hits.get(sensor.getName());
            if (nHits == null) {
                nHits = 0;
            }
            nHits++;
            hits.put(rawHit.getDetectorElement().getName(), nHits);

            // does this RawTrackerHit have a corresponding pulse fit?
            if (fittedRawTrackerHitMap.get(rawHit) != null) {
                // Get the hit amplitude
                double amplitude = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rawHit));

                aida.histogram1D(sensor.getName() + " amplitude", 250, 0., 5000.).fill(amplitude);
                // Get the t0 of the hit
                double t0 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(rawHit));
                aida.histogram1D(sensor.getName() + " t0", 250, -400., 100.).fill(t0);

                // Retrieve the channel ID of the raw hit
                int channel = rawHit.getIdentifierFieldValue("strip");
                aida.histogram1D(sensor.getName() + " channel", 640, 0, 640.).fill(channel);
                // Unfortunatley all the sensors and channels appear to be jumbled in the LCIO file.
//            System.out.println(hitNumber + " " + sensor.getName() + " " + channel);
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
                aida.histogram1D("deltaADC", 150, -5000., 7000.).fill(sign * deltaADC);
                aida.histogram2D("deltaADC vs amplitude", 150, -5000., 7000., 250, 0., 5000.).fill(sign * deltaADC, amplitude);
                aida.histogram2D("deltaADC vs time", 150, -5000., 7000., 250, -400., 100.).fill(sign * deltaADC, t0);
                aida.histogram2D("time vs amplitude", 250, -400., 100., 250, 0., 5000.).fill(t0, amplitude);

                aida.histogram1D("amplitude", 250, 0., 5000.).fill(amplitude);
                aida.histogram1D("t0", 250, -400., 100.).fill(t0);
                int deltaADCCut = 700;
                if (sign * deltaADC > deltaADCCut) {
                    aida.histogram1D("amplitude deltaADC> " + deltaADCCut, 250, 0., 5000.).fill(amplitude);
                    aida.histogram1D("t0 deltaADC> " + deltaADCCut, 250, -400., 100.).fill(t0);
                } else {
                    aida.histogram1D("amplitude deltaADC< " + deltaADCCut, 250, 0., 5000.).fill(amplitude);
                    aida.histogram1D("t0 deltaADC< " + deltaADCCut, 250, -400., 100.).fill(t0);
                }
                hitNumber++;
//            short[] adcs = rawHit.getADCValues();
//
//            double[] threshold = new double[6];
//
//            for (int i = 0; i < adcs.length; i++) {
//                adc.get(sensor.getName()).fill(i, adcs[i] - sensor.getPedestal(channel, i));
//                threshold[i] = adcs[i] - sensor.getPedestal(channel, i) - 2 * sensor.getNoise(channel, i);
//            }
//
//            boolean sample1 = threshold[0] > 0 && threshold[1] > 0 && threshold[2] > 0;
//            boolean sample2 = threshold[1] > 0 && threshold[2] > 0 && threshold[3] > 0;
//            boolean sample3 = threshold[2] > 0 && threshold[3] > 0 && threshold[4] > 0;
//            boolean sample4 = threshold[3] > 0 && threshold[4] > 0 && threshold[5] > 0;
//
//            if (sample1 || sample2 || sample3 || sample4) {
//                channelAboveThresh.get(sensor.getName()).fill(channel);
//            } else {
//                channelBelowThresh.get(sensor.getName()).fill(channel);
//            }
            } // end of check whether RawTrackerHit has a pulse fit.
        }
        for (HpsSiSensor sensor : sensors) {
            Integer nHits = hits.get(sensor.getName());
            if (nHits == null) {
                nRawHits.get(sensor.getName()).fill(0);
            } else {
                nRawHits.get(sensor.getName()).fill(nHits);
                if (nHits > nHitMax) {
                    for (RawTrackerHit rawHit : rawHits) {
                        if (!sensor.equals(rawHit.getDetectorElement())) {
                            continue;
                        }
                        int channel = rawHit.getIdentifierFieldValue("strip");
//                        short[] adcs = rawHit.getADCValues();
//                        for (int i = 0; i < adcs.length; i++) {
//                            adcMonster.get(sensor.getName()).fill(i, adcs[i] - sensor.getPedestal(channel, i));
//                        }
//                        channelMonster.get(sensor.getName()).fill(channel);
                    }
                }
            }
        }
    }

    private void analyzeFittedRawTrackerHits(List<LCRelation> fittedHits) {
        aida.tree().mkdirs("fittedHitsAnalysis");
        aida.tree().cd("fittedHitsAnalysis");
        // Map the fitted hits to their corresponding raw hits
        Map<RawTrackerHit, List<LCRelation>> fittedRawTrackerHitMap = new HashMap<>();
        // use a list because the waveforms can have more than one valid fit...
        for (LCRelation fittedHit : fittedHits) {
            RawTrackerHit rth = FittedRawTrackerHit.getRawTrackerHit(fittedHit);
            aida.histogram1D("all hit t0", 250, -400., 100.).fill(FittedRawTrackerHit.getT0(fittedHit));
            aida.histogram1D("all hit t0 -150:100", 250, -150., 100.).fill(FittedRawTrackerHit.getT0(fittedHit));
            if (fittedRawTrackerHitMap.containsKey(rth)) {
                // this fittedHit is the second fit
                double t1 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(rth).get(0));
                double a1 = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rth).get(0));
                double t2 = FittedRawTrackerHit.getT0(fittedHit);
                double a2 = FittedRawTrackerHit.getAmp(fittedHit);
                aida.histogram2D("t1 vs t2", 250, -400., 100., 250, -400., 100.).fill(t1, t2);
                if (t1 > t2) {
                    aida.histogram1D("early hit amplitude", 250, 0., 5000.).fill(t2);
                    aida.histogram1D("early hit t0", 250, -400., 100.).fill(t2);
                    aida.histogram1D("late hit amplitude", 250, 0., 5000.).fill(t1);
                    aida.histogram1D("late hit t0", 250, -400., 100.).fill(t1);
                } else {
                    aida.histogram1D("early hit amplitude", 250, 0., 5000.).fill(t1);
                    aida.histogram1D("early hit t0", 250, -400., 100.).fill(t1);
                    aida.histogram1D("late hit amplitude", 250, 0., 5000.).fill(t2);
                    aida.histogram1D("late hit t0", 250, -400., 100.).fill(t2);
                }
                if (abs(t1) < abs(t2)) {
                    aida.histogram1D("hit t0 closer to 0", 250, -400., 100.).fill(t1);
                    aida.histogram1D("hit t0 farther from 0", 250, -400., 100.).fill(t2);
                } else {
                    aida.histogram1D("hit t0 closer to 0", 250, -400., 100.).fill(t2);
                    aida.histogram1D("hit t0 farther from 0", 250, -400., 100.).fill(t1);

                }
                fittedRawTrackerHitMap.get(rth).add(fittedHit);
            } else {
                List<LCRelation> arrayList = new ArrayList<>();
                arrayList.add(fittedHit);
                fittedRawTrackerHitMap.put(rth, arrayList);
            }
        }
        aida.histogram1D("Number of RawTrackerHits in Event", 100, 0., 2000.).fill(fittedRawTrackerHitMap.size());
        aida.histogram1D("Number of FittedRawTrackerHits in Event", 100, 0., 2000.).fill(fittedHits.size());
        aida.histogram1D("Number of extra FittedRawTrackerHits in Event", 100, 0., 500.).fill(fittedHits.size() - fittedRawTrackerHitMap.size());

        aida.tree().cd("..");
    }

    private void analyzeSiTrackerHitStrip1D(List<SiTrackerHitStrip1D> stripHits) {
        aida.tree().mkdirs("stripHit1DAnalysis");
        aida.tree().cd("stripHit1DAnalysis");
        int nStripHits = stripHits.size();
        aida.histogram1D("SVT number of StripHits", 100, 0., 1000.).fill(nStripHits);
        for (TrackerHit thit : stripHits) {
            SiTrackerHitStrip1D hit = new SiTrackerHitStrip1D(thit);
            aida.histogram1D("StripHit cluster size", 10, 0.5, 10.5).fill(hit.getRawHits().size());
            aida.histogram1D("StripHit time", 250, -150., 100.).fill(hit.getTime());
        }
        aida.tree().cd("..");
    }

    private void analyzeHitsOnTracks(List<ReconstructedParticle> rpList) {
        aida.tree().mkdirs("hitsOnKFTracksAnalysis");
        aida.tree().cd("hitsOnKFTracksAnalysis");
        for (ReconstructedParticle rp : rpList) {
            if (!rp.getTracks().isEmpty()) {
                Track t = rp.getTracks().get(0);
                int nHitsOnTrack = t.getTrackerHits().size();
                aida.histogram1D("number of hits on track", 20, 0.5, 20.5).fill(nHitsOnTrack);
                makeTrackHitPlots(t, null);
                if (nHitsOnTrack >= 10) {
                    aida.tree().mkdirs("greaterThan10Hits");
                    aida.tree().cd("greaterThan10Hits");
                    makeTrackHitPlots(t, null);
                    aida.tree().cd("..");
                }
                if (!rp.getClusters().isEmpty()) {
                    aida.tree().mkdirs("tracksWithClusters");
                    aida.tree().cd("tracksWithClusters");
                    makeTrackHitPlots(t, rp.getClusters().get(0));
                    aida.tree().cd("..");
                }
            }
        }
        aida.tree().cd("..");
    }

    private void makeTrackHitPlots(Track t, Cluster c) {
        //let's also calculate the track time
        double[] hitTimes = new double[t.getTrackerHits().size()];
        int nHit = 0;
        List<TrackerHit> hits = t.getTrackerHits();
        for (TrackerHit hit : hits) {
            HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement();
            String sensorName = sensor.getName();
            int layer = sensor.getLayerNumber();
            String sensorType = layer < 3 ? "thin" : "thick";
            int clusSize = hit.getRawHits().size();
            double dEdx = 1000000 * hit.getdEdx();
            double time = hit.getTime();
            hitTimes[nHit++] = time;
            aida.histogram1D("hit dEdx " + sensorType, 100, 0., 20.).fill(dEdx);
            aida.histogram1D("hit time " + sensorType, 100, -20., 20.).fill(time);
            aida.histogram1D("hit cluster size " + sensorName, 10, 0.5, 10.5).fill(clusSize);
            aida.histogram1D("hit dEdx ", 100, 0., 20.).fill(dEdx);
            aida.histogram1D("hit time ", 100, -20., 20.).fill(time);
            aida.histogram1D("hit time " + sensorName, 100, -20., 20.).fill(time);

            aida.histogram1D("hit cluster size ", 10, 0.5, 10.5).fill(clusSize);

            aida.histogram1D("hit time  -150:100", 250, -150., 100.).fill(time);
            aida.histogram1D("hit time  -150:100" + sensorType, 250, -150., 100.).fill(time);
        }

        double trackTime = StatUtils.mean(hitTimes);
        double trackTimeRms = sqrt(StatUtils.variance(hitTimes));
        double deltaTimes = StatUtils.max(hitTimes) - StatUtils.min(hitTimes);
        double maxPlusDiff = StatUtils.max(hitTimes) - trackTime;
        double maxMinusDiff = trackTime - StatUtils.min(hitTimes);
        aida.histogram1D("Track hit mean time", 150, -100., 50.).fill(trackTime);
        aida.histogram1D("Track hit time rms (sqrt(variance))", 100, 0., 20.).fill(trackTimeRms);
        aida.histogram1D("Track time max plus hit diff", 100, 0., 30.).fill(maxPlusDiff);
        aida.histogram1D("Track time max minus hit diff", 100, 0., 30.).fill(maxMinusDiff);
        aida.histogram1D("Track hit deltaTime", 100, 0., 40.).fill(deltaTimes);
        aida.histogram2D("Track hit time max plus vs max minus", 100, 0., 30., 100, 0., 30.).fill(maxPlusDiff, maxMinusDiff);
        aida.histogram2D("Track time vs rms", 150, -100., 50., 100, 0., 20.).fill(trackTime, trackTimeRms);
        aida.histogram2D("Track time vs deltaTime", 150, -100., 50., 100, 0., 40.).fill(trackTime, deltaTimes);
        if (c != null) {
            double clusterTime = ClusterUtilities.getSeedHitTime(c);
            aida.histogram1D("cluster time", 100, 0., 100.).fill(clusterTime);
            aida.histogram1D("track time - cluster time", 100, -50., -30.).fill(trackTime - clusterTime);
            aida.histogram2D("track time vs cluster time", 100, -40., 40., 100, 0., 100.).fill(trackTime, clusterTime);
            aida.histogram1D("Ecal cluster energy", 100, 0., 6.).fill(c.getEnergy());
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
