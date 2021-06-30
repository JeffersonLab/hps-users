package org.hps.users.ngraf.calibration;

import java.util.List;
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
public class SvtBaselineCalibrationHistogramFiller extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _nYbins = 500;
    private double _yLo = 3500.;
    private double _yHi = 7000.;

    protected void process(EventHeader event) {
        int run = event.getRunNumber();
        List<RawTrackerHit> rawTrackerHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");
        setupSensors(event);
        for (RawTrackerHit hit : rawTrackerHits) {
            int channel = hit.getIdentifierFieldValue("strip");
            short adc0 = hit.getADCValues()[0];
            String sensorName = ((HpsSiSensor) hit.getDetectorElement()).getName();
            if (sensorName.contains("_L1") || sensorName.contains("_L2")) {
                aida.histogram2D(run + " " + sensorName + " channel id vs APV25 channel 0", 512, 0, 512., _nYbins, _yLo, _yHi).fill(channel, adc0);
            } else {
                aida.histogram2D(run + " " + sensorName + " channel id vs APV25 channel 0", 640, 0, 640., _nYbins, _yLo, _yHi).fill(channel, adc0);
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

    public void setNYbins(int i) {
        _nYbins = i;
    }

    public void setYLo(double d) {
        _yLo = d;
    }

    public void setYHi(double d) {
        _yHi = d;
    }
}
