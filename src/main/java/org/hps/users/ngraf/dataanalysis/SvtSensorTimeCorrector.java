package org.hps.users.ngraf.dataanalysis;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;

/**
 *
 * @author Norman A. Graf
 */
public class SvtSensorTimeCorrector {

    private Map<String, Double> timeCorrections = new HashMap<>();

    public SvtSensorTimeCorrector() {
//        timeCorrections.put("module_L1b_halfmodule_axial_sensor0", -6.31168681391);
//        timeCorrections.put("module_L1b_halfmodule_stereo_sensor0", -8.0071655711);
//        timeCorrections.put("module_L2b_halfmodule_axial_sensor0", -6.021908338);
//        timeCorrections.put("module_L2b_halfmodule_stereo_sensor0", -5.52646419087);
//        timeCorrections.put("module_L3b_halfmodule_axial_sensor0", -8.05310652042);
//        timeCorrections.put("module_L3b_halfmodule_stereo_sensor0", -7.57543395597);
//        timeCorrections.put("module_L4b_halfmodule_axial_sensor0", -8.71221376122);
//        timeCorrections.put("module_L4b_halfmodule_stereo_sensor0", -7.86551785118);
//        timeCorrections.put("module_L5b_halfmodule_axial_hole_sensor0", -8.09560445872);
////        timeCorrections.put("module_L5b_halfmodule_stereo_hole_sensor0", 0.0);
//        timeCorrections.put("module_L5b_halfmodule_axial_slot_sensor0", -6.84603583651);
//        timeCorrections.put("module_L5b_halfmodule_stereo_slot_sensor0", -8.34041703412);
//        timeCorrections.put("module_L6b_halfmodule_axial_hole_sensor0", -11.1326845999);
//        timeCorrections.put("module_L6b_halfmodule_stereo_hole_sensor0", -10.6374436607);
//        timeCorrections.put("module_L6b_halfmodule_axial_slot_sensor0", -10.3205591365);
//        timeCorrections.put("module_L6b_halfmodule_stereo_slot_sensor0", -9.20881029289);
//        timeCorrections.put("module_L7b_halfmodule_axial_hole_sensor0", -11.1615456135);
//        timeCorrections.put("module_L7b_halfmodule_stereo_hole_sensor0", -9.39704538043);
//        timeCorrections.put("module_L7b_halfmodule_axial_slot_sensor0", -9.39470277174);
//        timeCorrections.put("module_L7b_halfmodule_stereo_slot_sensor0", -8.57903819588);
//        timeCorrections.put("module_L1t_halfmodule_axial_sensor0", -7.53200287342);
//        timeCorrections.put("module_L1t_halfmodule_stereo_sensor0", -7.5161613747);
//        timeCorrections.put("module_L2t_halfmodule_axial_sensor0", -8.72637333333);
//        timeCorrections.put("module_L2t_halfmodule_stereo_sensor0", -8.38213950864);
//        timeCorrections.put("module_L3t_halfmodule_axial_sensor0", -9.84669834908);
//        timeCorrections.put("module_L3t_halfmodule_stereo_sensor0", -9.69494606311);
//        timeCorrections.put("module_L4t_halfmodule_axial_sensor0", -10.0568006322);
//        timeCorrections.put("module_L4t_halfmodule_stereo_sensor0", -11.2607220152);
//        timeCorrections.put("module_L5t_halfmodule_axial_hole_sensor0", -8.80756657068);
//        timeCorrections.put("module_L5t_halfmodule_stereo_hole_sensor0", -10.0282613082);
//        timeCorrections.put("module_L5t_halfmodule_axial_slot_sensor0", -8.87331912841);
//        timeCorrections.put("module_L5t_halfmodule_stereo_slot_sensor0", -8.42071159597);
//        timeCorrections.put("module_L6t_halfmodule_axial_hole_sensor0", -8.97056066007);
//        timeCorrections.put("module_L6t_halfmodule_stereo_hole_sensor0", -10.6115963871);
//        timeCorrections.put("module_L6t_halfmodule_axial_slot_sensor0", -9.19052690163);
//        timeCorrections.put("module_L6t_halfmodule_stereo_slot_sensor0", -9.815320983);
//        timeCorrections.put("module_L7t_halfmodule_axial_hole_sensor0", -10.5540336142);
//        timeCorrections.put("module_L7t_halfmodule_stereo_hole_sensor0", -10.9966381764);
//        timeCorrections.put("module_L7t_halfmodule_axial_slot_sensor0", -8.54410719796);
//        timeCorrections.put("module_L7t_halfmodule_stereo_slot_sensor0", -9.00695297064);
        timeCorrections.put("module_L1b_halfmodule_axial_sensor0", -45.3650908009);
        timeCorrections.put("module_L1b_halfmodule_stereo_sensor0", -47.0886955565);
        timeCorrections.put("module_L2b_halfmodule_axial_sensor0", -45.1140425957);
        timeCorrections.put("module_L2b_halfmodule_stereo_sensor0", -44.6117976223);
        timeCorrections.put("module_L3b_halfmodule_axial_sensor0", -46.9972926303);
        timeCorrections.put("module_L3b_halfmodule_stereo_sensor0", -46.4972758514);
        timeCorrections.put("module_L4b_halfmodule_axial_sensor0", -47.7460323919);
        timeCorrections.put("module_L4b_halfmodule_stereo_sensor0", -46.8006035338);
        timeCorrections.put("module_L5b_halfmodule_axial_hole_sensor0", -47.203608838);
        timeCorrections.put("module_L5b_halfmodule_axial_slot_sensor0", -46.0600708111);
        timeCorrections.put("module_L5b_halfmodule_stereo_slot_sensor0", -47.9001175183);
        timeCorrections.put("module_L6b_halfmodule_axial_hole_sensor0", -49.9797378607);
        timeCorrections.put("module_L6b_halfmodule_stereo_hole_sensor0", -49.5027781972);
        timeCorrections.put("module_L6b_halfmodule_axial_slot_sensor0", -49.1458941531);
        timeCorrections.put("module_L6b_halfmodule_stereo_slot_sensor0", -48.0826812988);
        timeCorrections.put("module_L7b_halfmodule_axial_hole_sensor0", -50.1135299177);
        timeCorrections.put("module_L7b_halfmodule_stereo_hole_sensor0", -48.2011297762);
        timeCorrections.put("module_L7b_halfmodule_axial_slot_sensor0", -48.2401946106);
        timeCorrections.put("module_L7b_halfmodule_stereo_slot_sensor0", -47.3917337673);
        timeCorrections.put("module_L1t_halfmodule_axial_sensor0", -46.3807592829);
        timeCorrections.put("module_L1t_halfmodule_stereo_sensor0", -46.4296785525);
        timeCorrections.put("module_L2t_halfmodule_axial_sensor0", -47.6107862045);
        timeCorrections.put("module_L2t_halfmodule_stereo_sensor0", -47.2649694096);
        timeCorrections.put("module_L3t_halfmodule_axial_sensor0", -48.5455507938);
        timeCorrections.put("module_L3t_halfmodule_stereo_sensor0", -48.3322018429);
        timeCorrections.put("module_L4t_halfmodule_axial_sensor0", -48.7006956478);
        timeCorrections.put("module_L4t_halfmodule_stereo_sensor0", -49.9815525769);
        timeCorrections.put("module_L5t_halfmodule_axial_hole_sensor0", -47.2412363932);
        timeCorrections.put("module_L5t_halfmodule_stereo_hole_sensor0", -48.4693565272);
        timeCorrections.put("module_L5t_halfmodule_axial_slot_sensor0", -47.901580325);
        timeCorrections.put("module_L5t_halfmodule_stereo_slot_sensor0", -47.1583942748);
        timeCorrections.put("module_L6t_halfmodule_axial_hole_sensor0", -47.2994007137);
        timeCorrections.put("module_L6t_halfmodule_stereo_hole_sensor0", -49.0185414994);
        timeCorrections.put("module_L6t_halfmodule_axial_slot_sensor0", -48.0730406695);
        timeCorrections.put("module_L6t_halfmodule_stereo_slot_sensor0", -48.6469552228);
        timeCorrections.put("module_L7t_halfmodule_axial_hole_sensor0", -49.0809792435);
        timeCorrections.put("module_L7t_halfmodule_stereo_hole_sensor0", -49.5656646857);
        timeCorrections.put("module_L7t_halfmodule_axial_slot_sensor0", -47.3897038371);
        timeCorrections.put("module_L7t_halfmodule_stereo_slot_sensor0", -47.7726522099);
    }

    public double correctTrackTime(Track track) {
        double trackTime = 0.;
        int nHits = 0;
        for (TrackerHit hit : track.getTrackerHits()) {
            List<RawTrackerHit> rthits = hit.getRawHits();
            RawTrackerHit rawHit = rthits.get(0);
            HpsSiSensor sensor = ((HpsSiSensor) rawHit.getDetectorElement());
            String sensorName = sensor.getName();
            trackTime += hit.getTime() - timeCorrections.get(sensorName);
            nHits++;
        }
        return trackTime / nHits;
    }

}
