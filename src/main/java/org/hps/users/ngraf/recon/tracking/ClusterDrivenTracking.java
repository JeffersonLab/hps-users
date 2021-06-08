package org.hps.users.ngraf.recon.tracking;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.circlefit.CircleFit;
import org.hps.recon.tracking.circlefit.TwoPointRadiusCircleFitter;
import org.lcsim.constants.Constants;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.Transform3D;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.Identifier;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.detector.tracker.silicon.SiTrackerIdentifierHelper;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;

/**
 * Track-finding driven by Ecal Clusters
 *
 * Using the ECal cluster position and energy and assuming the track came from
 * the origin, derive the track parameters from a circle fit.
 *
 * Use this trajectory to find intercepts with the SVT sensors to get a list of
 * possible hits on the track.
 *
 * Fit the list of 1D strip clusters.
 *
 *
 * @author Norman A. Graf
 */
public class ClusterDrivenTracking extends Driver {

    private List<String> _allSensorNames = new ArrayList<>();
    private List<String> _topSensorNames = new ArrayList<>();
    private List<String> _bottomSensorNames = new ArrayList<>();
    private Map<String, Hep3Vector> _sensorPositionMap = new HashMap<>();
    private Map<String, Hep3Vector> _sensorDirectionMap = new HashMap<>();
    private org.lcsim.geometry.FieldMap _fieldMap;
    private double[] _origin = {-7.5, 0.};
    private Hep3Vector _ip = new BasicHep3Vector(0., 0., -7.5);
    private double _bField;

    @Override
    public void detectorChanged(Detector det) {
        _fieldMap = det.getFieldMap();
        _bField = TrackUtils.getBField(det).magnitude();
        List<HpsSiSensor> sensors = det.getSubdetector("Tracker").getDetectorElement().findDescendants(HpsSiSensor.class);
        for (HpsSiSensor sensor : sensors) {
            String name = sensor.getName();
            _allSensorNames.add(name);
            if (name.contains("t_half")) {
                _topSensorNames.add(name);
            }
            if (name.contains("b_half")) {
                _bottomSensorNames.add(name);
            }
            System.out.println(name + " layer number: " + sensor.getLayerNumber());
            Hep3Vector position = sensor.getGeometry().getPosition();
            System.out.println("position " + position);
            Hep3Vector[] unitVectors = getUnitVectors(sensor);
            System.out.println("direction " + unitVectors[2]);
            _sensorPositionMap.put(name, position);
            _sensorDirectionMap.put(name, unitVectors[2]);
        }
    }

    @Override
    public void process(EventHeader event) {
        // get the list of calorimeter clusters 
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        for (Cluster c : ecalClusters) {
            // want to fit a circle in the x,z plane
            double[] ecalClusterPos = c.getPosition();
            double[] pos = {c.getPosition()[2], c.getPosition()[0]};

            double radius = c.getEnergy() / (_bField * Constants.fieldConversion);
            System.out.println("energy: " + c.getEnergy() + " bField " + _bField + " radius " + radius + " omega " + 1 / radius);
            CircleFit[] circles = TwoPointRadiusCircleFitter.findCircles(_origin, pos, radius);
            if (circles == null) {
                System.out.println("found no solutions");
            } else {
                for (int j = 0; j < circles.length; ++j) {
                    System.out.println(circles[j] + " dx/dz " + circles[j].tangentAtPoint(pos) + " at " + pos[0] + "," + pos[1]);
                }
            }
            Hep3Vector line = new BasicHep3Vector(ecalClusterPos[0], ecalClusterPos[1], ecalClusterPos[2]);
            // get straight-line intercept with sensors this should work well for identifying axial strips...
            List<String> sensorsToUse = ecalClusterPos[1] > 0 ? _topSensorNames : _bottomSensorNames;
            for (String sensorName : sensorsToUse) {
                Hep3Vector xcept = getLinePlaneIntercept(line, _ip, _sensorPositionMap.get(sensorName), _sensorDirectionMap.get(sensorName));
                System.out.println("intercept with " + sensorName + " at " + xcept);
            }
        }
        System.out.println("\n\n\n");
    }

    private Hep3Vector[] getUnitVectors(SiSensor sensor) {

        Hep3Vector unit_vecU = new BasicHep3Vector(-99, -99, -99);
        Hep3Vector unit_vecV = new BasicHep3Vector(-99, -99, -99);
        Hep3Vector unit_vecW = new BasicHep3Vector(-99, -99, -99);

        for (ChargeCarrier carrier : ChargeCarrier.values()) {
            if (sensor.hasElectrodesOnSide(carrier)) {
                int channel = 1;
                long cell_id = sensor.makeStripId(channel, carrier.charge()).getValue();
                IIdentifier id = new Identifier(cell_id);
                SiTrackerIdentifierHelper _sid_helper = (SiTrackerIdentifierHelper) sensor.getIdentifierHelper();
                SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(carrier);
                ITransform3D local_to_global = new Transform3D();// sensor.getGeometry().getLocalToGlobal();
                ITransform3D electrodes_to_global = electrodes.getLocalToGlobal();
                ITransform3D global_to_hit = local_to_global.inverse();
                ITransform3D electrodes_to_hit = Transform3D.multiply(global_to_hit, electrodes_to_global);

                unit_vecU = electrodes_to_hit.rotated(electrodes.getMeasuredCoordinate(0));
                unit_vecV = electrodes_to_hit.rotated(electrodes.getUnmeasuredCoordinate(0));
                unit_vecW = VecOp.cross(unit_vecU, unit_vecV);
            }
        }
        return new Hep3Vector[]{unit_vecU, unit_vecV, unit_vecW};
    }

    /**
     * Finds point of intercept between a generic straight line and a plane.
     *
     * @param l - vector pointing along the line
     * @param l0 - point on the line
     * @param p0 - point on the plane
     * @param n - normal vector of the plane.
     * @return point of intercept.
     */
    private static Hep3Vector getLinePlaneIntercept(Hep3Vector l, Hep3Vector l0, Hep3Vector p0, Hep3Vector n) {
        if (VecOp.dot(l, n) == 0) {
            throw new RuntimeException("This line and plane are parallel!");
        }
        final double d = VecOp.dot(VecOp.sub(p0, l0), n) / VecOp.dot(l, n);
        Hep3Vector p = VecOp.add(VecOp.mult(d, l), l0);
        return p;

    }
}
