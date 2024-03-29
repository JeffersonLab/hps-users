package org.hps.users.ngraf.alignment.straighttrack;

import Jama.Matrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MaterialSupervisor;
import org.hps.recon.tracking.MaterialSupervisor.ScatteringDetectorVolume;
import org.hps.recon.tracking.MaterialSupervisor.SiStripPlane;
import static org.hps.users.ngraf.alignment.straighttrack.FitTracks.GEN_ROTMAT;
import static org.hps.users.ngraf.alignment.straighttrack.FitTracks.PROD_ROT;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.detector.ILogicalVolume;
import org.lcsim.detector.IPhysicalVolume;
import org.lcsim.detector.IPhysicalVolumeContainer;
import org.lcsim.detector.material.IMaterial;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiTrackerModule;
import org.lcsim.geometry.Detector;

/**
 *
 * @author ngraf
 */
public class DetectorBuilder {

    Detector _det;
    Map<String, String[]> sensorNameMap = new HashMap<String, String[]>();
    Map<String, SiStripPlane> stripPlaneNameMap = new HashMap<String, SiStripPlane>();
    Map<String, List<DetectorPlane>> trackerMap = new HashMap<String, List<DetectorPlane>>();
    Map<String, DetectorPlane> planeMap = new HashMap<String, DetectorPlane>();

    MaterialSupervisor materialManager;

    String _detectorName;

    Boolean _drawDetector = true;

    // should be enums?
    public enum Tracker {
        TOPHOLE("topHole"),
        TOPSLOT("topSlot"),
        BOTTOMHOLE("bottomHole"),
        BOTTOMSLOT("bottomSlot");
        private String _name;

        Tracker(String name) {
            _name = name;
        }

        String trackerName() {
            return _name;
        }
    }

    Vector3D vX = Vector3D.PLUS_I;
    Vector3D vY = Vector3D.PLUS_J;

    // TODO fix this...
    double[] SIGS = {0.0055, 0.00};     //  detector resolutions in mm

    boolean _debug = false;

    String[] bottomSlotNames = {
        "module_L1b_halfmodule_stereo_sensor0",
        "module_L1b_halfmodule_axial_sensor0",
        "module_L2b_halfmodule_stereo_sensor0",
        "module_L2b_halfmodule_axial_sensor0",
        "module_L3b_halfmodule_stereo_sensor0",
        "module_L3b_halfmodule_axial_sensor0",
        "module_L4b_halfmodule_stereo_sensor0",
        "module_L4b_halfmodule_axial_sensor0",
        "module_L5b_halfmodule_stereo_slot_sensor0",
        "module_L5b_halfmodule_axial_slot_sensor0",
        "module_L6b_halfmodule_stereo_slot_sensor0",
        "module_L6b_halfmodule_axial_slot_sensor0",
        "module_L7b_halfmodule_stereo_slot_sensor0",
        "module_L7b_halfmodule_axial_slot_sensor0"};
    String[] bottomHoleNames = {
        "module_L1b_halfmodule_stereo_sensor0",
        "module_L1b_halfmodule_axial_sensor0",
        "module_L2b_halfmodule_stereo_sensor0",
        "module_L2b_halfmodule_axial_sensor0",
        "module_L3b_halfmodule_stereo_sensor0",
        "module_L3b_halfmodule_axial_sensor0",
        "module_L4b_halfmodule_stereo_sensor0",
        "module_L4b_halfmodule_axial_sensor0",
        "module_L5b_halfmodule_stereo_hole_sensor0",
        "module_L5b_halfmodule_axial_hole_sensor0",
        "module_L6b_halfmodule_stereo_hole_sensor0",
        "module_L6b_halfmodule_axial_hole_sensor0",
        "module_L7b_halfmodule_stereo_hole_sensor0",
        "module_L7b_halfmodule_axial_hole_sensor0"};

    String[] topSlotNames = {
        "module_L1t_halfmodule_axial_sensor0",
        "module_L1t_halfmodule_stereo_sensor0",
        "module_L2t_halfmodule_axial_sensor0",
        "module_L2t_halfmodule_stereo_sensor0",
        "module_L3t_halfmodule_axial_sensor0",
        "module_L3t_halfmodule_stereo_sensor0",
        "module_L4t_halfmodule_axial_sensor0",
        "module_L4t_halfmodule_stereo_sensor0",
        "module_L5t_halfmodule_axial_slot_sensor0",
        "module_L5t_halfmodule_stereo_slot_sensor0",
        "module_L6t_halfmodule_axial_slot_sensor0",
        "module_L6t_halfmodule_stereo_slot_sensor0",
        "module_L7t_halfmodule_axial_slot_sensor0",
        "module_L7t_halfmodule_stereo_slot_sensor0"};
    String[] topHoleNames = {
        "module_L1t_halfmodule_axial_sensor0",
        "module_L1t_halfmodule_stereo_sensor0",
        "module_L2t_halfmodule_axial_sensor0",
        "module_L2t_halfmodule_stereo_sensor0",
        "module_L3t_halfmodule_axial_sensor0",
        "module_L3t_halfmodule_stereo_sensor0",
        "module_L4t_halfmodule_axial_sensor0",
        "module_L4t_halfmodule_stereo_sensor0",
        "module_L5t_halfmodule_axial_hole_sensor0",
        "module_L5t_halfmodule_stereo_hole_sensor0",
        "module_L6t_halfmodule_axial_hole_sensor0",
        "module_L6t_halfmodule_stereo_hole_sensor0",
        "module_L7t_halfmodule_axial_hole_sensor0",
        "module_L7t_halfmodule_stereo_hole_sensor0"};

    public DetectorBuilder(String detectorName) {
        final DatabaseConditionsManager manager = DatabaseConditionsManager.getInstance();
        try {
            manager.setDetector(detectorName, 5772);
        } catch (ConditionsManager.ConditionsNotFoundException ex) {
            Logger.getLogger(DetectorBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }
        _det = manager.getCachedConditions(Detector.class, "compact.xml").getCachedData();
        System.out.println(_det.getName());
        setup();
        _detectorName = detectorName;
    }

    public DetectorBuilder(Detector detector) {
        _det = detector;
        System.out.println(_det.getName());
        setup();
        _detectorName = detector.getDetectorName();
    }

    public DetectorBuilder(Path path) {
        sensorNameMap.put("topHole", topHoleNames);
        sensorNameMap.put("bottomHole", bottomHoleNames);
        sensorNameMap.put("topSlot", topSlotNames);
        sensorNameMap.put("bottomSlot", bottomSlotNames);
        int oIndex = 2;
        int uIndex = 8;
        int vIndex = 14;
        int wIndex = 20;
        int aIndex = 26;
        List<String> lines = null;
        try {
            lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        } catch (IOException ex) {
            Logger.getLogger(DetectorBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }
        //List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
        // first line is the detector name
        _detectorName = lines.remove(0);
        // rest is plane data (get it?)
        for (String line : lines) {
            if (_debug) {
                System.out.println(line);
            }
            String[] s = line.split("\\s+");
            int id = Integer.valueOf(s[0]);
            String stripPlaneName = s[1];
            double[] hpsAngles = arrayFromStrings(s, aIndex);
            Hep3Vector origin = hep3vectorFromStrings(s, oIndex);
            Hep3Vector uDir = hep3vectorFromStrings(s, uIndex);
            Hep3Vector vDir = hep3vectorFromStrings(s, vIndex);
            Hep3Vector normal = hep3vectorFromStrings(s, wIndex);
            double width = Double.valueOf(s[32]);
            double height = Double.valueOf(s[33]);
            double zmin = Double.valueOf(s[34]);
            double zmax = Double.valueOf(s[35]);
            if (_debug) {
                System.out.println(" " + stripPlaneName);
                System.out.println("   origin: " + origin);
                System.out.println("   normal: " + normal);
                System.out.println("   uDir: " + uDir);
                System.out.println("   vDir: " + vDir);
                System.out.println("   Apache commons angles: " + Arrays.toString(hpsAngles));
                System.out.println("zmin " + zmin + " zmax " + zmax);
            }
            DetectorPlane dp = new DetectorPlane(id++, matrixFromAngles(hpsAngles), origin.v(), SIGS);
            dp.setName(stripPlaneName);
            dp.setUVWR(uDir, vDir, normal, origin);
            dp.setDimensions(width, height, zmin, zmax);
            dp.setAngles(hpsAngles);
            //planes.add(dp);
            if (_debug) {
                System.out.println("  " + dp);
            }
            planeMap.put(dp.name(), dp);
        }
        if (_debug) {
            System.out.println("Populated planemap with " + planeMap.size() + " planes");
        }
        for (Tracker t : Tracker.values()) {
            String s = t.trackerName();
            if (_debug) {
                System.out.println("building " + s + " : ");
            }
            String[] stripNames = sensorNameMap.get(s);
            List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
            int id = 1;
            for (String stripPlaneName : stripNames) {
                planes.add(planeMap.get(stripPlaneName));
            }
            trackerMap.put(s, planes);
            if (_debug) {
                System.out.println("TrackerMap " + s + " has " + planes.size() + " planes");
            }
        }

    }

    public DetectorBuilder(List<String> lines) {
        sensorNameMap.put("topHole", topHoleNames);
        sensorNameMap.put("bottomHole", bottomHoleNames);
        sensorNameMap.put("topSlot", topSlotNames);
        sensorNameMap.put("bottomSlot", bottomSlotNames);
        int oIndex = 2;
        int uIndex = 8;
        int vIndex = 14;
        int wIndex = 20;
        int aIndex = 26;
        //List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
        // first line is the detector name
        _detectorName = lines.remove(0);
        // rest is plane data (get it?)
        for (String line : lines) {
            if (_debug) {
                System.out.println(line);
            }
            String[] s = line.split("\\s+");
            int id = Integer.valueOf(s[0]);
            String stripPlaneName = s[1];
            double[] hpsAngles = arrayFromStrings(s, aIndex);
            Hep3Vector origin = hep3vectorFromStrings(s, oIndex);
            Hep3Vector uDir = hep3vectorFromStrings(s, uIndex);
            Hep3Vector vDir = hep3vectorFromStrings(s, vIndex);
            Hep3Vector normal = hep3vectorFromStrings(s, wIndex);
            double width = Double.valueOf(s[32]);
            double height = Double.valueOf(s[33]);
            double zmin = Double.valueOf(s[34]);
            double zmax = Double.valueOf(s[35]);
            if (_debug) {
                System.out.println(" " + stripPlaneName);
                System.out.println("   origin: " + origin);
                System.out.println("   normal: " + normal);
                System.out.println("   uDir: " + uDir);
                System.out.println("   vDir: " + vDir);
                System.out.println("   Apache commons angles: " + Arrays.toString(hpsAngles));
                System.out.println("zmin " + zmin + " zmax " + zmax);
            }
            DetectorPlane dp = new DetectorPlane(id++, matrixFromAngles(hpsAngles), origin.v(), SIGS);
            dp.setName(stripPlaneName);
            dp.setUVWR(uDir, vDir, normal, origin);
            dp.setDimensions(width, height, zmin, zmax);
            dp.setAngles(hpsAngles);
            //planes.add(dp);
            if (_debug) {
                System.out.println("  " + dp);
            }
            planeMap.put(dp.name(), dp);
        }
        if (_debug) {
            System.out.println("Populated planemap with " + planeMap.size() + " planes");
        }
        for (Tracker t : Tracker.values()) {
            String s = t.trackerName();
            if (_debug) {
                System.out.println("building " + s + " : ");
            }
            String[] stripNames = sensorNameMap.get(s);
            List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
            int id = 1;
            for (String stripPlaneName : stripNames) {
                planes.add(planeMap.get(stripPlaneName));
            }
            trackerMap.put(s, planes);
            if (_debug) {
                System.out.println("TrackerMap " + s + " has " + planes.size() + " planes");
            }
        }

    }

    private Matrix matrixFromAngles(double[] hpsAngles) {
        Matrix[] mats = new Matrix[3];
        double[][] RW = new double[3][9];
        for (int j = 0; j < 3; ++j) {
            double angl = hpsAngles[j];
            mats[j] = GEN_ROTMAT(angl, j);
            GEN_ROTMAT(angl, j, RW[j]);
            if (_debug) {
                System.out.println("Rotation matrix for axis " + (j + 1) + " angle " + angl);
            }
            if (_debug) {
                mats[j].print(6, 4);
            }
        }
        return PROD_ROT(mats[0], mats[1], mats[2]);
    }

    private double[] arrayFromStrings(String[] s, int index) {
        return new double[]{Double.valueOf(s[index + 2]), Double.valueOf(s[index + 3]), Double.valueOf(s[index + 4])};
    }

    private Hep3Vector hep3vectorFromStrings(String[] s, int index) {
        return new BasicHep3Vector(Double.valueOf(s[index + 2]), Double.valueOf(s[index + 3]), Double.valueOf(s[index + 4]));
    }

    private void setup() {
        sensorNameMap.put("topHole", topHoleNames);
        sensorNameMap.put("bottomHole", bottomHoleNames);
        sensorNameMap.put("topSlot", topSlotNames);
        sensorNameMap.put("bottomSlot", bottomSlotNames);
        materialManager = new MaterialSupervisor();
        materialManager.buildModel(_det);

        for (ScatteringDetectorVolume vol : materialManager.getMaterialVolumes()) {
            SiStripPlane plane = (SiStripPlane) vol;
            stripPlaneNameMap.put(plane.getName(), plane);
            if (_debug) {
                System.out.println(plane.getName());
            }
        }

        for (Tracker t : Tracker.values()) {
            String s = t.trackerName();
            if (_debug) {
                System.out.println("building " + s + " : ");
            }
            String[] stripNames = sensorNameMap.get(s);
            List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
            int id = 1;
            for (String stripPlaneName : stripNames) {

                SiStripPlane plane = stripPlaneNameMap.get(stripPlaneName);
                if (_debug) {
                    System.out.println(stripPlaneName);
                }
                // the following method is off by 100 microns for layers 1 and 2 and 160 microns for the rest
                // changed 20200524
//                Hep3Vector origin = CoordinateTransformations.transformVectorToDetector(plane.origin());
                Hep3Vector origin = plane.getSensor().getGeometry().getLocalToGlobal().getTranslation().getTranslationVector();

                Hep3Vector normal = CoordinateTransformations.transformVectorToDetector(plane.normal());
                // orient all normals to point along +ive z
                if (normal.z() < 0.) {
                    normal = VecOp.neg(normal);
                }
                Hep3Vector vDir = CoordinateTransformations.transformVectorToDetector(plane.getUnmeasuredCoordinate());
                Hep3Vector uDir = VecOp.cross(normal, vDir);

                // extract the rotation angles...
                // Note that these are "reversed" because the code assumes the u measurement direction is along x, not y
                // so we need an additional PI/2 rotation...
                Vector3D vYprime = new Vector3D(vDir.x(), vDir.y(), vDir.z());  // nominal x
                Vector3D vXprime = new Vector3D(uDir.x(), uDir.y(), uDir.z());   // nominal y
                // create a rotation matrix from this pair of vectors
                Rotation xyVecRot = new Rotation(vX, vY, vXprime, vYprime);
                double[] hpsAngles = xyVecRot.getAngles(RotationOrder.XYZ, RotationConvention.VECTOR_OPERATOR);
//                if (_debug) {
                Matrix[] mats = new Matrix[3];
                double[][] RW = new double[3][9];
                for (int j = 0; j < 3; ++j) {
                    double angl = hpsAngles[j];
                    mats[j] = GEN_ROTMAT(angl, j);
                    GEN_ROTMAT(angl, j, RW[j]);
                    if (_debug) {
                        System.out.println("Rotation matrix for axis " + (j + 1) + " angle " + angl);

                        mats[j].print(6, 4);
                    }
                }
                Matrix prodrot = PROD_ROT(mats[0], mats[1], mats[2]);
                // calculate zmin, zmax
                double width = plane.getUnmeasuredDimension() / 2.;
                double height = plane.getMeasuredDimension() / 2.;
                double[] bounds = findZBounds(origin, VecOp.mult(width, vDir), VecOp.mult(height, uDir));
                double zmin = bounds[0];
                double zmax = bounds[1];
                if (_debug) {
                    System.out.println(" " + stripPlaneName);
                    System.out.println("   origin: " + origin);
                    System.out.println("   normal: " + normal);
                    System.out.println("   uDir: " + uDir);
                    System.out.println("   vDir: " + vDir);
                    System.out.println("   Apache commons angles: " + Arrays.toString(hpsAngles));
                    System.out.println("zmin " + zmin + " zmax " + zmax);
                }
                DetectorPlane dp = new DetectorPlane(id++, prodrot, origin.v(), SIGS);
                dp.setName(stripPlaneName);
                dp.setUVWR(uDir, vDir, normal, origin);
                dp.setDimensions(width, height, zmin, zmax);
                dp.setAngles(hpsAngles);
                planes.add(dp);
                if (_debug) {
                    System.out.println("  " + dp);
                }
                planeMap.put(plane.getName(), dp);
            }
            trackerMap.put(s, planes);
            if (_debug) {
                System.out.println("");
            }
        }
        if (_debug) {
            System.out.println("Populated planemap with " + planeMap.size() + " planes");
        }
    }

    public static double[] findZBounds(Hep3Vector origin, Hep3Vector width, Hep3Vector height) {
        Hep3Vector[] corners = new Hep3Vector[4];
        double zmin = 999.;
        double zmax = -999;
        // o + w*vDir + h*uDir

        Hep3Vector edge = VecOp.add(origin, width);
        corners[0] = VecOp.add(edge, height);
        corners[1] = VecOp.sub(edge, height);
        edge = VecOp.sub(origin, width);
        corners[2] = VecOp.add(edge, height);
        corners[3] = VecOp.sub(edge, height);

        for (int i = 0; i < 4; ++i) {
            //System.out.println("corner " + i + " : " + corners[i]);
            if (corners[i].z() > zmax) {
                zmax = corners[i].z();
            }
            if (corners[i].z() < zmin) {
                zmin = corners[i].z();
            }
        }
        //System.out.println("zmin " + zmin + " zmax " + zmax);
        return new double[]{zmin, zmax};
    }

    public List<DetectorPlane> getTracker(String trackerName) {
        return trackerMap.get(trackerName);
    }

    public String[] getTrackerSensorNames(String trackerName) {
        return sensorNameMap.get(trackerName);
    }

    public void drawDetector() throws Exception {
        drawDetector(_det);
    }

    private void drawDetector(Detector det) throws Exception {
        //TODO provide this functionality working from the planemap
        boolean debug = false;
        List<MaterialSupervisor.ScatteringDetectorVolume> stripPlanes = materialManager.getMaterialVolumes();
        //TODO replace these lists with a helper class.
        List<Hep3Vector> oList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> uList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> vList = new ArrayList<Hep3Vector>();
        List<Double> measDim = new ArrayList<Double>();
        List<Double> unmeasDim = new ArrayList<Double>();
        List<Boolean> isAxial = new ArrayList<Boolean>();

        for (MaterialSupervisor.ScatteringDetectorVolume vol : stripPlanes) {
            MaterialSupervisor.SiStripPlane plane = (MaterialSupervisor.SiStripPlane) vol;
            if (debug) {
                System.out.println(plane.getName());
            }
            Hep3Vector oprime = CoordinateTransformations.transformVectorToDetector(plane.origin());
            Hep3Vector nprime = CoordinateTransformations.transformVectorToDetector(plane.normal());
            if (debug) {
                System.out.println(" origin: " + oprime);
                System.out.println(" normal: " + nprime);
                System.out.println(" Plane is: " + plane.getMeasuredDimension() + " x " + plane.getUnmeasuredDimension());
            }
            HpsSiSensor sensor = (HpsSiSensor) plane.getSensor();

//            if (debug) {
//                System.out.println(SvtUtils.getInstance().isAxial(sensor) ? "axial" : "stereo");
//            }
            Hep3Vector unmeasDir = CoordinateTransformations.transformVectorToDetector(plane.getUnmeasuredCoordinate());
            // by default, the HPS unmeasDir points along -x, so invert here...
            // what!?
//            unmeasDir = VecOp.neg(unmeasDir);
            Hep3Vector measDir = CoordinateTransformations.transformVectorToDetector(plane.getMeasuredCoordinate());
            //measDir points either up or down depending on orientation of the sensor. I want it always to point up, along y.
            // so let's check something here...
            if (measDir.y() < 0.) {
                measDir = VecOp.neg(measDir);
            }
            Hep3Vector tst = VecOp.cross(unmeasDir, measDir);
            // if pointing along the z axis, OK. if not, invert unmeasDir...
            if (tst.z() < 0.) {
                unmeasDir = VecOp.neg(unmeasDir);
            }
            if (debug) {
                System.out.println("measured coordinate:    " + measDir);
                System.out.println("unmeasured coordinate:   " + unmeasDir);
                System.out.println("thickness: " + plane.getThickness() + " in X0: " + plane.getThicknessInRL());
            }
            SiTrackerModule module = (SiTrackerModule) plane.getSensor().getGeometry().getDetectorElement().getParent();
            IPhysicalVolume parent = module.getGeometry().getPhysicalVolume();
            IPhysicalVolumeContainer daughters = parent.getLogicalVolume().getDaughters();
            if (debug) {
                System.out.printf("%s found %d daughters to SiTrackerModule\n", this.getClass().getSimpleName(), daughters.size());
            }
            for (IPhysicalVolume daughter : daughters) {
                ILogicalVolume logicalVolume = daughter.getLogicalVolume();
                IMaterial material = logicalVolume.getMaterial();
                //System.out.println(material);
                String name = material.getName();
                double X0 = 10.0 * material.getRadiationLength() / material.getDensity();
                double X0cm = material.getRadiationLengthWithDensity();
                if (debug) {
                    System.out.println("material: " + name + " with X0= " + X0 + " mm " + X0cm);
                }
                Box solid = (Box) logicalVolume.getSolid();
                if (debug) {
                    System.out.printf("%s x %f y %f z %f box\n", this.getClass().getSimpleName(), solid.getXHalfLength(), solid.getYHalfLength(), solid.getZHalfLength());
                }
                double halfThickness = solid.getZHalfLength();
            }
            oList.add(oprime);
            uList.add(measDir);
            vList.add(unmeasDir);
            measDim.add(plane.getMeasuredDimension());
            unmeasDim.add(plane.getUnmeasuredDimension());
            isAxial.add(sensor.isAxial());
        }

//        // now get ECal
//        List<List<Hep3Vector>> calcells = new ArrayList<List<Hep3Vector>>();
//        IIdentifierHelper helper = _api.getIdentifierHelper();
//        for (EcalCrystal crystal : _api.getCrystals()) {
//            IGeometryInfo geom = crystal.getGeometry();
//            Trd trd = (Trd) geom.getLogicalVolume().getSolid();
//            List<Point3D> vertices = trd.getVertices();
//            List<Hep3Vector> verts = new ArrayList<Hep3Vector>();
//            for (Point3D p : vertices) {
//                verts.add(geom.transformLocalToGlobal(p.getHep3Vector()));
//            }
//            calcells.add(verts);
//            //
////            IExpandedIdentifier expandedId = helper.createExpandedIdentifier();
////            expandedId.setValue(helper.getFieldIndex("ix"), crystal.getX());
////            expandedId.setValue(helper.getFieldIndex("iy"), crystal.getY());
////            expandedId.setValue(helper.getFieldIndex("system"), _api.getSystemID());
////            IIdentifier id = helper.pack(expandedId);
////            if (id.getValue() != crystal.getIdentifier().getValue()) {
////                throw new RuntimeException("Reencoded ID " + id.getValue() + " does not match crystal ID " + crystal.getIdentifier().getValue());
////            }
////
////            IDetectorElementStore deStore = DetectorElementStore.getInstance();
////            if (deStore.find(crystal.getIdentifier()).size() == 0) {
////                throw new RuntimeException("Failed to find crystal ID in store.");
////            }
////
////            if (deStore.find(id).size() == 0) {
////                throw new RuntimeException("Failed to find repacked ID in store.");
////            }
//        }
//        System.out.println("calcells has " + calcells.size() + " entries");
//        // test print...
//        List<Hep3Vector> vs0 = calcells.get(0);
//        System.out.println("cell 0 has " + vs0.size() + "points");
//        System.out.println(vs0);
//        List<Hep3Vector> vs1 = calcells.get(1);
//        System.out.println("cell 1 has " + vs1.size() + "points");
//        System.out.println(vs1);
        // now draw
        if (_drawDetector) {
            FileOutputStream fos = new FileOutputStream(_detectorName + "_" + myDate() + "_3Dgeom.asy");
            Writer w = new BufferedWriter(new OutputStreamWriter(fos, "Cp850"));

            //TODO extract the following into a separate method
            w.write("\n\n\n");

            w.write("size(600);\n");
            w.write("import three;\n");
            w.write("currentprojection=orthographic(\n");
            w.write("camera=(-1250,850,-1500),\n");
            w.write("up=(0.374061249845773,1.48257908496856,0.568508408053946),\n");
            w.write("target=(0,0,0),\n");
            w.write("zoom=1);\n");

            w.write("\n");
            w.write("pen fill=gray(0.6)+opacity(0.2);\n");
            w.write("pen edge=black;\n");
            w.write("pen axialEdge=red;\n");
            w.write("pen stereoEdge=blue;\n");
            w.write("// o is the origin point in the plane, u is the measurement direction, v is the unmeasured direction\n");
            w.write("// clockwise normal orientation corresponds to z = u x v\n");
            w.write("//\n");
            w.write("//     +----------------------+  ^\n");
            w.write("//     |                      |  |\n");
            w.write("//     |                      |  u\n");
            w.write("//     |           o          |  |\n");
            w.write("//     |                      |   \n");
            w.write("//     |                      |   \n");
            w.write("//     +----------------------+   \n");
            w.write("//               ------ v ---->   \n");
            w.write("path3 myPlane(triple o, triple u, triple v)\n");
            w.write("{\n");
            w.write("    path3 P;\n");
            w.write("    P=o+(u+v)--o-u+v--o-u-v--o+u-v--cycle;\n");
            w.write("    return P;\n");
            w.write("}\n");
            w.write(" void drawAxes(triple o, triple u, triple v)\n");
            w.write(" {\n");
            w.write("pen p = red;\n");
            w.write(" draw(Label(\"$u$\",1),(o--o+u),p,Arrow3);\n");
            w.write("p=green;\n");
            w.write(" draw(Label(\"$v$\",1),(o--o+v),p,Arrow3);\n");
            w.write("p=blue;\n");
            w.write(" draw(Label(\"$w$\",1),(o--o+10.*cross(unit(v),unit(u))),p,Arrow3);\n");
            w.write(" }\n");

            // ECal crystal frustum volume
            // use volumes for hit crystals in event data
            w.write("void drawvolume(triple p1, triple p2, triple p3, triple p4, triple p5, triple p6, triple p7, triple p8, pen surfpen)\n");
            w.write("{\n");
            w.write("surface[] surfs;\n");
            w.write("\n");
            w.write("surfs.push( surface(p1--p2--p4--p3--cycle) );\n");
            w.write("surfs.push( surface(p4--p8--p7--p3--cycle) );\n");
            w.write("surfs.push( surface(p3--p7--p5--p1--cycle) );\n");
            w.write("surfs.push( surface(p2--p1--p5--p6--cycle) );\n");
            w.write("surfs.push( surface(p2--p6--p8--p4--cycle) );\n");
            w.write("surfs.push( surface(p6--p5--p7--p8--cycle) );\n");
            w.write("\n");
            w.write("for (int i = 0; i < surfs.length; ++i) {\n");
            w.write("        draw(surfs[i], surfpen);\n");
            w.write("    }\n");
            w.write("\n");
            w.write("}\n");

            //ECal crystal frustum edges
            // Use edges to draw detector
            // this causes out-of-memory crashes.
            // keep this here, but don't use it for now.
            w.write("void drawEdges(triple p1, triple p2, triple p3, triple p4, triple p5, triple p6, triple p7, triple p8, pen edgepen)\n");
            w.write("{\n");
            w.write("draw(p1--p2--p3--p4--cycle, edgepen);\n");
            w.write("draw(p2--p6--p7--p3--cycle, edgepen);\n");
            w.write("draw(p3--p7--p8--p4--cycle, edgepen);\n");
            w.write("draw(p1--p4--p8--p5--cycle, edgepen);\n");
            w.write("draw(p1--p5--p6--p2--cycle, edgepen);\n");
            w.write("draw(p5--p8--p7--p6--cycle, edgepen);\n");
            w.write("}\n");

            DecimalFormat df = new DecimalFormat("###.######\n");

            w.write("draw(Label(\"$x$\",1),(O--100.*X),Arrow3(HookHead3));\n");
            w.write("draw(Label(\"$y$\",1),(O--100.*Y),Arrow3(HookHead3));\n");
            w.write("draw(Label(\"$z$\",1),(O--1000.*Z),Arrow3(HookHead3));\n");

            for (int i = 0; i < oList.size(); ++i) {
                Hep3Vector o = oList.get(i);
                Hep3Vector u = uList.get(i);
                Hep3Vector v = vList.get(i);
                double l = unmeasDim.get(i) / 2.;
                double h = measDim.get(i) / 2.;
                //redefine op as lower left corner, viz.
                // o' = o-(l*V+h*u)
                // Hep3Vector op = VecOp.sub(o, VecOp.add(VecOp.mult(l, v), VecOp.mult(h, u)));
                String planeType = isAxial.get(i) ? "axial" : "stereo";
                w.write("//" + planeType + "\n");
//                StringBuffer sb = new StringBuffer("draw(surface( myPlane( ");
//                //origin
//                sb.append("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), ");
//                //u coordinate
//                sb.append(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), ");
//                //v coordinate
//                sb.append(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " )  ) ), fill, " + planeType + "Edge); ");
//                // sensor axes
//                sb = new StringBuffer("drawAxes( ");
//                //origin
//                sb.append("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), ");
//                //u coordinate
//                sb.append(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), ");
//                //v coordinate
//                sb.append(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " ) ); ");
//
//                w.write(sb.toString()+"\n");
                w.write("draw(surface( myPlane( " + "\n");
                //origin
                w.write("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), " + "\n");
                //u coordinate
                w.write(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), " + "\n");
                //v coordinate
                w.write(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " )  ) ), fill, " + planeType + "Edge); " + "\n");
                // sensor axes
                w.write("drawAxes( " + "\n");
                //origin
                w.write("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), " + "\n");
                //u coordinate
                w.write(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), " + "\n");
                //v coordinate
                w.write(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " ) ); " + "\n");

            }
            // now the calorimeter...
//            StringBuffer sbcal = new StringBuffer();
//            for (List<Hep3Vector> vertices : calcells) {
//                addcell(sbcal, vertices);
//            }
//            w.write(sbcal.toString());
            w.flush();
            w.close();
        }
    }

    public static void main(String[] args) throws Exception {
        String detectorName = "HPS-PhysicsRun2019-v2-4pt5";
        if (args.length > 0) {
            detectorName = args[0];
        }
        System.out.println("Building "+detectorName);
        DetectorBuilder db = new DetectorBuilder(detectorName);

        List<DetectorPlane> planes = db.getTracker("topSlot");
        for (DetectorPlane p : planes) {
            System.out.println(p);
            System.out.println("");
        }

        db.drawDetector();
        db.archiveIt("test");
    }

    public Map<String, DetectorPlane> planeMap() {
        return planeMap;
    }

    public void archiveIt(String suffix) throws Exception {
        //TODO just work from the complete set of planes in planemap...
        FileOutputStream fos = new FileOutputStream(_detectorName + "_" + myDate() + "_" + suffix + ".txt");
        Writer w = new BufferedWriter(new OutputStreamWriter(fos, "Cp850"));
        w.write(_detectorName + "\n");
        for (Tracker t : Tracker.values()) {
            String s = t.trackerName();
            if (_debug) {
                System.out.println("Archiving " + _detectorName + " " + s + " : ");
            }
            List<DetectorPlane> planes = getTracker(s);
            for (DetectorPlane p : planes) {
                StringBuffer sb = new StringBuffer(p.id() + " " + p.name() + " ");
                sb.append("o ( " + p.origin().x() + " " + p.origin().y() + " " + p.origin().z() + " ) ");
                sb.append("u ( " + p.u().x() + " " + p.u().y() + " " + p.u().z() + " ) ");
                sb.append("v ( " + p.v().x() + " " + p.v().y() + " " + p.v().z() + " ) ");
                sb.append("w ( " + p.normal().x() + " " + p.normal().y() + " " + p.normal().z() + " ) ");
                sb.append("a ( " + p.rotationAngles()[0] + " " + p.rotationAngles()[1] + " " + p.rotationAngles()[2] + " ) ");
                sb.append(p.width() + " " + p.height() + " " + p.zmin() + " " + p.zmax() + " \n");
                w.write(sb.toString());
            }
        }
        w.flush();
        w.close();
    }

    static private String myDate() {
        Calendar cal = new GregorianCalendar();
        Date date = new Date();
        cal.setTime(date);
        DecimalFormat formatter = new DecimalFormat("00");
        String day = formatter.format(cal.get(Calendar.DAY_OF_MONTH));
        String month = formatter.format(cal.get(Calendar.MONTH) + 1);
        return cal.get(Calendar.YEAR) + month + day;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer(" Detector " + _detectorName);
        List<String> list = new ArrayList<>(planeMap.keySet());
        Collections.sort(list);
        for (String s : list) {
            System.out.println(planeMap.get(s));
        }
        return sb.toString();
    }
}
