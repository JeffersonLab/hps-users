package org.hps.users.ngraf.alignment;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Formatter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.FittedRawTrackerHit;
import org.hps.recon.tracking.MaterialSupervisor;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.ILogicalVolume;
import org.lcsim.detector.IPhysicalVolume;
import org.lcsim.detector.IPhysicalVolumeContainer;
import org.lcsim.detector.IRotation3D;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.ITranslation3D;
import org.lcsim.detector.RotationPassiveXYZ;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.material.IMaterial;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.detector.tracker.silicon.SiStrips;
import org.lcsim.detector.tracker.silicon.SiTrackerIdentifierHelper;
import org.lcsim.detector.tracker.silicon.SiTrackerModule;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.recon.tracking.digitization.sisim.TrackerHitType;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author ngraf
 */
public class StraightTrackSvtAlignmentDriver extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    RelationalTable hitToStrips;
    RelationalTable hitToRotated;

    boolean _useWeights = true;
    private double _oneClusterErr = 1 / Math.sqrt(12);
    private double _twoClusterErr = 1 / 5.;
    private double _threeClusterErr = 1 / 3.;
    private double _fourClusterErr = 1 / 2.;
    private double _fiveClusterErr = 1;

    // let's store some geometry here...
    Map<String, double[]> sensorAngles = new ConcurrentSkipListMap<String, double[]>();
    Map<String, double[]> sensorShifts = new ConcurrentSkipListMap<String, double[]>();
    Map<String, ITransform3D> localToGlobalMap = new ConcurrentSkipListMap<String, ITransform3D>();
    Map<String, ITransform3D> globalToLocalMap = new ConcurrentSkipListMap<String, ITransform3D>();

    Map<String, IRotation3D> sensorRotations = new ConcurrentSkipListMap<String, IRotation3D>();
    Map<String, ITranslation3D> sensorTranslations = new ConcurrentSkipListMap<String, ITranslation3D>();

    Map<String, Double> uLocal = new ConcurrentSkipListMap<String, Double>();
    Map<String, Double> uSigLocal = new ConcurrentSkipListMap<String, Double>();

    Map<Integer, String> layerModules = new ConcurrentSkipListMap<Integer, String>();
    Formatter topEvents;
    Formatter bottomEvents;
    Formatter topExtrap;

    @Override
    protected void detectorChanged(Detector detector) {

        System.out.println(detector.getName());

        boolean debug = true;
        MaterialSupervisor materialManager = new MaterialSupervisor();
//        MultipleScattering scattering = new MultipleScattering(materialManager);
        materialManager.buildModel(detector);
        List<MaterialSupervisor.ScatteringDetectorVolume> stripPlanes = materialManager.getMaterialVolumes();
        //TODO replace these lists with a helper class.
        List<String> names = new ArrayList<String>();
        List<Hep3Vector> oList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> uList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> vList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> nList = new ArrayList<Hep3Vector>();
        List<Double> measDim = new ArrayList<Double>();
        List<Double> unmeasDim = new ArrayList<Double>();
        List<Boolean> isAxial = new ArrayList<Boolean>();
        List<IPhysicalVolume> volumeList = new ArrayList<IPhysicalVolume>();

        for (MaterialSupervisor.ScatteringDetectorVolume vol : stripPlanes) {
            MaterialSupervisor.SiStripPlane plane = (MaterialSupervisor.SiStripPlane) vol;
            if (debug) {
                System.out.println(plane.getName());
            }
            names.add(plane.getName());
            Hep3Vector oprime = CoordinateTransformations.transformVectorToDetector(plane.origin());
            Hep3Vector nprime = CoordinateTransformations.transformVectorToDetector(plane.normal());
            if (debug) {
                System.out.println(" origin: " + oprime);
            }
            if (debug) {
                System.out.println(" normal: " + nprime);
            }
            if (debug) {
                System.out.println(" Plane is: " + plane.getMeasuredDimension() + " x " + plane.getUnmeasuredDimension());
            }
            HpsSiSensor sensor = (HpsSiSensor) plane.getSensor();

//            if (debug) {
//                System.out.println(SvtUtils.getInstance().isAxial(sensor) ? "axial" : "stereo");
//            }
            Hep3Vector measDir = CoordinateTransformations.transformVectorToDetector(plane.getMeasuredCoordinate());
            if (debug) {
                System.out.println("measured coordinate:    " + measDir);
            }
            Hep3Vector unmeasDir = CoordinateTransformations.transformVectorToDetector(plane.getUnmeasuredCoordinate());
            if (debug) {
                System.out.println("unmeasured coordinate:   " + unmeasDir);
            }
            if (debug) {
                System.out.println("thickness: " + plane.getThickness() + " in X0: " + plane.getThicknessInRL());
            }
            SiTrackerModule module = (SiTrackerModule) plane.getSensor().getGeometry().getDetectorElement().getParent();
            IPhysicalVolume parent = module.getGeometry().getPhysicalVolume();
            IPhysicalVolumeContainer daughters = parent.getLogicalVolume().getDaughters();
            if (debug) {
                System.out.printf(" found %d daughters to SiTrackerModule\n", daughters.size());
            }
            for (IPhysicalVolume daughter : daughters) {
                volumeList.add(daughter);
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
                    System.out.printf(" x %f y %f z %f box\n", solid.getXHalfLength(), solid.getYHalfLength(), solid.getZHalfLength());
                }
                double halfThickness = solid.getZHalfLength();
            }
            oList.add(oprime);
            nList.add(nprime);
            uList.add(measDir);
            vList.add(unmeasDir);
            measDim.add(plane.getMeasuredDimension());
            unmeasDim.add(plane.getUnmeasuredDimension());
            isAxial.add(sensor.isAxial());
        }
        DecimalFormat df = new DecimalFormat("###.######");
//        for (int i = 0; i < oList.size(); ++i) {
//            Hep3Vector o = oList.get(i);
//            Hep3Vector n = nList.get(i);
//            Hep3Vector u = uList.get(i);
//            Hep3Vector v = vList.get(i);
//            double l = unmeasDim.get(i) / 2.;
//            double h = measDim.get(i) / 2.;
////            String planeType = isAxial.get(i) ? "axial" : "stereo";
////            System.out.println("//" + planeType);
//            System.out.println(names.get(i));
//            System.out.println(df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " " + df.format(n.x()) + " , " + df.format(n.y()) + " , " + df.format(n.z()));
//
//        }

        // let's test the rotations...
        for (int i = 10; i < 12; ++i) {
            Hep3Vector X = new BasicHep3Vector(1., 0., 0.);
            Hep3Vector Y = new BasicHep3Vector(0., 1., 0.);
            Hep3Vector Z = new BasicHep3Vector(0., 0., 1.);
            Hep3Vector n = nList.get(i);
            Hep3Vector u = uList.get(i);
            Hep3Vector v = vList.get(i);

            Hep3Vector calcNormal = VecOp.cross(v, u);
            System.out.println("calculated normal v x u " + calcNormal);
            IPhysicalVolume volume = volumeList.get(i);

            Vector3D vX = Vector3D.PLUS_I;
            Vector3D vY = Vector3D.PLUS_J;
            Vector3D vZ = Vector3D.PLUS_K;

            Vector3D vXprime = new Vector3D(v.x(), v.y(), v.z());
            Vector3D vYprime = new Vector3D(u.x(), u.y(), u.z());
//            Vector3D vZprime = new Vector3D(n.x(), n.y(), n.z());
            Vector3D vZprime = new Vector3D(calcNormal.x(), calcNormal.y(), calcNormal.z());

            System.out.println(names.get(i));
            System.out.println("vX " + vX);
            System.out.println("vXprime " + vXprime);
            System.out.println("vY " + vY);
            System.out.println("vYprime " + vYprime);
            System.out.println("vZ " + vZ);
            System.out.println("vZprime " + vZprime);

            Rotation xyVecRot = new Rotation(vX, vY, vXprime, vYprime);
            double[][] xyVecRotMat = xyVecRot.getMatrix();
//            for (int ii = 0; ii < 3; ++ii) {
//                System.out.println(xyVecRotMat[ii][0] + " " + xyVecRotMat[ii][1] + " " + xyVecRotMat[ii][2]);
//            }
            System.out.println("");
            Rotation xzVecRot = new Rotation(vX, vZ, vXprime, vZprime);
            double[][] xzVecRotMat = xzVecRot.getMatrix();
//            for (int ii = 0; ii < 3; ++ii) {
//                System.out.println(xzVecRotMat[ii][0] + " " + xzVecRotMat[ii][1] + " " + xzVecRotMat[ii][2]);
//            }
            System.out.println("");
            Rotation yzVecRot = new Rotation(vY, vZ, vYprime, vZprime);
            double[][] yzVecRotMat = yzVecRot.getMatrix();
//            for (int ii = 0; ii < 3; ++ii) {
//                System.out.println(yzVecRotMat[ii][0] + " " + yzVecRotMat[ii][1] + " " + yzVecRotMat[ii][2]);
//            }
//            System.out.println("");
            IRotation3D volRot = volume.getRotation();
//            System.out.println(volRot);

            // so to within the signs of the rotations I am getting consistent answers.
            // Now the big question is how this relates to the usual alpha, beta, gamma rotations....
            //let's experiment
            //
            double alpha = 0.; //rotation about X axis
            double beta = 0.0305; // 30.5mr rotation about y
            double gammaAxial = 0.; // rotation about z axis for axial layer
            double gammaStereo = 0.050; // rotation about z axis for stereo layer

            double gamma = gammaAxial;
            if (names.get(i).contains("stereo")) {
                gamma = gammaStereo;
            }
//            Rotation testRot = new Rotation(RotationOrder.XYX, RotationConvention.VECTOR_OPERATOR, alpha, beta, gamma);
            Rotation testRot = new Rotation(RotationOrder.XYX, alpha, beta, gamma);
//            System.out.println("testing...");

            double[][] testRotMat = testRot.getMatrix();
            System.out.println("testRotMat");
            for (int ii = 0; ii < 3; ++ii) {
                System.out.println(testRotMat[ii][0] + " " + testRotMat[ii][1] + " " + testRotMat[ii][2]);
            }
            System.out.println("testRotMat angles");
//            double[] angles = testRot.getAngles(RotationOrder.XYX, RotationConvention.VECTOR_OPERATOR);
            double[] angles = testRot.getAngles(RotationOrder.XYX);
            System.out.println(Arrays.toString(angles));

            // so it appears that xzRot is the correct one to use.
            // test it.
            System.out.println("");
            System.out.println("xzVecRotMat");
            for (int ii = 0; ii < 3; ++ii) {
                System.out.println(xzVecRotMat[ii][0] + " " + xzVecRotMat[ii][1] + " " + xzVecRotMat[ii][2]);
            }
            System.out.println("xzVecRotMat angles");
//            double[] hpsAngles = xzVecRot.getAngles(RotationOrder.XYX, RotationConvention.VECTOR_OPERATOR);
            double[] hpsAngles = xzVecRot.getAngles(RotationOrder.XYX);
            System.out.println(Arrays.toString(hpsAngles));

            //hmmm, seems to work OK for the bottom axial, but not so much for the stereo
        }

    }

    protected void process(EventHeader event) {
        // only keep events with one and only one cluster
        List<Cluster> ecalClusters = event.get(Cluster.class, "EcalClustersCorr");
        if (ecalClusters.size() != 1) {
            return;
        }
        uLocal.clear();
        uSigLocal.clear();
        layerModules.clear();

        hitToStrips = TrackUtils.getHitToStripsTable(event);
        hitToRotated = TrackUtils.getHitToRotatedTable(event);

        setupSensors(event);
        List<TrackerHit> stripClusters = event.get(TrackerHit.class, "StripClusterer_SiTrackerHitStrip1D");
        // Get the list of fitted hits from the event
        List<LCRelation> fittedHits = event.get(LCRelation.class, "SVTFittedRawTrackerHits");
        // Map the fitted hits to their corresponding raw hits
        Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();
        for (LCRelation fittedHit : fittedHits) {
            fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
        }

        Map<Integer, double[]> globalPos = new ConcurrentSkipListMap<Integer, double[]>();
        Map<Integer, double[]> localPos = new ConcurrentSkipListMap<Integer, double[]>();

        Map<Integer, Hep3Vector> sensorOrigins = new ConcurrentSkipListMap<Integer, Hep3Vector>();
        Map<Integer, Hep3Vector> sensorNormals = new ConcurrentSkipListMap<Integer, Hep3Vector>();
        Map<Integer, String> sensorNames = new ConcurrentSkipListMap<Integer, String>();

        Cluster c = ecalClusters.get(0);
        double[] ecalClusterPos = c.getPosition();

        aida.cloud2D("Ecal cluster x vs y").fill(ecalClusterPos[0], ecalClusterPos[1]);

        List<Track> tracks = event.get(Track.class, "MatchedTracks");
        //System.out.println("found " + tracks.size() + " tracks ");
        aida.histogram1D("Number of Tracks", 10, 0., 10.).fill(tracks.size());
        for (Track t : tracks) {
            if (!isGoodTrack(t)) {
                continue;
            }
            plotTrack(t);

        }

        if (stripClusters.size() == 12) {
            aida.cloud2D("Ecal cluster x vs y 12-hit events").fill(ecalClusterPos[0], ecalClusterPos[1]);
            for (TrackerHit hit : stripClusters) {

                List rthList = hit.getRawHits();
                int size = rthList.size();
                double sense_pitch = 0;
                List<Double> signals = new ArrayList<Double>();
                List<Hep3Vector> positions = new ArrayList<Hep3Vector>();
                List<FittedRawTrackerHit> cluster = new ArrayList<FittedRawTrackerHit>();
                ITransform3D local_to_global = null;
                ITransform3D global_to_local = null;
                Hep3Vector sensorOrigin;
                Hep3Vector sensorNormal;
                for (int i = 0; i < size; ++i) {
                    RawTrackerHit rth = ((RawTrackerHit) rthList.get(i));
                    IIdentifier id = rth.getIdentifier();
                    SiSensor sensor = (SiSensor) rth.getDetectorElement();
                    local_to_global = sensor.getGeometry().getLocalToGlobal();
                    global_to_local = sensor.getGeometry().getGlobalToLocal();
                    SiTrackerIdentifierHelper _sid_helper = (SiTrackerIdentifierHelper) sensor.getIdentifierHelper();
                    ChargeCarrier carrier = ChargeCarrier.getCarrier(_sid_helper.getSideValue(id));
                    SiSensorElectrodes electrodes = ((SiSensor) rth.getDetectorElement()).getReadoutElectrodes(carrier);
                    sense_pitch = sensor.getSenseElectrodes(electrodes.getChargeCarrier()).getPitch(0);
                    Hep3Vector stripPosition = ((SiStrips) electrodes).getStripCenter(_sid_helper.getElectrodeValue(id));
                    double stripAmp = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rth));
                    signals.add(stripAmp);
                    positions.add(stripPosition);

                } // loop over strips in cluster
                Hep3Vector weightedPos = weightedAveragePosition(signals, positions);

                double measured_resolution;
                switch (size) {
                    case 1:
                        measured_resolution = sense_pitch * _oneClusterErr;
                        break;
                    case 2:
                        measured_resolution = sense_pitch * _twoClusterErr;
                        break;
                    case 3:
                        measured_resolution = sense_pitch * _threeClusterErr;
                        break;
                    case 4:
                        measured_resolution = sense_pitch * _fourClusterErr;
                        break;
                    default:
                        measured_resolution = sense_pitch * _fiveClusterErr;
                        break;
                }

                String moduleName = ((RawTrackerHit) rthList.get(0)).getDetectorElement().getName();
                int layer = TrackUtils.getLayer(hit);
                globalPos.put(layer, hit.getPosition());
                localPos.put(layer, weightedPos.v());
                layerModules.put(layer, moduleName);
                boolean isTop = moduleName.contains("t_halfmodule") ? true : false;
                uLocal.put(moduleName, weightedPos.x());
                uSigLocal.put(moduleName, 1. / (measured_resolution * measured_resolution));
                SiTrackerHitStrip1D stripHit = new SiTrackerHitStrip1D(hit);
                //
                sensorOrigin = getOrigin(stripHit);
                sensorNormal = getNormal(stripHit);

                sensorOrigins.put(layer, sensorOrigin);
                sensorNormals.put(layer, sensorNormal);
                sensorNames.put(layer, moduleName);
                localToGlobalMap.put(moduleName, local_to_global);
                globalToLocalMap.put(moduleName, global_to_local);

                //
                Hep3Vector uMeas = stripHit.getMeasuredCoordinate();
                Hep3Vector vMeas = stripHit.getUnmeasuredCoordinate();
                Hep3Vector pos = stripHit.getPositionAsVector();
                Hep3Vector calcNormal = VecOp.cross(vMeas, uMeas);
                Hep3Vector transformed_ltg = local_to_global.transformed(weightedPos);
                Hep3Vector transformed_gtl = global_to_local.transformed(pos);
                IRotation3D gtl_rot = global_to_local.getRotation();
                ITranslation3D gtl_trans = global_to_local.getTranslation();
                IRotation3D ltg_rot = local_to_global.getRotation();
                ITranslation3D ltg_trans = local_to_global.getTranslation();
                Vector3D vX = Vector3D.PLUS_I;
                Vector3D vY = Vector3D.PLUS_J;
                //Vector3D vZ = Vector3D.PLUS_K;

                Vector3D vXprime = new Vector3D(uMeas.x(), uMeas.y(), uMeas.z());
                Vector3D vYprime = new Vector3D(vMeas.x(), vMeas.y(), vMeas.z());
                // create a rotation matrix from this pair of vectors
                Rotation xyVecRot = new Rotation(vX, vY, vXprime, vYprime);
//            double[] hpsAngles = xyVecRot.getAngles(RotationOrder.XYZ, RotationConvention.VECTOR_OPERATOR);
                double[] hpsAngles = xyVecRot.getAngles(RotationOrder.XYZ);
                sensorAngles.put(moduleName, hpsAngles);
                sensorShifts.put(moduleName, ltg_trans.getTranslationVector().v());
                sensorRotations.put(moduleName, ltg_rot);
                sensorTranslations.put(moduleName, ltg_trans);

                boolean debug = false;
                if (debug) {
                    System.out.println("measured_resolution " + measured_resolution);
                    System.out.println(moduleName);
                    System.out.println("layer: " + layer);
                    System.out.println("u: " + uMeas);
                    System.out.println("v: " + vMeas);
                    System.out.println("calculated normal v x u " + calcNormal);
                    System.out.println("weighted pos " + weightedPos);
                    System.out.println("transformed ltg" + transformed_ltg);
                    System.out.println("pos: " + pos);
                    System.out.println("transformed gtl" + transformed_gtl);
                    System.out.println("gtl_rot " + gtl_rot);
                    System.out.println("gtl_trans " + gtl_trans);
                    System.out.println("ltg_rot " + ltg_rot);
                    System.out.println("ltg_trans " + ltg_trans);
//do some testing here...
                    Hep3Vector X = new BasicHep3Vector(1., 0., 0.); // this is local u
                    Hep3Vector Y = new BasicHep3Vector(0., 1., 0.); // this is local v
                    Hep3Vector Z = new BasicHep3Vector(0., 0., 1.); // this is local z

                    double[][] xyVecRotMat = xyVecRot.getMatrix();
                    System.out.println("Apache commons rotation:");
                    for (int ii = 0; ii < 3; ++ii) {
                        System.out.println(xyVecRotMat[ii][0] + " " + xyVecRotMat[ii][1] + " " + xyVecRotMat[ii][2]);
                    }
                    System.out.println("Apache commons angles");
                    System.out.println(Arrays.toString(hpsAngles));

                    double alpha = hpsAngles[0];
                    double beta = hpsAngles[1];
                    double gamma = hpsAngles[2];
                    IRotation3D tstRotPassive = new RotationPassiveXYZ(alpha, beta, gamma);
                    // this equals local_to_global
                    System.out.println("rotPassive " + tstRotPassive);
                    IRotation3D tstRotPassiveInv = tstRotPassive.inverse();
                    System.out.println("rotPassiveInv " + tstRotPassiveInv);
                } //debug output
            } // loop over strip clusters

//            for (Map.Entry<String, Double> entry : uLocal.entrySet()) {
//                aida.cloud1D(entry.getKey() + " u Local").fill(entry.getValue());
//                System.out.println(entry);
//            }
            for (int i = 0; i < 12; ++i) {
                String module = layerModules.get(i + 1);
                //System.out.println(module);
                aida.histogram1D(module + " u Local", 41, -20., 20.).fill(uLocal.get(module));
                //System.out.println("layer "+(i+1)+" "+module+uLocal.get(module));
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

    public Hep3Vector weightedAveragePosition(List<Double> signals, List<Hep3Vector> positions) {
        double total_weight = 0;
        Hep3Vector position = new BasicHep3Vector(0, 0, 0);
        for (int istrip = 0; istrip < signals.size(); istrip++) {
            double signal = signals.get(istrip);

            double weight = _useWeights ? signal : 1;
            total_weight += weight;
            position = VecOp.add(position, VecOp.mult(weight, positions.get(istrip)));
            /*if (_debug) {
                System.out.println(this.getClass().getSimpleName() + "strip " + istrip + ": signal " + signal + " position " + positions.get(istrip) + " -> total_position " + position.toString() + " ( total charge " + total_charge + ")");
            }*/
        }
        return VecOp.mult(1 / total_weight, position);
    }

    static Hep3Vector getOrigin(SiTrackerHitStrip1D stripCluster) {
        SiTrackerHitStrip1D local = stripCluster.getTransformedHit(TrackerHitType.CoordinateSystem.SENSOR);
        ITransform3D trans = local.getLocalToGlobal();
        return trans.transformed(new BasicHep3Vector(0, 0, 0));
    }

    static Hep3Vector getNormal(SiTrackerHitStrip1D s2) {
        Hep3Vector u2 = s2.getMeasuredCoordinate();
        Hep3Vector v2 = s2.getUnmeasuredCoordinate();
        return VecOp.cross(u2, v2);
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

    private void plotTrack(Track t) {
        Map<Integer, TrackerHit> orderedHits = new ConcurrentSkipListMap<>();
        // in principle, tracks with multi-strip hits are better measured...
        // 1st axial layer has greatest influence on theta, so require 2 strips in hit
        // TODO should I also require 2 strips in stereo layers?
        int tL1AxialNstrips = 0;
        int tL1StereoNstrips = 0;
        int tL2AxialNstrips = 0;
        int tL2StereoNstrips = 0;
        int tL1AxialStripNumber = 0;
        int tL1StereoStripNumber = 0;
        int tL2AxialStripNumber = 0;
        int tL2StereoStripNumber = 0;
        TrackState ts = t.getTrackStates().get(0);
        double d0 = ts.getD0();
        double z0 = ts.getZ0();
        double tanL = ts.getTanLambda();

        String half = isTopTrack(t) ? "top" : "bottom";

        for (TrackerHit hit : TrackUtils.getStripHits(t, hitToStrips, hitToRotated)) {
            int layer = TrackUtils.getLayer(hit);
            orderedHits.put(layer, hit);
//            System.out.println(layer);
            List rthList = hit.getRawHits();
            String moduleName = ((RawTrackerHit) rthList.get(0)).getDetectorElement().getName();
            if (moduleName.contains("module_L1")) {
                if (moduleName.contains("axial")) {
                    tL1AxialNstrips = rthList.size();
                    if (rthList.size() == 1) // look at single strip clusters
                    {
                        tL1AxialStripNumber = ((RawTrackerHit) hit.getRawHits().get(0)).getIdentifierFieldValue("strip");
                        aida.histogram1D(moduleName + "single strip cluster strip number", 640, 0., 640.).fill(tL1AxialStripNumber);
                    }
                }
                if (moduleName.contains("stereo")) {
                    tL1StereoNstrips = rthList.size();
                    if (rthList.size() == 1) // look at single strip clusters
                    {
                        tL1StereoStripNumber = ((RawTrackerHit) hit.getRawHits().get(0)).getIdentifierFieldValue("strip");
                        aida.histogram1D(moduleName + "single strip cluster strip number", 640, 0., 640.).fill(tL1StereoStripNumber);
                    }
                }
            }
            if (moduleName.contains("module_L2")) {
                if (moduleName.contains("axial")) {
                    tL2AxialNstrips = rthList.size();
                    if (rthList.size() == 1) // look at single strip clusters
                    {
                        tL2AxialStripNumber = ((RawTrackerHit) hit.getRawHits().get(0)).getIdentifierFieldValue("strip");
                        aida.histogram1D(moduleName + "single strip cluster strip number", 640, 0., 640.).fill(tL2AxialStripNumber);
                    }
                }
                if (moduleName.contains("stereo")) {
                    tL2StereoNstrips = rthList.size();
                    if (rthList.size() == 1) // look at single strip clusters
                    {
                        tL2StereoStripNumber = ((RawTrackerHit) hit.getRawHits().get(0)).getIdentifierFieldValue("strip");
                        aida.histogram1D(moduleName + "single strip cluster strip number", 640, 0., 640.).fill(tL2StereoStripNumber);
                    }
                }
            }
        } // end of hits on track
        //System.out.println(orderedHits.size());
    }

    private boolean isGoodTrack(Track t) {
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
            return true;
        }
        return false;
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
}
