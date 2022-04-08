package org.hps.users.ngraf.dataanalysis.fieldOff;

import Jama.Matrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import static java.lang.Math.signum;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.hps.users.ngraf.alignment.straighttrack.DetectorBuilder;
import org.hps.users.ngraf.alignment.straighttrack.DetectorPlane;
import org.hps.users.ngraf.alignment.straighttrack.FitTracks;
import org.hps.users.ngraf.alignment.straighttrack.Hit;
import org.hps.users.ngraf.alignment.straighttrack.ImpactPoint;
import org.hps.users.ngraf.alignment.straighttrack.StraightTrackAnalysisDriver;
import org.hps.users.ngraf.alignment.straighttrack.StraightTrackUtils;
import org.hps.users.ngraf.alignment.straighttrack.TrackFit;
import org.hps.users.ngraf.alignment.straighttrack.vertex.Track;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.solids.LineSegment3D;
import org.lcsim.detector.solids.Point3D;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.math.chisq.ChisqProb;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class FieldOffAnalysis extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;

    private DetectorBuilder _defaultDetector;
    private double _minClusterEnergy = 2.0;
    private double _maxChisq = 500.;
    private int _minHitsToFit = 8;
    // initial guess for (x,y,z) of track origin
    // TODO get estimate for x of beam on wire. Was x=-63 in 2016
    private double[] A0 = {0., 0., -2267.};
    // initial guess for the track direction
    private double[] B0 = {0., 0., 1.};

    //TODO check these numbers
    // z from blueprints
    // x from fit (was -63? in 2016)
    private double[] H02Wire = {-68.0, 0., -(672.71 - 583.44) * 25.4};

    // field-off data does not hit the first two layers
    // top is missing last layer
    List<String> topHoleSensorNamesToFit = Arrays.asList(
            "module_L3t_halfmodule_axial_sensor0",
            "module_L3t_halfmodule_stereo_sensor0",
            "module_L4t_halfmodule_axial_sensor0",
            "module_L4t_halfmodule_stereo_sensor0",
            "module_L5t_halfmodule_axial_hole_sensor0",
            "module_L5t_halfmodule_stereo_hole_sensor0",
            "module_L6t_halfmodule_axial_hole_sensor0",
            "module_L6t_halfmodule_stereo_hole_sensor0",
            "module_L7t_halfmodule_axial_hole_sensor0",
            "module_L7t_halfmodule_stereo_hole_sensor0");

    List<String> bottomHoleSensorNamesToFit = Arrays.asList(
            "module_L3b_halfmodule_stereo_sensor0",
            "module_L3b_halfmodule_axial_sensor0",
            "module_L4b_halfmodule_stereo_sensor0",
            "module_L4b_halfmodule_axial_sensor0",
            "module_L5b_halfmodule_stereo_hole_sensor0",
            "module_L5b_halfmodule_axial_hole_sensor0",
            "module_L6b_halfmodule_stereo_hole_sensor0",
            "module_L6b_halfmodule_axial_hole_sensor0",
            "module_L7b_halfmodule_stereo_hole_sensor0",
            "module_L7b_halfmodule_axial_hole_sensor0");

    List<String> topSlotSensorNamesToFit = Arrays.asList(
            "module_L3t_halfmodule_axial_sensor0",
            "module_L3t_halfmodule_stereo_sensor0",
            "module_L4t_halfmodule_axial_sensor0",
            "module_L4t_halfmodule_stereo_sensor0",
            "module_L5t_halfmodule_axial_slot_sensor0",
            "module_L5t_halfmodule_stereo_slot_sensor0",
            "module_L6t_halfmodule_axial_slot_sensor0",
            "module_L6t_halfmodule_stereo_slot_sensor0",
            "module_L7t_halfmodule_axial_slot_sensor0",
            "module_L7t_halfmodule_stereo_slot_sensor0");

    List<String> bottomSlotSensorNamesToFit = Arrays.asList(
            "module_L3b_halfmodule_stereo_sensor0",
            "module_L3b_halfmodule_axial_sensor0",
            "module_L4b_halfmodule_stereo_sensor0",
            "module_L4b_halfmodule_axial_sensor0",
            "module_L5b_halfmodule_stereo_slot_sensor0",
            "module_L5b_halfmodule_axial_slot_sensor0",
            "module_L6b_halfmodule_stereo_slot_sensor0",
            "module_L6b_halfmodule_axial_slot_sensor0",
            "module_L7b_halfmodule_stereo_slot_sensor0",
            "module_L7b_halfmodule_axial_slot_sensor0");

    protected void detectorChanged(Detector detector) {
        Path resourcePath = null;
        try {
            resourcePath = Paths.get(getClass().getClassLoader().getResource("org/hps/alignment/HPS_Run2021Pass1Top_20211110.txt").toURI());
        } catch (URISyntaxException ex) {
            Logger.getLogger(StraightTrackAnalysisDriver.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println(resourcePath);
        _defaultDetector = new DetectorBuilder(resourcePath);
    }

    public void process(EventHeader event) {
        boolean debug = false;
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        int eventNumber = event.getEventNumber();
        if (debug) {
            System.out.println(eventNumber + "");
        }
        if (event.hasCollection(Cluster.class, "EcalClustersCorr")) {
            List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
            aida.histogram1D("Number of Clusters in Event", 10, -0.5, 9.5).fill(clusters.size());
            Cluster clus = analyzeClusters(clusters);
            if (clus != null) {
                setupSensors(event);
                // run 14764 uses the HARP wire at z=-2270
                // run 14768 uses the collimator wire at z=-3080
                if (event.getRunNumber() == 14768) {
                    H02Wire[2] = -3080;
                    A0[2] = -3080.;
                }
                aida.histogram1D("Cluster Energy", 100, 0., 6.).fill(clus.getEnergy());
                if (debug) {
                    System.out.println(eventNumber + " found cluster with " + clus.getEnergy());
                }
                //OK, we have a good, high-energy cluster in tthe calorimeter...
                boolean isFiducial = TriggerModule.inFiducialRegion(clus);
                String fid = isFiducial ? "fiducial" : "";

                String topOrBottom = clus.getPosition()[1] > 0 ? "top " : "bottom ";
                boolean isTop = clus.getPosition()[1] > 0;
                List<String> sensorNames = isTop ? topHoleSensorNamesToFit : bottomHoleSensorNamesToFit;
//                List<String> sensorNames = isTop ? topSlotSensorNamesToFit : bottomSlotSensorNamesToFit;

                Point3D P0 = new Point3D(H02Wire[0], H02Wire[1], H02Wire[2]);
                double[] cPos = clus.getPosition();
                Point3D P1 = new Point3D(cPos[0], cPos[1], cPos[2]);
                List<SiTrackerHitStrip1D> stripClusters = event.get(SiTrackerHitStrip1D.class, "StripClusterer_SiTrackerHitStrip1D");
                int nStripHits = stripClusters.size();
                aida.histogram1D(topOrBottom + fid + "number of strip clusters", 100, 0., 200.).fill(nStripHits);
                // lets partition the strip clusters into each module
                Map<String, List<SiTrackerHitStrip1D>> hitsPerModuleMap = new HashMap<>();
                for (TrackerHit hit : stripClusters) {
                    List rthList = hit.getRawHits();
                    String moduleName = ((RawTrackerHit) rthList.get(0)).getDetectorElement().getName();
                    if (sensorNames.contains(moduleName)) {
                        if (!hitsPerModuleMap.containsKey(moduleName)) {
                            hitsPerModuleMap.put(moduleName, new ArrayList<SiTrackerHitStrip1D>());
                            hitsPerModuleMap.get(moduleName).add(new SiTrackerHitStrip1D(hit));
                        } else {
                            hitsPerModuleMap.get(moduleName).add(new SiTrackerHitStrip1D(hit));
                        }
                    }
                } // end of loop over strip clusters

                Map<String, SiTrackerHitStrip1D> hitsToFit = new LinkedHashMap<>();
                Map<String, DetectorPlane> detectorPlanesInFit = new LinkedHashMap<>();

                String trackingDetectorName = cPos[1] > 0 ? "topHole" : "bottomHole"; // work on slot later
                List<DetectorPlane> td = _defaultDetector.getTracker(trackingDetectorName);
                String[] trackerSensorNames = _defaultDetector.getTrackerSensorNames(trackingDetectorName);
                double maxDist = 5.0; //50.;//20.;
                for (DetectorPlane dp : td) {
                    String moduleName = trackerSensorNames[dp.id() - 1];
//                System.out.println(moduleName);
                    if (hitsPerModuleMap.containsKey(moduleName)) {
//                    System.out.println(moduleName + " has " + hitsPerModuleMap.get(moduleName).size() + " strip hit clusters");

                        // get the best hit in this layer associated with this cluster                     
                        SiTrackerHitStrip1D closest = null;
                        double d = 9999.;
                        for (SiTrackerHitStrip1D stripHit : hitsPerModuleMap.get(moduleName)) {
                            // are we in time?
                            double t = stripHit.getTime();
                            if (-20. < t && t < 0.) {
                                // calculate the intercept of the straight track with this sensor...
                                Hep3Vector intercept = StraightTrackUtils.linePlaneIntersect(P0, P1, dp.origin(), dp.normal());
                                // calculate the distance between this point and the strip
                                LineSegment3D stripLine = stripHit.getHitSegment();
                                double dist = stripLine.distanceTo(new Point3D(intercept));
//                        double d2 = VecOp.cross(stripLine.getDirection(),VecOp.sub(intercept,stripLine.getStartPoint())).magnitude();
//                        System.out.println("dist "+dist+" d2 "+d2);
                                if (abs(dist) < d) {
                                    d = dist;
                                    closest = stripHit;
                                }
                            }
                        }
                        aida.histogram1D(moduleName + fid + " distance to hit", 100, 0., 50.).fill(d);
                        // are we within a reasonable distance?
                        if (abs(d) < maxDist) {
                            hitsToFit.put(moduleName, closest);
                            detectorPlanesInFit.put(moduleName, dp);
                        }
                    }
//                System.out.println(dp.id() + " " + trackerSensorNames[dp.id()-1]);
//                System.out.println(dp);
                } // end of loop over detector planes
                // we now have a list of hits to fit.
                List<Hit> hits = new ArrayList<Hit>();
                List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
                //for now, assign an error based on the size of the strip cluster
                // unbiased residual plots indicate a working resolution of ~40 to 60um
                double[] fixedDu = {0., .012, .006};
                List<SiTrackerHitStrip1D> hitsInFit = new ArrayList<>();
                for (String s : hitsToFit.keySet()) {
                    SiTrackerHitStrip1D stripHit = hitsToFit.get(s);
//                System.out.println(s + " has a hit at " + stripHit.getPositionAsVector());
//TODO need to fix this local coordinate issue...
                    double u = ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement().getGeometry().getGlobalToLocal().transformed(new BasicHep3Vector(stripHit.getPosition())).x();
                    int size = stripHit.getRawHits().size();
                    double du;
                    if (size < 3) {
                        du = fixedDu[size];
                    } else {
                        du = .04;
                    }
                    double[] pos = stripHit.getPosition();
                    Hit hStrip = makeHit(detectorPlanesInFit.get(s), pos, du);
                    Hit hPlane = makeHit(detectorPlanesInFit.get(s), signum(u * hStrip.uvm()[0]) * u, du);
                    Hit defaultHit = makeHit(_defaultDetector.planeMap().get(s), pos, du);
//                    System.out.println(detectorPlanesInFit.get(s).name());
//                    System.out.println("strip position " + Arrays.toString(pos));
//                    System.out.println("u " + u);
//                    System.out.println("hStrip " + hStrip);
//                    System.out.println("hPlane " + hPlane);
//                    System.out.println("defaultHit " + defaultHit);
                    hits.add(hPlane);
                    planes.add(detectorPlanesInFit.get(s));
                    hitsInFit.add(stripHit);
                }
                aida.histogram1D(topOrBottom + fid + " number of hits to fit", 20, 0., 20.).fill(hits.size());
                if (debug) {
                    System.out.println(eventNumber + " has " + hits.size() + " SVT hits to fit");
                }
                if (hits.size() >= _minHitsToFit) {
                    analyzeHitsInFit(hitsInFit, clus);
                    aida.histogram1D(topOrBottom + fid + " number of hits in fit", 20, 0., 20.).fill(hits.size());
                    // fit the track!
                    TrackFit fit = FitTracks.STR_LINFIT(planes, hits, A0, B0);
                    if (debug) {
                        System.out.println(eventNumber + " track was fit with chisq/ndf " + fit.chisq() / fit.ndf());
                    }
                    aida.histogram1D(topOrBottom + fid + " fit chisq per ndf", 100, 0., 1000.).fill(fit.chisq() / fit.ndf());
//                    aida.cloud1D(topOrBottom + fid + " fit chisq per ndf cloud").fill(fit.chisq() / fit.ndf());

                    if (fit.chisq() / fit.ndf() < _maxChisq) {
                        // quick check of track predicted impact points...
//                    List<double[]> impactPoints = fit.impactPoints();
//                    for (double[] pos : impactPoints) {
//                        System.out.println(Arrays.toString(pos));
//                    }
                        if (debug) {
                            System.out.println(eventNumber + " refitting track");
                        }
                        // calculate unbiased residuals here
                        refitTrack(planes, hits, A0, B0, isTop);
                        // Note that track position parameters x & y are reported at the input z.
                        double[] pars = fit.pars();
                        double[] cov = fit.cov();

                        aida.histogram1D(topOrBottom + fid + " x at z= " + fit.zPosition(), 100, -100., 0.).fill(pars[0]);
                        aida.histogram1D(topOrBottom + fid + " y at z= " + fit.zPosition(), 100, -20., 20.).fill(pars[1]);
                        aida.histogram1D(topOrBottom + fid + " dXdZ at z= " + fit.zPosition(), 100, 0., 0.050).fill(pars[2]);
                        aida.histogram1D(topOrBottom + fid + " dYdZ at z= " + fit.zPosition(), 100, -0.050, 0.050).fill(pars[3]);
                        aida.histogram1D(topOrBottom + fid + " track fit chiSquared per ndf", 100, 0., _maxChisq).fill(fit.chisq() / fit.ndf());
                        double zMin = -3500.;
                        int nSteps = 250;//200;
                        int stepSize = 20;//10;
                        double zMax = zMin + nSteps * stepSize;
                        for (int i = 0; i < nSteps; ++i) {
                            double d = zMin + i * stepSize;
                            double[] fitPos = fit.predict(d);
                            aida.histogram2D(topOrBottom + fid + " z vs x", nSteps, zMin, zMax, 200, -100., 100.).fill(d, fitPos[0]);
                            aida.histogram2D(topOrBottom + fid + " z vs y", nSteps, zMin, zMax, 200, -100., 100.).fill(d, fitPos[1]);
                            aida.histogram2D(fid + " z vs x", nSteps, zMin, zMax, 200, -100., 100.).fill(d, fitPos[0]);
                            aida.histogram2D(fid + " z vs y", nSteps, zMin, zMax, 200, -100., 100.).fill(d, fitPos[1]);
                        }
                        double chisqProb = ChisqProb.gammp(fit.ndf(), fit.chisq());
                        aida.histogram1D(topOrBottom + fid + " track fit chiSquared probability", 100, 0., 1.).fill(chisqProb);
                        aida.histogram2D("Final Cal Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(clus.getPosition()[0], clus.getPosition()[1]);
                        aida.histogram1D("Final Cluster energy", 100, 0., 7.).fill(clus.getEnergy());

                        // let's check the impact point at the calorimeter...
                        double[] tAtEcal = fit.predict(cPos[2]);
                        aida.histogram2D(fid + " Track at Ecal x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(tAtEcal[0], tAtEcal[1]);
                        aida.histogram1D(topOrBottom + fid + " dx at ECal", 100, -20., 20.).fill(cPos[0] - tAtEcal[0]);
                        aida.histogram1D(topOrBottom + fid + " dy at Ecal", 100, -10., 10.).fill(cPos[1] - tAtEcal[1]);
                        // Have a good candidate...
                        skipEvent = false;
//                        if (beamConstrain) {
//                            aida.tree().mkdirs("beam-constrained");
//                            aida.tree().cd("beam-constrained");
//                            planes.add(xPlaneAtWire);
//                            hits.add(beamAtWire);
//                            planes.add(yPlaneAtWire);
//                            hits.add(beamAtWire);
//
//                            aida.histogram1D(topOrBottom + fid + " number of hits in fit", 20, 0., 20.).fill(hits.size());
//                            // fit the track!
//                            TrackFit fitbc = FitTracks.STR_LINFIT(planes, hits, A0, B0);
//                            // calculate unbiased residuals here
//                            refitTrack(planes, hits, A0, B0, isTop);
//                            // Note that track position parameters x & y are reported at the input z.
//                            double[] parsbc = fitbc.pars();
//                            double[] covbc = fitbc.cov();
//
//                            aida.histogram1D(topOrBottom + fid + " x at z=-2267", 100, -100., 0.).fill(parsbc[0]);
//                            aida.histogram1D(topOrBottom + fid + " y at z=-2267", 100, -20., 20.).fill(parsbc[1]);
//                            aida.histogram1D(topOrBottom + fid + " dXdZ at z=-2267", 100, 0., 0.050).fill(parsbc[2]);
//                            aida.histogram1D(topOrBottom + fid + " dYdZ at z=-2267", 100, -0.050, 0.050).fill(parsbc[3]);
//                            aida.histogram1D(topOrBottom + fid + " track fit chiSquared per ndf", 100, 0., 100.).fill(fitbc.chisq() / fitbc.ndf());
//                            double chisqProbbc = ChisqProb.gammp(fitbc.ndf(), fitbc.chisq());
//                            aida.histogram1D(topOrBottom + fid + " track fit chiSquared probability", 100, 0., 1.).fill(chisqProbbc);
//
//                            // let's check the impact point at the calorimeter...
//                            double[] tAtEcalbc = fitbc.predict(cPos[2]);
//                            aida.histogram2D(fid + " Track at Ecal x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(tAtEcalbc[0], tAtEcalbc[1]);
//                            aida.histogram1D(topOrBottom + fid + " dx at ECal", 100, -20., 20.).fill(cPos[0] - tAtEcalbc[0]);
//                            aida.histogram1D(topOrBottom + fid + " dy at Ecal", 100, -10., 10.).fill(cPos[1] - tAtEcalbc[1]);
//                            aida.tree().cd("..");
//                        }
                        // let's apply a few cuts here to enable us to skim events...
                        // beam spot x at wire is -63
//                        if (abs(pars[0] + 63) < 20) {
//                            // beam spot y at wire is 0
//                            if (abs(pars[1]) < 15) {
//                                // keep this event
//                                skipEvent = false;
//                                // good clean event let's use it for alignment
//                                if (isTop) {
////                                    if (topPlanes == null) {
////                                        topPlanes = new ArrayList<DetectorPlane>();
////                                        topPlanes.addAll(planes);
////                                    }
////                                    if (alignit) {
////                                        topEventsToAlign.add(hits);
////                                    }
//////                                    topTracks.add(new Track(fit));
//                                    topAndBottomTracks.add(new Track(fit));
//                                } else {
////                                    if (bottomPlanes == null) {
////                                        bottomPlanes = new ArrayList<DetectorPlane>();
////                                        bottomPlanes.addAll(planes);
////                                    }
////                                    if (alignit) {
////                                        bottomEventsToAlign.add(hits);
////                                    }
////                                    bottomTracks.add(new Track(fit));
//                                    topAndBottomTracks.add(new Track(fit));
//                                }
//                            }
//                        }
                    }
                }
            } // end of processing of good cluster
        } // end of check on cluster collection
        // did we find a good candidate?
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }

    /**
     * Given a DetectorPlane and a global position, return a hit in local
     * coordinates
     *
     * @param p
     * @param pos
     * @return
     */
    public Hit makeHit(DetectorPlane p, double[] pos, double du) {
        Matrix R = p.rot();
        double[] r0 = p.r0();
        Matrix diff = new Matrix(3, 1);
        for (int i = 0; i < 3; ++i) {
            diff.set(i, 0, pos[i] - r0[i]);
        }
//        diff.print(6, 4);
//        System.out.println("pos " + Arrays.toString(pos));
//        System.out.println("r0  " + Arrays.toString(r0));
        Matrix local = R.times(diff);
//        local.print(6, 4);
        double[] u = new double[2];  // 2dim for u and v measurement 
        double[] wt = new double[3]; // lower diag cov matrix
        double[] sigs = p.sigs();
        u[0] = local.get(0, 0);
        wt[0] = 1 / (du * du); //(sigs[0] * sigs[0]);

        return new Hit(u, wt);
    }

    public Hit makeHit(DetectorPlane p, double uPos, double du) {
        double[] u = new double[2];  // 2dim for u and v measurement 
        double[] wt = new double[3]; // lower diag cov matrix
        double[] sigs = p.sigs();
        u[0] = uPos;
        wt[0] = 1 / (du * du); //(sigs[0] * sigs[0]);

        return new Hit(u, wt);

    }

    private void setupSensors(EventHeader event) {
        List<RawTrackerHit> rawTrackerHits = null;
        if (event.hasCollection(RawTrackerHit.class,
                "SVTRawTrackerHits")) {
            rawTrackerHits = event.get(RawTrackerHit.class,
                    "SVTRawTrackerHits");
        }
        if (event.hasCollection(RawTrackerHit.class,
                "RawTrackerHitMaker_RawTrackerHits")) {
            rawTrackerHits = event.get(RawTrackerHit.class,
                    "RawTrackerHitMaker_RawTrackerHits");
        }
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

    private Cluster analyzeClusters(List<Cluster> clusters) {
        aida.tree().mkdirs("Cluster analysis");
        aida.tree().cd("Cluster analysis");
        double maxEnergy = 0;
        Cluster clus = null;
        for (Cluster cluster : clusters) {
            if (cluster.getEnergy() > maxEnergy) {
                maxEnergy = cluster.getEnergy();
                if (maxEnergy > _minClusterEnergy) {
                    clus = cluster;
                }
            }
            double[] cPos = cluster.getPosition();
            CalorimeterHit seedHit = ClusterUtilities.findSeedHit(cluster);
            boolean clusterIsFiducial = TriggerModule.inFiducialRegion(cluster);
            aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cPos[0], cPos[1]);
            int ix = seedHit.getIdentifierFieldValue("ix");
            int iy = seedHit.getIdentifierFieldValue("iy");
            aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
            String half = cluster.getPosition()[1] > 0. ? "Top " : "Bottom ";
            aida.histogram1D(half + "Cluster Energy", 100, 0., 5.).fill(cluster.getEnergy());
            aida.histogram1D("Cluster Energy", 100, 0., 5.0).fill(cluster.getEnergy());
        }
        aida.histogram1D("Max Cluster Energy", 100, 0., 6.).fill(maxEnergy);
        aida.tree().cd("..");
        return clus;
    }

    /**
     * Refit a track to obtain the unbiased hit residuals one plane at a time
     *
     * @param planes The planes used in the original fit
     * @param hits The hits on the track
     * @param A0 Initial guess for (x,y,z) of track
     * @param B0 Initial guess for the track direction
     * @param isTop true if track is in the top SVT
     */
    public void refitTrack(List<DetectorPlane> planes, List<Hit> hits, double[] A0, double[] B0, boolean isTop) {
        String topOrBottom = isTop ? "top " : "bottom ";
        String path = "Track Refit ";
        aida.tree().mkdirs(path);
        aida.tree().cd(path);
        // refit this track dropping one different hit each time...
        int nHitsOnTrack = planes.size();
        // loop over all the hits
        for (int i = 0; i < nHitsOnTrack; ++i) {
            List<DetectorPlane> newPlanes = new ArrayList<>();
            List<Hit> newHits = new ArrayList<>();
            // remove each hit and refit track without this hit
            for (int j = 0; j < nHitsOnTrack; ++j) {
                if (j != i) {
                    newPlanes.add(planes.get(j));
                    newHits.add(hits.get(j));
                }
            }
            // refit without the one hit
            TrackFit fit = FitTracks.STR_LINFIT(newPlanes, newHits, A0, B0);
            // get the predicted impact point for the missing hit...

            // get the hit...
            DetectorPlane missingPlane = planes.get(i);
            Hit missingHit = hits.get(i);
            // get the unbiased residual for the missing hit
            double[] resid = unbiasedResidual(fit, missingPlane, missingHit, A0, B0);
            aida.histogram1D(topOrBottom + "unbiased residual " + missingPlane.id(), 100, -1.0, 1.0).fill(resid[1]);
//            aida.histogram2D(topOrBottom + "unbiased x vs residual " + missingPlane.id(), 300, -200, 100, 100, -1.0, 1.0).fill(resid[0], resid[1]);
//            aida.histogram1D(topOrBottom + "unbiased x " + missingPlane.id(), 300, -200, 100).fill(resid[0]);
        }
        aida.tree().cd("..");
    }

    /**
     * Calculate the unbiased residual for a hit not included in the fit. Note
     * that because our axial/stereo pairs move in concert it might be necessary
     * to remove two hits from the fit...
     *
     * @param fit The track fit excluding a hit
     * @param dp The DetectorPlane for the excluded hit
     * @param h The excluded hit
     * @param A0 The TrackFit position
     * @param B0 The TrackFit direction
     * @return
     */
    public double[] unbiasedResidual(TrackFit fit, DetectorPlane dp, Hit h, double[] A0, double[] B0) {
        double resid = 9999.;
        boolean debug = false;
        Matrix rot = dp.rot();
        Matrix[] uvwg = new Matrix[3];
        double[][] UVW = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        double[] BUVW = new double[3];
        double[] PAR = fit.pars();
        double[] A = {PAR[0], PAR[1], A0[2]};
        double[] B = {PAR[2], PAR[3], B0[2]};
        Matrix b = new Matrix(B, 1);
        for (int j = 0; j < 3; ++j) {
//                    if (debug) {
//                        System.out.println("  CALLING VMATR");
//                    }
            Matrix uvw = new Matrix(UVW[j], 3);
            if (debug) {
                System.out.println("  UVW(" + (j + 1) + ") " + uvw.get(0, 0) + " " + uvw.get(1, 0) + " " + uvw.get(2, 0));
            }
            uvwg[j] = rot.transpose().times(uvw);
            if (debug) {
                System.out.println("UVWG(" + (j + 1) + ") " + uvwg[j].get(0, 0) + " " + uvwg[j].get(1, 0) + " " + uvwg[j].get(2, 0) + " ");
            }
//                    System.out.println("j "+j);
//                    System.out.println("b");
//                    b.print(6,4);
            BUVW[j] = b.times(uvwg[j]).get(0, 0);
        }
        if (debug) {
            System.out.println("   BUVW " + BUVW[0] + " " + BUVW[1] + " " + BUVW[2]);
        }
        ImpactPoint ip = FitTracks.GET_IMPACT(A, B, rot, dp.r0(), uvwg[2], BUVW[2]);
//        System.out.println(ip);
//        System.out.println(h);
        return new double[]{ip.q()[1], h.uvm()[0] - ip.q()[0]};
    }

    private void analyzeHitsInFit(List<SiTrackerHitStrip1D> hitsInFit, Cluster clus) {
        aida.tree().mkdirs("hits in fit analysis");
        aida.tree().cd("hits in fit analysis");
        double tClus = ClusterUtilities.findSeedHit(clus).getTime();
        String torb = clus.getPosition()[1] > 0 ? "top" : "bottom";
        aida.histogram1D("cluster time " + torb, 100, 30., 50.).fill(tClus);
        for (SiTrackerHitStrip1D hit : hitsInFit) {
            String sensorName = hit.getSensor().getName();
            double t = hit.getTime();
            aida.histogram1D(sensorName + " hit time", 100, -20., 0.).fill(t);
            aida.histogram1D(sensorName + " hit time - cluster time", 100, -60., -30.).fill(t - tClus);
//            aida.cloud1D(sensorName + " hit time cloud").fill(t);
//            aida.cloud1D(sensorName + " hit time - cluster time cloud").fill(t - tClus);
        }
        aida.tree().cd("..");
    }

    public void setMinHitsToFit(int i) {
        _minHitsToFit = i;
    }

    public void setMaxChisq(double d) {
        _maxChisq = d;
    }
}
