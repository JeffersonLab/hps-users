package org.hps.users.ngraf.alignment.straighttrack;

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
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.hps.record.triggerbank.TriggerModule;
import org.hps.users.ngraf.alignment.straighttrack.vertex.Track;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.solids.LineSegment3D;
import org.lcsim.detector.solids.Point3D;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.SimTrackerHit;

import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.math.chisq.ChisqProb;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author ngraf
 */
public class StraightTrackAnalysisDriver extends Driver {

    DetectorBuilder _alignedDetector;
    DetectorBuilder _defaultDetector;
    String _detectorName;
    double _minClusterEnergy = 3.5;

    private AIDA aida = AIDA.defaultInstance();
    // z from blueprints
    // x from fit (was -63? in 2016)
    private double[] H02Wire = {-68.0, 0., -(672.71 - 583.44) * 25.4};
    // initial guess for (x,y,z) of track origin
    // TODO get estimate for x of beam on wire. Was x=-63 in 2016
    double[] A0 = {0., 0., -2267.};
    // initial guess for the track direction
    double[] B0 = {0., 0., 1.};

    // try some vertexing here...
    List<Track> topAndBottomTracks = new ArrayList<Track>();
    int nTracksToVertex = 50;
    double target_x = -68.0;
    double target_y = 0.;
    double target_z = -(672.71 - 583.44) * 25.4; // -2267.458

    double maxChisq = 200.;

    Map<String, SimTrackerHit> simTrackerHitByModule = new TreeMap<String, SimTrackerHit>();
    Map<String, SiTrackerHitStrip1D> stripTrackerHitByModule = new TreeMap<String, SiTrackerHitStrip1D>();
    Map<String, Hit> stripHitByModule = new TreeMap<String, Hit>();
    Map<String, Integer> moduleFlipMap = new HashMap<>();

//    Writer topEventsWriter = null;
//    Writer bottomEventsWriter = null;
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
            "module_L6t_halfmodule_stereo_hole_sensor0");

    //runs 10099 and 10101 have all 10 good layers
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

    protected void detectorChanged(Detector detector) {
        Path resourcePath = null;
        try {
            resourcePath = Paths.get(getClass().getClassLoader().getResource("org/hps/analysis/alignment/HPS-PhysicsRun2019-v1-4pt5_20200205_topAlignment250000EventsIteration_9.txt").toURI());
        } catch (URISyntaxException ex) {
            Logger.getLogger(StraightTrackAnalysisDriver.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println(resourcePath);
        _alignedDetector = new DetectorBuilder(resourcePath);
        _defaultDetector = new DetectorBuilder(detector);
    }

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        setupSensors(event);
        // Will use calorimeter clusters to define the road in which we search for hits.
        List<Cluster> clusters = null;
        if (event.hasCollection(Cluster.class, "EcalClustersCorr")) {
            clusters = event.get(Cluster.class, "EcalClustersCorr");
        } else if (event.hasCollection(Cluster.class, "EcalClusters")) {
            clusters = event.get(Cluster.class, "EcalClusters");
        }
        aida.histogram1D("number of clusters", 10, 0., 10.).fill(clusters.size());
//        Cluster c = null;
        for (Cluster cluster : clusters) {
            aida.histogram1D("All Cluster energy", 100, 0., 7.).fill(cluster.getEnergy());
            aida.histogram2D("All Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            boolean isFiducial = TriggerModule.inFiducialRegion(cluster);
            String fid = isFiducial ? "fiducial" : "";
            aida.histogram1D("All Cluster " + fid + " energy", 100, 0., 7.).fill(cluster.getEnergy());
            aida.histogram2D("All Cluster " + fid + " x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            if (cluster.getPosition()[1] > 0) {
                aida.histogram1D("All Cluster " + fid + " x top", 320, -270.0, 370.0).fill(cluster.getPosition()[0]);
            } else {
                aida.histogram1D("All Cluster " + fid + " x bottom", 320, -270.0, 370.0).fill(cluster.getPosition()[0]);
            }
        }
        for (Cluster cluster : clusters) {
//            aida.histogram1D("All Cluster energy", 100, 0., 7.).fill(cluster.getEnergy());
            if (cluster.getEnergy() > _minClusterEnergy) {
                boolean isFiducial = TriggerModule.inFiducialRegion(cluster);
                String fid = isFiducial ? "fiducial" : "";
                Cluster c = cluster;
                String topOrBottom = c.getPosition()[1] > 0 ? "top " : "bottom ";
                boolean isTop = c.getPosition()[1] > 0;
                List<String> sensorNames = isTop ? topHoleSensorNamesToFit : bottomHoleSensorNamesToFit;
                aida.histogram2D("Cal " + fid + " Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                // calculate some slopes and intercepts...
                aida.histogram1D(topOrBottom + fid + "Cluster z", 10, 1443., 1445.).fill(cluster.getPosition()[2]);
                aida.histogram1D(topOrBottom + fid + "Cluster energy", 100, 0., 7.).fill(cluster.getEnergy());

//            }
//        }
//        if (c != null) {
                //OK, we have a good, high-energy cluster in tthe calorimeter...
                Point3D P0 = new Point3D(H02Wire[0], H02Wire[1], H02Wire[2]);
                double[] cPos = c.getPosition();
                Point3D P1 = new Point3D(cPos[0], cPos[1], cPos[2]);
                //let's get some SVT strip clusters...

//                // Get the list of fitted hits from the event
//                List<LCRelation> fittedHits = event.get(LCRelation.class,
//                        "SVTFittedRawTrackerHits");
//                // Map the fitted hits to their corresponding raw hits
//                Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();
//                for (LCRelation fittedHit : fittedHits) {
//                    fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
//
//                }
                List<SiTrackerHitStrip1D> stripClusters = event.get(SiTrackerHitStrip1D.class,
                        "StripClusterer_SiTrackerHitStrip1D");
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
                }
//            for (String s : hitsPerModuleMap.keySet()) {
//                System.out.println(s + " has " + hitsPerModuleMap.get(s).size() + " strip hit clusters");
//            }
                Map<String, SiTrackerHitStrip1D> hitsToFit = new LinkedHashMap<>();
                Map<String, DetectorPlane> detectorPlanesInFit = new LinkedHashMap<>();

                String trackingDetectorName = cPos[1] > 0 ? "topHole" : "bottomHole"; // work on slot later
                List<DetectorPlane> td = _alignedDetector.getTracker(trackingDetectorName);
                String[] trackerSensorNames = _alignedDetector.getTrackerSensorNames(trackingDetectorName);
                double maxDist = 5.;
                for (DetectorPlane dp : td) {
                    String moduleName = trackerSensorNames[dp.id() - 1];
//                System.out.println(moduleName);
                    if (hitsPerModuleMap.containsKey(moduleName)) {
//                    System.out.println(moduleName + " has " + hitsPerModuleMap.get(moduleName).size() + " strip hit clusters");

                        // get the best hit in this layer associated with this cluster                     
                        SiTrackerHitStrip1D closest = null;
                        double d = 9999.;
                        for (SiTrackerHitStrip1D stripHit : hitsPerModuleMap.get(moduleName)) {
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
                        // are we within a reasonable distance?
                        if (abs(d) < maxDist) {
                            hitsToFit.put(moduleName, closest);
                            detectorPlanesInFit.put(moduleName, dp);
//                            aida.histogram1D(moduleName + fid + " distance to hit", 100, -maxDist, maxDist).fill(d);
                        }
                    }
//                System.out.println(dp.id() + " " + trackerSensorNames[dp.id()-1]);
//                System.out.println(dp);
                }
                // we now have a list of hits to fit.
                List<Hit> hits = new ArrayList<Hit>();
                List<DetectorPlane> planes = new ArrayList<DetectorPlane>();
                //for now, assign an error based on the size of the strip cluster
                // unbiased residual plots indicate a working resolution of ~40 to 60um
                double[] fixedDu = {0., .012, .006};
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
                    System.out.println(detectorPlanesInFit.get(s).name());
                    System.out.println("strip position " + Arrays.toString(pos));
                    System.out.println("u " + u);
                    System.out.println("hStrip " + hStrip);
                    System.out.println("hPlane " + hPlane);
                    System.out.println("defaultHit " + defaultHit);
                    hits.add(hPlane);
                    planes.add(detectorPlanesInFit.get(s));
                }
                // require at least 8 hits for fit in top, 10 in bottom (for early runs), 9 in bottom in later runs due to missing layer 4
                int minHitsToFit = isTop ? 8 : 10;  // only 10101 has a working layer 4...
                if (hits.size() >= minHitsToFit) {
                    aida.histogram1D(topOrBottom + fid + " number of hits in fit", 20, 0., 20.).fill(hits.size());
                    // fit the track!
                    TrackFit fit = FitTracks.STR_LINFIT(planes, hits, A0, B0);
                    if (fit.chisq() / fit.ndf() < maxChisq) {
                        // quick check of track predicted impact points...
//                    List<double[]> impactPoints = fit.impactPoints();
//                    for (double[] pos : impactPoints) {
//                        System.out.println(Arrays.toString(pos));
//                    }
                        // calculate unbiased residuals here
                        refitTrack(planes, hits, A0, B0, isTop);
                        // Note that track position parameters x & y are reported at the input z.
                        double[] pars = fit.pars();
                        double[] cov = fit.cov();

                        aida.histogram1D(topOrBottom + fid + " x at z=-2267", 100, -100., 0.).fill(pars[0]);
                        aida.histogram1D(topOrBottom + fid + " y at z=-2267", 100, -20., 20.).fill(pars[1]);
                        aida.histogram1D(topOrBottom + fid + " dXdZ at z=-2267", 100, 0., 0.050).fill(pars[2]);
                        aida.histogram1D(topOrBottom + fid + " dYdZ at z=-2267", 100, -0.050, 0.050).fill(pars[3]);
                        aida.histogram1D(topOrBottom + fid + " track fit chiSquared per ndf", 100, 0., maxChisq).fill(fit.chisq() / fit.ndf());
                        double chisqProb = ChisqProb.gammp(fit.ndf(), fit.chisq());
                        aida.histogram1D(topOrBottom + fid + " track fit chiSquared probability", 100, 0., 1.).fill(chisqProb);
                        aida.histogram2D("Final Cal Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(c.getPosition()[0], c.getPosition()[1]);
                        aida.histogram1D("Final Cluster energy", 100, 0., 7.).fill(c.getEnergy());

                        // let's check the impact point at the calorimeter...
                        double[] tAtEcal = fit.predict(cPos[2]);
                        aida.histogram2D(fid + " Track at Ecal x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(tAtEcal[0], tAtEcal[1]);
                        aida.histogram1D(topOrBottom + fid + " dx at ECal", 100, -20., 20.).fill(cPos[0] - tAtEcal[0]);
                        aida.histogram1D(topOrBottom + fid + " dy at Ecal", 100, -10., 10.).fill(cPos[1] - tAtEcal[1]);
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
                        if (abs(pars[0] + 63) < 20) {
                            // beam spot y at wire is 0
                            if (abs(pars[1]) < 15) {
                                // keep this event
                                skipEvent = false;
                                // good clean event let's use it for alignment
                                if (isTop) {
//                                    if (topPlanes == null) {
//                                        topPlanes = new ArrayList<DetectorPlane>();
//                                        topPlanes.addAll(planes);
//                                    }
//                                    if (alignit) {
//                                        topEventsToAlign.add(hits);
//                                    }
////                                    topTracks.add(new Track(fit));
                                    topAndBottomTracks.add(new Track(fit));
                                } else {
//                                    if (bottomPlanes == null) {
//                                        bottomPlanes = new ArrayList<DetectorPlane>();
//                                        bottomPlanes.addAll(planes);
//                                    }
//                                    if (alignit) {
//                                        bottomEventsToAlign.add(hits);
//                                    }
//                                    bottomTracks.add(new Track(fit));
                                    topAndBottomTracks.add(new Track(fit));
                                }
                            }
                        }
                    }
                }
            } // end of loop over clusters with energy greater that min cluster energy
        } // end of loop over clusters        
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
}
