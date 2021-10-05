package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.List;
import java.util.Map;
import org.hps.recon.tracking.TrackData;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.Vertex;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.fourvec.Lorentz4Vector;
import org.lcsim.util.fourvec.Momentum4Vector;

/**
 *
 * @author Norman A. Graf
 */
public class SimpleMollerAnalysis2021 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix _beamAxisRotation = new BasicHep3Matrix();

    private double _deltaTrackTimeCut = 10;
    private int _minNhitsOnTrack = 11;

    protected void detectorChanged(Detector detector) {
        _beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    @Override
    protected void process(EventHeader event) {
        List<Vertex> vertices = event.get(Vertex.class, "UnconstrainedMollerVertices_KF");
        if (!vertices.isEmpty()) {
            RelationalTable trackToData = getKFTrackDataRelations(event);
            aida.histogram1D("Number of Moller vertices", 10, -0.5, 9.5).fill(vertices.size());
            for (Vertex v : vertices) {
                ReconstructedParticle rp = v.getAssociatedParticle();
                double invMass = rp.getMass();
                List<ReconstructedParticle> parts = rp.getParticles();
                ReconstructedParticle rp1 = parts.get(0);
                ReconstructedParticle rp2 = parts.get(1);
                double p1 = rp1.getMomentum().magnitude();
                double p2 = rp2.getMomentum().magnitude();
                double psum = p1 + p2;
                Track t1 = rp1.getTracks().get(0);
                Track t2 = rp2.getTracks().get(0);
                double track1_time = ((GenericObject) trackToData.from(t1)).getFloatVal(0);
                double track2_time = ((GenericObject) trackToData.from(t2)).getFloatVal(0);
                double deltaTrackTime = track1_time - track2_time;
                int track1_nHits = t1.getTrackerHits().size();
                int track2_nHits = t2.getTrackerHits().size();
                aida.histogram1D("track delta time", 100, -20., 20.).fill(deltaTrackTime);
                aida.histogram1D("track1 nHits", 20, -0.5, 19.5).fill(track1_nHits);
                aida.histogram1D("track2 nHits", 20, -0.5, 19.5).fill(track2_nHits);
                aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                aida.histogram1D("invariant mass", 100, 0., 0.2).fill(invMass);
                aida.histogram2D("psum vs invariant mass", 100, 0., 7., 100, 0., 0.2).fill(psum, invMass);
                // add some selection cuts
                // timing and nhits cut
                if (track1_nHits >= _minNhitsOnTrack && track2_nHits >= _minNhitsOnTrack && abs(deltaTrackTime) < _deltaTrackTimeCut) {
                    aida.tree().mkdirs("after cuts");
                    aida.tree().cd("after cuts");
                    aida.histogram1D("track delta time", 100, -20., 20.).fill(deltaTrackTime);
                    aida.histogram1D("track1 nHits", 20, -0.5, 19.5).fill(track1_nHits);
                    aida.histogram1D("track2 nHits", 20, -0.5, 19.5).fill(track2_nHits);
                    aida.histogram1D("psum", 100, 0., 7.).fill(psum);
                    aida.histogram1D("invariant mass", 100, 0., 0.2).fill(invMass);
                    aida.histogram2D("psum vs invariant mass", 100, 0., 7., 100, 0., 0.2).fill(psum, invMass);
                    reAnalyzeMoller(v);
                    aida.tree().cd("..");
                }
            }
        }
    }

    private void reAnalyzeMoller(Vertex v) {
        aida.tree().mkdirs("momentum-corrected");
        aida.tree().cd("momentum-corrected");
        double emass = 0.000511;
        double emass2 = emass * emass;
        String[] names = {"X", "Y", "Z"};
        double[] p1 = new double[3];
        double[] p2 = new double[3];
        double[] pV = new double[3];
        Hep3Vector vertexPosition = v.getPosition();
        Map<String, Double> vals = v.getParameters();
        // System.out.println(vals);
        p1[0] = vals.get("p1X");
        p1[1] = vals.get("p1Y");
        p1[2] = vals.get("p1Z");
        p2[0] = vals.get("p2X");
        p2[1] = vals.get("p2Y");
        p2[2] = vals.get("p2Z");
//        double v0p = vals.get("V0P");
//        pV[0] = vals.get("V0Px");
//        pV[1] = vals.get("V0Py");
//        pV[2] = vals.get("V0Pz");
//        double e1 = 0;
//        double e2 = 0.;
//        for (int i = 0; i < 3; ++i) {
//            e1 += p1[i] * p1[i];
//            e2 += p2[i] * p2[i];
//        }
//        e1 = sqrt(e1 + emass2);
//        e2 = sqrt(e2 + emass2);
//        Momentum4Vector evec1 = new Momentum4Vector(p1[0], p1[1], p1[2], e1);
//        Momentum4Vector evec2 = new Momentum4Vector(p2[0], p2[1], p2[2], e2);
//        Lorentz4Vector eesum = evec1.plus(evec2);
//        double eemass = eesum.mass();
//        double invMass = vals.get("invMass");
//        aida.histogram1D("vertex invariant mass", 100, 0., 0.2).fill(invMass);
//        aida.histogram1D("vertex invariant mass e+e-", 100, 0., 0.2).fill(eemass);
//
//        double invMass2 = invMass(p1, p2);
//        aida.histogram1D("vertex invariant mass e+e- recalculated", 100, 0., 0.2).fill(invMass2);

        double[] p1corr = correctMomentum(p1);
        double[] p2corr = correctMomentum(p2);
        double invMasscorrected = invMass(p1corr, p2corr);
        aida.histogram1D("vertex invariant mass e+e- corrected", 100, 0., 0.2).fill(invMasscorrected);

        aida.histogram1D("vtx px corrected", 100, -0.3, 0.3).fill(p1corr[0] + p2corr[0]);
        aida.histogram1D("vtx py corrected", 100, -0.1, 0.1).fill(p1corr[1] + p2corr[1]);
        aida.histogram1D("vtx pz corrected", 100, 0., 5.).fill(p1corr[2] + p2corr[2]);

        //
        Hep3Vector p1corrmom = new BasicHep3Vector(p1corr);
        Hep3Vector p2corrmom = new BasicHep3Vector(p2corr);
        Hep3Vector p1corrmomrot = VecOp.mult(_beamAxisRotation, p1corrmom);
        Hep3Vector p2corrmomrot = VecOp.mult(_beamAxisRotation, p2corrmom);
        double theta1 = Math.acos(p1corrmomrot.z() / p1corrmomrot.magnitude());
        double theta2 = Math.acos(p2corrmomrot.z() / p2corrmomrot.magnitude());
        double thetasum = theta1 + theta2;

        aida.histogram1D("theta sum", 100, 0.02, 0.06).fill(thetasum);
        aida.histogram1D("p1", 100, 1., 3.).fill(p1corrmomrot.magnitude());
        aida.histogram1D("p2", 100, 1., 3.).fill(p2corrmomrot.magnitude());
        aida.histogram1D("theta1", 100, 0.01, 0.05).fill(theta1);
        aida.histogram1D("theta2", 100, 0.01, 0.05).fill(theta2);
        aida.histogram2D("p1 vs theta1", 100, 1., 3., 100, 0.01, 0.05).fill(p1corrmomrot.magnitude(), theta1);
        aida.histogram2D("p2 vs theta2", 100, 1., 3., 100, 0.01, 0.05).fill(p2corrmomrot.magnitude(), theta2);
        aida.histogram2D("p1 vs p2", 100, 1., 3., 100, 1., 3.).fill(p1corrmomrot.magnitude(), p2corrmomrot.magnitude());
        aida.histogram2D("theta1 vs theta2", 100, 0.01, 0.05, 100, 0.01, 0.05).fill(theta1, theta2);
        // let's look at top/bottom
        double psum = p1corrmom.magnitude() + p2corrmomrot.magnitude();
        aida.histogram1D("psum corrected", 100, 0., 7.).fill(psum);
        aida.histogram2D("psum corrected vs invariant mass corrected ", 100, 0., 7., 100, 0., 0.2).fill(psum, invMasscorrected);
                    

//        if (clus.getPosition()[1] > 0.) {
//            aida.histogram1D("matched track momentum top", 50, 1., 3.).fill(p);
//            aida.histogram1D("unmatched track momentum bottom", 50, 1., 3.).fill(p2corrmom.magnitude());
//            aida.histogram1D("theta1 top", 100, 0.01, 0.05).fill(theta1);
//            aida.histogram1D("theta2 bottom", 100, 0.01, 0.05).fill(theta2);
//
//            aida.histogram2D("p1 vs theta1 top", 100, 1., 3., 100, 0.01, 0.05).fill(p1corrmomrot.magnitude(), theta1);
//            aida.histogram2D("p2 vs theta2 bottom", 100, 1., 3., 100, 0.01, 0.05).fill(p2corrmomrot.magnitude(), theta2);
//        } else {
//            aida.histogram1D("matched track momentum bottom", 50, 1., 3.).fill(p);
//            aida.histogram1D("unmatched track momentum top", 50, 1., 3.).fill(p2corrmom.magnitude());
//            aida.histogram1D("theta1 bottom", 100, 0.01, 0.05).fill(theta1);
//            aida.histogram2D("p1 vs theta1 bottom", 100, 1., 3., 100, 0.01, 0.05).fill(p1corrmomrot.magnitude(), theta1);
//            aida.histogram2D("p2 vs theta2 top", 100, 1., 3., 100, 0.01, 0.05).fill(p2corrmomrot.magnitude(), theta2);
//        }
        aida.tree().cd("..");
    }

    private double invMass(double[] p1, double[] p2) {
        double emass = 0.000511;
        double emass2 = emass * emass;
        double e1 = 0;
        double e2 = 0.;
        for (int i = 0; i < 3; ++i) {
            e1 += p1[i] * p1[i];
            e2 += p2[i] * p2[i];
        }
        e1 = sqrt(e1 + emass2);
        e2 = sqrt(e2 + emass2);
        Momentum4Vector evec1 = new Momentum4Vector(p1[0], p1[1], p1[2], e1);
        Momentum4Vector evec2 = new Momentum4Vector(p2[0], p2[1], p2[2], e2);
        Lorentz4Vector eesum = evec1.plus(evec2);
        return eesum.mass();
    }

    private double[] correctMomentum(double[] p) {
        //20210923 using fiducial electron clusters from WABs
        double[] parTop = {1.0894, -0.069986};
        double[] parBottom = {1.0376, -0.10535};
        double[] pout = new double[3];
        double pmom = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        boolean isTop = p[1] > 0 ? true : false;
        double correctedTrackMomentum = 0.;
        if (isTop) {
            correctedTrackMomentum = (parTop[0] + parTop[1] * pmom) * pmom;
        } else {
            correctedTrackMomentum = (parBottom[0] + parBottom[1] * pmom) * pmom;
        }
        // correction for the ECal energy scale 3.74/3.5
        double eCorr = 1.07;
        correctedTrackMomentum *= eCorr;
        for (int i = 0; i < 3; ++i) {
            pout[i] = p[i] * correctedTrackMomentum / pmom;
        }
        return pout;
    }

    public RelationalTable getKFTrackDataRelations(EventHeader event) {
        List<TrackData> TrackData;
        RelationalTable trackToData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackRelations;
        TrackData trackdata;
        TrackData = event.get(TrackData.class, "KFTrackData");
        trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
        for (LCRelation relation : trackRelations) {
            if (relation != null && relation.getTo() != null) {
                trackToData.add(relation.getFrom(), relation.getTo());
            }
        }
        return trackToData;
    }

    public void setDeltaTrackTimeCut(double d) {
        _deltaTrackTimeCut = d;
    }

    public void setMinNhitsOnTrack(int i) {
        _minNhitsOnTrack = i;
    }

}
