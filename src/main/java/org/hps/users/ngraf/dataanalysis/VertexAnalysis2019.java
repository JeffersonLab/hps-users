package org.hps.users.ngraf.dataanalysis;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.util.List;
import org.hps.recon.tracking.TrackType;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 *
 */
public class VertexAnalysis2019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    public void process(EventHeader event) {
        analyzeV0(event);
    }

    private void analyzeV0(EventHeader event) {
        String[] v0Dirs = {"UnconstrainedV0Candidates", "UnconstrainedV0Candidates_KF"};
        for (String dir : v0Dirs) {
            if (event.hasCollection(ReconstructedParticle.class, dir)) {
                aida.tree().mkdirs(dir);
                aida.tree().cd(dir);
                List<ReconstructedParticle> v0List = event.get(ReconstructedParticle.class, dir);
                aida.histogram1D("Number of V0s in event", 10, -0.5, 9.5).fill(v0List.size());
                for (ReconstructedParticle v0 : v0List) {
                    String trackType = "SeedTrack ";
                    int minNhits = 5;
                    if (TrackType.isGBL(v0.getType())) {
                        trackType = "GBL ";
                    }
                    if (v0.getType() == 1) {
                        trackType = "Kalman ";
                        minNhits = 10;
                    }
                    Vertex uncVert = v0.getStartVertex();
                    Hep3Vector pVtxRot = VecOp.mult(beamAxisRotation, v0.getMomentum());
                    Hep3Vector vtxPosRot = VecOp.mult(beamAxisRotation, uncVert.getPosition());
                    double theta = Math.acos(pVtxRot.z() / pVtxRot.magnitude());
                    double phi = Math.atan2(pVtxRot.y(), pVtxRot.x());

                    // this always has 2 tracks.
                    List<ReconstructedParticle> trks = v0.getParticles();
                    ReconstructedParticle ele = trks.get(0);
                    ReconstructedParticle pos = trks.get(1);
                    int eNhits = ele.getTracks().get(0).getTrackerHits().size();
                    int pNhits = pos.getTracks().get(0).getTrackerHits().size();
                    double eMom = ele.getMomentum().magnitude();
                    double pMom = pos.getMomentum().magnitude();
                    aida.histogram1D("electron track nHits", 20, 0., 20.).fill(eNhits);
                    aida.histogram1D("positron track nHits", 20, 0., 20.).fill(pNhits);
                    aida.histogram1D("electron momentum", 100, 0., 6.0).fill(eMom);
                    aida.histogram1D("positron momentum", 100, 0., 6.0).fill(pMom);

                    if (eNhits >= minNhits && pNhits >= minNhits) {
                        aida.histogram1D("v0 x", 50, -5., 5.).fill(vtxPosRot.x());
                        aida.histogram1D("v0 y", 50, -2., 2.).fill(vtxPosRot.y());
                        aida.histogram1D("v0 z", 50, -25., 0.).fill(vtxPosRot.z());
                        aida.histogram1D("v0 x ele " + eNhits + " pos " + pNhits + " hits on track", 50, -5., 5.).fill(vtxPosRot.x());
                        aida.histogram1D("v0 y ele " + eNhits + " pos " + pNhits + " hits on track", 50, -2., 2.).fill(vtxPosRot.y());
                        aida.histogram1D("v0 z ele " + eNhits + " pos " + pNhits + " hits on track", 50, -25., 0.).fill(vtxPosRot.z());
                        aida.histogram1D("v0 energy", 100, 0., 10.).fill(v0.getEnergy());
                        aida.histogram1D("v0 mass", 50, 0., 0.5).fill(v0.getMass());
                        aida.histogram2D("v0 mass vs Z vertex", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
                        aida.histogram2D("v0 mass vs Z vertex ele " + eNhits + " pos " + pNhits + " hits on track", 50, 0., 0.5, 100, -20., 20.).fill(v0.getMass(), vtxPosRot.z());
                        aida.profile1D("v0 mass vs Z vertex profile", 50, 0.05, 0.25).fill(v0.getMass(), vtxPosRot.z());
                        if (ele.getClusters().isEmpty()) {
                            aida.histogram1D("psum no electron ECal Cluster", 100, 0., 6.0).fill(eMom + pMom);
                            aida.histogram1D("psum no electron ECal Cluster ele " + eNhits + " pos " + pNhits + " hits on track", 100, 0., 6.0).fill(eMom + pMom);
                        }
                        aida.histogram1D("psum", 100, 0., 6.0).fill(eMom + pMom);
                        aida.histogram1D("psum ele " + eNhits + " pos " + pNhits + " hits on track", 100, 0., 6.0).fill(eMom + pMom);
                        aida.histogram2D("electron vs positron momentum", 100, 0., 6.0, 100, 0., 6.).fill(eMom, pMom);
                        if (ele.getClusters().size() > 0 && pos.getClusters().size() > 0) {
                            aida.histogram1D("psum both ECal Clusters", 100, 0., 6.0).fill(eMom + pMom);
                            aida.histogram1D("esum both ECal Clusters", 100, 0., 6.0).fill(ele.getClusters().get(0).getEnergy() + pos.getClusters().get(0).getEnergy());
                        }
                    }
                }
                aida.tree().cd("..");
            }
        }
    }
}
