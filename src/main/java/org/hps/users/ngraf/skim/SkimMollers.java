package org.hps.users.ngraf.skim;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.util.List;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimMollers extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;

    private String[] mollerCollectionNames = {"UnconstrainedMollerCandidates", "UnconstrainedMollerCandidates_KF"};
    private String[] mollerVertexCollectionNames = {"UnconstrainedMollerVertices", "UnconstrainedMollerVertices_KF"};

    protected void detectorChanged(Detector detector) {
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    protected void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;

        for (String mollerCollectionName : mollerCollectionNames) {
            if (event.hasCollection(ReconstructedParticle.class, mollerCollectionName)) {
                List<ReconstructedParticle> mollers = event.get(ReconstructedParticle.class, mollerCollectionName);
                if (!mollers.isEmpty()) {
                    skipEvent = false;
                }
            }
        }

        for (String mollerCollectionName : mollerVertexCollectionNames) {
            if (event.hasCollection(Vertex.class, mollerCollectionName)) {
                aida.tree().mkdirs("Moller analysis " + mollerCollectionName);
                aida.tree().cd("Moller analysis " + mollerCollectionName);
                List<Vertex> vertices = event.get(Vertex.class, mollerCollectionName);
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                for (Vertex v : vertices) {
                    ReconstructedParticle rp = v.getAssociatedParticle();
                    List<ReconstructedParticle> parts = rp.getParticles();
                    ReconstructedParticle rp1 = parts.get(0);
                    ReconstructedParticle rp2 = parts.get(1);
                    double p1 = rp1.getMomentum().magnitude();
                    double p2 = rp2.getMomentum().magnitude();
                    int nclusters = 0;
                    Cluster c1 = null;
                    Cluster c2 = null;
                    if (!rp1.getClusters().isEmpty()) {
                        nclusters++;
                        c1 = rp1.getClusters().get(0);
                        CalorimeterHit seed = c1.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double energy = c1.getEnergy();
                        aida.histogram2D("cluster 1 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram1D("cluster 1 energy", 100, 0., 3.0).fill(energy);
                    }
                    if (!rp2.getClusters().isEmpty()) {
                        nclusters++;
                        c2 = rp2.getClusters().get(0);
                        CalorimeterHit seed = c2.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double energy = c2.getEnergy();
                        aida.histogram2D("cluster 2 ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram1D("cluster 2 energy", 100, 0., 3.0).fill(energy);
                    }

                    //rotate into physics frame of reference
                    Hep3Vector rprot = VecOp.mult(beamAxisRotation, rp.getMomentum());

                    Hep3Vector p1rot = VecOp.mult(beamAxisRotation, rp1.getMomentum());
                    Hep3Vector p2rot = VecOp.mult(beamAxisRotation, rp2.getMomentum());
                    double theta1 = Math.acos(p1rot.z() / p1rot.magnitude());
                    double theta2 = Math.acos(p2rot.z() / p2rot.magnitude());
                    double thetasum = theta1 + theta2;

                    aida.histogram1D("invariant mass", 100, 0.02, 0.1).fill(rp.getMass());
                    aida.histogram1D("px", 100, -0.1, 0.1).fill(rprot.x());
                    aida.histogram1D("py", 100, -0.05, 0.05).fill(rprot.y());
                    aida.histogram1D("pz", 100, 2.5, 5.0).fill(rprot.z());
                    aida.histogram1D("p", 100, 2.5, 5.0).fill(rprot.magnitude());
                    aida.histogram1D("number of clusters", 3, -0.5, 2.5).fill(nclusters);

                    aida.histogram1D("theta sum", 100, 0.02, 0.06).fill(thetasum);
                    aida.histogram2D("p1 vs theta1", 100, 1.2, 3., 100, 0.012, 0.024).fill(p1, theta1);
                    aida.histogram2D("p2 vs theta2", 100, 1.2, 3., 100, 0.012, 0.024).fill(p2, theta2);
                    aida.histogram2D("p1 vs p2", 100, 1.2, 3., 100, 1.2, 3.).fill(p1, p2);
                    aida.histogram2D("theta1 vs theta2", 100, 0.012, 0.024, 100, 0.012, 0.024).fill(theta1, theta2);

                    aida.tree().mkdirs(nclusters + "clusters");
                    aida.tree().cd(nclusters + "clusters");
                    aida.histogram1D("invariant mass", 100, 0.02, 0.1).fill(rp.getMass());
                    aida.histogram1D("px", 100, -0.1, 0.1).fill(rprot.x());
                    aida.histogram1D("py", 100, -0.05, 0.05).fill(rprot.y());
                    aida.histogram1D("pz", 100, 2.5, 5.0).fill(rprot.z());
                    aida.histogram1D("p", 100, 2.5, 5.0).fill(rprot.magnitude());
                    aida.histogram1D("number of clusters", 3, -0.5, 2.5).fill(nclusters);

                    aida.histogram1D("theta sum", 100, 0.02, 0.06).fill(thetasum);
                    aida.histogram2D("p1 vs theta1", 100, 1.2, 3., 100, 0.012, 0.024).fill(p1, theta1);
                    aida.histogram2D("p2 vs theta2", 100, 1.2, 3., 100, 0.012, 0.024).fill(p2, theta2);
                    aida.histogram2D("p1 vs p2", 100, 1.2, 3., 100, 1.2, 3.).fill(p1, p2);
                    aida.histogram2D("theta1 vs theta2", 100, 0.012, 0.024, 100, 0.012, 0.024).fill(theta1, theta2);

                    aida.tree().cd("..");
                }
                aida.tree().cd("..");
            }
        }

        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }
}
