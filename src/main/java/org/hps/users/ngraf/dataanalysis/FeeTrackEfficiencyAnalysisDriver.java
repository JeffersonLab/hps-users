package org.hps.users.ngraf.dataanalysis;

import java.util.List;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class FeeTrackEfficiencyAnalysisDriver extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles", "FinalStateParticles_KF"};

    protected void process(EventHeader event) {
        List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
        int nClusters = clusters.size();
        aida.histogram1D("number of Clusters", 5, -0.5, 4.5).fill(nClusters);
        for (Cluster cluster : clusters) {
            aida.histogram2D("Cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            if (cluster.getPosition()[1] > 0.) {
                aida.histogram1D("Top cluster energy ", 100, 3.5, 5.5).fill(cluster.getEnergy());
            } else {
                aida.histogram1D("Bottom cluster energy ", 100, 3.5, 5.5).fill(cluster.getEnergy());
            }
        }
    }
}
