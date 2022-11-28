package org.hps.users.ngraf.skim.golden;

import static java.lang.Math.abs;
import java.util.List;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimGoldenFee2016 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed;

    private int _maxNClusters = 1;
    private double _minClusterEnergy = 2.0;
    private double _minSeedHitEnergy = 0.5;
    private int _minHitsOnTrack = 11;
    private boolean _requireFiducialFee = true;

    String[] ReconstructedParticleCollectionNames = {"FinalStateParticles_KF", "OtherElectrons"};

    protected void process(EventHeader event) {
        boolean skipEvent = true;
        boolean isGoldenFee = false;
        _numberOfEventsProcessed++;
        for (String rpCollectionName : ReconstructedParticleCollectionNames) {
            List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, rpCollectionName);
            if (rpList.size() <= 2) {
                for (ReconstructedParticle rp : rpList) {
                    int pdgId = rp.getParticleIDUsed().getPDG();
                    if (pdgId == 11) {
                        if (!rp.getClusters().isEmpty()) {
                            Cluster c = rp.getClusters().get(0);
                            Track t = rp.getTracks().get(0);
                            String torb = c.getPosition()[1] > 0 ? " top " : " bottom ";
                            if (rp.getEnergy() > _minClusterEnergy && rp.getMomentum().magnitude() > _minClusterEnergy) {
                                if (t.getTrackerHits().size() > _minHitsOnTrack) {
                                    isGoldenFee = true;
                                    aida.histogram1D("track momentum" + torb, 100, 0., 3.).fill(rp.getMomentum().magnitude());
                                    aida.histogram1D("cluster energy" + torb, 100, 0., 3.).fill(c.getEnergy());
                                    aida.histogram1D("track nHits" + torb, 20, 0., 20.).fill(t.getTrackerHits().size());
                                    CalorimeterHit seed = c.getCalorimeterHits().get(0);
                                    int ix = seed.getIdentifierFieldValue("ix");
                                    int iy = seed.getIdentifierFieldValue("iy");
                                    aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                                    aida.histogram2D("event cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(c.getPosition()[0], c.getPosition()[1]);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (isGoldenFee) {
            skipEvent = false;
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    @Override
    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " golden FEEs of " + _numberOfEventsProcessed + "  events processed");
    }

    public void setMaxNClusters(int i) {
        _maxNClusters = i;
    }

    public void setMinHitsOnTrack(int i) {
        _minHitsOnTrack = i;
    }

    public void setMinClusterEnergy(double d) {
        _minClusterEnergy = d;
    }

    public void setMinSeedHitEnergy(double d) {
        _minSeedHitEnergy = d;
    }

    public void setRequireFiducialFee(boolean b) {
        _requireFiducialFee = b;
    }
}
