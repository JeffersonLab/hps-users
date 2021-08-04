package org.hps.users.ngraf.dataanalysis;

import java.util.List;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TSData2019;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

public class Singles0EcalClusterAnalysis2019 extends Driver {

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    private AIDA aida = AIDA.defaultInstance();

    public void process(EventHeader event) {
        boolean skipEvent = true;
        _numberOfEventsProcessed++;
        boolean triggered = false;
        if (event.hasCollection(GenericObject.class, "TSBank")) {
            List<GenericObject> triggerList = event.get(GenericObject.class, "TSBank");
            for (GenericObject data : triggerList) {
                if (AbstractIntData.getTag(data) == TSData2019.BANK_TAG) {
                    TSData2019 triggerData = new TSData2019(data);
                    triggered = triggerData.isSingle0BotTrigger() || triggerData.isSingle0TopTrigger();
                }
            }
        }
        aida.histogram1D("Triggered singles0", 2, -0.5, 1.5).fill(triggered ? 1 : 0);
        if (triggered) {
            skipEvent = false;
            if (event.hasCollection(Cluster.class, "EcalClustersCorr")) {
                List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
                aida.histogram1D("Number of Clusters in triggered Event", 10, 0., 10.).fill(clusters.size());
                if (clusters.size() == 1) {
                    aida.tree().mkdirs("Singles0 single cluster");
                    aida.tree().cd("Singles0 single cluster");
                    for (Cluster cluster : clusters) {
                        CalorimeterHit seed = cluster.getCalorimeterHits().get(0);
                        int ix = seed.getIdentifierFieldValue("ix");
                        int iy = seed.getIdentifierFieldValue("iy");
                        double[] cPos = cluster.getPosition();
                        double energy = cluster.getEnergy();
                        aida.histogram2D("cluster ix vs iy", 47, -23.5, 23.5, 11, -5.5, 5.5).fill(ix, iy);
                        aida.histogram2D("cluster ix vs energy", 47, -23.5, 23.5, 100, 0., 5.0).fill(ix, energy);
                        aida.histogram2D("cluster x vs y", 320, -270.0, 370.0, 90, -90.0, 90.0).fill(cPos[0], cPos[1]);
                        String half = cluster.getPosition()[1] > 0. ? "Top " : "Bottom ";
                        aida.histogram1D(half + "Cluster Energy", 100, 0., 5.).fill(cluster.getEnergy());
                        aida.histogram1D("Cluster Energy", 100, 0., 5.0).fill(cluster.getEnergy());
                        //Mollers?
                        if (energy > 0.6 && ix < -10) {
                            aida.histogram1D("Cluster Energy signal region", 100, 0., 5.0).fill(cluster.getEnergy());
                        }
                    }
                    aida.tree().cd("..");
                }
            }
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsSelected++;
        }
    }

    protected void endOfData() {
        System.out.println("Selected " + _numberOfEventsSelected + " of " + _numberOfEventsProcessed + "  events processed");
    }
}
