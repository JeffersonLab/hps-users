package org.hps.users.ngraf.dataanalysis;

import java.util.BitSet;
import java.util.List;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TSData2019;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A> Graf
 */
public class TriggerAnalysis extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    public void process(EventHeader event) {
        boolean triggered = false;
        if (event.hasCollection(GenericObject.class, "TSBank")) {
            aida.histogram1D("triggers", 6, -0.5, 5.5).fill(0);
            List<GenericObject> triggerList = event.get(GenericObject.class, "TSBank");
            for (GenericObject data : triggerList) {
                if (AbstractIntData.getTag(data) == TSData2019.BANK_TAG) {
                    TSData2019 triggerData = new TSData2019(data);
                    BitSet bs = triggerData.getBitSet();
                    for (int i = 0; i < bs.length(); ++i) {
                        if (bs.get(i)) {
                            aida.histogram1D("trigger bits", 32, -0.5, 31.5).fill(i);
                        }
                    }
                    // pulser
                    if (triggerData.isPulserTrigger()) {
                        aida.histogram1D("triggers", 6, -0.5, 5.5).fill(1);
                    }
                    // faraday cup
                    if (triggerData.isFaradayCupTrigger()) {
                        aida.histogram1D("triggers", 6, -0.5, 5.5).fill(2);
                    }
                    //FEE
                    if (triggerData.isFEEBotTrigger() || triggerData.isFEETopTrigger()) {
                        aida.histogram1D("triggers", 6, -0.5, 5.5).fill(3);
                    }
//                    //Moller
//                    if (triggerData.isFEEBotTrigger() || triggerData.isFEETopTrigger()) {
//                        aida.histogram1D("triggers", 6, -0.5, 5.5).fill(1);
//                    }
                }
            }
        }
    }
}
