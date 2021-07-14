package org.hps.users.ngraf.testing;

import java.util.List;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.util.Driver;

/**
 *
 * @author Norman A. Graf
 */
public class CheckHelicalTrackHitRelations extends Driver {

    @Override
    protected void process(EventHeader event) {
        boolean OK = true;
        List<LCRelation> hitrelations = event.get(LCRelation.class, "HelicalTrackHitRelations");
        for (LCRelation relation : hitrelations) {
            if (relation != null && relation.getFrom() != null && relation.getTo() != null) {
                System.out.println(" from " + relation.getFrom().getClass().getName());
                System.out.println("   to " + relation.getTo().getClass().getName());
            } else {
                OK = false;
            }
        }
        System.out.println("HelicalTrackHitRelations " + (OK ? " OK " : " NOT OK "));
    }
}
