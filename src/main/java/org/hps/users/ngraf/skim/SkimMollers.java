package org.hps.users.ngraf.skim;

import java.util.List;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SkimMollers extends Driver {

    private AIDA aida = AIDA.defaultInstance();

    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;

    private String[] mollerCollectionNames = {"UnconstrainedMollerCandidates", "UnconstrainedMollerCandidates_KF"};

    protected void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;

        for (String mollerCollectionName : mollerCollectionNames) {
            List<ReconstructedParticle> mollers = event.get(ReconstructedParticle.class, mollerCollectionName);
            if (!mollers.isEmpty()) {
                skipEvent = false;
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
