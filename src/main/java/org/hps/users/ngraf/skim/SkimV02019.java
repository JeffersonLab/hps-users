package org.hps.users.ngraf.skim;

import java.util.List;
import org.lcsim.event.EventHeader;
import org.lcsim.event.Vertex;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Skim off events containing Vertex object(s) Loop over both GBL and KF
 * collections
 *
 * @author Norman A. Graf
 */
public class SkimV02019 extends Driver {

    private AIDA aida = AIDA.defaultInstance();
    private int _numberOfEventsSelected;
    private int _numberOfEventsProcessed = 0;
    boolean skipEvent = true;
    String[] vertexCollectionNames = {"UnconstrainedV0Vertices", "UnconstrainedV0Vertices_KF"};

    private double _trackDeltaTimeCut = 20.;

    public void process(EventHeader event) {
        skipEvent = true;
        _numberOfEventsProcessed++;
        for (String vertexCollectionName : vertexCollectionNames) {
            aida.tree().mkdirs(vertexCollectionName);
            aida.tree().cd(vertexCollectionName);
            if (event.hasCollection(Vertex.class, vertexCollectionName)) {
                List<Vertex> vertices = event.get(Vertex.class, vertexCollectionName);
                aida.histogram1D("number of vertices in event", 10, -0.5, 9.5).fill(vertices.size());
                if (vertices.size() > 0) {
                    skipEvent = false;
                }
            } // end of loop over check on collection
            aida.tree().cd("..");
        } // end of loop over vertex collections 
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

    public void setTrackDeltaTimeCut(double d) {
        _trackDeltaTimeCut = d;
    }
}
