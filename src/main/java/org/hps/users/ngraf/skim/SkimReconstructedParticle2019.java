package org.hps.users.ngraf.skim;

import java.util.List;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;

/**
 * Utility class to skim off events with ReconstructedParticles
 *
 * @author Norman A. Graf
 */
public class SkimReconstructedParticle2019 extends Driver {

    private boolean _writeRunAndEventNumbers = false;
    private double _minMomentumCut = 0.5;
    private double _maxMomentumCut = 10.;
    private int _minNumberOfHitsOnTrack = 6;
    private int _maxNumberOfReconstructedParticles = 1;
    private int _pdgId = 11;
    private boolean _requireCluster = true;
    private int _numberOfEventsWritten = 0;
    private String _reconstructedParticleCollectionName = "FinalStateParticles";

    /**
     *
     * @param event
     */
    @Override
    protected void process(EventHeader event) {
        boolean skipEvent = true;
        if (event.hasCollection(ReconstructedParticle.class, _reconstructedParticleCollectionName)) {
            List<ReconstructedParticle> rps = event.get(ReconstructedParticle.class, _reconstructedParticleCollectionName);
            if (rps.size() <= _maxNumberOfReconstructedParticles) {
                for (ReconstructedParticle rp : rps) {
                    if (rp.getParticleIDUsed().getPDG() == _pdgId) {
                        double mom = rp.getMomentum().magnitude();
                        if (mom >= _minMomentumCut && mom <= _maxMomentumCut) {
                            skipEvent = false;
                            if (_requireCluster && rp.getClusters().isEmpty()) {
                                skipEvent = true;
                            }
                            if (!rp.getTracks().isEmpty()) {
                                if (rp.getTracks().get(0).getTrackerHits().size() < _minNumberOfHitsOnTrack) {
                                    skipEvent = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            if (_writeRunAndEventNumbers) {
                System.out.println(event.getRunNumber() + " " + event.getEventNumber());
            }
            _numberOfEventsWritten++;
        }
    }

    public void setRequireCluster(boolean b) {
        _requireCluster = b;
    }

    /**
     * Tracks having fewer than the number of hits will be rejected.
     *
     * @param i
     */
    public void setMinNumberOfHitsOnTrack(int i) {
        _minNumberOfHitsOnTrack = i;
    }

    /**
     * Events having more than the number of ReconstructedParticles will be
     * rejected.
     *
     * @param i
     */
    public void setMaxNumberOfReconstructedParticles(int i) {
        _maxNumberOfReconstructedParticles = i;
    }

    /**
     * Write out run and event numbers of events passing the cuts if desired
     *
     * @param b
     */
    public void setWriteRunAndEventNumbers(boolean b) {
        _writeRunAndEventNumbers = b;
    }

    public void setMinMomentumCut(double d) {
        _minMomentumCut = d;
    }

    public void setMaxMomentumCut(double d) {
        _maxMomentumCut = d;
    }

    public void setPdgId(int i) {
        _pdgId = i;
    }

    public void setReconstructedParticleCollectionName(String s) {
        _reconstructedParticleCollectionName = s;
    }

    @Override
    protected void endOfData() {
        System.out.println("Wrote " + _numberOfEventsWritten + " events");
    }

}
