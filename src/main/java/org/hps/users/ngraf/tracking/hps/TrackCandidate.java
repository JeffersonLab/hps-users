/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking.hps;

import java.util.List;

/**
 *
 * @author ngraf
 */
public class TrackCandidate {
    List<Hit> _hits;
    double[] _pars;
    public TrackCandidate(List<Hit> hits, double[] pars)
    {
        _hits = hits;
        _pars = pars;
    }
    
    public List<Hit> hits()
    {
        return _hits;
    }
    
    public double[] pars()
    {
        return _pars;
    }
}
