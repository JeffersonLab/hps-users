/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking;

import java.util.List;

/**
 *
 * DOC PAR = parameters. 
 * DOC COV = covariance matrix of PAR.
 * DOC RX = array of impact points. 
 * DOC CHI = chi squared of the fit. 
 * DOC NDF = number of degrees of freedom. 
 * DOC NIT = number of iterations.
 *
 * @author ngraf
 */
public class TrackFit {
    private double[] _par;
    private double[] _cov;
    private List<double[]> _rx;
    private double _chi;
    private int _ndf;
    private int _nit;
    
    public TrackFit( double[] par, double[] cov, List<double[]> rx, double chi, int ndf, int nit)
    {
        _par = par;
        _cov = cov;
        _rx = rx;
        _ndf = ndf;
        _nit = nit;
    }
    
    public List<double[]> impactPoints()
    {
        return _rx;
    }
    
    public double[] pars()
    {
        return _par;
    }
    
    public double[] cov()
    {
        return _cov;
    }
}
