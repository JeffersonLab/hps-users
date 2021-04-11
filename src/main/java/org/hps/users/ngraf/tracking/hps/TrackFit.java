package org.hps.users.ngraf.tracking.hps;

import java.util.Arrays;
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
 * @author Norman Graf
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
        _chi = chi;
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
    
    public double chisq()
    {
        return _chi;
    }
    public int ndf()
    {
        return _ndf;
    }
    public int niterations()
    {
        return _nit;
    }
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer("TrackFit:");
        sb.append("chisq "+_chi+ " ndf "+_ndf+" nit "+_nit);
        sb.append(" pars: "+Arrays.toString(_par));
        return sb.toString();
    }
}
