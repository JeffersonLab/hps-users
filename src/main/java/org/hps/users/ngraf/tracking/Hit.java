/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking;

/**
 *
 * @author ngraf
 */
public class Hit {
    
    private double[] _uvm; //measurements (local coordinates)
    private double[] _wt; //weights
    
    public Hit(double[] u, double[] w)
    {
        _uvm = u;
        _wt = w;
    }
    
    public double[] uvm()
    {
        return _uvm;
    }
    
    public double[] wt()
    {
        return _wt;
    }
}
