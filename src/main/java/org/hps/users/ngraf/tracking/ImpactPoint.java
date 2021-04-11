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
public class ImpactPoint {

    private double _ti; //parameter at impact point
    private double[] _q;//local  u,v,w coordinates
    private double[] _r;//global x,y,z coordinates

    public ImpactPoint(double ti, double[] q, double[] r) {
        _ti = ti;
        _q = q;
        _r = r;
    }

    public double ti() {
        return _ti;
    }

    public double[] q() {
        return _q;
    }

    public double[] r() {
        return _r;
    }
}
