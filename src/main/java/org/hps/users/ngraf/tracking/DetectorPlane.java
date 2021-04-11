package org.hps.users.ngraf.tracking;

import Jama.Matrix;

/**
 *
 * @author ngraf
 */
public class DetectorPlane {
    private int _id;
    private Matrix _rot; // rotation matrices local --> global
    private double[] _r0; //detector origin in global coordinates (mm)
    private double[] _sigs; //detector resolutions in v and w (mm)
    private Offset _offset; // pssible offsets in position and rotation
    
    
    
    public DetectorPlane(int id, Matrix m, double[] pos, double[] s)
    {
        this(id,m,pos,s,null);
    }
    public DetectorPlane(int id, Matrix m, double[] pos, double[] s, Offset off)
    {
        _id = id;
        _rot = m;
        _r0 = pos;
        _sigs = s;
        _offset = off;
    }
    
    public DetectorPlane(DetectorPlane p, Offset off)
    {
        _id = p.id();
        _rot = p.rot();
        _r0 = p.r0();
        _sigs = p.sigs();
        _offset = off;
    }
    
    public Matrix rot()
    {
        return _rot;
    }
    public int id()
    {
        return _id;
    }
    public double[] r0()
    {
        return _r0;
    }
    
    public double[] sigs()
    {
        return _sigs;
    }
    
    public Offset offset()
    {
        return _offset;
    }
}
