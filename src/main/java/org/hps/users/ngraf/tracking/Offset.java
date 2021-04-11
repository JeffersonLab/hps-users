
package org.hps.users.ngraf.tracking;

import Jama.Matrix;

/**
 *
 * @author ngraf
 */
public class Offset {
// *DOC          MASK = flags to indicate offset parameters
//*DOC          DANG = offsets in orientation
//*DOC          DUVW = offsets in location (local coordinates)
//*DOC          R0C  = new offsets in global
//*DOC          ROTC = new rotation  
    
    private int _id;
    private Matrix _rot; // rotation matrices local --> global
    private double[] _r0; //detector origin in global coordinates (mm)
    private int[] _mask; // whether to try to align this parameter
                         // 0-2 offset in location
                         // 3-5 offset in orientation (tilt)
    private double[] _angles; // smeared angles (tilts)
    private double[] _offsets; // smeared offsets
public Offset(int id, Matrix r, double[] pos, double[] angles, double[] offsets, int[] m)
    {
        _id = id;
        _rot = r;
        _r0 = pos;
        _mask = m;
        _angles = angles;
        _offsets = offsets;
    }
    
    public Matrix rot()
    {
        return _rot;
    }
    
    public double[] r0()
    {
        return _r0;
    }
    
    public int[] mask()
    {
        return _mask;
    }    
    
    public double[] angles()
    {
        return _angles;
    }
    
    public double[] offsets()
    {
        return _offsets;
    }
}
