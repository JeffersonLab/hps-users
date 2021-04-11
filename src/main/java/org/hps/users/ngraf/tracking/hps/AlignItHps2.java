package org.hps.users.ngraf.tracking.hps;

import Jama.Matrix;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import static org.hps.users.ngraf.tracking.hps.FitTracks.readEvents;
import org.lcsim.math.chisq.ChisqProb;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman Graf
 */
public class AlignItHps2 {

    public static void main(String[] args) throws Exception {
        boolean alignit = true;
        String hitFileName = "misalignedHits.txt";
        //TODO read this from args.

        AIDA aida = AIDA.defaultInstance();

        // generate the ideal detector
        List<DetectorPlane> planes = FitTracks.GENER_DET();
        for (DetectorPlane plane : planes) {
            System.out.println(plane);
        }
        // fit the hits with this ideal, unaligned detector
        fitHits(hitFileName, planes, aida, "unaligned");

        if (alignit) {
            // Generate offset rotation matrices and translation vectors
            GEN_OFFSETS(planes);
            System.out.println("SIMULATED MIS-ALIGNMENTS ARE:");
            int NN = planes.size();
            int[] nprs = new int[NN];
            for (int i = 0; i < NN; ++i) {
                DetectorPlane p = planes.get(i);
                Offset o = p.offset();
                System.out.println(" PLANE " + i + " OFFSETS: " + Arrays.toString(o.offsets()) + " TILTS: " + Arrays.toString(o.angles()));
                int[] mask = o.mask();
                int doit = 0;
                for (int k = 0; k < 6; ++k) {
                    doit = doit + mask[k];
                }
                nprs[i] = doit;
            }
            int NITER = 7;
            // a map of the Alignment objects keyed by detector plane
            Map<Integer, Alignment> alignmentMap = new HashMap<Integer, Alignment>();
            double[] PAR = new double[6]; // local offsets and tilt angles
            double[] COV = new double[21]; // covariance matrix
            double[] A0 = {0., 0., 0.}; // initial guess for (x,y,z) of track
            double[] B0 = {0., 0., 1.}; // initial guess for the track direction
            for (int ITER = 0; ITER < NITER; ++ITER) { // iterate the alignment
                //           System.out.println("Iteration " + (ITER + 1));
                // read in some events generated using the misaligned detector...
                List<String> lines = readEvents(hitFileName);
                int NTIME = lines.size();
                for (int i = 0; i < NTIME; ++i) { // loop over events
                    double[] parin = new double[4];  // generated track parameters
                    List<Hit> hits = FitTracks.nextEvent(lines.get(i), parin);
                    // reconstruct using the ideal detector... 
//                System.out.println("TRACK FIT: "+(i+1)+" "+(ITER+1)+" A0 "+Arrays.toString(A0)+" B0 "+Arrays.toString(B0));
                    TrackFit fit = FitTracks.STR_LINFIT(planes, hits, A0, B0);
//                System.out.println("TRACK FIT: "+(i+1)+" "+(ITER+1)+" "+Arrays.toString(fit.pars()));
                    List<double[]> rx = fit.impactPoints();
                    //find alignment parameters
                    double[] B = new double[3];  // track direction at impact point.
                    B[0] = fit.pars()[2];  // dx/dz
                    B[1] = fit.pars()[3];  // dy/dz
                    B[2] = 1.;             // z
                    for (int j = 0; j < NN; ++j) { //loop over all the detector planes
                        int NPR = nprs[j];
                        if (NPR > 0) {  // found a detector plane which we want to align
                            if (debug()) {
                                System.out.println("J " + (j + 1) + " NPR " + NPR);
                            }
                            DetectorPlane p = planes.get(j);
                            Offset o = p.offset();
                            int[] mask = o.mask();
                            Alignment a = null;
                            if (alignmentMap.containsKey(j)) {
                                a = alignmentMap.get(j);
                            } else {
                                //best guess for plane rotation and offset is ideal position
                                //this will be updated with each iteration
                                a = new Alignment(j, mask, p.rotArray(), p.r0());
                                alignmentMap.put(j, a);
                            }
                            // need to calculate track impact point RR at this plane
                            Hit hit = hits.get(j);
                            double[] QM = hit.uvm();
                            double[] W = hit.wt();
                            double[] RX = rx.get(j);
                            //System.out.println((i+1)+" "+(j+1)+" "+Arrays.toString(RX));
                            a.accumulate(RX, B, QM, W);
                            double[] ROT = new double[9];
                            double[] R0 = new double[3];
                            if (i == NTIME - 1) {
                                a.solve(PAR, COV, ROT, R0);
                                // this updates ROT and R0, which needs to be fed to the accumulator.
//                            System.out.println("Iteration " + (ITER+1));
//                            System.out.println("Update Plane "+j);
//                            System.out.println("R0 " + Arrays.toString(R0));
//                            System.out.println("ROT " + Arrays.toString(ROT));
                                p.setUpdatedPosition(R0);
                                p.setUpdatedRotation(ROT);
                            }
                        }
                    } //loop over planes
                } // loop over events
            }// loop over iterations
        }
        // fit with aligned planes
        fitHits(hitFileName, planes, aida, "aligned");
        aida.saveAs("AlignItHps2.aida");
    }

    public static void fitHits(String filename, List<DetectorPlane> planes, AIDA aida, String dir) {
        // Testing the track fit code   
        double[] A0 = {0., 0., 0.}; //-2337.1810}; // initial guess for (x,y,z) of track
        double[] B0 = {0., 0., 1.}; // initial guess for the track direction
        List<String> lines = FitTracks.readEvents(filename);
        int NTIME = lines.size();
        System.out.println("found " + NTIME + " events in file");
        for (int i = 0; i < NTIME; ++i) { // loop over events
            double[] parin = new double[4];  // generated track parameters
            List<Hit> hits = FitTracks.nextEvent(lines.get(i), parin);
            // reconstruct using the ideal detector... 
            TrackFit fit = FitTracks.STR_LINFIT(planes, hits, A0, B0);
            double[] pars = fit.pars();
            double[] cov = fit.cov();
            if (debug()) {
                System.out.println("fit: " + fit);
                System.out.println("pars " + Arrays.toString(pars));
                System.out.println("parin " + Arrays.toString(parin));
                System.out.println("fit cov " + Arrays.toString(cov));
            }
//            System.out.println(Arrays.toString(cov));
            aida.cloud1D(dir + "/fit chisq per ndf").fill(fit.chisq() / fit.ndf());
            aida.cloud1D(dir + "/fit ndf ").fill(fit.ndf());
            double chisqProb = ChisqProb.gammp(fit.ndf(), fit.chisq());
            aida.cloud1D(dir + "/fit chisq prob ").fill(chisqProb);
            aida.cloud1D(dir + "/fit niterations").fill(fit.niterations());
            aida.cloud1D(dir + "/x meas-pred").fill(pars[0] - parin[0]);
            aida.cloud1D(dir + "/y meas-pred").fill(pars[1] - parin[1]);
            aida.cloud1D(dir + "/dxdz meas-pred").fill(pars[2] - parin[2]);
            aida.cloud1D(dir + "/dydz meas-pred").fill(pars[3] - parin[3]);
            aida.cloud1D(dir + "/x meas-pred pull").fill((pars[0] - parin[0]) / sqrt(cov[0]));
            aida.cloud1D(dir + "/y meas-pred pull").fill((pars[1] - parin[1]) / sqrt(cov[2]));
            aida.cloud1D(dir + "/dxdz meas-pred pull").fill((pars[2] - parin[2]) / sqrt(cov[5]));
            aida.cloud1D(dir + "/dydz meas-pred pull").fill((pars[3] - parin[3]) / sqrt(cov[9]));

        }

    }

    static void GEN_OFFSETS(List<DetectorPlane> planes) {

        // idea is to loop over the existing list of planes and add offsets to
        // each. that keeps the ideal and misaligned detectors together
        //
//*DOC  Input:
//*DOC          NN   = number of detectors
//*DOC          R0   = un-offset locations of detector origins
//*DOC          ROT  = un-tilted rot mat
//*DOC  Output:
//*DOC          MASK = flags to indicate offset parameters
//*DOC          DANG = offsets in orientation
//*DOC          DUVW = offsets in location (local coordinates)
//*DOC          R0C  = new offsets in global
//*DOC          ROTC = new rotation      
        // a bit hinky here since we're hard-coding some of this.
        // can generalize this later, though.
        //offsets in orientation (tilts)
        int NN = planes.size();
        List<double[]> dq = new ArrayList<double[]>();
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{-0.1, -0.05, 0.0});
        dq.add(new double[]{0., -0., -0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        dq.add(new double[]{0., 0., 0.});
        //offsets in location      
        List<double[]> da = new ArrayList<double[]>();
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{-0.001, 0.003, 0.002});
        da.add(new double[]{0.0, 0.0, 0.0});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        da.add(new double[]{0., 0., 0.});
        // loop over the existing planes...
        for (int i = 0; i < NN; ++i) {
            if (debug()) {
                System.out.println("Generating offsets for detector " + (i + 1));
            }
            DetectorPlane p = planes.get(i);
            int[] MASK = new int[6]; // keep track of which parameters have been offset
            // translation offset (shifts)
            double[] q = dq.get(i);
            // rotation offsets (tilts)
            double[] a = da.get(i);
            double[][] RW = new double[3][9];
            for (int j = 0; j < 3; ++j) {
                if (q[j] != 0.) {
                    MASK[j] = 1;  // displaced
                }
                if (a[j] != 0.) {
                    MASK[j + 3] = 1; //rotated
                }
                double angl = a[j];//Math.toRadians(a[j]);
                GEN_ROTMAT(angl, j, RW[j]);
            }
            // combined rotation
            double[] DROT = new double[9];
            PROD_ROT(RW[2], RW[1], RW[0], DROT); //delta rotation
            double[] ROT = p.rotArray();
            double[] ROTC = new double[9];
            MXMPY(ROT, DROT, ROTC, 3, 3, 3); // corr rot
            //translations
            double[] DR0 = new double[3];

            VMATR(q, ROTC, DR0, 3, 3);
            double[] R0 = p.r0();
            double[] R0C = VADD(R0, DR0, 3);

            //create an Offset
            Offset o = new Offset(i, ROTC, R0C, a, q, MASK);
            // add this offset to the existing detector plane...
            p.addOffset(o);
        }
    }

    static double[][] UCOPY(double[][] source) {
        double[][] destination = new double[source.length][];
        for (int i = 0; i < source.length; ++i) {
            // allocating space for each row of destination array
            destination[i] = new double[source[i].length];
            System.arraycopy(source[i], 0, destination[i], 0, destination[i].length);
        }
        // displaying destination array
        if (debug()) {
            System.out.println(Arrays.deepToString(destination));
        }
        return destination;
    }

    static double[] TRSA(double[] S, double[] A, int N) {
//CNG CALL TRSA (S,A,C,M,N) SA -> C
//CNG A and C are M X N rectangular matrices
//CNG S is an M X M symmetrix matrix
        Matrix s = new Matrix(N, N);
        s.set(0, 0, S[0]);
        s.set(1, 0, S[1]);
        s.set(1, 1, S[2]);
        s.set(2, 0, S[3]);
        s.set(2, 1, S[4]);
        s.set(2, 2, S[5]);
        s.set(3, 0, S[6]);
        s.set(3, 1, S[7]);
        s.set(3, 2, S[8]);
        s.set(3, 3, S[9]);

        s.set(0, 1, s.get(1, 0));
        s.set(0, 2, s.get(2, 0));
        s.set(0, 3, s.get(3, 0));
        s.set(1, 2, s.get(2, 1));
        s.set(1, 3, s.get(3, 1));
        s.set(2, 3, s.get(3, 2));

        if (debug()) {
            s.print(6, 4);
        }
        Matrix a = new Matrix(4, 1);
        a.set(0, 0, A[0]);
        a.set(1, 0, A[1]);
        a.set(2, 0, A[2]);
        a.set(3, 0, A[3]);
        if (debug()) {
            a.print(6, 4);
        }

        Matrix d = s.times(a);
        if (debug()) {
            d.print(6, 4);
        }

        double[] r = new double[4];
        r[0] = d.get(0, 0);
        r[1] = d.get(1, 0);
        r[2] = d.get(2, 0);
        r[3] = d.get(3, 0);

        return r;
    }

    static double[] TRSINV(double[] A, int M) {
        // invert symmetric matrix given by lower diagonal element array...
        Matrix a = new Matrix(M, M);
        a.set(0, 0, A[0]);
        a.set(1, 0, A[1]);
        a.set(1, 1, A[2]);
        a.set(2, 0, A[3]);
        a.set(2, 1, A[4]);
        a.set(2, 2, A[5]);
        a.set(3, 0, A[6]);
        a.set(3, 1, A[7]);
        a.set(3, 2, A[8]);
        a.set(3, 3, A[9]);

        a.set(0, 1, a.get(1, 0));
        a.set(0, 2, a.get(2, 0));
        a.set(0, 3, a.get(3, 0));
        a.set(1, 2, a.get(2, 1));
        a.set(1, 3, a.get(3, 1));
        a.set(2, 3, a.get(3, 2));

        if (debug()) {
            a.print(6, 4);
        }

        Matrix b = a.inverse();
        if (debug()) {
            b.print(6, 4);
        }

        double[] r = new double[10];
        r[0] = b.get(0, 0);
        r[1] = b.get(1, 0);
        r[2] = b.get(1, 1);
        r[3] = b.get(2, 0);
        r[4] = b.get(2, 1);
        r[5] = b.get(2, 2);
        r[6] = b.get(3, 0);
        r[7] = b.get(3, 1);
        r[8] = b.get(3, 2);
        r[9] = b.get(3, 3);
        return r;

    }

    static double[] MXMPY(double[][] A, double[] B, int I, int J, int K) {
        //MXMPY(A,B,C,NI,NJ,NK) (Aij)(Bjk) -> (Cik)
        // called with I=4, J=2, K=1
        double[] C = new double[4];
        C[0] = A[0][0] * B[0] + A[1][0] * B[1];
        C[1] = A[2][0] * B[0] + A[3][0] * B[1];
        C[2] = A[0][1] * B[0] + A[1][1] * B[1];
        C[3] = A[2][1] * B[0] + A[3][1] * B[1];
        return C;
    }

    static double[][] TRATS(double[][] B, double[] R, int N, int M) {
        // s is symmetric matrix packed as lower diagonal
        // R = B'S
        // N=2, M=4 
        Matrix b = new Matrix(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[i][j]);
            }
        }
        Matrix s = new Matrix(N, N); // symmetric
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) {
                s.set(i, j, R[i * N + j]);
                s.set(j, i, R[i * N + j]);
            }
        }
        if (debug()) {
            s.print(6, 4);
        }
        if (debug()) {
            b.print(6, 4);
        }
        Matrix r = b.transpose().times(s);
        if (debug()) {
            r.print(6, 4);
        }
        return r.getArray();
    }

    static double[] TRATSA(double[][] B, double[] S, int N, int M) {
        // s is symmetric matrix packed as lower(?) diagonal
        // R = B'SB
        // N=2, M=4
        Matrix b = new Matrix(N, M);
        Matrix s = new Matrix(N, N); // symmetric
        s.set(0, 0, S[0]);
        if (debug()) {
            s.print(6, 4);
        }
        if (debug()) {
            b.print(6, 4);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[i][j]);
            }
        }
        if (debug()) {
            s.print(6, 4);
        }
        if (debug()) {
            b.print(6, 4);
        }
        Matrix r = b.transpose().times(s.times(b));
        if (debug()) {
            r.print(6, 4);
        }
        //OK, now get the lower diagonal
        //matrix is now MxM
        double[] R = new double[10];
        int n = 0;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j <= i; ++j) {
                R[n++] = r.get(i, j);
            }
        }
        return R;
    }

    static double TRASAT(double[] a, double[] s, int N) {
        // s is a symmetric (lower diagonal?) matrix
        // a is a vector
        double r = 0.;
        if (N == 1) {
            r = a[0] * s[0] * a[0];
        } else if (N == 2) {
            r = a[0] * (a[0] * s[0] + a[1] * s[1]) + a[1] * (a[0] * s[1] + a[1] * s[2]);
        }
        return r;
    }

    static double[] asArray(Matrix a) {
        double[] r = new double[a.getRowDimension() * a.getColumnDimension()];
        int c = 0;
        for (int j = 0; j < a.getColumnDimension(); ++j) {
            for (int i = 0; i < a.getRowDimension(); ++i) {
                r[c++] = a.get(i, j);
            }
        }
        return r;
    }

    static void myMXMPY(double[][] A, double[] B, double[] C, int I, int J, int K) {
        int rdim = A.length;
        int cdim = A[0].length;
        double[] a = new double[rdim * cdim];
        int c = 0;
        for (int j = 0; j < cdim; ++j) {
            for (int i = 0; i < rdim; ++i) {
                a[c++] = A[i][j];
            }
        }
        myMXMPY(a, B, C, I, J, K);
    }

    static void myMXMPY(double[] A, double[] B, double[] C, int I, int J, int K) {
        int IIA = 1;
        int IOA = J;
        int IIB = K;
        int IOB = 1;
        int IA = 0;
        int IC = 0;
        for (int L = 0; L < I; ++L) {
            int IB = 0;
            for (int M = 0; M < K; ++M) {
                C[IC] = 0.;
                if (J > 0) {
                    int JA = IA;
                    int JB = IB;
                    for (int N = 0; N < J; ++N) {
                        C[IC] = C[IC] + A[JA] * B[JB];
                        JA = JA + IIA;
                        JB = JB + IIB;
                    } //ENDDO
                    IB = IB + IOB;
                    IC = IC + 1;
                }
            }
            IA = IA + IOA;
        }
    }

    static void GEN_ROTMAT(double ANG, int IAX, double[] R) {
        double[][] r = new double[3][3];
        if (debug()) {
            System.out.println("ANG " + ANG + " IAX " + (IAX + 1));
        }
        r[IAX][IAX] = 1;
        double C = cos(ANG);
        double S = sin(ANG);
        switch (IAX) {
            case 0:
                r[1][1] = C;
                r[1][2] = S;
                r[2][1] = -S;
                r[2][2] = C;
                break;
            case 1:
                r[0][0] = C;
                r[0][2] = -S;
                r[2][0] = S;
                r[2][2] = C;
                break;
            case 2:
                r[0][0] = C;
                r[0][1] = S;
                r[1][0] = -S;
                r[1][1] = C;
                break;
            default:
                break;
        }
        int n = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R[n++] = r[i][j];
            }
        }
    }

    static void PROD_ROT(double[] RA, double[] RB, double[] RC, double[] RTOT) {
        double[] RAPU = new double[9];
        MXMPY(RA, RB, RAPU, 3, 3, 3);
        MXMPY(RAPU, RC, RTOT, 3, 3, 3);
    }

    static void MXMPY(double[] A, double[] B, double[] C, int I, int J, int K) {
        int IIA = 1;
        int IOA = J;
        int IIB = K;
        int IOB = 1;
        int IA = 0;
        int IC = 0;
        for (int L = 0; L < I; ++L) {
            int IB = 0;
            for (int M = 0; M < K; ++M) {
                C[IC] = 0.;
                if (J > 0) {
                    int JA = IA;
                    int JB = IB;
                    for (int N = 0; N < J; ++N) {
                        C[IC] = C[IC] + A[JA] * B[JB];
                        JA = JA + IIA;
                        JB = JB + IIB;
                    }
                    IB = IB + IOB;
                    IC = IC + 1;
                }
            }
            IA = IA + IOA;
        }
    }

    static void VMATR(double[] A, double[] G, double[] X, int M, int N) {
        // A,X One-dimensional arrays of length X
        // G Two-dimensional array of dimension (M,N), stored row-wise
        // GA = X
        for (int i = 0; i < M; ++i) {
            X[i] = 0;
            for (int j = 0; j < N; ++j) {
                X[i] = X[i] + G[i + j * M] * A[j];
            }
        }
    }

    static double[] VADD(double[] a, double[] b, int NELEMF) {
        double[] c = new double[NELEMF];
        for (int i = 0; i < NELEMF; ++i) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    //TODO understand this...
    // WTF!? this is wrong, but seems to work for the alignment!!!
    static Matrix GEN_ROTMAT(double ANG, int IAX) {
        Matrix r = Matrix.identity(3, 3);
        if (debug()) {
            System.out.println("ANG " + ANG + " IAX " + (IAX + 1));
        }
        double C = cos(ANG);
        double S = sin(ANG);
        switch (IAX) {
            case 0:
                r.set(1, 1, C);
                r.set(1, 2, S);
                r.set(2, 1, -S);
                r.set(2, 2, C);
                break;
            case 1:
                r.set(0, 0, C);
                r.set(0, 2, -S);
                r.set(2, 0, S);
                r.set(2, 2, C);
                break;
            case 2:
                r.set(0, 0, C);
                r.set(0, 1, S);
                r.set(1, 0, -S);
                r.set(1, 1, C);
                break;
            default:
                break;
        }
        return r;
    }

    static void GEN_ROTMAT(double ANG, int IAX, double[][] ROTM) {
// *DOC  SUBROUTINE GEN_ROTMAT(ANG,IAX,ROTM)
//*DOC  Action: Generate rotation matrix around axis IAX
//*DOC
//*DOC  Input: ANG = rotation angle (rad) around axis IAX
//*DOC         IAX = axis number (1, 2 or 3) of rotation
//*DOC
//*DOC  Output: ROTM = rotation matrix  
        ROTM[IAX][IAX] = 1.;
        double C = cos(ANG);
        double S = sin(ANG);
        if (IAX == 0) {
            ROTM[1][1] = C;
            ROTM[2][1] = S;
            ROTM[1][2] = -S;
            ROTM[2][2] = C;
        } else if (IAX == 1) {
            ROTM[0][0] = C;
            ROTM[2][0] = -S;
            ROTM[0][2] = S;
            ROTM[2][2] = C;
        } else if (IAX == 2) {
            ROTM[0][0] = C;
            ROTM[1][0] = S;
            ROTM[0][1] = -S;
            ROTM[1][1] = C;
        }
    }

    static Matrix PROD_ROT(Matrix rx, Matrix ry, Matrix rz) {
        return rx.times(ry.times(rz));
    }

    static List<Hit> GENER_EVT(List<DetectorPlane> dets, double[] par, double z) {
        // PAR = parameters a1,a2,b1,b2 (generated)
        // x, y, dx/dz, dy/dz
//CNG   FLAT BETWEEN -25 AND 25
//      PAR(1) = 25.*(2*RNDM(I+1)-1.)
//CNG   FLAT BETWEEN -60 AND 60
//      PAR(2) = 60.*(2*RNDM(I+3)-1.)
//      PAR(3) = (PAR(1)-25.*(2*RNDM(I+2)-1.))/500.
//      PAR(4) = (PAR(2)-60.*(2*RNDM(I+4)-1.))/500.  
        Random ran = new Random();
        par[0] = 25. * (2. * ran.nextDouble() - 1.);
        par[1] = 60. * (2. * ran.nextDouble() - 1.);
        par[2] = (par[0] - 25. * (2. * ran.nextDouble() - 1.)) / 500.;
        par[3] = (par[1] - 60. * (2. * ran.nextDouble() - 1.)) / 500.;

        // should be randomly generating track parameters, but for now lets used a fixed set
        List<Hit> hitlist = new ArrayList<Hit>();
        int N = dets.size();
        double[] W = {0., 0., 1.};
        //initial guess is along zhat
        double[] A = {par[0], par[1], z};
        double[] B = {par[2], par[3], 1.0000000};
        Matrix w = new Matrix(W, 3);
        Matrix a = new Matrix(A, 3);
        Matrix b = new Matrix(B, 3);
        for (int i = 0; i < N; ++i) {
            Hit hit;
            DetectorPlane dp = dets.get(i);
            double[] sigs = dp.sigs();
            Matrix ROT = dp.rot();
            double[] R0 = dp.r0();
            if (debug()) {
                System.out.println("CALLING VMATR(W, ROT(1," + (i + 1) + "), WG, 3, 3) ");
            }
            Matrix wg = ROT.times(w);
            if (debug()) {
                System.out.println(" WG " + wg.get(0, 0) + " " + wg.get(1, 0) + " " + wg.get(2, 0));
            }
            Matrix bwg = b.transpose().times(wg);
            if (debug()) {
                System.out.println(" BWG " + bwg.get(0, 0));
            }
            ImpactPoint ip = GET_IMPACT(A, B, ROT, R0, wg, bwg.get(0, 0));
            //TODO will have to modify when I restrict code to 1D strips
            double[] u = new double[2];  // 2dim for u and v measurement 
            double[] wt = new double[3]; // lower diag cov matrix
            for (int j = 0; j < 2; ++j) // do fluctuations
            {
                double coor = ip.q()[j]; // local u coordinate of impact point
                double smear = sigs[j] * ran.nextGaussian(); // fluctuation
                coor += smear;
                u[j] = coor;
                int k = j * (j + 1); // should be 0 and 2
                wt[1] = 0.; // explicitly list here in case we ever have cov between u and v
                if (sigs[j] > 0.) {
                    wt[k] = 1 / (sigs[j] * sigs[j]);
                }
                if (debug()) {
                    System.out.println("MEASUREMENT Q(" + (j + 1) + ") " + ip.q()[j]);
                }
                if (debug()) {
                    System.out.println("SMEARED MEASUREMENT " + (j + 1) + ") " + coor);
                }
            }
            hit = new Hit(u, wt);
            if (debug()) {
                System.out.println("SMEARED UVM(JJ," + (i + 1) + ") " + hit.uvm()[0] + " " + hit.uvm()[1]);
            }
            if (debug()) {
                System.out.println("WT(JJ, " + (i + 1) + ") " + hit.wt()[0] + " " + hit.wt()[1] + " " + hit.wt()[2]);
            }
            hitlist.add(hit);
        }
        return hitlist;
    }

    static ImpactPoint GET_IMPACT(double[] A, double[] B, Matrix ROT, double[] R0, Matrix wg, double BWG) {
        if (debug()) {
            System.out.println("GET_IMPACT A " + Arrays.toString(A));
            System.out.println("GET_IMPACT B " + Arrays.toString(B));
            System.out.println("GET_IMPACT R0 " + Arrays.toString(R0));
        }
        double[] AMR = VSUB(A, R0);
        Matrix amr = new Matrix(AMR, 3);
        if (debug()) {
            System.out.println("GET_IMPACT AMR " + AMR[0] + " " + AMR[1] + " " + AMR[2]);
        }
        double ti = -amr.transpose().times(wg).get(0, 0) / BWG;
        if (debug()) {
            System.out.println("GET_IMPACT PARAMETER AT IMPACT TI " + ti);
        }
        double[] APU = VLINE(AMR, 1., B, ti);
        if (debug()) {
            System.out.println("GET_IMPACT APU " + APU[0] + " " + APU[1] + " " + APU[2]);
        }
        double[] q = VMATL(ROT, APU);
        if (debug()) {
            System.out.println("GET_IMPACT LOCAL U,V,W Q " + q[0] + " " + q[1] + " " + q[2]);
        }
        double[] r = VADD(APU, R0);
        if (debug()) {
            System.out.println("GET_IMPACT GLOBAL X,Y,Z R " + r[0] + " " + r[1] + " " + r[2]);
        }
        return new ImpactPoint(ti, q, r);
    }

    static double[] VMATL(Matrix ROT, double[] APU) {
        double[] c = new double[3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                c[i] += ROT.get(j, i) * APU[j];
            }
        }
        return c;
    }

    static double[] VLINE(double[] A, double F1, double[] B, double F2) {
        double[] c = new double[A.length];
        for (int i = 0; i < A.length; ++i) {
            c[i] = A[i] * F1 + B[i] * F2;
        }
        return c;
    }

    static double[] VSUB(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = a[i] - b[i];
        }
        return c;
    }

    static double[] VADD(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    static boolean debug() {
        return false;
    }

    static boolean debug2() {
        return false;
    }
}
