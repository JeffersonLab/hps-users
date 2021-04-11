package org.hps.users.ngraf.tracking;

import Jama.Matrix;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author ngraf
 */
public class LineFit {

    public static void main(String[] args) throws Exception {

        AIDA aida = AIDA.defaultInstance();
        int NN = 8;
        System.out.println("Generating detector");
        List<DetectorPlane> planes = GENER_DET(NN);
        // Testing the track fit code   
        int NTIMES = 10000;
        if (debug()) {
            System.out.println("Generating " + NTIMES + " events");
        }
        for (int i = 0; i < NTIMES; ++i) {
            //Generate tracks and hits in ideal detectors
            // contains the generated track parameters
            double[] parin = new double[4];

            List<Hit> hits = GENER_EVT(planes, parin);
            if (debug()) {
                System.out.println("*** END OF EVENT " + (i + 1) + " ***");
            }
            if (debug()) {
                System.out.println("FITTING TRACK");
            }
            TrackFit tf = STR_LINFIT(planes, hits);
            double[] par = tf.pars();
            double[] cov = tf.cov();
            for (int jj = 0; jj < 4; ++jj) {
                double error = sqrt(cov[jj + (jj * (jj + 1) / 2)]);//sqrt(cov[jj * (jj + 1) / 2]); // check this
                double pull = (par[jj] - parin[jj]) / error;
                aida.cloud1D(" parin " + jj + " ").fill(parin[jj]);
                aida.cloud1D(" par " + jj + " - parin " + jj).fill(par[jj] - parin[jj]);
                aida.cloud1D(" par " + jj + " pull").fill(pull);
                if (debug()) {
                    System.out.println(" PAR(" + (jj + 1) + ") " + par[jj] + " PARIN(" + (jj + 1) + ") " + parin[jj]);
                }
            }
        }

        aida.saveAs("LineFit.aida");
    }

    static List<DetectorPlane> GENER_DET(int NN) {
        double[] ANGLS = {0., 0., 10., 10., 0., 95., 10., 20., 30., 0., 20., 90., 7., 0., 0., 0., 7., 80., 3., 0., 0., 0., 0., 93.};
        double[] XCOOR = {0., 0., 31., 0., 0., 64., 0., 0., 210., 0., 0., 422., 0., 0., 372., 0., 0., 406., 0., 0., 539., 0., 0., 573.};
        double[] SIGS = {0.010, 0.200};     //  detector resolutions in mm

        List<double[]> angles = new ArrayList<double[]>();
        angles.add(new double[]{0., 0., 10.});
        angles.add(new double[]{10., 0., 95.});
        angles.add(new double[]{10., 20., 30.});
        angles.add(new double[]{0., 20., 90.});
        angles.add(new double[]{7., 0., 0.});
        angles.add(new double[]{0., 7., 80.});
        angles.add(new double[]{3., 0., 0.});
        angles.add(new double[]{0., 0., 93.});

        List<double[]> r0 = new ArrayList<double[]>();
        r0.add(new double[]{0., 0., 31.});
        r0.add(new double[]{0., 0., 64.});
        r0.add(new double[]{0., 0., 210.});
        r0.add(new double[]{0., 0., 422.});
        r0.add(new double[]{0., 0., 372.});
        r0.add(new double[]{0., 0., 406.});
        r0.add(new double[]{0., 0., 539.});
        r0.add(new double[]{0., 0., 573.});

        List<Matrix> rot = new ArrayList<Matrix>();

        List<DetectorPlane> planes = new ArrayList<DetectorPlane>();

        if (debug()) {
            System.out.println("GENERATING DETECTOR WITH " + NN + " DETECTOR PLANES");
        }
        for (int i = 0; i < NN; ++i) {
            if (debug()) {
                System.out.println("Generating rotation matrices for detector " + (i + 1));
            }
            double[] a = angles.get(i);
            Matrix[] mats = new Matrix[3];
            for (int j = 0; j < 3; ++j) {
                double angl = Math.toRadians(a[j]);
                mats[j] = GEN_ROTMAT(angl, j);
                if (debug()) {
                    System.out.println("Rotation matrix for axis " + (j + 1) + " angle " + angl);
                }
                if (debug()) {
                    mats[j].print(6, 4);
                }
            }
            Matrix prodrot = PROD_ROT(mats[0], mats[1], mats[2]);
            if (debug()) {
                System.out.println("Combined rotation matrix is");
            }
            if (debug()) {
                prodrot.print(6, 4);
            }
            rot.add(prodrot);
            if (debug()) {
                System.out.println("Translation vector is");
            }
            double[] off = r0.get(i);
            if (debug()) {
                System.out.println(off[0] + " " + off[1] + " " + off[2]);
            }
            if (debug()) {
                System.out.println("Sigs: " + SIGS[0] + " " + SIGS[1]);
            }

            planes.add(new DetectorPlane(i, prodrot, off, SIGS));
        }
        return planes;
    }

    static TrackFit STR_LINFIT(List<DetectorPlane> planes, List<Hit> hits) {

        double[][] UVW = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        if (debug()) {
            System.out.println("\n\n\nIN STR_LINFIT\n\n\n");
        }
        double CHI = 1.E10;
        double CHI0 = 0.;
        int NDF = -4;
        int NIT = 0;
        double[] PAR = new double[4];
        double[] COV = new double[10];
        List<ImpactPoint> rx = new ArrayList<ImpactPoint>();
        // LOOP comes here...
        do {
            rx.clear();
            NIT++;
            // initial guess is track along the z axis
            double[] A = {0., 0., 0.};
            double[] B = {0., 0., 1.};

            A[0] = PAR[0];
            A[1] = PAR[1];
            B[0] = PAR[2];
            B[1] = PAR[3];

            Matrix b = new Matrix(B, 1);

            double[] BUVW = new double[3];

            CHI0 = CHI;
            CHI = 0.;
            double DCHI = 0;
            NDF = -4;

            double[] WMT = new double[10];
            double[] WVT = new double[4];
//
// See: http://cmsdoc.cern.ch/documents/99/note99_041.ps.Z
//
            int N = planes.size();
            for (int i = 0; i < N; ++i) {
                if (debug()) {
                    System.out.println("DETECTOR " + (i + 1));
                }
                DetectorPlane dp = planes.get(i);
                Matrix rot = dp.rot();
                Hit h = hits.get(i);

                Matrix[] uvwg = new Matrix[3];
                for (int j = 0; j < 3; ++j) {
                    if (debug()) {
                        System.out.println("  CALLING VMATR");
                    }
                    Matrix uvw = new Matrix(UVW[j], 3);
                    if (debug()) {
                        System.out.println("  UVW(" + (j + 1) + ") " + uvw.get(0, 0) + " " + uvw.get(1, 0) + " " + uvw.get(2, 0));
                    }
                    uvwg[j] = rot.times(uvw);
                    if (debug()) {
                        System.out.println("UVWG(" + (j + 1) + ") " + uvwg[j].get(0, 0) + " " + uvwg[j].get(1, 0) + " " + uvwg[j].get(2, 0) + " ");
                    }
                    BUVW[j] = b.times(uvwg[j]).get(0, 0);
                }
                if (debug()) {
                    System.out.println("   BUVW " + BUVW[0] + " " + BUVW[1] + " " + BUVW[2]);
                }
                ImpactPoint ip = GET_IMPACT(A, B, rot, dp.r0(), uvwg[2], BUVW[2]);
                rx.add(ip);
                int NM = 0;
                double[] wt = h.wt();
                double[] eps = new double[2];
                if (wt[0] > 0.) {
                    NM = NM + 1; // precise measurement
                }
                if (wt[2] > 0.) {
                    NM = NM + 1; // also coarse
                }
                double[] q = ip.q();    //impact point in local coordinates
                double[] uvm = h.uvm();  // hit in local coordinates
                double[][] DER = new double[4][2];
                double ti = ip.ti();
                if (debug()) {
                    System.out.println("uvwg[0]");
                }
                if (debug()) {
                    uvwg[0].print(6, 4);
                }
                if (debug()) {
                    System.out.println("uvwg[1]");
                }
                if (debug()) {
                    uvwg[1].print(6, 4);
                }
                if (debug()) {
                    System.out.println("uvwg[2]");
                }
                if (debug()) {
                    uvwg[2].print(6, 4);
                }

                for (int j = 0; j < NM; ++j) {
                    eps[j] = q[j] - uvm[j]; //predicted minus measured
                    DER[0][j] = uvwg[j].get(0, 0) - uvwg[2].get(0, 0) * BUVW[j] / BUVW[2];
                    DER[1][j] = uvwg[j].get(1, 0) - uvwg[2].get(1, 0) * BUVW[j] / BUVW[2];
                    DER[2][j] = ti * DER[0][j];
                    DER[3][j] = ti * DER[1][j];
                    if (debug()) {
                        System.out.println("j " + (j + 1));
                    }
                    if (debug()) {
                        System.out.println(uvwg[j].get(0, 0) + " " + uvwg[2].get(0, 0) + " " + BUVW[j] + " " + BUVW[2]);
                    }
                    if (debug()) {
                        System.out.println(uvwg[j].get(1, 0) + " " + uvwg[2].get(1, 0) + " " + BUVW[j] + " " + BUVW[2]);
                    }
                    if (debug()) {
                        System.out.println(DER[0][j]);
                    }
                    if (debug()) {
                        System.out.println(DER[1][j]);
                    }
                    if (debug()) {
                        System.out.println(DER[2][j]);
                    }
                    if (debug()) {
                        System.out.println(DER[3][j]);
                    }
//                DER[0][j] = UVWG(1, J) - UVWG(1, 3) * BUVW(J) / BUVW(3);
//                            uvwg[0].get(0,0)
//                DER[1][j] = UVWG(2, J) - UVWG(2, 3) * BUVW(J) / BUVW(3);
//                DER[2][j] = TI * DER(1, J);
//                DER[3][j] = TI * DER(2, J);
                    NDF = NDF + 1;
                }

                if (NM > 0) {
                    if (debug()) {
                        System.out.println("CALLING TRASAT");
                    }
                    DCHI = TRASAT(eps, wt, NM);
                    if (debug()) {
                        System.out.println(" EPS " + eps[0] + " " + eps[1]);
                    }
                    if (debug()) {
                        System.out.println(" WT(" + (i + 1) + ") " + wt[0] + " " + wt[1] + " " + wt[2]);
                    }
                    if (debug()) {
                        System.out.println(" DCHI " + DCHI);
                    }
                    CHI = CHI + DCHI;
                    if (debug()) {
                        System.out.println("CALLING TRATSA");
                    }
                    double[] WP = TRATSA(DER, wt, NM, 4);
                    for (int jj = 0; jj < 2; ++jj) {
                        if (debug()) {
                            System.out.println(" DER(" + (jj + 1) + ") " + DER[0][jj] + " " + DER[1][jj] + " " + DER[2][jj] + " " + DER[3][jj] + " ");
                        }
                    }
                    if (debug()) {
                        System.out.println("wt(" + wt[0] + " " + wt[1] + " " + wt[2]);
                    }
                    if (debug()) {
                        System.out.println(" WP " + Arrays.toString(WP));
                    }
                    if (debug()) {
                        System.out.println("CALLING TRATS");
                    }
                    double[][] C = TRATS(DER, wt, NM, 4);
                    for (int jj = 0; jj < NM; ++jj) {
                        if (debug()) {
                            System.out.println(" C(" + (jj + 1) + ") " + C[0][jj] + " " + C[1][jj] + " " + C[2][jj] + " " + C[3][jj]);
                        }
                    }
                    if (debug()) {
                        System.out.println(" EPS " + eps[0] + " " + eps[1]);
                    }
                    if (debug()) {
                        System.out.println("CALLING MXMPY(C, EPS, 4, 2, 1)");
                    }
                    double[] WV = new double[4]; //MXMPY(C, eps, 4, 2, 1);
                    if (debug()) {
                        System.out.println(" WV " + WV[0] + " " + WV[1] + " " + WV[2] + " " + WV[3]);
                    }
                    myMXMPY(C, eps, WV, 4, 2, 1);
                    if (debug()) {
                        System.out.println(" WV " + WV[0] + " " + WV[1] + " " + WV[2] + " " + WV[3]);
                    }
                    if (debug()) {
                        System.out.println("CALLING VADD");
                    }
                    WMT = VADD(WMT, WP);
                    if (debug()) {
                        System.out.println(" WMT " + Arrays.toString(WMT));
                    }
                    if (debug()) {
                        System.out.println("CALLING VADD");
                    }
                    WVT = VADD(WVT, WV);
                    if (debug()) {
                        System.out.println("WVT(" + WVT[0] + " " + WVT[1] + " " + WVT[2] + " " + WVT[3]);
                    }
                }
            } // end of loop over planes
            if (NDF > 0) {
                if (debug()) {
                    System.out.println("CALLING TRSINV");
                }
                COV = TRSINV(WMT, 4);
                if (debug()) {
                    System.out.println(" COV " + Arrays.toString(COV));
                }
                if (debug()) {
                    System.out.println("CALLING TRSA");
                }
                double[] DPAR = TRSA(COV, WVT, 4); // parameters delta
                if (debug()) {
                    System.out.println("DPAR " + Arrays.toString(DPAR));
                }
                if (debug()) {
                    System.out.println("CALLING VSUB");
                }
                if (debug()) {
                    System.out.println(" PAR " + Arrays.toString(PAR));
                }
                PAR = VSUB(PAR, DPAR);
                if (debug()) {
                    System.out.println("PAR-DPAR " + Arrays.toString(PAR));
                }
                if (debug()) {
                    System.out.println(" CHI0 " + CHI0 + " CHI " + CHI + " CHI0-CHI " + (CHI0 - CHI));
                }
            }
            // end of while loop
        } while (CHI0 - CHI > 1);

        if (debug()) {
            System.out.println("PAR " + Arrays.toString(PAR));
        }
        if (debug()) {
            System.out.println("COV " + Arrays.toString(COV));
        }
        int count = 0;
        List<double[]> impactPoints = new ArrayList<double[]>();
        for (ImpactPoint ip : rx) {
            if (debug()) {
                System.out.println("RX( " + (count + 1) + " )" + Arrays.toString(ip.r()));
            }
            impactPoints.add(ip.r());
            count++;
        }
        if (debug()) {
            System.out.println("CHI " + CHI);
        }
        if (debug()) {
            System.out.println("NDF " + NDF);
        }
        if (debug()) {
            System.out.println("NIT " + NIT);
        }
        return new TrackFit(PAR, COV, impactPoints, CHI, NDF, NIT);
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

    static double[][] TRATS(double[][] B, double[] S, int N, int M) {
        // s is symmetric matrix packed as lower diagonal
        // R = B'S
        // N=2, M=4 
        Matrix b = new Matrix(N, M);
        Matrix s = new Matrix(N, N); // symmetric
        s.set(0, 0, S[0]);
        s.set(0, 1, S[1]);
        s.set(1, 0, S[1]);
        s.set(1, 1, S[2]);
        if (debug()) {
            s.print(6, 4);
        }
        if (debug()) {
            b.print(6, 4);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[j][i]);
            }
        }
        // if(debug())System.out.println("s");
//        s.print(6, 4);
        //if(debug())System.out.println("b");
//        b.print(6, 4);
//TODO figure out what is going on here...
        Matrix r = s.times(b);
        // if(debug())System.out.println("r");
//        r.print(6, 4);
        //for now, fill in by hand...
        double[][] C = new double[4][2];
        C[0][0] = r.get(0, 0);
        C[1][0] = r.get(1, 0);
        C[2][0] = r.get(0, 1);
        C[3][0] = r.get(1, 1);
        C[0][1] = r.get(0, 2);
        C[1][1] = r.get(1, 2);
        C[2][1] = r.get(0, 3);
        C[3][1] = r.get(1, 3);
        return C;
    }

    static double[] TRATSA(double[][] B, double[] S, int N, int M) {
        // s is symmetric matrix packed as lower(?) diagonal
        // R = B'SB
        // N=2, M=4
        Matrix b = new Matrix(N, M);
        Matrix s = new Matrix(N, N); // symmetric
        s.set(0, 0, S[0]);
        s.set(0, 1, S[1]);
        s.set(1, 0, S[1]);
        s.set(1, 1, S[2]);
        if (debug()) {
            s.print(6, 4);
        }
        if (debug()) {
            b.print(6, 4);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[j][i]);
            }
        }
        // if(debug())System.out.println("s");
        if (debug()) {
            s.print(6, 4);
        }
        //if(debug())System.out.println("b");
        if (debug()) {
            b.print(6, 4);
        }

        Matrix r = b.transpose().times(s.times(b));
        // if(debug())System.out.println("r");
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

    static Matrix GEN_ROTMAT(double ANG, int IAX) {
        Matrix r = Matrix.identity(3, 3);
        if (debug()) {
            System.out.println("ANG " + ANG + " IAX " + (IAX + 1));
        }
        double C = cos(ANG);
        double S = sin(ANG);
        if (IAX == 0) {
            r.set(1, 1, C);
            r.set(1, 2, -S);
            r.set(2, 1, S);
            r.set(2, 2, C);
        } else if (IAX == 1) {
            r.set(0, 0, C);
            r.set(0, 2, -S);
            r.set(2, 0, S);
            r.set(2, 2, C);
        } else if (IAX == 2) {
            r.set(0, 0, C);
            r.set(0, 1, -S);
            r.set(1, 0, S);
            r.set(1, 1, C);
        }
        return r;
    }

    static Matrix PROD_ROT(Matrix rx, Matrix ry, Matrix rz) {
        return rx.times(ry.times(rz));
    }

    static List<Hit> GENER_EVT(List<DetectorPlane> dets, double[] par) {
        // PAR = parameters a1,a2,b1,b2 (generated)
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
//         GENER_EVT A   -15.073771       47.742912       0.0000000    
//         GENER_EVT B   7.09767174E-03  7.03085661E-02   1.0000000   
        List<Hit> hitlist = new ArrayList<Hit>();
        int N = dets.size();
        double[] W = {0., 0., 1.};
        //initial guess is along zhat
        double[] A = {par[0], par[1], 0.0000000};
        double[] B = {par[2], par[3], 1.0000000};
        Matrix w = new Matrix(W, 3);
        Matrix a = new Matrix(A, 3);
        Matrix b = new Matrix(B, 3);

//        boolean use_fixed_hits = false;
//
//        if (use_fixed_hits) {
////            double[] W = {0., 0., 1.};
////            double[] A = {-15.073771, 47.742912, 0.0000000};
////            double[] B = {7.09767174E-03, 7.03085661E-02, 1.0000000};
//
//            w = new Matrix(W, 3);
//            a = new Matrix(A, 3);
//            b = new Matrix(B, 3);
//
//            par[0] = A[0];
//            par[1] = A[1];
//            par[2] = B[0];
//            par[3] = B[1];
//
//            if(debug())System.out.println("GENER_EVT A " + A[0] + " " + A[1] + " " + A[2]);
//            if(debug())System.out.println("GENER_EVT B " + B[0] + " " + B[1] + " " + B[2]);
//
//            //List<Hit> hits = new ArrayList<Hit>();
//            hitlist = fixedHits();
//        }
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
//            if (use_fixed_hits) {
//                hit = hitlist.get(i);
//            }
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

//    static List<Hit> fixedHits() {
//        double[] wt = {10000.000, 0.0000000, 24.999998};
//        List<double[]> hits = new ArrayList<double[]>();
//        hits.add(new double[]{-5.9537010, 51.489368});
//        hits.add(new double[]{54.770451, 9.8119898});
//        hits.add(new double[]{19.061852, 61.684181});
//        hits.add(new double[]{77.115250, 13.132210});
//        hits.add(new double[]{-12.350903, 74.988289});
//        hits.add(new double[]{72.895874, 25.362803});
//        hits.add(new double[]{-11.212708, 86.184692});
//        hits.add(new double[]{88.477890, 6.4977794});
//        List<Hit> hitlist = new ArrayList<Hit>();
//        for (int i = 0; i < 8; ++i) {
//            hitlist.add(new Hit(hits.get(i), wt));
//        }
//        return hitlist;
//    }
    static ImpactPoint GET_IMPACT(double[] A, double[] B, Matrix ROT, double[] R0, Matrix wg, double BWG) {
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

}
