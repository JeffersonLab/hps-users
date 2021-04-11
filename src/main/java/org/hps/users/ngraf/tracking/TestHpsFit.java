package org.hps.users.ngraf.tracking;

import Jama.Matrix;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Norman Graf
 */
public class TestHpsFit {

    public static void main(String[] args) {

        int NN = 8;
        double[] ANGLS = {0., 0., 10., 10., 0., 95., 10., 20., 30., 0., 20., 90., 7., 0., 0., 0., 7., 80., 3., 0., 0., 0., 0., 93.};
        double[] XCOOR = {0., 0., 31., 0., 0., 64., 0., 0., 210., 0., 0., 422., 0., 0., 372., 0., 0., 406., 0., 0., 539., 0., 0., 573.};
        double[] SIGS = {0.010, 0.200};     //  detector resolutions in mm

        List<double[]> angles = new ArrayList<double[]>();
        //Top
        angles.add(new double[]{3.1406969851087094, -0.03177473972250014, -1.5707096750665268});
        angles.add(new double[]{-8.763843967516362E-4, 0.030857430620294213, -1.6707969967061551});
        angles.add(new double[]{-3.1371868210063276, -0.031982073016205884, -1.5708002765469602});
        angles.add(new double[]{-0.008253950195711953, 0.031206900985050606, -1.670332814679196});
        angles.add(new double[]{3.140449458676607, -0.02612204537702204, -1.5708809505162877});
        angles.add(new double[]{-0.005357612619242289, 0.025781003210562113, -1.6704142005768545});
        angles.add(new double[]{-3.135680720223676, -0.03153474000599398, -1.570835921868696});
        angles.add(new double[]{0.004079561821680232, 0.034697901407559295, -1.6205748135045286});

        angles.add(new double[]{-3.1358886943959, -0.030816539774356996, -1.5711458236970879});
        angles.add(new double[]{0.006434838235404464, 0.032235412982587405, -1.6204426401435137});
        angles.add(new double[]{3.1410568732583375, -0.030229266015029423, -1.5709559555513355});
        angles.add(new double[]{-2.1930870918929454E-4, 0.031152570497038165, -1.6201278501297305});

        List<double[]> r0 = new ArrayList<double[]>();
        r0.add(new double[]{3.083251888032933, 20.6939715764943, 87.90907758988494});
        r0.add(new double[]{3.3897589906322203, 20.723713416767932, 96.09347236561139});
        r0.add(new double[]{6.209507497482065, 22.251100764548063, 187.979198401254});
        r0.add(new double[]{6.448495884160572, 22.27056067869049, 195.97435895116467});
        r0.add(new double[]{9.266018127938018, 23.77017055031802, 287.706424933812});
        r0.add(new double[]{9.524206517587587, 23.79411639949488, 295.78541118932907});
        r0.add(new double[]{-35.45499313137205, 26.704165491252628, 489.9250042201181});
        r0.add(new double[]{-35.177966825203924, 29.23243054764656, 496.79862110072156});

        r0.add(new double[]{-29.40714364121509, 29.778421220699943, 689.7885600858451});
        r0.add(new double[]{-29.05999883738231, 32.30110797572307, 697.1004544791742});
        r0.add(new double[]{-23.24754288591039, 32.76458780556233, 889.4575319549996});
        r0.add(new double[]{-22.89425800691335, 35.2745708242259, 896.9124291322955});

        List<Matrix> rot = new ArrayList<Matrix>();

        List<DetectorPlane> planes = new ArrayList<DetectorPlane>();

        System.out.println("GENERATING DETECTOR WITH " + NN + " DETECTOR PLANES");
        for (int i = 0; i < NN; ++i) {
            System.out.println("Generating rotation matrices for detector " + (i + 1));
            double[] a = angles.get(i);
            Matrix[] mats = new Matrix[3];
            for (int j = 0; j < 3; ++j) {
                double angl = a[j];
                mats[j] = GEN_ROTMAT(angl, j);
                System.out.println("Rotation matrix for axis " + (j + 1) + " angle " + angl);
                mats[j].print(6, 4);
            }
            Matrix prodrot = PROD_ROT(mats[0], mats[1], mats[2]);
            System.out.println("Combined rotation matrix is");
            prodrot.print(6, 4);
            rot.add(prodrot);
            System.out.println("Translation vector is");
            double[] off = r0.get(i);
            System.out.println(off[0] + " " + off[1] + " " + off[2]);
            System.out.println("Sigs: " + SIGS[0] + " " + SIGS[1]);

            planes.add(new DetectorPlane(i, prodrot, off, SIGS));
        }

        // Testing the track fit code   
        int NTIMES = 1;
        System.out.println("Generating " + NTIMES + " events");
        for (int i = 0; i < NTIMES; ++i) {
            //Generate tracks and hits in ideal detectors
            // contains the generated track parameters
            double[] parin = new double[4];

            List<Hit> hits = GENER_EVT(planes, parin);
            System.out.println("*** END OF EVENT " + (i + 1) + " ***");
            System.out.println("FITTING TRACK");
            TrackFit tf = STR_LINFIT(planes, hits);
            double[] par = tf.pars();
            double[] cov = tf.cov();
            for (int jj = 0; jj < 4; ++jj) {
                double error = sqrt(cov[jj * (jj + 1) / 2]);
                double pull = (par[jj] - parin[jj]) / error;
                System.out.println(" PAR(" + (jj + 1) + ") " + par[jj] + " PARIN(" + (jj + 1) + ") " + parin[jj]);
            }
        }
    }

    static TrackFit STR_LINFIT(List<DetectorPlane> planes, List<Hit> hits) {

        double[][] UVW = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        System.out.println("\n\n\nIN STR_LINFIT\n\n\n");
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
                System.out.println("DETECTOR " + (i + 1));
                DetectorPlane dp = planes.get(i);
                Matrix rot = dp.rot();
                Hit h = hits.get(i);

                Matrix[] uvwg = new Matrix[3];
                for (int j = 0; j < 3; ++j) {
                    System.out.println("  CALLING VMATR");
                    Matrix uvw = new Matrix(UVW[j], 3);
                    System.out.println("  UVW(" + (j + 1) + ") " + uvw.get(0, 0) + " " + uvw.get(1, 0) + " " + uvw.get(2, 0));
                    uvwg[j] = rot.times(uvw);
                    System.out.println("UVWG(" + (j + 1) + ") " + uvwg[j].get(0, 0) + " " + uvwg[j].get(1, 0) + " " + uvwg[j].get(2, 0) + " ");
                    BUVW[j] = b.times(uvwg[j]).get(0, 0);
                }
                System.out.println("   BUVW " + BUVW[0] + " " + BUVW[1] + " " + BUVW[2]);
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
                System.out.println("uvwg[0]");
                uvwg[0].print(6, 4);
                System.out.println("uvwg[1]");
                uvwg[1].print(6, 4);
                System.out.println("uvwg[2]");
                uvwg[2].print(6, 4);

                for (int j = 0; j < NM; ++j) {
                    eps[j] = q[j] - uvm[j]; //predicted minus measured
                    DER[0][j] = uvwg[j].get(0, 0) - uvwg[2].get(0, 0) * BUVW[j] / BUVW[2];
                    DER[1][j] = uvwg[j].get(1, 0) - uvwg[2].get(1, 0) * BUVW[j] / BUVW[2];
                    DER[2][j] = ti * DER[0][j];
                    DER[3][j] = ti * DER[1][j];
                    System.out.println("j " + (j + 1));
                    System.out.println(uvwg[j].get(0, 0) + " " + uvwg[2].get(0, 0) + " " + BUVW[j] + " " + BUVW[2]);
                    System.out.println(uvwg[j].get(1, 0) + " " + uvwg[2].get(1, 0) + " " + BUVW[j] + " " + BUVW[2]);
                    System.out.println(DER[0][j]);
                    System.out.println(DER[1][j]);
                    System.out.println(DER[2][j]);
                    System.out.println(DER[3][j]);
//                DER[0][j] = UVWG(1, J) - UVWG(1, 3) * BUVW(J) / BUVW(3);
//                            uvwg[0].get(0,0)
//                DER[1][j] = UVWG(2, J) - UVWG(2, 3) * BUVW(J) / BUVW(3);
//                DER[2][j] = TI * DER(1, J);
//                DER[3][j] = TI * DER(2, J);
                    NDF = NDF + 1;
                }

                if (NM > 0) {
                    System.out.println("CALLING TRASAT");
                    DCHI = TRASAT(eps, wt, NM);
                    System.out.println(" EPS " + eps[0] + " " + eps[1]);
                    System.out.println(" WT(" + (i + 1) + ") " + wt[0] + " " + wt[1] + " " + wt[2]);
                    System.out.println(" DCHI " + DCHI);
                    CHI = CHI + DCHI;
                    System.out.println("CALLING TRATSA");
                    double[] WP = TRATSA(DER, wt, NM, 4);
                    for (int jj = 0; jj < 2; ++jj) {
                        System.out.println(" DER(" + (jj + 1) + ") " + DER[0][jj] + " " + DER[1][jj] + " " + DER[2][jj] + " " + DER[3][jj] + " ");
                    }
                    System.out.println("wt(" + wt[0] + " " + wt[1] + " " + wt[2]);
                    System.out.println(" WP " + Arrays.toString(WP));
                    System.out.println("CALLING TRATS");
                    double[][] C = TRATS(DER, wt, NM, 4);
                    for (int jj = 0; jj < NM; ++jj) {
                        System.out.println(" C(" + (jj + 1) + ") " + C[0][jj] + " " + C[1][jj] + " " + C[2][jj] + " " + C[3][jj]);
                    }
                    System.out.println(" EPS " + eps[0] + " " + eps[1]);
                    System.out.println("CALLING MXMPY");
                    double[] WV = new double[4]; //MXMPY(C, eps, 4, 2, 1);
                    System.out.println(" WV " + WV[0] + " " + WV[1] + " " + WV[2] + " " + WV[3]);
                    myMXMPY(C, eps, WV, 4, 2, 1);
                    System.out.println(" WV " + WV[0] + " " + WV[1] + " " + WV[2] + " " + WV[3]);
                    System.out.println("CALLING VADD");
                    WMT = VADD(WMT, WP);
                    System.out.println(" WMT " + Arrays.toString(WMT));
                    System.out.println("CALLING VADD");
                    WVT = VADD(WVT, WV);
                    System.out.println("WVT(" + WVT[0] + " " + WVT[1] + " " + WVT[2] + " " + WVT[3]);
                }
            } // end of loop over planes
            if (NDF > 0) {
                System.out.println("CALLING TRSINV");
                COV = TRSINV(WMT, 4);
                System.out.println(" COV " + Arrays.toString(COV));
                System.out.println("CALLING TRSA");
                double[] DPAR = TRSA(COV, WVT, 4); // parameters delta
                System.out.println("DPAR " + Arrays.toString(DPAR));
                System.out.println("CALLING VSUB");
                System.out.println(" PAR " + Arrays.toString(PAR));
                PAR = VSUB(PAR, DPAR);
                System.out.println("PAR-DPAR " + Arrays.toString(PAR));
                System.out.println(" CHI0 " + CHI0 + " CHI " + CHI + " CHI0-CHI " + (CHI0 - CHI));
            }
            // end of while loop
        } while (CHI0 - CHI > 1);

        System.out.println("PAR " + Arrays.toString(PAR));
        System.out.println("COV " + Arrays.toString(COV));
        int count = 0;
        List<double[]> impactPoints = new ArrayList<double[]>();
        for (ImpactPoint ip : rx) {
            System.out.println("RX( " + (count + 1) + " )" + Arrays.toString(ip.r()));
            impactPoints.add(ip.r());
            count++;
        }
        System.out.println("CHI " + CHI);
        System.out.println("NDF " + NDF);
        System.out.println("NIT " + NIT);
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

        s.print(6, 4);
        Matrix a = new Matrix(4, 1);
        a.set(0, 0, A[0]);
        a.set(1, 0, A[1]);
        a.set(2, 0, A[2]);
        a.set(3, 0, A[3]);
        a.print(6, 4);

        Matrix d = s.times(a);
        d.print(6, 4);

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

        a.print(6, 4);

        Matrix b = a.inverse();
        b.print(6, 4);

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
        s.print(6, 4);
        b.print(6, 4);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[j][i]);
            }
        }
        // System.out.println("s");
//        s.print(6, 4);
        //System.out.println("b");
//        b.print(6, 4);
//TODO figure out what is going on here...
        Matrix r = s.times(b);
        // System.out.println("r");
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
        s.print(6, 4);
        b.print(6, 4);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                b.set(i, j, B[j][i]);
            }
        }
        // System.out.println("s");
        s.print(6, 4);
        //System.out.println("b");
        b.print(6, 4);

        Matrix r = b.transpose().times(s.times(b));
        // System.out.println("r");
        r.print(6, 4);
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

    static Matrix GEN_ROTMAT(double ANG, int IAX) {
        Matrix r = Matrix.identity(3, 3);
        System.out.println("ANG " + ANG + " IAX " + (IAX + 1));
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
        // should be randomly generating track parameters, but for now lets used a fixed set
//         GENER_EVT A   -15.073771       47.742912       0.0000000    
//         GENER_EVT B   7.09767174E-03  7.03085661E-02   1.0000000   
        double[] W = {0., 0., 1.};
        double[] A = {-15.073771, 47.742912, 0.0000000};
        double[] B = {7.09767174E-03, 7.03085661E-02, 1.0000000};

        Matrix w = new Matrix(W, 3);
        Matrix a = new Matrix(A, 3);
        Matrix b = new Matrix(B, 3);

        par[0] = A[0];
        par[1] = A[1];
        par[2] = B[0];
        par[3] = B[1];

        System.out.println("GENER_EVT A " + A[0] + " " + A[1] + " " + A[2]);
        System.out.println("GENER_EVT B " + B[0] + " " + B[1] + " " + B[2]);

        //List<Hit> hits = new ArrayList<Hit>();
        int N = dets.size();
        List<Hit> hitlist = fixedHits();
        for (int i = 0; i < N; ++i) {
            Hit hit = hitlist.get(i);
            DetectorPlane dp = dets.get(i);
            Matrix ROT = dp.rot();
            double[] R0 = dp.r0();
            System.out.println("CALLING VMATR(W, ROT(1," + (i + 1) + "), WG, 3, 3) ");
            Matrix wg = ROT.times(w);
            System.out.println(" WG " + wg.get(0, 0) + " " + wg.get(1, 0) + " " + wg.get(2, 0));
            Matrix bwg = b.transpose().times(wg);
            System.out.println(" BWG " + bwg.get(0, 0));
            ImpactPoint ip = GET_IMPACT(A, B, ROT, R0, wg, bwg.get(0, 0));
            for (int j = 0; j < 2; ++j) // do fluctuations
            {
                System.out.println("MEASUREMENT Q(" + (j + 1) + ") " + ip.q()[j]);
                System.out.println("SMEARED MEASUREMENT " + (j + 1) + ") " + hit.uvm()[j]);
            }
            System.out.println("SMEARED UVM(JJ," + (i + 1) + ") " + hit.uvm()[0] + " " + hit.uvm()[1]);
            System.out.println("WT(JJ, " + (i + 1) + ") " + hit.wt()[0] + " " + hit.wt()[1] + " " + hit.wt()[2]);
        }
        return hitlist;
    }

    static List<Hit> fixedHits() {
        double[] wt = {10000.000, 0.0000000, 24.999998};
        List<double[]> hits = new ArrayList<double[]>();
        hits.add(new double[]{-5.9537010, 51.489368});
        hits.add(new double[]{54.770451, 9.8119898});
        hits.add(new double[]{19.061852, 61.684181});
        hits.add(new double[]{77.115250, 13.132210});
        hits.add(new double[]{-12.350903, 74.988289});
        hits.add(new double[]{72.895874, 25.362803});
        hits.add(new double[]{-11.212708, 86.184692});
        hits.add(new double[]{88.477890, 6.4977794});
        List<Hit> hitlist = new ArrayList<Hit>();
        for (int i = 0; i < 8; ++i) {
            hitlist.add(new Hit(hits.get(i), wt));
        }
        return hitlist;
    }

    static ImpactPoint GET_IMPACT(double[] A, double[] B, Matrix ROT, double[] R0, Matrix wg, double BWG) {
        double[] AMR = VSUB(A, R0);
        Matrix amr = new Matrix(AMR, 3);
        System.out.println("GET_IMPACT AMR " + AMR[0] + " " + AMR[1] + " " + AMR[2]);
        double ti = -amr.transpose().times(wg).get(0, 0) / BWG;
        System.out.println("GET_IMPACT PARAMETER AT IMPACT TI " + ti);
        double[] APU = VLINE(AMR, 1., B, ti);
        System.out.println("GET_IMPACT APU " + APU[0] + " " + APU[1] + " " + APU[2]);
        double[] q = VMATL(ROT, APU);
        System.out.println("GET_IMPACT LOCAL U,V,W Q " + q[0] + " " + q[1] + " " + q[2]);
        double[] r = VADD(APU, R0);
        System.out.println("GET_IMPACT GLOBAL X,Y,Z R " + r[0] + " " + r[1] + " " + r[2]);
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

    static void mxmad(int n_, double[] a, double[] b, double[] c, int i, int j, int k) {
        /* Local variables */
        int l, m, n, ia, ic, ib, ja, jb, iia, iib, ioa, iob;

        /* Parameter adjustments */
//?        --a;  --b;  --c;                                     
        /* Function Body */
 /*                      MXMAD MXMAD1 MXMAD2 MXMAD3 MXMPY MXMPY1 MXMPY2 MXMPY3 MXMUB MXMUB1 MXMUB2 MXMUB3 */
 /*  const int iandj1[] = {21,   22,    23,    24,   11,    12,    13,    14,    31,   32,   33,    34 }; */
        int[] iandj1 = {2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3};
        int[] iandj2 = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
        int n1 = iandj1[n_];
        int n2 = iandj2[n_];
        if (i == 0 || k == 0) {
            throw new RuntimeException();
        }

        switch (n2) {
            case 1:
                iia = 1;
                ioa = j;
                iib = k;
                iob = 1;
                break;
            case 2:
                iia = 1;
                ioa = j;
                iib = 1;
                iob = j;
                break;
            case 3:
                iia = i;
                ioa = 1;
                iib = k;
                iob = 1;
                break;
            case 4:
                iia = i;
                ioa = 1;
                iib = 1;
                iob = j;
                break;
            default:
                iia = ioa = iib = iob = 0;
                assert (iob != 0);
        };

        ia = 1;
        ic = 1;
        for (l = 1; l <= i; ++l) {
            ib = 1;
            for (m = 1; m <= k; ++m, ++ic) {
                switch (n1) {
                    case 1:
                        c[ic] = 0.;
                        break;
                    case 3:
                        c[ic] = -c[ic];
                        break;
                };
                if (j == 0) {
                    continue;
                }
                ja = ia;
                jb = ib;
                double cic = c[ic];
                for (n = 1; n <= j; ++n, ja += iia, jb += iib) {
                    cic += a[ja] * b[jb];
                }
                c[ic] = cic;
                ib += iob;
            }
            ia += ioa;
        }
    }
}
