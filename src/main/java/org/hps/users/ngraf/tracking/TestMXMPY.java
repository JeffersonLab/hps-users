/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking;

import java.util.Arrays;
import Jama.Matrix;

/**
 *
 * @author ngraf
 */
public class TestMXMPY {

    public static void main(String[] args) {
//        int r = 4;
//        int c = 2;
//        int count = 0;
//        Matrix a = new Matrix(r, c);
//        for (int i = 0; i < r; ++i) {
//            for (int j = 0; j < c; ++j) {
//                a.set(i, j, count++);
//            }
//        }
//        a.print(6, 4);
//        double[] smatrix = asArray(a);
//        System.out.println(Arrays.toString(smatrix));
//
//        
//        double[] C = {1, 2, 1, 2, 1, 3, 1, 4};
//        Matrix matc = new Matrix(4,2);
//        matc.set(0,0,1);
//        matc.set(0,1,2);
//        matc.set(1,0,1);
//        matc.set(1,1,2);
//        matc.set(2,0,1);
//        matc.set(2,1,3);
//        matc.set(3,0,1);
//        matc.set(3,1,4);
//        matc.print(6,4);
//        Matrix cmatT = matc.transpose();
//        cmatT.print(6,4);
//        double[] eps = {1.,1.};
//        eps[0]=1.;
//        eps[1]=1.;
//        double[] wv = new double[4];
//        
//        double[]cmatTarray = asArray(cmatT);
//        
//        System.out.println("eps "+Arrays.toString(eps));
//        System.out.println(" C "+Arrays.toString(C));
//        mxmad(4,C, eps, wv, 3,1,1);
//        System.out.println("wv "+Arrays.toString(wv));
//        mxmad(4,cmatTarray, eps, wv, 3,1,1);
//        System.out.println("wv "+Arrays.toString(wv));
        double[][] C = new double[4][2];
        C[0][0] = 9848.0771;
        C[1][0] = -4.3412046;
        C[2][0] = 1736.4819;
        C[3][0] = 24.620192;

        C[0][1] = 305290.41;
        C[1][1] = -134.57733;
        C[2][1] = 53830.938;
        C[3][1] = 763.22595;
        double[] EPS = {-5.43117523E-03, 0.25400162};

        Matrix c = new Matrix(C, 4, 2);
        Matrix eps = new Matrix(2, 1);
        eps.set(0, 0, EPS[0]);
        eps.set(1, 0, EPS[1]);
        c.print(6, 4);
        eps.print(6, 4);
        Matrix wv = c.times(eps);
        wv.print(6, 4);
        Matrix tst = eps.transpose().times(c.transpose());
        tst.print(6, 4);
        System.out.println("CALLING MXMPY(C, EPS, 4, 2, 1)");
        double[] WV = MXMPY(C, EPS, 4, 2, 1);
        System.out.println("WV " + Arrays.toString(WV));
        System.out.println(" WV should be: -54.589306      -3.1775694      -1692.2686      -98.504639 ");

        myMXMPY(asArray(c), EPS, WV, 4, 2, 1);
        System.out.println("WV " + Arrays.toString(WV));
        myMXMPY(C, EPS, WV, 4, 2, 1);
        System.out.println("WV " + Arrays.toString(WV));
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
//        if (i == 0 || k == 0) {
//            throw new RuntimeException();
//        }

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
