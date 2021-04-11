/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.alignment;

import java.util.Arrays;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

/**
 *
 * @author ngraf
 */
public class FailedTestAlignDebug {

    public static void main(String[] args) throws Exception {

        int NN = 8;
        double[][] ROT = new double[NN][9]; //rotation matrices local -> global
        double[][] R0 = new double[NN][3]; // detector origin in global coordinates (mm)
        double[] SIGS = new double[2]; //detector resolutions in u and v (mm)
        GENER_DET(NN, ROT, R0, SIGS);
    }

    static void GENER_DET(int NN, double[][] ROT, double[][] R0, double[] SIGS) {

        double[][] ANGLS = {{0., 0., 10.}, {10., 0., 95.}, {10., 20., 30.}, {0., 20., 90.}, {7., 0., 0.}, {0., 7., 80.}, {3., 0., 0.}, {0., 0., 93.}};
        double[][] XCOOR = {{0., 0., 31.}, {0., 0., 64.}, {0., 0., 210.}, {0., 0., 422.}, {0., 0., 372.}, {0., 0., 406.}, {0., 0., 539.}, {0., 0., 573.}};

        double[][][] RAX = new double[3][3][3];
        UCOPY(XCOOR, R0);
        SIGS[0] = 0.010;  //detector u resolution in mm
        SIGS[1] = 0.200;  //detector v resolution in mm
        for (int i = 0; i < NN; ++i) {
            System.out.println("generating rotation matrices for detector " + (i+1));
            for (int j = 0; j < 3; ++j) {
                double ANGL = Math.toRadians(ANGLS[i][j]);
                GEN_ROTMAT(ANGL, j, RAX[j]);
                System.out.println("rotation matrix for axis " + (j+1) + " angle " + ANGL);
                //System.out.println(Arrays.deepToString(RAX));
                printmat(RAX[j]);
            }
            double[][] PROT = new double[3][3];
            PROD_ROT(RAX[2], RAX[1], RAX[0], PROT);
        }
    }
    
    static void PROD_ROT(double[][] RA, double[][] RB, double[][] RC, double[][]RTOT)
    {
       double[][] RAPU = new double[3][3];
       MXMPY(RA, RB, RAPU, 3, 3, 3);
       System.out.println("PROD_ROT RAPU rotation matrix:");
        printmat(RAPU);
       MXMPY(RAPU, RC, RTOT,3,3,3);
        System.out.println("PROD_ROT combined rotation matrix:");
        printmat(RTOT);
    }

    static void MXMPY(double[][] A, double[][] B, double[][] C, int I, int J, int K) {
        int rdim = A.length;
        int cdim = A[0].length;
        double[] a = new double[rdim * cdim];
        double[] b = new double[rdim * cdim];
        double[] c = new double[rdim * cdim];
        int count = 0;
        for (int j = 0; j < cdim; ++j) {
            for (int i = 0; i < rdim; ++i) {
                a[count] = A[i][j];
                b[count] = B[i][j];
                count++;
            }
        }
        MXMPY(a, b, c, I, J, K);
        count = 0;
        for(int i=0; i<rdim; ++i){
            for(int j=0; j<cdim; ++j)
            {
                C[i][j]=c[count++];
            }
        }
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
                    } //ENDDO
                    IB = IB + IOB;
                    IC = IC + 1;
                }
            }
            IA = IA + IOA;
        }
    }
    
    static void GEN_ROTMAT(double ANG, int IAX, double[][] ROTM) {
//Generate rotation matrix around axis IAX
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
            ROTM[2][0] = S;
            ROTM[0][2] = -S;
            ROTM[2][2] = C;
        } else if (IAX == 2) {
            ROTM[0][0] = C;
            ROTM[1][0] = S;
            ROTM[0][1] = -S;
            ROTM[0][0] = C;
        }
    }

    static void UCOPY(double[][] source, double[][] destination) {
        for (int i = 0; i < source.length; ++i) {
            System.arraycopy(source[i], 0, destination[i], 0, destination[i].length);
        }
        // displaying destination array
        System.out.println(Arrays.deepToString(destination));
    }

    static void printmat(double[][] m) {
        for (int i = 0; i < m.length; ++i) {
            System.out.println(Arrays.toString(m[i]));
        }
    }
}
