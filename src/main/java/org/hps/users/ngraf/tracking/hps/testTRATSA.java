/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking.hps;

import Jama.Matrix;
import java.util.Arrays;
import static org.hps.users.ngraf.tracking.hps.LineFit.TRATS;

/**
 *
 * @author ngraf
 */
public class testTRATSA {

    public static void main(String[] args) {

        double[][] DER = new double[4][2];
        DER[0][0] = -4.85610478E-02;
        DER[1][0] = -0.99878895;
        DER[2][0] = -157.12111;
        DER[3][0] = -3231.6194;
        for (int jj = 0; jj < 2; ++jj) {

            System.out.println(" DER(" + (jj + 1) + ") " + DER[0][jj] + " " + DER[1][jj] + " " + DER[2][jj] + " " + DER[3][jj] + " ");

        }
        double[] wt = {27777.777, 0.0000000, 0.0000000};
        System.out.println(" wt = " + wt[0] + " " + wt[1] + " " + wt[2]);
        int NM = 1;
        double[] WP = TRATSA(DER, wt, NM, 4);
        System.out.println("WP " + Arrays.toString(WP));

        double[][] C = TRATS(DER, wt, NM, 4);
        for (int jj = 0; jj < NM; ++jj) {
            if (debug()) {
                System.out.println(" C(" + (jj + 1) + ") " + C[0][jj] + " " + C[1][jj] + " " + C[2][jj] + " " + C[3][jj]);
            }
        }
    }

    static double[] TRATSA(double[][] B, double[] S, int N, int M) {
        // s is symmetric matrix packed as lower(?) diagonal
        // R = B'SB
        // N=2, M=4
        Matrix b = new Matrix(N, M);
        Matrix s = new Matrix(N, N); // symmetric
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                s.set(i, j, S[i + j]);
//                s.set(0, 0, S[0]);
//                s.set(0, 1, S[1]);
//                s.set(1, 0, S[1]);
//                s.set(1, 1, S[2]);
            }
        }
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

    static double[][] TRATS(double[][] B, double[] S, int N, int M) {
        // s is symmetric matrix packed as lower diagonal
        // R = B'S
        // N=2, M=4 
        Matrix b = new Matrix(N, M);
        Matrix s = new Matrix(N, N); // symmetric
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                s.set(i, j, S[i + j]);
            }
        }
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
        // System.out.println("s");
//        s.print(6, 4);
        //System.out.println("b");
//        b.print(6, 4);
//TODO figure out what is going on here...
        Matrix r = s.times(b);
        System.out.println("r");
        r.print(6, 4);
        //for now, fill in by hand...
        double[][] C = new double[4][2];

               if(N==2)
        {
        C[0][0] = r.get(0, 0);
        C[1][0] = r.get(1, 0);
        C[2][0] = r.get(0, 1);
        C[3][0] = r.get(1, 1);
        C[0][1] = r.get(0, 2);
        C[1][1] = r.get(1, 2);
        C[2][1] = r.get(0, 3);
        C[3][1] = r.get(1, 3);
        }
        if(N==1)
        {
        C[0][0] = r.get(0, 0);
        C[1][0] = r.get(0, 1);
        C[2][0] = r.get(0, 2);
        C[3][0] = r.get(0, 3);
        }
        return C;
    }

    public static boolean debug() {
        return true;
    }

}
