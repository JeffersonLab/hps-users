/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.ngraf.tracking;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author ngraf
 */
public class testRotDer {

    public static void main(String[] args) {
        double[] ROT = {0.81379771, 0.44096959, 0.37852234, -0.46984631, 0.88256413, -1.80282891E-02, -0.34202015, -0.16317593, 0.92541653};
        List<double[]> DROT = ROT_DER(ROT);
        for(int i=0; i<DROT.size(); ++i)
        {
        System.out.println("DROT " + i+" "+ Arrays.toString(DROT.get(i)));
        }
    }

    static List<double[]> ROT_DER(double[] ROT) {
        double[][] DA = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0},
            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
        List<double[]> l = new ArrayList<double[]>();
        for (int i = 0; i < 3; ++i) {
            double[] c = new double[9];
            System.out.println("ROT "+Arrays.toString(ROT));
            System.out.println("DA("+i+") "+Arrays.toString(DA[i]));
            
            myMXMPY(ROT, DA[i], c, 3, 3, 3);
            System.out.println("DROT "+Arrays.toString(c));
            l.add(c);
        }
        return l;
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
}
