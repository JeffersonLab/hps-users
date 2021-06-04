package org.hps.users.ngraf.plotutils;

import java.util.Random;

import org.junit.Test;
import org.lcsim.util.aida.AIDA;

import junit.framework.TestCase;

/**
 * Test of generating and writing comparison plots
 */
public class ComparePlotsTest extends TestCase {

    /**
     * Override AIDA setup from hps-java so graphics can be generated and shown,
     * which is generally not enabled for test cases.
     */
    static {
        System.getProperties().setProperty("hep.aida.IAnalysisFactory", "hep.aida.ref.AnalysisFactory");
        System.getProperties().setProperty("java.awt.headless", "false");
    }

    /**
     * Test of main method, of class ComparePlots.
     */
    @Test
    public void testMain() throws Exception {
        AIDA aida = AIDA.defaultInstance();
        Random r = new Random();
        for (int i = 0; i < 10000; i++) {
            aida.histogram1D("Histogram 1D", 50, -3, 3).fill(r.nextGaussian());
        }
        aida.saveAs("reference.aida");
        for (int i = 0; i < 10; i++) {
            aida.histogram1D("Histogram 1D", 50, -3, 3).fill(r.nextGaussian());
        }
        aida.saveAs("test.aida");
        String[] args = {"reference.aida", "test.aida", "png"};
        ComparePlots.main(args);
    }
}
