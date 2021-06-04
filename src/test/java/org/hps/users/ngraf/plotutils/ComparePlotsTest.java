package org.hps.users.ngraf.plotutils;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class ComparePlotsTest {

    /**
     * Test of main method, of class ComparePlots.
     */
    @Test
    public void testMain() throws Exception {

        AIDA aida = AIDA.defaultInstance();
        IAnalysisFactory af = aida.analysisFactory();
        ITree tree = af.createTreeFactory().create();
        IHistogramFactory hf = af.createHistogramFactory(tree);

        Random r = new Random();
        for (int i = 0; i < 10000; i++) {
            aida.histogram1D("Histogram 1D", 50, -3, 3).fill(r.nextGaussian());
        }
        aida.saveAs("reference.aida");
        for (int i = 0; i < 10; i++) {
            aida.histogram1D("Histogram 1D", 50, -3, 3).fill(r.nextGaussian());
        }
        aida.saveAs("test.aida");
        String[] args = {"reference.aida", "test.aida","png"};
        ComparePlots.main(args);
    }
}
