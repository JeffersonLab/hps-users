package org.hps.users.ngraf.plotutils;

/**
 *
 * @author Norman A. Graf
 */
import hep.aida.IAnalysisFactory;
import hep.aida.ICloud1D;
import hep.aida.ICloud2D;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.IPlotterFactory;
import hep.aida.ITree;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.hps.monitoring.plotting.ExportPdf;

import org.lcsim.util.aida.AIDA;

public class AidaTest {

    public static void main(String args[]) throws Exception {
        AIDA aida = AIDA.defaultInstance();
        IAnalysisFactory analysisFactory = aida.analysisFactory();
        System.out.println("analysisFactory: " + analysisFactory.getClass().getCanonicalName());
        IPlotterFactory plotterFactory = analysisFactory.createPlotterFactory();
        System.out.println("plotterFactory: " + plotterFactory.getClass().getCanonicalName());
        IPlotter plotter = plotterFactory.create();
        System.out.println("unnamed plotter: " + plotter.getClass().getCanonicalName());
        plotter.createRegion();
        System.out.println("region: " + plotter.region(0).getClass().getCanonicalName());
        plotterFactory = analysisFactory.createPlotterFactory("dummy");
        System.out.println("named plotterFactory: " + plotterFactory.getClass().getCanonicalName());
        plotter = plotterFactory.create("dummy");
        System.out.println("named plotter: " + plotter.getClass().getCanonicalName());

        //
        IAnalysisFactory af = IAnalysisFactory.create();
        ITree tree = af.createTreeFactory().create();
        IHistogramFactory hf = af.createHistogramFactory(tree);

        tree.mkdir("/Histograms");
        tree.cd("/Histograms");

        IHistogram1D h1 = hf.createHistogram1D("Histogram 1D", 50, -3, 3);
        IHistogram2D h2 = hf.createHistogram2D("Histogram 2D", 40, -3, 3, 40, -3, 3);

        tree.mkdir("/Clouds");
        tree.cd("/Clouds");

        ICloud1D c1 = hf.createCloud1D("Cloud 1D");
        ICloud2D c2 = hf.createCloud2D("Cloud 2D");

        plotter = af.createPlotterFactory().create("CreateAndPlotHistograms.java plot");

        String[] params = plotter.availableParameters();
        for (String s : params) {
            System.out.println("param " + s + " " + plotter.parameterValue(s));

            String[] paramVals = plotter.availableParameterOptions(s);
            for (String v : paramVals) {
                System.out.println("val: " + v + " " + plotter.parameterValue(v));
            }
        }
        plotter.setParameter("plotterWidth", "792");
        plotter.setParameter("plotterHeight", "459");
        
        System.out.println("plotter width " + plotter.parameterValue("plotterWidth"));
        plotter.show();
        plotter.createRegions(2, 2);

        plotter.region(0).plot(h1);
        plotter.region(1).plot(h2);
        plotter.region(2).plot(c1);
        plotter.region(3).plot(c2);

        Random r = new Random();

        for (int i = 0; i < 100000; i++) {
            h1.fill(r.nextGaussian());
            h2.fill(r.nextGaussian(), r.nextGaussian());
            c1.fill(r.nextGaussian());
            c2.fill(r.nextGaussian(), r.nextGaussian());
        }

        String fileName = "AidaTest.pdf";
        List<String> runData = new ArrayList<>();
        runData.add("Hello WOrld!");
        List<IPlotter> plotters = new ArrayList<>();
        plotters.add(plotter);

//        ExportPdf.write(plotters, fileName, runData);
    }
}
