package org.hps.users.ngraf.aidaanalysis;

import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IFunctionFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.ITree;
import hep.aida.ref.plotter.PlotterRegion;
import java.io.File;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class WabTrackEfficiencyAnalyzer {

    public static void main(String[] args) throws Exception {
        String tmp = "D://work/hps/analysis/physrun2019/wabAnalysis/20210607_TrackingStudies/hps_010031_TrackingStudies.aida";
        if (args.length > 0) {
            tmp = args[0];
        }
        File inputFile = new File(tmp);
        AIDA aida = AIDA.defaultInstance();
        IAnalysisFactory af = aida.analysisFactory();
        IHistogramFactory hf = aida.histogramFactory();

        ITree tree = af.createTreeFactory().create(inputFile.getAbsolutePath());
        String[] objectTypes = tree.listObjectTypes("/", true);
        String[] objectNames = tree.listObjectNames("/", true);
        for (int pathIndex = 0; pathIndex < objectNames.length; pathIndex++) {
            if (objectTypes[pathIndex].startsWith("IHistogram")) {
                String histogramName = objectNames[pathIndex];
                System.out.println(histogramName);
            }
        }

        IHistogram1D eg = (IHistogram1D) tree.find("/FinalStateParticles WAB Tracking Analysis/Cluster esum GBL eg");
        IHistogram1D ge = (IHistogram1D) tree.find("/FinalStateParticles WAB Tracking Analysis/Cluster esum GBL ge");
        IHistogram1D gg = (IHistogram1D) tree.find("/FinalStateParticles WAB Tracking Analysis/Cluster esum GBL gg");

        IHistogram1D egge = hf.add("eg + ge", eg, ge);

        // get the contents of eg+ge
        double maxBinHeight = egge.maxBinHeight();
        double maxBinHeightGG = gg.maxBinHeight();
        double scaleFactor = maxBinHeightGG / maxBinHeight;

        int eggeNumberOfBins = eg.axis().bins();
        System.out.println("hist dimension " + eggeNumberOfBins);

//        double[] eggeScaledContents = new double[egge.dimension()];
//        double[] eggeScaledValue = new double[egge.dimension()];
        IHistogram1D eggeScaled = hf.createCopy("Scaled eg+ge", egge);
        eggeScaled.setTitle("Scaled eg+ge");
        eggeScaled.reset();
        for (int i = 1; i < eggeNumberOfBins; ++i) {
//            eggeScaledContents[i - 1] = egge.binEntries(i);
            for (int j = 0; j < egge.binEntries(i) * scaleFactor; ++j) {
                eggeScaled.fill(egge.binMean(i - 1), 1.);
            }
        }

        IHistogram1D background = hf.subtract("Background: gg - eg+ge scaled", gg, eggeScaled);

        IPlotter plotter = af.createPlotterFactory().create();
        plotter.createRegions(2, 3);
        plotter.region(0).plot(eg);
        plotter.region(1).plot(ge);

        plotter.region(2).plot(gg);
        plotter.region(3).plot(egge);
        plotter.region(4).plot(eggeScaled);
        plotter.region(5).plot(gg);
        plotter.region(5).plot(eggeScaled);

        //        String[] params = plotter.region(0).style().xAxisStyle().availableParameters();
        //        System.out.println(Arrays.toString(params));
        // [isVisible, labelVertical, label, upperLimit, lowerLimit, type, allowZeroSuppression, scale]
        String parameter = "lowerLimit";

        IPlotter plotter2 = af.createPlotterFactory().create();
        plotter2.createRegions(2, 2);
        plotter2.region(0).style().xAxisStyle().setParameter("lowerLimit", "3.4");
        plotter2.region(0).style().xAxisStyle().setParameter("upperLimit", "5.5");

        plotter2.region(1).style().xAxisStyle().setParameter(parameter, "3.4");
        plotter2.region(1).style().xAxisStyle().setParameter("upperLimit", "5.5");

        plotter2.region(2).style().xAxisStyle().setParameter(parameter, "3.4");
        plotter2.region(2).style().xAxisStyle().setParameter("upperLimit", "5.5");

        plotter2.region(3).style().xAxisStyle().setParameter(parameter, "3.4");
        plotter2.region(3).style().xAxisStyle().setParameter("upperLimit", "5.5");

        plotter2.region(0).plot(gg);
        plotter2.region(1).plot(eggeScaled);
        plotter2.region(2).plot(gg);
        plotter2.region(2).plot(eggeScaled);
        plotter2.region(3).plot(background);

        // fit the esum plot
        IFitFactory ff = af.createFitFactory();
        IFunctionFactory functionFactory = af.createFunctionFactory(tree);
        IFitter jminuit = ff.createFitter("Chi2", "jminuit");
        IFunction gauss = functionFactory.createFunctionByName("gauss", "G");

        gauss.setParameter("amplitude", egge.maxBinHeight());
        gauss.setParameter("mean", 4.5);//egge.mean());
        gauss.setParameter("sigma", .1);//egge.rms());

        IFitResult jgaussResult = jminuit.fit(egge, gauss, "range = \"(4.45,5.0)\"");

        IPlotter plotter3 = af.createPlotterFactory().create();
        plotter3.region(0).style().xAxisStyle().setParameter(parameter, "3.4");
        plotter3.region(0).style().xAxisStyle().setParameter("upperLimit", "5.5");
        plotter3.region(0).plot(egge);
        plotter3.region(0).plot(jgaussResult.fittedFunction(), "range = \"(4.45,5.0)\"");
        ((PlotterRegion) plotter3.region(0)).getPlot().getStats().setVisible(true);
        // show the results
        plotter.show();
        plotter2.show();
        plotter3.show();
    }

}
