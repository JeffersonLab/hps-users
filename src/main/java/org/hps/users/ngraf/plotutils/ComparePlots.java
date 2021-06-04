package org.hps.users.ngraf.plotutils;

import hep.aida.IAnalysisFactory;
import hep.aida.IBaseHistogram;
import hep.aida.IDataStyle;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.ITree;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class ComparePlots {

    public static void main(String[] args) throws Exception {
        String testPlots = args[0];
        File testPlotsInputFile = new File(testPlots);
        String refPlots = args[1];
        File refPlotsInputFile = new File(refPlots);

        boolean writeToFile = false;
        String fileType = "png";
        if (args.length > 2) {
            writeToFile = true;
            fileType = args[2];
        }
        AIDA aida = AIDA.defaultInstance();
        IAnalysisFactory analysisFactory = aida.analysisFactory();
        IHistogramFactory hf = aida.histogramFactory();

        ITree testSourceTree = analysisFactory.createTreeFactory().create(testPlotsInputFile.getAbsolutePath());
        ITree refSourceTree = analysisFactory.createTreeFactory().create(refPlotsInputFile.getAbsolutePath());

        List<IPlotter> plotters = new ArrayList<>();

        // get the list of histograms in this file 
        String[] objectTypes = testSourceTree.listObjectTypes("/", true);
        String[] objectNames = testSourceTree.listObjectNames("/", true);
        for (int pathIndex = 0; pathIndex < objectNames.length; pathIndex++) {
            if (objectTypes[pathIndex].startsWith("IHistogram")) {
                String histogramName = objectNames[pathIndex];
                String[] tokens = histogramName.split("/");
                String plotterTitle = tokens[1];
                if (tokens.length > 2) {
                    for (int i = 2; i < tokens.length; ++i) {
                        plotterTitle += " ";
                        plotterTitle += tokens[i];
                    }
                }
                System.out.println("processing: " + histogramName);
                IBaseHistogram test = (IBaseHistogram) testSourceTree.find(histogramName);
                IBaseHistogram ref = (IBaseHistogram) refSourceTree.find(histogramName);
                IPlotter plotter = analysisFactory.createPlotterFactory().create();
                IDataStyle dataStyle = plotter.style().dataStyle();
                dataStyle.fillStyle().setVisible(false);
                dataStyle.errorBarStyle().setVisible(false);
                if (objectTypes[pathIndex].equals("IHistogram2D")) {
                    plotter.region(0).plot(ref);
                    plotter.region(0).plot(test);
                    plotter.setTitle(plotterTitle);
                    if (!writeToFile) {
                        plotter.show();
                    }
                    plotters.add(plotter);
                } else if (objectTypes[pathIndex].equals("IHistogram1D")) {
                    plotter.destroyRegions();
                    plotter.createRegion(0, .75, 1, .25);
                    plotter.createRegion(0, 0, 1, .75);
                    plotter.region(1).plot(ref);
                    plotter.region(1).plot(test);
                    plotter.region(0).plot(hf.subtract("test-reference", (IHistogram1D) test, (IHistogram1D) ref));
//                    plotter.region(0).style().statisticsBoxStyle().setVisible(false);
                    plotter.setTitle(plotterTitle);
                    if (!writeToFile) {
                        plotter.show();
                    }
                    plotters.add(plotter);
                }
            }
        }
        System.out.println("Added " + plotters.size() + " plotters");

        if (writeToFile) {
            for (IPlotter p : plotters) {
                p.writeToFile(p.title() + "." + fileType);
            }
        }
    }
}
