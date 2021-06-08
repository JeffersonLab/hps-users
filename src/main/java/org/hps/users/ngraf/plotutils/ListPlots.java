package org.hps.users.ngraf.plotutils;

import hep.aida.IAnalysisFactory;
import hep.aida.IBaseHistogram;
import hep.aida.IHistogramFactory;
import hep.aida.IManagedObject;
import hep.aida.ITree;
import java.io.File;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class ListPlots {

    public static void main(String[] args) throws Exception {
        String tmp = args[0];
        File inputFile = new File(tmp);

        AIDA aida = AIDA.defaultInstance();
        IAnalysisFactory af = aida.analysisFactory();
        IHistogramFactory hf = aida.histogramFactory();
        ITree myTree = aida.tree();

        ITree srcTree = af.createTreeFactory().create(inputFile.getAbsolutePath());
        // get the list of histograms in this file 
        String[] objectTypes = srcTree.listObjectTypes("/", true);
        String[] objectNames = srcTree.listObjectNames("/", true);
        for (int pathIndex = 0; pathIndex < objectNames.length; pathIndex++) {
            if (objectTypes[pathIndex].startsWith("IHistogram")) {
                String histogramName = objectNames[pathIndex];
                System.out.println("processing: " + histogramName);
                IBaseHistogram src = (IBaseHistogram) srcTree.find(histogramName);
                String path = histogramName.substring(0, histogramName.lastIndexOf('/'));
//                System.out.println(path);
                IManagedObject object = srcTree.find(histogramName);
                String filePath = srcTree.findPath(object);
//                System.out.println(filePath);
            }

        }
    }

}
