package org.hps.users.ngraf.calibration;

import hep.aida.IAnalysisFactory;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IFunctionFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.IPlotterStyle;
import hep.aida.ITree;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SvtOnlineBaselineCalibrationFitter {

    private String _histogramFileName;
    private AIDA _aida = AIDA.defaultInstance();
    private IAnalysisFactory _af;
    private IHistogramFactory _hf;
    private ITree _tree;
    private IHistogramFactory _outputhf;
    private ITree _outputTree;
    private IFunctionFactory _functionFactory;
    private IPlotterStyle _functionStyle;
    private IFitter _jminuit;
    // standard style for histogram data
    IPlotterStyle _dataStyle;

    private Map<String, IHistogram2D> _histograms = new TreeMap<>();
    private List<String> _failedSlices = new ArrayList<>();
    private int _successfulSlices;
    private int _runNumber;

    boolean writeEmOut = false;
    boolean showem = false;

    private int _minSliceEntries = 1000;

    public SvtOnlineBaselineCalibrationFitter(String histogramFileName) {
        _runNumber = Integer.parseInt(histogramFileName.substring(histogramFileName.indexOf("hpssvt_") + 7, histogramFileName.indexOf("hpssvt_") + 7 + 6));
        _histogramFileName = histogramFileName;
        File inputFile = new File(_histogramFileName);
        _af = _aida.analysisFactory();
        _hf = _aida.histogramFactory();
        try {
            _tree = _af.createTreeFactory().create(inputFile.getAbsolutePath());
        } catch (IllegalArgumentException | IOException ex) {
            Logger.getLogger(SvtOnlineBaselineCalibrationFitter.class.getName()).log(Level.SEVERE, "File " + _histogramFileName + " not found", ex);
        }
        try {
            _outputTree = _af.createTreeFactory().create("failedFits_" + _runNumber + ".aida", "xml", false, true);
        } catch (IOException | IllegalArgumentException ex) {
            Logger.getLogger(SvtOnlineBaselineCalibrationFitter.class.getName()).log(Level.SEVERE, "Could not create failedFits.aida", ex);

        }
        _outputhf = _af.createHistogramFactory(_outputTree);
        _functionFactory = _af.createFunctionFactory(_tree);
        _jminuit = _af.createFitFactory().createFitter("Chi2", "jminuit");
        String[] objectTypes = _tree.listObjectTypes("/", true);
        String[] objectNames = _tree.listObjectNames("/", true);
        for (int pathIndex = 0; pathIndex < objectNames.length; pathIndex++) {
            if (objectTypes[pathIndex].startsWith("IHistogram2D")) {
                String histogramName = objectNames[pathIndex];
                _histograms.put(((IHistogram2D) _tree.find(objectNames[pathIndex])).title(), (IHistogram2D) _tree.find(objectNames[pathIndex]));
            }
        }
        for (String s : _histograms.keySet()) {
            System.out.println(s + " " + _histograms.get(s).title());
        }

        _functionStyle = _af.createPlotterFactory().createPlotterStyle();
        _functionStyle.dataStyle().outlineStyle().setColor("RED");
        _functionStyle.dataStyle().lineStyle().setThickness(7);

        _dataStyle = _af.createPlotterFactory().createPlotterStyle();
        _dataStyle.dataStyle().lineStyle().setVisible(true);
        _dataStyle.dataStyle().lineStyle().setColor("BLACK");
        _dataStyle.dataStyle().fillStyle().setColor("RED");
        _dataStyle.dataStyle().fillStyle().setVisible(false);

    }

    public static void main(String[] args) {
        String fileName = "D://work/hps/analysis/physrun2019/unblinded/svtBaselineCalibration/hpssvt_010629_onlinebaseline.aida";
        if (args.length > 0) {
            fileName = args[0];
        }
        SvtOnlineBaselineCalibrationFitter svt = new SvtOnlineBaselineCalibrationFitter(fileName);
        long startTime = System.nanoTime();

        svt.analyzeAll();

        long endTime = System.nanoTime();
        long duration = (endTime - startTime) / 1000000;  //divide by 1000000 to get milliseconds.

        svt.finish();
        //
        System.out.println("processing took " + duration + " ms");
    }

    private void analyzeAll() {
        for (String s : _histograms.keySet()) {
//            System.out.println("\n \n \n \n analyzing " + s);
            analyzeOneHistogram(s);
        }
    }

    private void analyzeOneHistogram(String histogramName) {
        analyzeOneHistogram(_histograms.get(histogramName));
    }

    private void analyzeOneHistogram(IHistogram2D hist) {
        IPlotter plotter = null;
        if (writeEmOut || showem) {
            plotter = _af.createPlotterFactory().create(hist.title() + " slices");
            plotter.setParameter("plotterWidth", "1600");
            plotter.setParameter("plotterHeight", "1000");
            plotter.createRegions(4, 4);
            plotter.style().dataStyle().errorBarStyle().setVisible(false);
            plotter.style().legendBoxStyle().setVisible(false);
        }
        if (showem && !writeEmOut) {
            plotter.show();
        }
        String histo = hist.title();
        String[] parts = histo.split(" ");
        String[] sensor = parts[1].split("_");
        String sensorName = sensor[1] + "_" + sensor[3];
        int sensorNumber = Integer.parseInt((sensor[1].substring(1, 2)));
        if (sensorNumber > 4) {
            sensorName = sensorName + "_" + sensor[4];
        }
        int nXbins = hist.xAxis().bins();
        if (sensorNumber < 3) {
            nXbins = 512;
        }
//        System.out.println("histo " + histo + " has " + nXbins + " x bins");
        int startBin = 0;
        int stopBin = nXbins;
        for (int j = startBin; j < stopBin; ++j) {
            _successfulSlices++;
            IHistogram1D slice = _hf.sliceY(parts[0] + " " + parts[1] + " slice " + j, hist, j);
            int regionIndex = j % 16;
            int plotterIndex = j / 16;
            boolean analyzeIt = true;
            // don't process if we don't have enough entries in the slice
//            System.out.println("slice: " + j + " entries: " + slice.entries());
            boolean goodFit = analyzeOneSlice(slice, j, plotter, regionIndex);
            //save out bad slices, but not empty ones
            if (!goodFit && slice.entries() != 0) {
                _outputhf.createCopy(slice.title() + " bad fit", slice);
                //let's keep a record of it for our summary
                _failedSlices.add(slice.title());
            }
//            if (slice.entries() < _minSliceEntries) {
            //plot it anyway so each page always has the same layout
//            if (writeEmOut || showem) {
//                plotter.region(regionIndex).plot(slice);
//            }
//            }
            if (regionIndex == 15) {
                if (writeEmOut) {
                    try {
                        String plotFileName = "Run_" + parts[0] + "_" + sensorName + "_" + plotterIndex + ".png";
                        System.out.println("Writing " + plotFileName);
                        plotter.writeToFile(plotFileName);
                    } catch (IOException ex) {
                        System.out.println("Can't open plotter file");
                    }

                }
                if (showem || writeEmOut) {
                    plotter.clearRegions();
                }
            }
        }
    }

    private boolean analyzeOneSlice(IHistogram1D slice, int sliceNumber, IPlotter globalPlotter, int regionIndex) {
        boolean fitGood = false;
        // number of bins in this slice
        int nBins = slice.axis().bins();
        //accumulate data for the gaussian fit...
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = 0; i < nBins - 1; ++i) {
            if (slice.binHeight(i) > 0.) {
                // this is a good bin
//                System.out.println("bin mean " + slice.binMean(i) + " entries " + slice.binEntries(i));
                obs.add(slice.binError(i), slice.binMean(i), slice.binEntries(i));
            }
        }
        List<WeightedObservedPoint> pointsToFit = obs.toList();
        double[] parameters = new double[3];
        if (pointsToFit.size() > 5) {
            // fit...
            // robust linear parameter estimator
            double[] pars = fitFast(obs);
            // check that width of the gaussian is not too wide
            if (!Double.isNaN(pars[1]) && pars[1] < 300.) {
//            System.out.println("fitting " + pointsToFit.size() + " points");
                parameters = GaussianCurveFitter.create().fit(obs.toList());
                fitGood = true;
                if (parameters[2] > 300.) {
                    fitGood = false;
                }
            }

            if (globalPlotter != null) {
                globalPlotter.region(regionIndex).style().statisticsBoxStyle().setVisible(true);
                globalPlotter.region(regionIndex).style().statisticsBoxStyle().setVisibileStatistics("000");

                globalPlotter.region(regionIndex).style().legendBoxStyle().setVisible(false);
                globalPlotter.region(regionIndex).plot(slice, _dataStyle);
                if (fitGood) {
                    // overlay the gaussian fit...
                    IFunction gauss = _functionFactory.createFunctionByName("gauss", "G");
                    gauss.setParameter("amplitude", parameters[0]);
                    gauss.setParameter("mean", parameters[1]);
                    gauss.setParameter("sigma", parameters[2]);
                    globalPlotter.region(regionIndex).plot(gauss, _functionStyle);
                }
            }
        }
        return fitGood;
    }
    //TODO accumulate the sums as I process each slice

    private double[] fitFast(WeightedObservedPoints obs) {
        double sumX = 0.;
        double sumX2 = 0.;
        double sumX3 = 0.;
        double sumX4 = 0.;
        double sumlny = 0.;
        double sumxlny = 0.;
        double sumx2lny = 0.;
        int N = 0;
        for (WeightedObservedPoint p : obs.toList()) {
            double x = p.getX();
            double y = p.getY();
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            double lny = log(y);
            double xlny = x * log(y);
            double x2lny = x * xlny;
            sumX += x;
            sumX2 += x2;
            sumX3 += x3;
            sumX4 += x4;
            sumlny += lny;
            sumxlny += xlny;
            sumx2lny += x2lny;
            N++;
        }
        RealMatrix coefficients
                = new Array2DRowRealMatrix(new double[][]{{N, sumX, sumX2}, {sumX, sumX2, sumX3}, {sumX2, sumX3, sumX4}},
                false);
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        RealVector constants = new ArrayRealVector(new double[]{sumlny, sumxlny, sumx2lny}, false);
        RealVector solution = null;
        try {
            solution = solver.solve(constants);
        } catch (Exception ex) {
        }
        if (solution != null) {
            double mu = -solution.getEntry(1) / (2 * solution.getEntry(2));
            double sigma = sqrt(-1 / (2 * solution.getEntry(2)));
            return new double[]{mu, sigma};
        }
        return new double[]{Double.NaN, Double.NaN};
    }

    private void finish() {
        try {
            _outputTree.commit();
            _outputTree.close();
        } catch (IOException ex) {
            System.out.println("Can't save output histogram file");
        }
        System.out.println("Failed to fit " + _failedSlices.size() + " slices out of " + _successfulSlices + " total");

    }
}
