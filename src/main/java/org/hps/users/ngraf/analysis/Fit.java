package org.hps.users.ngraf.analysis;

import hep.aida.IAnalysisFactory;
import hep.aida.IAnnotation;
import hep.aida.IFunction;
import hep.aida.IFunctionFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.ITree;
import hep.aida.ITreeFactory;
import hep.aida.ref.Annotation;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.Random;
import static org.apache.commons.math3.special.Erf.erf;

public class Fit {

    public static void main(String[] args) {
        // Create factories
        IAnalysisFactory analysisFactory = IAnalysisFactory.create();
        ITreeFactory treeFactory = analysisFactory.createTreeFactory();
        ITree tree = treeFactory.create();
        IPlotter plotter = analysisFactory.createPlotterFactory().create("Fit.java Plot");
        IFunctionFactory functionFactory = analysisFactory.createFunctionFactory(tree);
        IFunction gauss = functionFactory.createFunctionFromScript("gauss", 1, "background+a*exp(-(x[0]-mean)*(x[0]-mean)/sigma/sigma)", "a,mean,sigma,background", "A Gaussian");
        IHistogramFactory histogramFactory = analysisFactory.createHistogramFactory(tree);
        IHistogram1D h1 = histogramFactory.createHistogram1D("Histogram 1D", 50, -10., 10.);

        Random r = new Random();

        IFunction cb = new CrystalBall();
        //                 a n xb sig
        double[] cbpars = {2., 3., 0., 1.2};
        cb.setParameters(cbpars);
        double[] x = {-3., -2.5, -2. -1., 0., 1., 2., 2.5, 3.};
        double[] xin = new double[1];
        for(int i = 0; i<x.length; ++i)
        {
            xin[0] = x[i];
            System.out.println(x[i]+" "+cb.value(xin));
        }
        
        int nTries = 10000;
        for (int i = 0; i < 100000; i++) {
            h1.fill(r.nextGaussian(), 1./nTries);
        }
  
        plotter.region(0).plot(h1);
        gauss.setParameter("a",h1.maxBinHeight());
        gauss.setParameter("mean",h1.mean());
      gauss.setParameter("sigma",h1.rms());
 
//      IFitter jminuit = fitFactory.createFitter("Chi2","jminuit");
//
//      IFitResult jminuitResult = jminuit.fit(h1,gauss);
//
//      plotter.region(0).plot(jminuitResult.fittedFunction());
      plotter.show();
      plotter.region(0).plot(cb);

    }
}

class CrystalBall implements IFunction {

    double a;
    double n;
    double xb;
    double sig;

    String[] vars = {"a", "n", "xb", "sig"};

    @Override
    public String title() {
        return "Crystal Ball Function";
    }

    @Override
    public void setTitle(String string) throws IllegalArgumentException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double value(double[] doubles) {

        if (a < 0) {
            a = -a;
        }
        if (n < 0) {
            n = -n;
        }
        double aa = abs(a);
        double x = doubles[0];
        double A = pow(n / aa, n) * exp(-aa * aa / 2.);

        double B = n / aa - aa;

        double C = (n / aa)*(1 / (n - 1.)) * exp(-aa * aa / 2.);

        double D = sqrt(PI / 2.) * (1. + erf(aa / sqrt(2.)));

        double N = 1. / (sig * (C + D));

        double total = 0.;

        if ((x - xb) / sig > -a) {
            total = N * exp(-(x - xb) * (x - xb) / (2. * sig * sig));
        } else {
            total = N * A * pow(B - (x - xb) / sig, -n);
        }
        return total;
    }

    @Override
    public int dimension() {
        return 1;
    }

    @Override
    public boolean isEqual(IFunction i) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] gradient(double[] doubles) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean providesGradient() {
        return false;
    }

    @Override
    public String variableName(int i) {
        return vars[i];
    }

    @Override
    public String[] variableNames() {
        return vars;
    }

    @Override
    public void setParameters(double[] doubles) throws IllegalArgumentException {
        a = doubles[0];
        n = doubles[1];
        xb = doubles[2];
        sig = doubles[3];
    }

    @Override
    public double[] parameters() {
        return new double[]{a, n, xb, sig};
    }

    @Override
    public int numberOfParameters() {
        return 4;
    }

    @Override
    public String[] parameterNames() {
        return vars;
    }

    @Override
    public void setParameter(String string, double d) throws IllegalArgumentException {

        switch (string) {
            case "a":
                a = d;
                break;
            case "n":
                n = d;
                break;
            case "xb":
                xb = d;
                break;
            case "sig":
                sig = d;
                break;
            default:
                throw new RuntimeException("variable " + string + " not recognized");
        }
    }

    @Override
    public double parameter(String string) {
        switch (string) {
            case "a":
                return a;
            case "n":
                return n;
            case "xb":
                return xb;
            case "sig":
                return sig;
            default:
                throw new RuntimeException("variable " + string + " not recognized");
        }
    }

    @Override
    public int indexOfParameter(String string) {
        switch (string) {
            case "a":
                return 0;
            case "n":
                return 1;
            case "xb":
                return 2;
            case "sig":
                return 3;
            default:
                throw new RuntimeException("variable " + string + " not recognized");
        }
    }

    @Override
    public IAnnotation annotation() {
        return new Annotation();
    }

    @Override
    public String codeletString() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String normalizationParameter() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
