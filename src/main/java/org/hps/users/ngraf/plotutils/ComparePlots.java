package org.hps.users.ngraf.plotutils;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import org.lcsim.util.aida.AIDA;

import com.itextpdf.text.BadElementException;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Element;
import com.itextpdf.text.PageSize;
import com.itextpdf.text.pdf.PdfWriter;

import hep.aida.IAnalysisFactory;
import hep.aida.IBaseHistogram;
import hep.aida.IDataStyle;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.ITree;

/**
 *
 * @author Norman A. Graf
 *
 * To run:
 * java -cp target/hps-users-5.1-SNAPSHOT-bin.jar org.hps.users.ngraf.plotutils.ComparePlots reference.aida test.aida
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

        /**
         * Write generated plots to PDF file --JM
         */
        if (writeToFile) {
            String pdfName = "plots.pdf";
            System.out.println("Writing plots to: " + pdfName);
            writePdf(plotters, pdfName, fileType);
        }
    }

    /**
     * Write plotter images to a PDF, one per page
     *
     * @param plotters The list of AIDA plotters
     * @param fileName The name of the output PDF file
     * @param imgType The type of images to generate (e.g. "png")
     * @throws IOException If there is an IO error
     */
    static void writePdf(List<IPlotter> plotters, String fileName, String imgType) throws IOException {

        // Generate a list of buffered images from the plotters
        List<BufferedImage> images = new ArrayList<BufferedImage>();
        String tempFileName = "temp." + imgType;
        File tempFile = new File(tempFileName);
        tempFile.delete();
        for (IPlotter p : plotters) {
            p.writeToFile(tempFileName, imgType);
            BufferedImage img = ImageIO.read(tempFile);
            images.add(img);
            tempFile.delete();
        }

        // Open the PDF document and the writer
        Document document = new Document(PageSize.LETTER.rotate(), 50, 50, 50, 50);
        PdfWriter writer = null;
        try {
            writer = PdfWriter.getInstance(document, new FileOutputStream(fileName));
        } catch (DocumentException e) {
            throw new IOException(e);
        }
        document.open();

        // Write the generated images to a PDF, one per page
        for (BufferedImage img : images) {
            document.newPage();

            // Write image into the document.
            com.itextpdf.text.Image iTextImage = null;
            try {
                iTextImage = com.itextpdf.text.Image.getInstance(writer, img, 1f);
            } catch (BadElementException e) {
                throw new IOException(e);
            }
            iTextImage.scaleAbsolute(document.getPageSize().getWidth(), (float) 0.75 * document.getPageSize().getHeight());
            iTextImage.setAlignment(Element.ALIGN_CENTER);
            try {
                document.add(iTextImage);
            } catch (DocumentException e) {
                throw new IOException(e);
            }
        }

        document.close();
    }

}
