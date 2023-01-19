package org.hps.users.ngraf.plotutils;

import com.itextpdf.io.font.constants.StandardFonts;
import com.itextpdf.kernel.colors.ColorConstants;
import com.itextpdf.kernel.font.PdfFontFactory;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfReader;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.kernel.utils.PdfMerger;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.extgstate.PdfExtGState;

import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.properties.TextAlignment;
import com.itextpdf.layout.properties.VerticalAlignment;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

/**
 *
 * @author Norman A. Graf
 */
public class CombinePdfs {

    public static final String SRC1 = "tst.pdf";
    public static final String SRC2 = "tst2.pdf";
    public static final String DEST = "combined_documents.pdf";

    public static void main(String args[]) throws IOException {
        File file = new File(DEST);
        new CombinePdfs().createAndAddToPdf("helloWorld.pdf");
        new CombinePdfs().createPdf(DEST);

        new CombinePdfs().manipulatePdf(DEST, "combined_annotated.pdf");
    }

    public void createAndAddToPdf(String dest) throws IOException {
        PdfWriter writer = new PdfWriter(dest);
        PdfDocument pdf = new PdfDocument(writer);
        System.out.println("height: " + pdf.getDefaultPageSize().getHeight() + " width: " + pdf.getDefaultPageSize().getWidth());

        Document document = new Document(pdf);
        document.add(new Paragraph("Hello World!"));
        document.showTextAligned(new Paragraph(String.format("page %s of %s", 1, 1)),
                200, 200, 1, TextAlignment.RIGHT, VerticalAlignment.TOP, 0);
        document.close();
        //Initialize PDF document with output intent
//        PdfDocument pdf = new PdfDocument(new PdfWriter(dest));
//        PdfMerger merger = new PdfMerger(pdf);
//
//        //Add pages from the first document
//        PdfDocument firstSourcePdf = new PdfDocument(new PdfReader(SRC1));
//        merger.merge(firstSourcePdf, 1, firstSourcePdf.getNumberOfPages());
//
//        //Add pages from the second pdf document
//        PdfDocument secondSourcePdf = new PdfDocument(new PdfReader(SRC2));
//        merger.merge(secondSourcePdf, 1, secondSourcePdf.getNumberOfPages());
//
//        firstSourcePdf.close();
//        secondSourcePdf.close();
//        pdf.close();
    }

    public void createPdf(String dest) throws IOException {
        //Initialize PDF document with output intent
        PdfDocument pdf = new PdfDocument(new PdfWriter(dest));
        PdfMerger merger = new PdfMerger(pdf);

        //Add pages from the first document
        PdfDocument firstSourcePdf = new PdfDocument(new PdfReader(SRC1));
        merger.merge(firstSourcePdf, 1, firstSourcePdf.getNumberOfPages());

        //Add pages from the second pdf document
        PdfDocument secondSourcePdf = new PdfDocument(new PdfReader(SRC2));
        merger.merge(secondSourcePdf, 1, secondSourcePdf.getNumberOfPages());

        firstSourcePdf.close();
        secondSourcePdf.close();
        pdf.close();
    }

    public void manipulatePdf(String src, String dest) throws IOException {

        //Initialize PDF document
        PdfDocument pdfDoc = new PdfDocument(new PdfReader(src), new PdfWriter(dest));

        Document document = new Document(pdfDoc);
        Rectangle pageSize;
        PdfCanvas canvas;
        int n = pdfDoc.getNumberOfPages();
        for (int i = 1; i <= n; i++) {
            PdfPage page = pdfDoc.getPage(i);
            pageSize = page.getPageSize();
            canvas = new PdfCanvas(page);
            //Draw header text
            canvas.beginText().setFontAndSize(PdfFontFactory.createFont(StandardFonts.HELVETICA), 7)
                    .moveText(pageSize.getWidth() / 2 - 24, pageSize.getHeight() - 10)
                    .showText("Histogram Comparisons "+myDate())
                    .endText();
            //Draw footer line
            canvas.setStrokeColor(ColorConstants.BLACK)
                    .setLineWidth(.2f)
                    .moveTo(pageSize.getWidth() / 2 - 30, 20)
                    .lineTo(pageSize.getWidth() / 2 + 30, 20).stroke();
            //Draw page number
            canvas.beginText().setFontAndSize(PdfFontFactory.createFont(StandardFonts.HELVETICA), 7)
                    .moveText(pageSize.getWidth() / 2 - 7, 10)
                    .showText(String.valueOf(i))
                    .showText(" of ")
                    .showText(String.valueOf(n))
                    .endText();
            //Draw watermark
            Paragraph p = new Paragraph("CONFIDENTIAL").setFontSize(60);
            canvas.saveState();
            PdfExtGState gs1 = new PdfExtGState().setFillOpacity(0.2f);
            canvas.setExtGState(gs1);
            document.showTextAligned(p,
                    pageSize.getWidth() / 2, pageSize.getHeight() / 2,
                    pdfDoc.getPageNumber(page),
                    TextAlignment.CENTER, VerticalAlignment.MIDDLE, 45);
            canvas.restoreState();
        }
        pdfDoc.close();
    }

    static private String myDate() {
        Calendar cal = new GregorianCalendar();
        Date date = new Date();
        cal.setTime(date);
        DecimalFormat formatter = new DecimalFormat("00");
        String day = formatter.format(cal.get(Calendar.DAY_OF_MONTH));
        String month = formatter.format(cal.get(Calendar.MONTH) + 1);
        return cal.get(Calendar.YEAR) + month + day;
    }
}
