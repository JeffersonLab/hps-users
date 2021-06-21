package org.hps.users.ngraf.mc;

import hep.io.stdhep.StdhepBeginRun;
import hep.io.stdhep.StdhepEndRun;
import hep.io.stdhep.StdhepExtendedEvent;
import hep.io.stdhep.StdhepWriter;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 *
 * @author Norman A. Graf
 */

/**
 * Reads in text file in DIPSI format and creates events in stdhep file format.
 */
/* File Format: 
 /* 0 event weight 
 /* 1 rho px 
 /* 2 rho py 
 /* 3 rho pz 
 /* 4 rho energy
 /* 5 rho mass
 /* 6 pi+ px 
 /* 7 pi+ py 
 /* 8 pi+ pz 
 /* 9 pi+ energy
 /* 10 pi+ mass
 /* 11 pi- px 
 /* 12 pi- py 
 /* 13 pi- pz 
 /* 14 pi- energy
 /* 15 pi- mass
 */
public class DipsiToStdhepConverter {

    public static void main(String args[]) throws IOException {

        // remind user of correct usage
        if (args.length < 1) {
            usage();
        }
        if (args.length == 1 && args[0].equals("-h")) {
            usage();
        }
        String fileName = args[0];

        FileInputStream fin1 = new FileInputStream(fileName);

        double zOffset = 0.;
        if (args.length > 1) {
            zOffset = Double.parseDouble(args[1]);
            System.out.println("Offsetting z vertex by " + zOffset + " mm.");
        }
        File outputDir = new File(".");
        if (args.length > 2) {
            outputDir = new File(args[2]);
        }
        // check if output directory exists

        if (!outputDir.exists()) {
            System.out.println("\n\n  Directory " + outputDir + " does not exist!");
            System.exit(1);
        }

        int dot = fileName.lastIndexOf(".txt");
        String stdhepFileName = fileName.substring(0, dot);

        stdhepFileName += ".stdhep";

        String outputFile = outputDir + "/" + stdhepFileName;
        System.out.println(outputFile);
        StdhepWriter w = null;
        try {
            w = new StdhepWriter(outputFile, "Imported Stdhep Events v1.0", "From file " + fileName, 10);
            w.setCompatibilityMode(false);
        } catch (java.io.IOException ex) {
            System.err.println("Error opening file: " + outputFile);
            ex.printStackTrace();
            System.exit(1);
        }

        /*--------------------------------------------------------*/
 /* NEVHEP          - event number (or some special meaning*/
 /*                    (see documentation for details)     */
 /* NHEP            - actual number of entries in current  */
 /*                    event.                              */
 /* ISTHEP[IHEP]    - status code for IHEP'th entry - see  */
 /*                    documentation for details           */
 /* IDHEP [IHEP]    - IHEP'th particle identifier according*/
 /*                    to PDG.                             */
 /* JMOHEP[IHEP][0] - pointer to position of 1st mother    */
 /* JMOHEP[IHEP][1] - pointer to position of 2nd mother    */
 /* JDAHEP[IHEP][0] - pointer to position of 1st daughter  */
 /* JDAHEP[IHEP][1] - pointer to position of 2nd daughter  */
 /* PHEP  [IHEP][0] - X momentum [Gev/c]                   */
 /* PHEP  [IHEP][1] - Y momentum [Gev/c]                   */
 /* PHEP  [IHEP][2] - Z momentum [Gev/c]                   */
 /* PHEP  [IHEP][3] - Energy [Gev]                         */
 /* PHEP  [IHEP][4] - Mass[Gev/c^2]                        */
 /* VHEP  [IHEP][0] - X vertex [mm]                        */
 /* VHEP  [IHEP][1] - Y vertex [mm]                        */
 /* VHEP  [IHEP][2] - Z vertex [mm]                        */
 /* VHEP  [IHEP][3] - production time [mm/c]               */
 /*========================================================*/
 /* IstHep convention: */
 /* 0      - final state particle if JdaHEP=0 */
 /*          intermediate particle if JdaHEP>0 */
 /*          (NEUGEN extension; was "null") */
 /* 1      - final state particle */
 /* 2      - intermediate state */
 /* 3      - documentation line */
 /* 4-10   - reserved for future */
 /* 11-200 - reserved for specific model use */
 /* 201+   - reserved for users */
        int _nmax = 2;
        int _eventnum = 0;
        int _nhep = 0;
        int[] _isthep = new int[_nmax];
        int[] _idhep = new int[_nmax];
        int[] _jmohep = new int[2 * _nmax];
        int[] _jdahep = new int[2 * _nmax];
        double[] _phep = new double[5 * _nmax];
        double[] _vhep = new double[4 * _nmax];

        //Dummy values...
        int nevtreq = 1;
        int nevtgen = 1;
        int nevtwrt = 1;
        float stdecom = 2.F;
        float stdxsec = 2.F;
        double stdseed1 = 137.;
        double stdseed2 = 138.;

        //hepev4
        double[] scale = new double[5];
        double[] spin = new double[_nmax * 3];
        int[] colorflow = new int[_nmax * 2];
        int idrupl = 0;

        double[] mom = new double[3];
        double[] pos = new double[3];

        String thisLine;

        // write a begin run record
        w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt, stdecom, stdxsec, stdseed1, stdseed2));
        // now loop over contents of this file...
        FileInputStream fin2 = new FileInputStream(fileName);
        try {
            BufferedReader myInput = new BufferedReader(new InputStreamReader(fin2));
            double[] values = new double[16];
            int nEvent = 0;

            while ((thisLine = myInput.readLine()) != null) {
                // tokenize the string and convert to double values
                StringTokenizer st = new java.util.StringTokenizer(thisLine, " ");
                _nhep = 0;
                int j = 0;
                int numTokens = st.countTokens();
                while (st.hasMoreElements()) {
                    values[j++] = Double.parseDouble(st.nextToken());
                }
//                System.out.println(thisLine);
                double weight = values[0];
                // now populate the HEPEVT "common block"
                // only add the pi+ and pi- to the event...
                //pi+
                _isthep[_nhep] = 1; // final state particle
                _idhep[_nhep] = 211;
                _phep[0 + 5 * _nhep] = values[6]; //px
                _phep[1 + 5 * _nhep] = values[7]; //py
                _phep[2 + 5 * _nhep] = -values[8]; //pz
                _phep[3 + 5 * _nhep] = values[9]; //E
                _phep[4 + 5 * _nhep] = values[10]; // mass
                _vhep[0 + 4 * _nhep] = 0.; //x
                _vhep[1 + 4 * _nhep] = 0.; //y
                _vhep[2 + 4 * _nhep] = zOffset; //z
                // increment the number of particles in this event
                _nhep++;
                //pi-
                _isthep[_nhep] = 1; // final state particle
                _idhep[_nhep] = -211;
                _phep[0 + 5 * _nhep] = values[11]; //px
                _phep[1 + 5 * _nhep] = values[12]; //py
                _phep[2 + 5 * _nhep] = -values[13]; //pz
                _phep[3 + 5 * _nhep] = values[14]; //E
                _phep[4 + 5 * _nhep] = values[15]; // mass
                _vhep[0 + 4 * _nhep] = 0.; //x
                _vhep[1 + 4 * _nhep] = 0.; //y
                _vhep[2 + 4 * _nhep] = zOffset; //z
                // increment the number of particles in this event
                _nhep++;

                //
                // Create an StdhepExtendedEvent top include event weight and write it out...
                // 20210619 seems not to work?
                //
                StdhepExtendedEvent event = new StdhepExtendedEvent(_eventnum++, _nhep, _isthep, _idhep, _jmohep, _jdahep, _phep, _vhep, weight, 0., 0., scale, spin, colorflow, idrupl);
//                System.out.println(event);
//                System.out.println(event.getEventWeight());
                w.writeRecord(event);
                nEvent++;
            }
            // write an end run record
            w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom, stdxsec, stdseed1, stdseed2));
            // close the file
            try {
                System.out.println(" Closing file: " + outputFile + " with " + nEvent + " events");
                w.close();
            } catch (java.io.IOException ex) {
                System.err.println("Error closing file: " + outputFile);
                ex.printStackTrace();
                System.exit(1);
            }
            fin2.close();
        } catch (EOFException ex) {
            ex.printStackTrace();
        }
    }

    public static void usage() {
        System.out.println("DipsiToStdhepConvertor: \n  an application to read in text files in DIPSI format and convert to stdhep format.\n");
        System.out.println("Usage: \n\n >> java DipsiToStdhepConvertor InputFile.txt <zOffSet[mm]> <output directory> \n");
        System.out.println(" Where: \n InputFile.txt    is an input text file in DIPSI format to process ");
        System.out.println(" Writes to the current working directory unless output directory is specified");
        System.out.println("\n e.g. >> java DipsiToStdhepConvertor DipsiEvents.txt -5 \n");
        System.out.println("  will convert the events in DipsiEvents.txt to DipsiEvents.stdhep offset by -5mm in z in the same directory");
        System.exit(0);
    }
}

