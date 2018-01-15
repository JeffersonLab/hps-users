package org.hps.users.spaul;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.util.Driver;

public class PrintHitCount extends Driver{
    PrintWriter pw;
    @Override
    public void process(EventHeader header){
        
        int eventNumber = header.getEventNumber();
        int hitCount = header.get(RawTrackerHit.class,"SVTRawTrackerHits").size();
        long timestamp = header.getTimeStamp();
        pw.println(eventNumber + "," + hitCount + "," + timestamp);
    }
    public void setOutputFileName(String file){
        try {
            pw = new PrintWriter(file);
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
}
