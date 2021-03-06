/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.users.luca;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.ecal.EcalChannel;
import org.hps.conditions.ecal.EcalChannel.EcalChannelCollection;
import org.hps.conditions.ecal.EcalChannelConstants;
import org.hps.conditions.ecal.EcalConditions;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.Subdetector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
/**
 *
 * @author Luca
 */
public class FeeCalibHistCreator extends Driver {
    
    double energyThreshold=0.5;
    protected String clusterCollectionName = "GTPEcalClusters";
    //create the writer to write the gains on txt
    private FileWriter writer;
    String outputFileName = "gains.txt";
    
    //things you need to get the gains
    private EcalConditions ecalConditions = null;
    final private String ecalName = "Ecal";
    private Subdetector ecal;
    private EcalChannelCollection channels= null;
    double[] gain=new double[442];
    
    //Declaration of histograms array
    AIDA aida = AIDA.defaultInstance();
    ArrayList<IHistogram2D> GTPHists = new ArrayList<IHistogram2D>(442);
   
    
    
    public void setEnergyThreshold (double threshold){
    this.energyThreshold=threshold;
       }
  
  public void setOutputFileName(String outputFileName){
  this.outputFileName = outputFileName+".txt";
  }
  
  @Override
    public void detectorChanged(Detector detector) {
        // Get the Subdetector.
        ecal = detector.getSubdetector(ecalName);
        
      /*  // ECAL combined conditions object.
        ecalConditions = ConditionsManager.defaultInstance()
                .getCachedConditions(EcalConditions.class, TableConstants.ECAL_CONDITIONS).getCachedData();*/
                       
        ecalConditions = DatabaseConditionsManager.getInstance().getEcalConditions();
         channels = ecalConditions.getChannelCollection();
    }  

   @Override   
public void startOfData(){
   //initialize histograms  
      for(int t=0; t<442; t++){
      String cristallo=String.valueOf(t+1);  
      String GTPhistName="GTPCluster-Seeds Energy " + cristallo;
      IHistogram2D GTPclushisto=aida.histogram2D(GTPhistName, 150, 0.0,1.5,150,0.0,1.5);
      GTPHists.add(GTPclushisto);
      }
   //initialize writer
      try{
         //initialize the writers
         writer=new FileWriter(outputFileName);
         //Clear the files
         writer.write("");
      }
      catch(IOException e ){
      System.err.println("Error initializing output file for gains.");
      }
    
}



    @Override
    public void endOfData(){
    
        for(int t=0;t<442;t++){
      int id =t+1;
      try{  writer.append(id + " " + gain[t] + "\n");}
      catch(IOException e){
      System.err.println("Error closing output file.");}
      }
      //close the file writer.  
      try{
      writer.close();
      }
      catch(IOException e){
      System.err.println("Error closing utput file.");
      }
} 
    @Override
    public void process (EventHeader event){
        /*
          if (event.hasCollection(GenericObject.class, "TriggerBank")) {
            List<GenericObject> triggerList = event.get(GenericObject.class, "TriggerBank");
            for (GenericObject data : triggerList)
                if (AbstractIntData.getTag(data) == TIData.BANK_TAG) {
                    TIData triggerData = new TIData(data);
                    if (!triggerData.isSingle1Trigger())//only process singles0 triggers...

                        return;
                }
        } else //if (debug)
            System.out.println(this.getClass().getSimpleName() + ":  No trigger bank found...running over all trigger types");
      */
        
        //here it writes the GTP clusters info
        if(event.hasCollection(Cluster.class,"EcalClustersGTP"))
        {List<Cluster> clusters= event.get(Cluster.class,"EcalClustersGTP");
         for(Cluster cluster : clusters){
           int idBack;
           
           
           idBack=getDBID(cluster);
           EcalChannel channel = channels.findGeometric(cluster.getCalorimeterHits().get(0).getCellID());
           EcalChannelConstants channelConstants = ecalConditions.getChannelConstants(channel);
           gain[idBack-1]=channelConstants.getGain().getGain();
           
           //riempio gli istogrammi
           if(cluster.getEnergy()>energyThreshold){
           GTPHists.get(idBack -1).fill(cluster.getEnergy(),cluster.getCalorimeterHits().get(0).getCorrectedEnergy());
          
           //System.out.println("Clu = " + cluster.getEnergy() + " Seed = " + cluster.getCalorimeterHits().get(0).getCorrectedEnergy() + "\n");
           }         
         }
        }
       
       
    }
    
    
  public int getDBID ( Cluster cluster ){
    int xx=  cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
    int yy=cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
    int XOFFSET=23;
    int YOFFSET=5;
    int ix = xx<0 ? xx+XOFFSET : xx+XOFFSET-1;
    int iy = yy<0 ? yy+YOFFSET : yy+YOFFSET-1;
    int dbid = ix + 2*XOFFSET*(YOFFSET*2-iy-1) + 1;
    if      (yy ==  1 && xx>-10){ dbid-=9;}
    else if (yy == -1 && xx<-10) {dbid-=9;}
    else if (yy < 0){dbid-=18;}
   return dbid;
}
  
public int getDBID ( CalorimeterHit hit ){
    int xx=  hit.getIdentifierFieldValue("ix");
    int yy=  hit.getIdentifierFieldValue("iy");
    int XOFFSET=23;
    int YOFFSET=5;
    int ix = xx<0 ? xx+XOFFSET : xx+XOFFSET-1;
    int iy = yy<0 ? yy+YOFFSET : yy+YOFFSET-1;
    int dbid = ix + 2*XOFFSET*(YOFFSET*2-iy-1) + 1;
    if      (yy ==  1 && xx>-10){ dbid-=9;}
    else if (yy == -1 && xx<-10) {dbid-=9;}
    else if (yy < 0){dbid-=18;}
   return dbid;
}  
    
}