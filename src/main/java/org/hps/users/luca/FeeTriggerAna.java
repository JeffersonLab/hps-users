/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.hps.users.luca;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.hps.readout.ecal.TriggerDriver;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;


/**
 * 
 * @author Luca Colaneri 
 */
public class FeeTriggerAna extends Driver {
    int posx, posy;
    int radius=2;
    int Clustercount=0;
    int clusterWindow=50;
    int TotalCluster=0;
    double timeDifference;
    double energyThreshold=0;
    private LinkedList<ArrayList<Cluster>> clusterBuffer;
    protected String clusterCollectionName = "EcalClusters";
    
 //AIDA aida = AIDA.defaultInstance();
// IHistogram1D clusterEne=aida.histogram1D("Clusters energy with Luca's trigger",300, 0, 3);
// ArrayList<IHistogram1D> SeedHistograms = new ArrayList<IHistogram1D>(442);
  //  ArrayList<IHistogram1D> ClustHistograms = new ArrayList<IHistogram1D>(442);
 //    ArrayList<IHistogram1D> HitHistograms = new ArrayList<IHistogram1D>(442);
    private FileWriter writer;
  //  private FileWriter writer2;
    String outputFileName = "LucaTriggerFEE.txt";
 //   String outputFileName2 = "LucaTriggerHits.txt";

   
 
   
   public void setRadius (int radius){
         this.radius=radius;
                 }
    
    public void setEnergyThreshold (double threshold){
    this.energyThreshold=threshold;
    }
   
    public void setClusterCollectionName(String clusterCollectionName) {
        this.clusterCollectionName = clusterCollectionName;
    }
 
   public void setOutputFileName(String outputFileName){
this.outputFileName = outputFileName;
}
///   public void setOutputFileName2(String outputFileName2){
//this.outputFileName2 = outputFileName2;
   //}
   public void settimeDifference(double time){
   this.timeDifference=time;
   
   }
  /*
   *
   *
   *
   */
   
   @Override   
public void startOfData(){

    //initialize the clusterbuffer
    clusterBuffer= new LinkedList<ArrayList<Cluster>>();
 //populate the clusterbuffer with (2*clusterWindow + 1)
 // empty events, representing the fact that the first few events will not have any events in the past portion of the buffer
    int bufferSize=(2*clusterWindow)+1;
    for(int i = 0;i<bufferSize; i++){
    clusterBuffer.add(new ArrayList<Cluster>(0));
    }
    
    
    
    
    try{
    //initialize the writers
    writer=new FileWriter(outputFileName);
   // writer2=new FileWriter(outputFileName2);
    //Clear the files
    writer.write("");
   // writer2.write("");
    
     //initialize histograms  
  /*  for(int t=0; t<442; t++){
      String cristallo=String.valueOf(t);  
      String seedhistogram="SeedHit_" + String.valueOf(t);
      String Clushistogram="Clusters_" + String.valueOf(t);
      String HitHistogram="Hits_" + String.valueOf(t);
      
      IHistogram1D seedhisto=aida.histogram1D(seedhistogram, 150, 0.0,3.0);
      IHistogram1D clushisto=aida.histogram1D(Clushistogram, 150, 0.0,3.0);
      IHistogram1D hitshisto=aida.histogram1D(HitHistogram,150,0.0,3.0);
    SeedHistograms.add(seedhisto);
    ClustHistograms.add(clushisto);
    HitHistograms.add(hitshisto);
    }*/
    
}
catch(IOException e ){
System.err.println("Error initializing output file for event display.");
}
}
  
@Override
public void endOfData(){
System.out.println("Ho contato" + TotalCluster + " clusters di cui " + Clustercount + "isolati\n");
    
    try{
//close the file writer.
    writer.close();
    //writer2.close();
    }
catch(IOException e){
    System.err.println("Error closing utput file for event display.");
}
} 
   
 @Override  
 public void process (EventHeader event){
   
          
           
           
     
     //get the clusters from the event
    if(TriggerDriver.triggerBit()){ //if they have triggered!
      
     if(event.hasCollection(Cluster.class, "EcalClusters")) {
         
        List<Cluster> clusterList =event.get(Cluster.class,clusterCollectionName );    
             
     //put the clusters in the arraylist
     
     ArrayList<Cluster> clusterSet = new ArrayList<Cluster>(); 
     for(Cluster cluster : clusterList){
      //   clusterEne.fill(cluster.getEnergy());
         TotalCluster++;
        // clusterSet.add(cluster);
     int id;
     id=getCrystal(cluster);
     try{
      writer.append(id + " " + cluster.getEnergy()+ " " + cluster.getSize() + " " + cluster.getCalorimeterHits().get(0).getRawEnergy() + " " + cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix")+" " +cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"));
     /*for(CalorimeterHit hit : cluster.getCalorimeterHits())
     {writer.append(hit.getRawEnergy()+ " ");
       }*/
     writer.append("\n");
     
     }
     catch(IOException e ){System.err.println("Error writing to output for event display");}  
     
     }
     //remove the last event from cluster buffer and add the new one
     //clusterBuffer.removeLast();
    // clusterBuffer.addFirst(clusterSet);
    //Run the sorting algorithm;
    // ClusterAnalyzer();
     
     }
     
      
     
    }// questa parentesi va scommentata se si scommenta l'if del trigger
     
}

 
 /**
  * For each crystal, looks for clusters that hit that clystar, if it is an isolated cluster, it's put in goodclusterqueue
  */
 public void ClusterAnalyzer(){
 //get the cluster list at the current time in the buffer
ArrayList<Cluster> currentClusters = clusterBuffer.get(clusterWindow+1);


 ///cerca i cluster nella posizione che ci interessa poi chiama la funzione che decide se sono "isolati"
   //System.out.println("Sta partendo il for sulla Queue \n");
 for(int y=-5;y<6;y++){
     for(int x=-23;x<24;x++){
      posx=x;
      posy=y;
         
         //ciclo for nel set di currentCluster, ovvero il set nel mezzo del buffer
    for(Cluster cluster : currentClusters){ 
    if((cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix")== posx) && (cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy")==posy )&& (cluster.getEnergy() > energyThreshold)){
        
           if(ClusterChecker(cluster)){
            int id;
            Clustercount++;
           id=getCrystal(cluster);
           try{
     writer.append(id + " " + cluster.getEnergy()+ " " + cluster.getSize() + " " + cluster.getCalorimeterHits().get(0).getRawEnergy() + " " + cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix")+" " +cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"));
     /*for(CalorimeterHit hit : cluster.getCalorimeterHits())
     {writer.append(hit.getRawEnergy()+ " ");
       }*/
     writer.append("\n");
  //  SeedHistograms.get(id-1).fill(cluster.getCalorimeterHits().get(0).getRawEnergy());
  //   ClustHistograms.get(id-1).fill(cluster.getEnergy());
     }
     
   catch(IOException e ){System.err.println("Error writing to output for event display");}   
           
           }
      }
     
     
    }
 
 
 
 }
 }
 
 
 
 

 }
 /**
  * Check if the cluster is isolaterd checking if there are clusters near it in time and in space in the buffer
  * @param cluster
  * @return 
  */
 
public boolean ClusterChecker (Cluster cluster){
//System.out.println("Sono nel clustercheck! \n");
    
boolean check=true;
  
    //ciclo sulle liste del buffer
loops:
     for(ArrayList<Cluster> currentList : clusterBuffer){
     //ciclo sui cluster della lista corrente
         for(Cluster currentcluster : currentList){
           if(currentcluster!= cluster){
             //if there is a cluster in the buffer that is in the considered radius in a time window lower than expected, the loop is brocken and the analyzed cluster is not good
         if(!((currentcluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix") < posx-radius || currentcluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix")> posx+radius)&& (currentcluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy")< posy-radius || currentcluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy")> posy+radius))&& Math.abs(cluster.getCalorimeterHits().get(0).getTime()-currentcluster.getCalorimeterHits().get(0).getTime())<timeDifference){
         check=false;
         break loops;
         }
           }
           
        
         
         }
      
     
     }
        
        
   
return check;

}
      
 
 
 public int getCrystal (Cluster cluster){
 int x,y,id=0;
 x= (-1)*cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
 y= cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
 
 if(y==5){
 if(x<0)
 {id=x+24;}
 else id= x+23;
 }
 
 else if(y==4)
 {if(x<0){
  id=x+70;}
 else id=x+69;}
 
 else if(y==3)
 {if(x<0){
  id=x+116;}
 else id=x+115;}
 
 else if(y==2)
 {if(x<0){
  id=x+162;}
 else id=x+161;}
 
 else if(y==1)
 {x=-x;
     if(x>0){
  id=-x+208;}
 else if(x==-1){id=208;}
 else if(x<-1) id=-x+198;}
 
  else if(y==-1)
 {x=-x;
     if(x>0){
  id=-x+245;}
 else if(x==-1 )id=245;
 else if(x<-1){id=-x+235;}}
 
 
 else if(y==-2)
 {if(x<0){
  id=x+282;}
 else id=x+281;}
 
  else if(y==-3)
 {if(x<0){
  id=x+328;}
 else id=x+327;}
 
 else if(y==-4)
 {if(x<0){
  id=x+374;}
 else id=x+373;}
 
 else if(y==-5)
 {if(x<0){
  id=x+420;}
 else id=x+419;}
 
 return id;
 
 }
 
 public int getCrystal (CalorimeterHit hit){
 int x,y,id=0;
 x= (-1)*hit.getIdentifierFieldValue("ix");
 y= hit.getIdentifierFieldValue("iy");
 
 if(y==5){
 if(x<0)
 {id=x+24;}
 else id= x+23;
 }
 
 else if(y==4)
 {if(x<0){
  id=x+70;}
 else id=x+69;}
 
 else if(y==3)
 {if(x<0){
  id=x+116;}
 else id=x+115;}
 
 else if(y==2)
 {if(x<0){
  id=x+162;}
 else id=x+161;}
 
 else if(y==1)
 {x=-x;
     if(x>0){
  id=-x+208;}
 else if(x==-1){id=208;}
 else if(x<-1) id=-x+198;}
 
  else if(y==-1)
 {x=-x;
     if(x>0){
  id=-x+245;}
 else if(x==-1 )id=245;
 else if(x<-1){id=-x+235;}}
 
 
 else if(y==-2)
 {if(x<0){
  id=x+282;}
 else id=x+281;}
 
  else if(y==-3)
 {if(x<0){
  id=x+328;}
 else id=x+327;}
 
 else if(y==-4)
 {if(x<0){
  id=x+374;}
 else id=x+373;}
 
 else if(y==-5)
 {if(x<0){
  id=x+420;}
 else id=x+419;}
 
 return id;
 
 }
 
 } //chiusura driver  
    
    
    
