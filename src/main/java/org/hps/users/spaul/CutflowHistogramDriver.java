package org.hps.users.spaul;

import java.util.ArrayList;
import java.util.List;

import org.hps.recon.particle.ReconParticleDriver;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;


public class CutflowHistogramDriver extends Driver{
    AIDA aida = AIDA.defaultInstance();

    private double goodnessPidThreshold = 6.1;
    private double feeThreshold = 1.76;
    private double pTotThreshold = 2.90;
    private double radThreshold = 1.51;//1.750;
    private double trackChi2Threshold = 70;
    private double trackClusterTimeDiffThresholdMean = 55;
    private double trackClusterTimeDiffThresholdAbs = 4.4;
    private int maxSharedHits = 3;
    public void setGoodnessPidThreshold(double goodnessPidThreshold) {
        this.goodnessPidThreshold = goodnessPidThreshold;
    }

    public void setFeeThreshold(double feeThreshold) {
        this.feeThreshold = feeThreshold;
    }

    public void setTrackChi2Threshold(double trackChi2Threshold) {
        this.trackChi2Threshold = trackChi2Threshold;
    }

    public void setMaxSharedHits(int maxSharedHits) {
        this.maxSharedHits = maxSharedHits;
    }



    public void setPosD0Threshold(double posD0Threshold) {
        this.posD0Threshold = posD0Threshold;
    }

    public void setClusterTimeDiffThreshold(double clusterTimeDiffThreshold) {
        this.clusterTimeDiffThreshold = clusterTimeDiffThreshold;
    }




    public void setRadThreshold(double val){
        radThreshold = val;
    }

    public void setPzTotThreshold(double val){
        pTotThreshold = val;
    }

    public void setMeanClusterTrackDt(double val){
        trackClusterTimeDiffThresholdMean = val;
    }
    public void setAbsClusterTrackDt(double val){
        trackClusterTimeDiffThresholdAbs = val;
    }
    /* boolean corrD0 = false;
    //correct the track D0 by subtracting 5 mm * phi0.  
    public void setTweakD0(boolean val){
        corrD0 = val;
    }*/

    private double posD0Threshold = 1.07;
    private double clusterTimeDiffThreshold = 2;

    int nMassBins = 6000;
    double maxMass = .3;

    int nCDTbins = 100;
    double maxCDT = 10;


    String levels[] = {"__pre", "__nsigma_cut", "__fee_cut", "__p_tot_max_cut", "__chi2_cut", "__tc_dt_cut", 
            "__pos_l1", "__pos_d0_cut", "__cl_dt_cut", "__p_tot_min_cut"
    };

    int nLevels = levels.length;

    IHistogram1D h_nPassingCuts = aida.histogram1D("n_events_passing_cuts", nLevels, 0, nLevels);
    IHistogram1D h_nsigma_ele[] = new IHistogram1D[nLevels];
    IHistogram1D h_nsigma_pos[] = new IHistogram1D[nLevels];
    IHistogram1D h_mass[] = new IHistogram1D[nLevels];
    IHistogram1D h_mass_tweak[] = new IHistogram1D[nLevels];
    IHistogram1D h_cdt[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleP[] = new IHistogram1D[nLevels];
    IHistogram1D h_posP[] = new IHistogram1D[nLevels];
    IHistogram1D h_totP[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleChisq[] = new IHistogram1D[nLevels];
    IHistogram1D h_posChisq[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleD0[] = new IHistogram1D[nLevels];
    IHistogram1D h_posD0[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleClTrkDt[] = new IHistogram1D[nLevels];
    IHistogram1D h_posClTrkDt[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleL1[] = new IHistogram1D[nLevels]; 
    IHistogram1D h_posL1[] = new IHistogram1D[nLevels];
    IHistogram1D h_eleL2[] = new IHistogram1D[nLevels];
    IHistogram1D h_posL2[] = new IHistogram1D[nLevels];
    IHistogram1D h_tarChisq[] = new IHistogram1D[nLevels];

    IHistogram1D h_nV0perEvent[] = new IHistogram1D[nLevels];
    IHistogram1D h_nSharedElectron[] = new IHistogram1D[nLevels];
    IHistogram1D h_nSharedPositron[] = new IHistogram1D[nLevels];

    @Override
    public void detectorChanged(Detector d ){
        setupHistograms();
    }

    void setupHistograms(){



        for(int i = 0; i< nLevels; i++){
            h_mass[i] = aida.histogram1D("mass" + levels[i], 6000, 0, .3);
            h_mass_tweak[i] = aida.histogram1D("mass_tweak" + levels[i], 6000, 0, .3);
            h_nsigma_ele[i] = aida.histogram1D("electron_nsigma" + levels[i], 100, 0, 20);
            h_nsigma_pos[i] = aida.histogram1D("positron_nsigma" + levels[i], 100, 0, 20);
            h_cdt[i] = aida.histogram1D("cluster_dt" + levels[i], 100, -10, 10);
            h_totP[i] = aida.histogram1D("total_p" + levels[i], 100, 0, 4.0);
            h_eleP[i] = aida.histogram1D("electron_p" + levels[i], 100, 0, 3.0);
            h_posP[i] = aida.histogram1D("positron_p" + levels[i], 100, 0, 3.0);
            h_eleChisq[i] = aida.histogram1D("electron_chi2" + levels[i], 100, 0, 100);
            h_posChisq[i] = aida.histogram1D("positron_chi2" + levels[i], 100, 0, 100);
            h_eleD0[i] = aida.histogram1D("electron_d0" + levels[i], 200, -10, 10);
            h_posD0[i] = aida.histogram1D("positron_d0" + levels[i], 200, -10, 10);
            h_eleClTrkDt[i] = aida.histogram1D("electron_cluster_track_dt" + levels[i], 100, -50, 50);
            h_posClTrkDt[i] = aida.histogram1D("positron_cluster_track_dt" + levels[i], 100, -50, 50);
            h_eleL1[i] = aida.histogram1D("electron_l1" + levels[i], 2, 0, 2);
            h_posL1[i] = aida.histogram1D("positron_l1" + levels[i], 2, 0, 2);
            h_eleL2[i] = aida.histogram1D("electron_l2" + levels[i], 2, 0, 2);
            h_posL2[i] = aida.histogram1D("positron_l2" + levels[i], 2, 0, 2);
            h_tarChisq[i] = aida.histogram1D("tar_chi2" + levels[i], 200, 0, 200);
            h_nV0perEvent[i] = aida.histogram1D("nV0" + levels[i], 20, 0, 20);
            h_nSharedPositron[i] = aida.histogram1D("nV0_with_shared_positron" + levels[i], 20, 0, 20);
            h_nSharedElectron[i] = aida.histogram1D("nV0_with_shared_electron" + levels[i], 20, 0, 20);
        }

    }

    /*
    // no cuts except event flags, GBL, clusters in opposite volumes, and tc match chi2 < 10.    
    // if tracks have shared hits, choose the track with the better fit chi2. (may change this criteria)

    IHistogram1D massPreliminary = aida.histogram1D("mass_pre", nMassBins, 0, maxMass); 
    IHistogram1D cdtPreliminary = aida.histogram1D("cluster_dt_pre", nCDTbins, -maxCDT, maxCDT); 
    // remove fees
    IHistogram1D massNoFee =  aida.histogram1D("mass_fee_cut", nMassBins, 0, maxMass); 
    IHistogram1D cdtNoFee = aida.histogram1D("cluster_dt_fee_cut", nCDTbins, -maxCDT, maxCDT); 
    IHistogram1D electronPz = aida.histogram1D("electron_pz", 100, 0, 3.0);
    IHistogram1D positronPz = aida.histogram1D("positron_pz", 100, 0, 3.0);

    //after a total pz cut
    IHistogram1D totPz = aida.histogram1D("total_pz", 100, 0, 4.0);
    IHistogram1D massPzCut =  aida.histogram1D("mass_pz_cut", nMassBins, 0, maxMass);
    IHistogram1D cdtPzCut  = aida.histogram1D("cluster_dt_pz_cut", nCDTbins, -maxCDT, maxCDT);

    //after a chi2 cut
    IHistogram1D eleTrackChi2 = aida.histogram1D("electron_track_chi2", 200, 0, 100);
    IHistogram1D posTrackChi2 = aida.histogram1D("positron_track_chi2", 200, 0, 100);
    IHistogram1D massTrackChi2Cut =  aida.histogram1D("mass_track_chi2_cut", nMassBins, 0, maxMass);
    IHistogram1D cdtTrackChi2Cut  = aida.histogram1D("cluster_dt_track_chi2_cut", nCDTbins, -maxCDT, maxCDT);

    //after a posD0 cut
    IHistogram1D eleD0 = aida.histogram1D("electron_d0", 200, -10, 10);
    IHistogram1D posD0 = aida.histogram1D("positron_d0", 200, -10, 10);
    IHistogram1D massPosD0Cut =  aida.histogram1D("mass_pos_d0_cut", nMassBins, 0, maxMass);
    IHistogram1D cdtPosD0Cut = aida.histogram1D("cluster_dt_pos_d0_cut", nCDTbins, -maxCDT, maxCDT);

    //track cluster time cut
    IHistogram1D clusterTrackDt = aida.histogram1D("cluster_track_dt", 200, 0, 100);
    IHistogram1D massClusterTrackDtCut =  aida.histogram1D("mass_tc_dt_cut", nMassBins, 0, maxMass);
    IHistogram1D cdtClusterTrackDtCut  = aida.histogram1D("cluster_dt_tc_dt_cut", nCDTbins, -maxCDT, maxCDT);

    //L1 cut:
    IHistogram1D eleL1 = aida.histogram1D("electron_l1", 5, 0, 5);
    IHistogram1D posL1 = aida.histogram1D("positron_l1", 5, 0, 5);
    IHistogram1D eleL2 = aida.histogram1D("electron_l2", 5, 0, 5);
    IHistogram1D posL2 = aida.histogram1D("positron_l2", 5, 0, 5);
    IHistogram1D massL1Cut =  aida.histogram1D("mass_l1l2_cut", nMassBins, 0, maxMass);
    IHistogram1D cdtL1Cut  = aida.histogram1D("cluster_dt_l1l2_cut", nCDTbins, -maxCDT, maxCDT);

    //px total
    IHistogram1D totPx = aida.histogram1D("total_px", 200, -.2, .2);

    //pt asymmetry
    IHistogram1D ptAsymmetry = aida.histogram1D("pt_asymmetry", 200, -1, 1);
    //pt asymmetry
    IHistogram1D pzAsymmetry = aida.histogram1D("pz_asymmetry", 200, -1, 1);
    IHistogram1D ptpzAsymmetry = aida.histogram1D("pt_over_pz_asymmetry", 200, -1, 1);


    // all cuts
    IHistogram1D massFinal =  aida.histogram1D("mass_final", nMassBins, 0, maxMass);
    IHistogram1D cdtFinal  = aida.histogram1D("cluster_dt_final", nCDTbins, -maxCDT, maxCDT);
    //after a total pz cut
    IHistogram1D totPzFinal = aida.histogram1D("total_pz_final", 100, 0, 4.0);

    IHistogram2D massVsPzFinal = aida.histogram2D("mass_vs_tot_pz_final", nMassBins, 0, maxMass, 100, 0, 4);

    IHistogram1D massRad = aida.histogram1D("mass_rad", nMassBins, 0, maxMass);
    IHistogram1D cdtRad = aida.histogram1D("cluster_dt_rad", nCDTbins, -maxCDT, maxCDT);*/

    /**
     * check if two tracks share hits.  If so, return the track that is "worse"
     * than the other one, by some criteria (for now, just use track fit chi2)
     * @param t1
     * @param t2
     * @return
     */
    Track getWorseTrack(Track t1, Track t2){
        if(t1 == t2)
            return null;
        if(TrackUtils.numberOfSharedHits(t1, t2) <= maxSharedHits)
            return null;
        //for now, just use track fit chi2.  
        if(t1.getChi2() > t2.getChi2())
            return t1;
        return t2;
    }

    void preliminaryCleanup(List<ReconstructedParticle> v0s){
        //first get rid of v0s made from matched tracks.  (gbl only).  
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            if(v1.getType()<32)
                trash.add(v1);
        }
        v0s.removeAll(trash);
        trash.clear();

        //then, remove v0s where there is a missing cluster, 
        //or both clusters are not on opposite sides of the ECal

        for(ReconstructedParticle v1 : v0s){
            //do they have clusters?  
            if(v1.getParticles().get(0).getClusters().size() == 0 
                    || v1.getParticles().get(1).getClusters().size() == 0 )
                trash.add(v1);
            // on opposite sides of detector?
            else if(v1.getParticles().get(0).getClusters().get(0).getPosition()[1]
                    *v1.getParticles().get(1).getClusters().get(0).getPosition()[1] > 0)
                trash.add(v1);

        }
        v0s.removeAll(trash);
    }

    public void cleanupDuplicates(List<ReconstructedParticle> v0s){

        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            Track e1 = v1.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0);
            Track p1 = v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0);
            for(ReconstructedParticle v2 : v0s){
                Track e2 = v2.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0);
                Track p2 = v2.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0);
                Track worse = getWorseTrack(e1, e2);
                if(worse == e1)
                    trash.add(v1);
                else if(worse == e2)
                    trash.add(v2);
                worse = getWorseTrack(p1, p2);
                if(worse == p1)
                    trash.add(v1);
                else if(worse == p2)
                    trash.add(v2);

            }
        }
        v0s.removeAll(trash);
    }

    public void feeCut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            if(v1.getParticles().get(ReconParticleDriver.ELECTRON).getMomentum().magnitude()>feeThreshold)
                trash.add(v1);
        }
        v0s.removeAll(trash);
    }
    /**
     *  pz max cut.  
     * @param v0s
     */
    public void pMaxCut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            if(v1.getMomentum().magnitude()>pTotThreshold)
                trash.add(v1);
        }
        v0s.removeAll(trash);
    }
    public void d0Cut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            double posD0 = TrackUtils.getDoca(v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0));

            if(posD0>posD0Threshold)
                trash.add(v1);

        }
        v0s.removeAll(trash);
    }

    public void trackClusterMatchCut(List<ReconstructedParticle> v0s){

        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            if(v1.getParticles().get(0).getGoodnessOfPID() > goodnessPidThreshold 
                    || v1.getParticles().get(1).getGoodnessOfPID() > goodnessPidThreshold)
                trash.add(v1);

        }
        v0s.removeAll(trash);
    }

    public void clusterTrackDTCut(List<ReconstructedParticle> v0s, EventHeader event) {
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            double trackEleTime = TrackUtils.getTrackTime(v1.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0),TrackUtils.getHitToStripsTable(event),TrackUtils.getHitToRotatedTable(event));
            double trackPosTime = TrackUtils.getTrackTime(v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0),TrackUtils.getHitToStripsTable(event),TrackUtils.getHitToRotatedTable(event));

            double clusterEleTime =  v1.getParticles().get(ReconParticleDriver.ELECTRON).getClusters().get(0).getCalorimeterHits().get(0).getTime();
            double clusterPosTime =  v1.getParticles().get(ReconParticleDriver.POSITRON).getClusters().get(0).getCalorimeterHits().get(0).getTime();

            if(clusterEleTime - trackEleTime > trackClusterTimeDiffThresholdMean + trackClusterTimeDiffThresholdAbs 
                    || clusterEleTime - trackEleTime < trackClusterTimeDiffThresholdMean - trackClusterTimeDiffThresholdAbs)
                trash.add(v1);
            if(clusterPosTime - trackPosTime > trackClusterTimeDiffThresholdMean + trackClusterTimeDiffThresholdAbs 
                    || clusterPosTime - trackPosTime < trackClusterTimeDiffThresholdMean - trackClusterTimeDiffThresholdAbs)
                trash.add(v1);
        }
        v0s.removeAll(trash);
    }

    public void trackChi2Cut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){
            if(v1.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0).getChi2()>trackChi2Threshold)
                trash.add(v1);
            if(v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0).getChi2()>trackChi2Threshold)
                trash.add(v1);
        }
        v0s.removeAll(trash);
    }

    public void clusterDtCut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){

            if(Math.abs(getClusterTimeDiff(v1))>clusterTimeDiffThreshold)
                trash.add(v1);
        }
        v0s.removeAll(trash);
    }

    public void L1L2Cut(List<ReconstructedParticle> v0s){
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){

            if(hasL1(v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0)) == 0)
                trash.add(v1);
            /*if(hasL1(v1.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0)) == 0)
                trash.add(v1);
            if(hasL2(v1.getParticles().get(ReconParticleDriver.POSITRON).getTracks().get(0)) == 0)
                trash.add(v1);
            if(hasL2(v1.getParticles().get(ReconParticleDriver.ELECTRON).getTracks().get(0)) == 0)
                trash.add(v1);*/
        }
        v0s.removeAll(trash);
    }

    private void radCut(List<ReconstructedParticle> v0s) {
        List<ReconstructedParticle> trash = new ArrayList<ReconstructedParticle>();
        for(ReconstructedParticle v1 : v0s){

            if(v1.getMomentum().magnitude() < radThreshold)
                trash.add(v1);

        }
        v0s.removeAll(trash);
    }


    @Override
    public void process(EventHeader event){

        List<ReconstructedParticle> v0s = event.get(ReconstructedParticle.class, "TargetConstrainedV0Candidates");

        preliminaryCleanup(v0s);
        cleanupDuplicates(v0s);
        fillHistograms(0, v0s, event);

        trackClusterMatchCut(v0s);
        fillHistograms(1, v0s, event);

        feeCut(v0s);
        fillHistograms(2, v0s, event);

        pMaxCut(v0s);
        fillHistograms(3, v0s, event);

        trackChi2Cut(v0s);
        fillHistograms(4, v0s, event);

        clusterTrackDTCut(v0s, event);
        fillHistograms(5, v0s, event);

        //now do L1 cut.  
        L1L2Cut(v0s);
        fillHistograms(6, v0s, event);

        //now look at the positron d0        
        d0Cut(v0s);
        fillHistograms(7, v0s, event);

        //cluster dt.  
        clusterDtCut(v0s);
        fillHistograms(8, v0s, event);

        radCut(v0s);
        fillHistograms(9, v0s, event);


    }



    double getMass(ReconstructedParticle v0, boolean useBeamEnergyConstraint){
        double mass = v0.getMass();
        if(useBeamEnergyConstraint && v0.getMomentum().magnitude()>2.306)
            mass*=  2.306/v0.getMomentum().magnitude();
        return mass;
    }


    void fillHistograms(int level, List<ReconstructedParticle> v0s, EventHeader event){
        for(ReconstructedParticle v0 : v0s){
            ReconstructedParticle electron = v0.getParticles().get(ReconParticleDriver.ELECTRON);
            ReconstructedParticle positron = v0.getParticles().get(ReconParticleDriver.POSITRON);

            Track eleTrack = electron.getTracks().get(0);
            Track posTrack = positron.getTracks().get(0);
            double trackEleTime = TrackUtils.getTrackTime(eleTrack,TrackUtils.getHitToStripsTable(event),TrackUtils.getHitToRotatedTable(event));
            double trackPosTime = TrackUtils.getTrackTime(posTrack,TrackUtils.getHitToStripsTable(event),TrackUtils.getHitToRotatedTable(event));
            double clusterEleTime =  electron.getClusters().get(0).getCalorimeterHits().get(0).getTime();
            double clusterPosTime =  positron.getClusters().get(0).getCalorimeterHits().get(0).getTime();

            h_eleClTrkDt[level].fill(clusterEleTime-trackEleTime-55);
            h_posClTrkDt[level].fill(clusterPosTime-trackPosTime-55);
            h_mass[level].fill(getMass(v0, false));
            h_mass_tweak[level].fill(getMass(v0, true));
            h_cdt[level].fill(getClusterTimeDiff(v0));

            h_nsigma_ele[level].fill(electron.getGoodnessOfPID());
            h_nsigma_pos[level].fill(positron.getGoodnessOfPID());

            h_totP[level].fill(v0.getMomentum().magnitude());
            h_eleP[level].fill(electron.getMomentum().magnitude());
            h_posP[level].fill(positron.getMomentum().magnitude());
            h_eleChisq[level].fill(eleTrack.getChi2());
            h_posChisq[level].fill(posTrack.getChi2());


            h_eleD0[level].fill(TrackUtils.getDoca(eleTrack));
            h_posD0[level].fill(TrackUtils.getDoca(posTrack));

            h_tarChisq[level].fill(v0.getStartVertex().getChi2());

            h_eleL1[level].fill(hasL1(eleTrack));
            h_posL1[level].fill(hasL1(posTrack));
            h_eleL2[level].fill(hasL2(eleTrack));
            h_posL2[level].fill(hasL2(posTrack));
            h_nPassingCuts.fill(level);
        }

        h_nV0perEvent[level].fill(v0s.size());

        int nSharedElectron = getNSharedElectron(v0s);
        h_nSharedElectron[level].fill(nSharedElectron);

        int nSharedPositron = getNSharedPositron(v0s, level == nLevels -1);
        h_nSharedPositron[level].fill(nSharedPositron);
    }

    private int getNSharedElectron(List<ReconstructedParticle> v0s) {
        int nSharedV0 = 0;
        for(ReconstructedParticle v1 : v0s){
            for(ReconstructedParticle v2 : v0s){
                if(v1 == v2)
                    continue;
                Track e1 = v1.getParticles().get(0).getTracks().get(0);
                Track e2 = v2.getParticles().get(0).getTracks().get(0);
                if(e1 == e2){
                    nSharedV0++;
                }
            }
        }

        if(nSharedV0 != 0){

        }
        return nSharedV0;
    }
    private int getNSharedPositron(List<ReconstructedParticle> v0s, boolean debug) {
        int nSharedV0 = 0;
        for(ReconstructedParticle v1 : v0s){
            for(ReconstructedParticle v2 : v0s){
                if(v1 == v2)
                    continue;
                Track e1 = v1.getParticles().get(0).getTracks().get(0);
                Track e2 = v2.getParticles().get(0).getTracks().get(0);
                Track p1 = v1.getParticles().get(1).getTracks().get(0);
                Track p2 = v2.getParticles().get(1).getTracks().get(0);
                if(p2 == p1){
                    nSharedV0++;
                    if(debug && e1.getChi2() > e2.getChi2()){
                        System.out.println("\nShared positron");
                        System.out.println(v1.getParticles().get(1).getMomentum().magnitude());
                        System.out.println("electron 1");
                        System.out.println(v1.getParticles().get(0).getMomentum().magnitude());
                        System.out.println(v1.getMass());
                        System.out.println("electron 2");
                        System.out.println(v2.getParticles().get(0).getMomentum().magnitude());
                        System.out.println(v2.getMass());
                        System.out.println("vtx chi2 1");
                        System.out.println(v1.getStartVertex().getChi2());
                        System.out.println("vtx chi2 2");
                        System.out.println(v2.getStartVertex().getChi2());
                    }
                }
            }
        }


        return nSharedV0;
    }

    /*private double getD0(Track track){
        double d0 =  TrackUtils.getDoca(track); 
        //make correction due to target being at -5 mm 
        if(corrD0) d0 -= TrackUtils.getPhi0(track)*5;
        return d0;

    }*/
    private int hasL1(Track track) {
        for(org.lcsim.event.TrackerHit hit : track.getTrackerHits()){
            if(TrackUtils.getLayer(hit) == 1)
                return 1;
        }
        return 0;
    }

    private int hasL2(Track track) {    
        for(TrackerHit hit : track.getTrackerHits()){
            if(TrackUtils.getLayer(hit) == 3)
                return 1;
        }
        return 0;
    }
    private double getClusterTimeDiff(ReconstructedParticle v0) {

        return v0.getParticles().get(0).getClusters().get(0).getCalorimeterHits().get(0).getTime()
                - v0.getParticles().get(1).getClusters().get(0).getCalorimeterHits().get(0).getTime();
    }
}
