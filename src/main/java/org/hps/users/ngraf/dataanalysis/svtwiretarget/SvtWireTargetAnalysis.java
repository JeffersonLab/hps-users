package org.hps.users.ngraf.dataanalysis.svtwiretarget;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import static java.lang.Math.abs;
import java.util.ArrayList;
import java.util.List;
import org.hps.recon.vertexing.BilliorTrack;
import org.hps.recon.vertexing.BilliorVertex;
import org.hps.recon.vertexing.BilliorVertexer;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author Norman A. Graf
 */
public class SvtWireTargetAnalysis extends Driver {

    private AIDA aida = AIDA.defaultInstance();
//    private List<Track> accumulatedTracks = new ArrayList<Track>();
//    private List<BilliorTrack> accumulatedBTracks = new ArrayList<BilliorTrack>();

//    private List<Track> accumulatedTracksTop = new ArrayList<Track>();
    private List<BilliorTrack> accumulatedBTracksTop = new ArrayList<BilliorTrack>();

//    private List<Track> accumulatedTracksBot = new ArrayList<Track>();
    private List<BilliorTrack> accumulatedBTracksBot = new ArrayList<BilliorTrack>();

    private String reconParticleCollectionName = "FinalStateParticles_KF";
    private String trackCollectionName = "KalmanFullTracks";
    protected double[] beamSize = {0.001, 0.130, 0.050}; // rough estimate from harp scans during engineering run
    private double[] beamPositionToUse = new double[3];
    private boolean _storeCovTrkMomList = false;
    private boolean debug = false;
    protected double bField;
    private String vtxFolder = "MultiEventVtx/";
    private int _ntrks = 50;
    private FieldMap fieldMap;
    BilliorVertexer vtxFitter;

    public void setNtrks(int ntrks) {
        _ntrks = ntrks;
    }

    @Override
    protected void detectorChanged(Detector detector) {

        // Set the magnetic field parameters to the appropriate values.
        // we want the field near the target...
        Hep3Vector ip = new BasicHep3Vector(0., 0., 0.0);
        fieldMap = detector.getFieldMap();
        bField = fieldMap.getField(ip).y();
        for (int i = -10; i < 490; ++i) {
            double by = fieldMap.getField(new BasicHep3Vector(0., 0., i)).y();
//            System.out.println(i + " : " + by);
            aida.histogram1D("Fieldmap By at z", 500, -10., 490.).fill(i, -by);
        }

        vtxFitter = new BilliorVertexer(bField);

        List<String> volumes = new ArrayList<String>();
        volumes.add("_top");
        volumes.add("_bottom");
//        volumes.add("");

        for (String vol : volumes) {
            aida.histogram1D(vtxFolder + "vtx_x" + vol, 100, -3, 3);
            aida.histogram1D(vtxFolder + "vtx_y" + vol, 100, -1, 1);
            aida.histogram1D(vtxFolder + "vtx_z" + vol, 100, 0., 50.);
            aida.histogram2D(vtxFolder + "vtx_x_y" + vol, 200, -3, 3, 200, -1, 1);
            aida.histogram1D(vtxFolder + "vtx_chi2" + vol, 200, 0, 8000);
            aida.histogram1D(vtxFolder + "vtx_ntrks" + vol, 200, 0, 200);

            aida.histogram1D("track sign top", 3, -1., 2.);
            aida.histogram1D("track sign bottom", 3, -1., 2.);

            aida.histogram1D("nhits track top electron", 15, -.5, 14.5);
            aida.histogram1D("track momentum top electron", 100, 0., 7.);
            aida.histogram1D("cluster energy top electron", 100, 0., 7.);
            aida.histogram1D("E over P top electron", 100, 0., 2.);
            aida.histogram1D("nhits track top positron", 15, -.5, 14.5);
            aida.histogram1D("track momentum top positron", 100, 0., 7.);
            aida.histogram1D("cluster energy top positron", 100, 0., 7.);
            aida.histogram1D("E over P top positron", 100, 0., 2.);

            aida.histogram1D("nhits track bottom electron", 15, -.5, 14.5);
            aida.histogram1D("track momentum bottom electron", 100, 0., 7.);
            aida.histogram1D("cluster energy bottom electron", 100, 0., 7.);
            aida.histogram1D("E over P bottom electron", 100, 0., 2.);
            aida.histogram1D("nhits track bottom positron", 15, -.5, 14.5);
            aida.histogram1D("track momentum bottom positron", 100, 0., 7.);
            aida.histogram1D("cluster energy bottom positron", 100, 0., 7.);
            aida.histogram1D("E over P bottom positron", 100, 0., 2.);
        }
    }

    @Override
    protected void process(EventHeader event) {
        List<ReconstructedParticle> rpList = event.get(ReconstructedParticle.class, "FinalStateParticles_KF");
        for (ReconstructedParticle rp : rpList) {
            int pdgId = rp.getParticleIDUsed().getPDG();
            if (abs(pdgId) == 11) {
                Track track = rp.getTracks().get(0);
                int charge = -1;
                if (pdgId == -11) {
                    charge = 1;
                }
                String pType = "";
                switch (charge) {
                    case -1:
                        pType = " electron";
                        break;
                    case 1:
                        pType = " positron";
                        break;
                }
//            accumulatedTracks.add(track);
//            accumulatedBTracks.add(new BilliorTrack(track));

//TODO add some quality cuts on the track, e.g. nHits, associated cluster, charge, etc.
                if (!rp.getClusters().isEmpty()) {
                    double eclus = rp.getClusters().get(0).getEnergy();
                    double trackMom = rp.getMomentum().magnitude();
                    double EoverP = eclus / trackMom;
                    int nHits = track.getTrackerHits().size();
                    if (nHits > 9) {
                        if (track.getTrackStates().get(0).getTanLambda() > 0) {
//                accumulatedTracksTop.add(track);
                            accumulatedBTracksTop.add(new BilliorTrack(track));
                            aida.histogram1D("track sign top").fill(charge);
                            aida.histogram1D("nhits track top" + pType).fill(nHits);
                            aida.histogram1D("track momentum top" + pType).fill(trackMom);
                            aida.histogram1D("cluster energy top" + pType).fill(eclus);
                            aida.histogram1D("E over P top" + pType).fill(EoverP);
                        } else {
//                accumulatedTracksBot.add(track);
                            accumulatedBTracksBot.add(new BilliorTrack(track));
                            aida.histogram1D("track sign bottom").fill(charge);
                            aida.histogram1D("nhits track bottom" + pType).fill(nHits);
                            aida.histogram1D("track momentum bottom" + pType).fill(trackMom);
                            aida.histogram1D("cluster energy bottom" + pType).fill(eclus);
                            aida.histogram1D("E over P bottom" + pType).fill(EoverP);
                        }
                    } // end of check on nHits
                } // end of check on associated ECal cluster
            } // end of check on electron or positron
            if (accumulatedBTracksBot.size() >= _ntrks) {
//            System.out.println("Fitting " + accumulatedBTracksBot.size() + " bottom tracks");
                String vol = "_bottom";
                BilliorVertex vtx = vtxFitter.fitVertex(accumulatedBTracksBot);
                aida.histogram1D(vtxFolder + "vtx_x" + vol).fill(vtx.getPosition().x());
                aida.histogram1D(vtxFolder + "vtx_y" + vol).fill(vtx.getPosition().y());
                aida.histogram1D(vtxFolder + "vtx_z" + vol).fill(vtx.getPosition().z());
                aida.histogram2D(vtxFolder + "vtx_x_y" + vol).fill(vtx.getPosition().x(), vtx.getPosition().y());
                aida.histogram1D(vtxFolder + "vtx_ntrks" + vol).fill(accumulatedBTracksBot.size());
                aida.histogram1D(vtxFolder + "vtx_chi2" + vol).fill(vtx.getChi2());
                accumulatedBTracksBot.clear();
            }
            if (accumulatedBTracksTop.size() >= _ntrks) {
//            System.out.println("Fitting " + accumulatedBTracksTop.size() + " top tracks");
                String vol = "_top";
                BilliorVertex vtx = vtxFitter.fitVertex(accumulatedBTracksTop);
                aida.histogram1D(vtxFolder + "vtx_x" + vol).fill(vtx.getPosition().x());
                aida.histogram1D(vtxFolder + "vtx_y" + vol).fill(vtx.getPosition().y());
                aida.histogram1D(vtxFolder + "vtx_z" + vol).fill(vtx.getPosition().z());
                aida.histogram2D(vtxFolder + "vtx_x_y" + vol).fill(vtx.getPosition().x(), vtx.getPosition().y());
                aida.histogram1D(vtxFolder + "vtx_ntrks" + vol).fill(accumulatedBTracksTop.size());
                aida.histogram1D(vtxFolder + "vtx_chi2" + vol).fill(vtx.getChi2());
                accumulatedBTracksTop.clear();
            }
        }
    }

    @Override
    protected void endOfData() {

//        System.out.println("Size of accumulated Billior tracks: " + accumulatedBTracks.size());
//        System.out.println("Fitting Vertex");
//
//        beamPositionToUse[0] = 0.;
//        beamPositionToUse[0] = 0.;
//        beamPositionToUse[0] = -10;
//
//        BilliorVertexer vtxFitter = new BilliorVertexer(bField);
//        vtxFitter.setBeamSize(beamSize);
//        vtxFitter.setBeamPosition(beamPositionToUse);
//        vtxFitter.setStoreCovTrkMomList(_storeCovTrkMomList);
//        vtxFitter.setDebug(debug);
//        vtxFitter.doBeamSpotConstraint(false);
//
//        //Randomize the tracks ? Better to do something more reproducible
////        Collections.shuffle(accumulatedBTracks);
////        FitMultiVtx(accumulatedBTracks, vtxFitter, "");
//        FitMultiVtx(accumulatedBTracksTop, vtxFitter, "_top");
//        FitMultiVtx(accumulatedBTracksBot, vtxFitter, "_bottom");
    }

    private void FitMultiVtx(List<BilliorTrack> accTrks, BilliorVertexer vtxFitter, String vol) {

        int n_chunks = accTrks.size() / _ntrks;
        int n_rest = accTrks.size() % _ntrks;
        System.out.println("n_chunks = " + n_chunks);
        System.out.println("n_rest = " + n_rest);
        System.out.printf("size  = %d \n", (n_chunks * _ntrks + n_rest));

        if (accTrks.size() < 2) {
            return;
        }

        for (int i_chunk = 0; i_chunk < n_chunks; i_chunk++) {
            List<BilliorTrack> tracksForFit = accTrks.subList(i_chunk * _ntrks, (i_chunk + 1) * _ntrks);
            BilliorVertex vtx = vtxFitter.fitVertex(tracksForFit);
            aida.histogram1D(vtxFolder + "vtx_x" + vol).fill(vtx.getPosition().x());
            aida.histogram1D(vtxFolder + "vtx_y" + vol).fill(vtx.getPosition().y());
            aida.histogram1D(vtxFolder + "vtx_z" + vol).fill(vtx.getPosition().z());
            aida.histogram2D(vtxFolder + "vtx_x_y" + vol).fill(vtx.getPosition().x(), vtx.getPosition().y());
            aida.histogram1D(vtxFolder + "vtx_chi2" + vol).fill(vtx.getChi2());
            aida.histogram1D(vtxFolder + "vtx_ntrks" + vol).fill(_ntrks);
            //System.out.printf("vtx  %.5f %.5f %.5f\n",vtx.getPosition().x(), vtx.getPosition().y(),vtx.getPosition().z());
        }

        if (n_rest > 1) {

            List<BilliorTrack> tracksForFit = accTrks.subList(n_chunks * _ntrks, n_chunks * _ntrks + n_rest);
            BilliorVertex vtx = vtxFitter.fitVertex(tracksForFit);
            aida.histogram1D(vtxFolder + "vtx_x" + vol).fill(vtx.getPosition().x());
            aida.histogram1D(vtxFolder + "vtx_y" + vol).fill(vtx.getPosition().y());
            aida.histogram1D(vtxFolder + "vtx_z" + vol).fill(vtx.getPosition().z());
            aida.histogram2D(vtxFolder + "vtx_x_y" + vol).fill(vtx.getPosition().x(), vtx.getPosition().y());
            aida.histogram1D(vtxFolder + "vtx_ntrks" + vol).fill(n_rest);
            aida.histogram1D(vtxFolder + "vtx_chi2" + vol).fill(vtx.getChi2());
        }

        return;
    }

}
