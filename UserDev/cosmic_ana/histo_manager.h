#ifndef HISTO_MANAGER_H
#define HISTO_MANAGER_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"

namespace histogram {

void gen_histograms();
void draw_save();

extern TH1D * h_nue_fv_cuts;
extern TH1D * h_nue_fv_top_cuts;
extern TH1D * h_numu_fv_cuts;
extern TH2D * h_nue_like_daughters;
extern TH2D * h_nue_like_daughters_cuts;
extern TH2D * h_nue_like_daughters_cuts_logz;
extern TH2D * h_nue_like_daughters_logz;
extern TH1D * h_nue_like_trk_daughters;
extern TH2D * h_numu_like_daughters;

extern TH2D * h_nue_like_shwr_daughters_xy;
extern TH2D * h_nue_like_shwr_daughters_yz;
extern TH2D * h_nue_like_trk_daughters_xy;
extern TH2D * h_nue_like_trk_daughters_yz;
extern TH2D * h_nue_like_vtx_xy;
extern TH2D * h_nue_like_vtx_yz;

extern TH2D * h_nue_like_shwr_lrgDist_vtx_xy;
extern TH2D * h_nue_like_shwr_lrgDist_vtx_zy;
extern TH1D * h_nue_like_shwr_lrgDist_dist_to_track;
extern TH1D * h_nue_like_shwr_lrgDist_num_trks;
extern TH1D * h_nue_like_shwr_lrgDist_dist_to_cosmic;

extern TH2D * h_numu_like_shwr_daughters_xy;
extern TH2D * h_numu_like_shwr_daughters_yz;
extern TH2D * h_numu_like_trk_daughters_xy;
extern TH2D * h_numu_like_trk_daughters_yz;
extern TH2D * h_numu_like_vtx_xy;
extern TH2D * h_numu_like_vtx_yz;

extern TH1D * h_nue_cosmic_closest;
extern TH1D * h_nue_shwr_cosmic_closest;
extern TH1D * h_nue_shwr_vtx_dist;

extern TH1D * h_nue_shwr_E;
extern TH2D * h_nue_shwr_cosmic_closest_vs_E;
extern TH2D * h_nue_shwr_cosmic_closest_vs_y;
extern TH2D * h_nue_shwr_cosmic_closest_vs_E_zoom;
extern TH2D * h_nue_shwr_cosmic_closest_vs_y_zoom;

extern TH1D * h_cosmic_trk_length;
extern TH1D * h_nue_trk_length;

extern TH1D * h_nue_trk_closest;
extern TH1D * h_nue_trk_closest_zoom;
extern TH1D * h_nue_shwr_trk_closest;

extern TH1D * h_num_trks_nearby;

extern TH2D * h_nue_shwr_cut_vtx_xy;
extern TH2D * h_nue_shwr_cut_vtx_zy;

extern TH1D * h_num_nue_per_event;

extern TH1D * h_cylinder_vol;

extern TH2D * h_shwr_direction_xy;
extern TH2D * h_shwr_direction_zy;
extern TH2D * h_shwr_direction_cut_xy;
extern TH2D * h_shwr_direction_cut_zy;
extern TH2D * h_shwr_theta_phi;
extern TH2D * h_shwr_cut_theta_phi;

extern TH2D * h_shwr_direction_y_vs_nearest_cosmic;

extern TH1D * h_shwr_length;
extern TH1D * h_shwr_end_width;
extern TH1D * h_shwr_open_angle;

extern TCanvas * c1;
extern TCanvas * c1b;
extern TCanvas * c1c;
extern TCanvas * c1d;
extern TCanvas * c1e;
extern TCanvas * c2;
extern TCanvas * c3;
extern TCanvas * c3b;
extern TCanvas * c4;
extern TCanvas * c5;
extern TCanvas * c6;
extern TCanvas * c7;
extern TCanvas * c8;
extern TCanvas * c9;
extern TCanvas * c10;
extern TCanvas * c11;
extern TCanvas * c12;
extern TCanvas * c13;
extern TCanvas * c14;
extern TCanvas * c15;
extern TCanvas * c16;
extern TCanvas * c17;
extern TCanvas * c17b;
extern TCanvas * c18;
extern TCanvas * c19;
extern TCanvas * c19b;
extern TCanvas * c19c;
extern TCanvas * c19d;
extern TCanvas * c19e;
extern TCanvas * c20a;
extern TCanvas * c20b;
extern TCanvas * c21a;
extern TCanvas * c21b;
extern TCanvas * c21c;
extern TCanvas * c22;
extern TCanvas * c23a;
extern TCanvas * c23b;
extern TCanvas * c24;
extern TCanvas * c25a;
extern TCanvas * c25b;
extern TCanvas * c25c;
extern TCanvas * c25d;
extern TCanvas * c25e;
extern TCanvas * c26;
extern TCanvas * c27a;
extern TCanvas * c27b;
extern TCanvas * c27c;
extern TCanvas * c27d;
extern TCanvas * c27e;
extern TCanvas * c27f;
extern TCanvas * c27g;
extern TCanvas * c28a;
extern TCanvas * c28b;
extern TCanvas * c28c;

}

#endif
