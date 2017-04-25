#ifndef HISTO_MANAGER_H
#define HISTO_MANAGER_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"

#include <iostream>

//namespace histogram {

class h_manager {


public:

void gen_histograms(TDirectory * dir);
void draw_save();
void gen_histograms_fv(TDirectory * dir);
void draw_save_fv();

TH1D * h_nue_fv_cuts;
TH1D * h_nue_fv_top_cuts;
TH1D * h_numu_fv_cuts;
TH2D * h_nue_like_daughters;
TH2D * h_nue_like_daughters_cuts;
TH2D * h_nue_like_daughters_cuts_logz;
TH2D * h_nue_like_daughters_logz;
TH1D * h_nue_like_trk_daughters;
TH2D * h_numu_like_daughters;

TH2D * h_nue_like_shwr_daughters_xy;
TH2D * h_nue_like_shwr_daughters_yz;
TH2D * h_nue_like_trk_daughters_xy;
TH2D * h_nue_like_trk_daughters_yz;
TH2D * h_nue_like_vtx_xy;
TH2D * h_nue_like_vtx_yz;

TH2D * h_nue_like_shwr_lrgDist_vtx_xy;
TH2D * h_nue_like_shwr_lrgDist_vtx_zy;
TH1D * h_nue_like_shwr_lrgDist_dist_to_track;
TH1D * h_nue_like_shwr_lrgDist_num_trks;
TH1D * h_nue_like_shwr_lrgDist_dist_to_cosmic;

TH2D * h_numu_like_shwr_daughters_xy;
TH2D * h_numu_like_shwr_daughters_yz;
TH2D * h_numu_like_trk_daughters_xy;
TH2D * h_numu_like_trk_daughters_yz;
TH2D * h_numu_like_vtx_xy;
TH2D * h_numu_like_vtx_yz;

TH1D * h_nue_cosmic_closest;
TH1D * h_nue_shwr_cosmic_closest;
TH1D * h_nue_shwr_vtx_dist;

TH1D * h_nue_shwr_E;
TH2D * h_nue_shwr_cosmic_closest_vs_E;
TH2D * h_nue_shwr_cosmic_closest_vs_y;
TH2D * h_nue_shwr_cosmic_closest_vs_E_zoom;
TH2D * h_nue_shwr_cosmic_closest_vs_y_zoom;

TH1D * h_cosmic_trk_length;
TH1D * h_nue_trk_length;

TH1D * h_nue_trk_closest;
TH1D * h_nue_trk_closest_zoom;
TH1D * h_nue_shwr_trk_closest;

TH1D * h_num_trks_nearby;

TH2D * h_nue_shwr_cut_vtx_xy;
TH2D * h_nue_shwr_cut_vtx_zy;

TH1D * h_num_nue_per_event;

TH1D * h_cylinder_vol;

TH2D * h_shwr_direction_xy;
TH2D * h_shwr_direction_zy;
TH2D * h_shwr_direction_cut_xy;
TH2D * h_shwr_direction_cut_zy;
TH2D * h_shwr_theta_phi;
TH2D * h_shwr_cut_theta_phi;

TH2D * h_shwr_direction_y_vs_nearest_cosmic;

TH1D * h_shwr_length;
TH1D * h_shwr_end_width;
TH1D * h_shwr_open_angle;
TH2D * h_shwr_length_width;

TH1D * h_nue_total_energy;

TH1D * h_shwr_to_boundary;
TH1D * h_trk_intersect;
TH1D * h_trk_to_boundary;

TH1D * h_nue_shwr_max_vtx_to_vtx_dist;
TH1D * h_nue_shwr_vtx_to_vtx_dist;

TH1D * h_nu_vtx_x2;
TH1D * h_nu_vtx_y2;
TH1D * h_nu_vtx_z2;

TCanvas * c1;
TCanvas * c1b;
TCanvas * c1c;
TCanvas * c1d;
TCanvas * c1e;
TCanvas * c2;
TCanvas * c3;
TCanvas * c3b;
TCanvas * c4;
TCanvas * c5;
TCanvas * c6;
TCanvas * c7;
TCanvas * c8;
TCanvas * c9;
TCanvas * c10;
TCanvas * c11;
TCanvas * c12;
TCanvas * c13;
TCanvas * c14;
TCanvas * c15;
TCanvas * c16;
TCanvas * c17;
TCanvas * c17b;
TCanvas * c18;
TCanvas * c19;
TCanvas * c19b;
TCanvas * c19c;
TCanvas * c19d;
TCanvas * c19e;
TCanvas * c20a;
TCanvas * c20b;
TCanvas * c21a;
TCanvas * c21b;
TCanvas * c21c;
TCanvas * c22;
TCanvas * c23a;
TCanvas * c23b;
TCanvas * c24;
TCanvas * c25a;
TCanvas * c25b;
TCanvas * c25c;
TCanvas * c25d;
TCanvas * c25e;
TCanvas * c26;
TCanvas * c27a;
TCanvas * c27b;
TCanvas * c27c;
TCanvas * c27d;
TCanvas * c27e;
TCanvas * c27f;
TCanvas * c27g;
TCanvas * c28a;
TCanvas * c28b;
TCanvas * c28c;
TCanvas * c28d;
TCanvas * c29;
TCanvas * c30a;
TCanvas * c30b;
TCanvas * c30c;
TCanvas * c31a;
TCanvas * c31b;
TCanvas * c32a;
TCanvas * c32b;
TCanvas * c32c;


};

#endif
