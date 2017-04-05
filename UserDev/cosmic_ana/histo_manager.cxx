#include "histo_manager.h"


void h_manager::gen_histograms(TDirectory * dir) {
//max fv cut for plotting (in cm)

	std::cout << "Generating Histograms" << std::endl;

	TH1::AddDirectory(kFALSE);

	h_nue_fv_cuts = new TH1D("h_nue_fv_cuts", "h_nue_fv_cuts", 50, 0, 50);
	h_nue_fv_top_cuts = new TH1D("h_nue_fv_top_cuts", "h_nue_fv_top_cuts", 50, 0, 50);
	h_numu_fv_cuts = new TH1D("h_numu_fv_cuts", "h_numu_fv_cuts", 50, 0, 50);

	h_nue_like_daughters = new TH2D("h_nue_like_daughters", "h_nue-like_daughters", 6, 0, 6, 6, 0, 6);
	h_nue_like_daughters_cuts = new TH2D ("h_nue_like_daughters_cuts", "h_nue_like_daughters_cuts", 6, 0, 6, 6, 0, 6);
	h_nue_like_daughters_cuts_logz = new TH2D ("h_nue_like_daughters_cuts_logz", "h_nue_like_daughters_cuts_logz", 6, 0, 6, 6, 0, 6);
	h_nue_like_daughters_logz = new TH2D ("h_nue_like_daughters_logz", "h_nue_like_daughters_logz", 6, 0, 6, 6, 0, 6);
	h_nue_like_trk_daughters = new TH1D ("h_nue_like_trk_daughters", "h_nue-like_trk_daughters", 6, 0, 6);
	h_numu_like_daughters = new TH2D("h_numu_like_daughters", "h_numu-like_daughters", 6, 0, 6, 6, 0, 6);

	h_nue_like_shwr_daughters_xy = new TH2D("h_nue_like_shwr_daughters_xy", "h_nue_like_shwr_daughters_xy", 52, 0, 260, 60, -120, 120);
	h_nue_like_shwr_daughters_yz = new TH2D("h_nue_like_shwr_daughters_yz", "h_nue_like_shwr_daughters_yz", 50, 0, 1050, 60, -120, 120);
	h_nue_like_trk_daughters_xy = new TH2D("h_nue_like_trk_daughters_xy", "h_nue_like_trk_daughters_xy", 52, 0, 260, 60, -120, 120);
	h_nue_like_trk_daughters_yz = new TH2D("h_nue_like_trk_daughters_yz", "h_nue_like_trk_daughters_yz", 50, 0, 1050, 60, -120, 120);
	h_nue_like_vtx_xy = new TH2D("h_nue_like_vtx_xy", "h_nue_like_vtx_xy", 52, 0, 260, 60, -120, 120);
	h_nue_like_vtx_yz = new TH2D("h_nue_like_vtx_yz", "h_nue_like_vtx_yz", 50, 0, 1050, 60, -120, 120);

	h_numu_like_shwr_daughters_xy = new TH2D("h_numu_like_shwr_daughters_xy", "h_numu_like_shwr_daughters_xy", 50, 0, 260, 60, -120, 120);
	h_numu_like_shwr_daughters_yz = new TH2D("h_numu_like_shwr_daughters_yz", "h_numu_like_shwr_daughters_yz", 50, 0, 1050, 60, -120, 120);
	h_numu_like_trk_daughters_xy = new TH2D("h_numu_like_trk_daughters_xy", "h_numu_like_trk_daughters_xy", 50, 0, 260, 60, -120, 120);
	h_numu_like_trk_daughters_yz = new TH2D("h_numu_like_trk_daughters_yz", "h_numu_like_trk_daughters_yz", 50, 0, 1050, 60, -120, 120);
	h_numu_like_vtx_xy = new TH2D("h_numu_like_vtx_xy", "h_numu_like_vtx_xy", 50, 0, 260, 60, -120, 120);
	h_numu_like_vtx_yz = new TH2D("h_numu_like_vtx_yz", "h_numu_like_vtx_yz", 50, 0, 1050, 60, -120, 120);

	h_nue_like_shwr_lrgDist_vtx_xy = new TH2D("h_nue_like_shwr_lrgDist_vtx_xy", "h_nue_like_shwr_lrgDist_vtx_xy", 50, 0, 260, 60, -120, 120);
	h_nue_like_shwr_lrgDist_vtx_zy = new TH2D("h_nue_like_shwr_lrgDist_vtx_zy", "h_nue_like_shwr_lrgDist_vtx_zy", 50, 0, 1050, 60, -120, 120);
	h_nue_like_shwr_lrgDist_dist_to_track = new TH1D("h_nue_like_shwr_lrgDist_dist_to_track", "h_nue_like_shwr_lrgDist_dist_to_track", 50, 0, 100);
	h_nue_like_shwr_lrgDist_num_trks = new TH1D("h_nue_like_shwr_lrgDist_num_trks", "h_nue_like_shwr_lrgDist_num_trks", 10, 0, 10);
	h_nue_like_shwr_lrgDist_dist_to_cosmic = new TH1D("h_nue_like_shwr_lrgDist_dist_to_cosmic", "h_nue_like_shwr_lrgDist_dist_to_cosmic", 60, 0, 120);

	h_nue_cosmic_closest = new TH1D("h_nue_cosmic_closest", "h_nue_cosmic_closest", 60, 0, 60);
	h_nue_shwr_cosmic_closest = new TH1D("h_nue_shwr_cosmic_closest", "h_nue_shwr_cosmic_closest", 60, 0, 60);
	h_nue_shwr_vtx_dist = new TH1D("h_nue_shwr_vtx_dist", "h_nue_shwr_vtx_dist", 60, 0, 120);

	h_nue_shwr_E = new TH1D("h_nue_shwr_E", "h_nue_shwr_E", 100, 0, 2);
	h_nue_shwr_cosmic_closest_vs_E = new TH2D("h_nue_shwr_cosmic_closest_vs_E", "h_nue_shwr_cosmic_closest_vs_E", 30, 0, 120, 30, 0, 2);
	h_nue_shwr_cosmic_closest_vs_y = new TH2D("h_nue_shwr_cosmic_closest_vs_y", "h_nue_shwr_cosmic_closest_vs_y", 30, 0, 120, 30, -120, 120);
	h_nue_shwr_cosmic_closest_vs_E_zoom = new TH2D("h_nue_shwr_cosmic_closest_vs_E_zoom", "h_nue_shwr_cosmic_closest_vs_E_zoom", 30, 0, 20, 30, 0, 2);
	h_nue_shwr_cosmic_closest_vs_y_zoom = new TH2D("h_nue_shwr_cosmic_closest_vs_y_zoom", "h_nue_shwr_cosmic_closest_vs_y_zoom", 30, 0, 20, 30, -120, 120);

	h_cosmic_trk_length = new TH1D ("h_cosmic_trk_length", "h_cosmic_trk_length", 50, 0, 100);
	h_nue_trk_length = new TH1D("h_nue_trk_length", "h_nue_trk_length", 50, 0, 100);

	h_nue_trk_closest = new TH1D("h_nue_trk_closest", "h_nue_trk_closest", 60, 0, 60);
	h_nue_trk_closest_zoom = new TH1D ("h_nue_trk_closest_zoom", "h_nue_trk_closest_zoom", 60, 0, 20);
	h_nue_shwr_trk_closest = new TH1D("h_nue_shwr_trk_closest", "h_nue_shwr_trk_closest", 60, 0, 60);

	h_num_trks_nearby = new TH1D("h_num_trks_nearby", "h_num_trks", 10, 0, 10);

	h_nue_shwr_cut_vtx_xy = new TH2D("h_nue_shwr_cut_vtx_xy", "h_nue_shwr_cut_vtx_xy", 50, 0, 260, 60, -120, 120);
	h_nue_shwr_cut_vtx_zy = new TH2D("h_nue_shwr_cut_vtx_zy", "h_nue_shwr_cut_vtx_zy", 50, 0, 1050, 60, -120, 120);

	h_num_nue_per_event = new TH1D("h_num_nue_per_event", "h_num_nue_per_event", 10, 0, 10);

	h_cylinder_vol = new TH1D("h_cylinder_vol", "h_cylinder_vol", 50, 0, 100);

	h_shwr_direction_xy = new TH2D("h_shwr_direction_xy", "h_shwr_direction_xy", 100, -1, 1, 100, -1, 1);
	h_shwr_direction_zy = new TH2D("h_shwr_direction_zy", "h_shwr_direction_zy", 100, -1, 1, 100, -1, 1);
	h_shwr_direction_cut_xy = new TH2D("h_shwr_direction_cut_xy", "h_shwr_direction_cut_xy", 100, -1, 1, 100, -1, 1);
	h_shwr_direction_cut_zy = new TH2D("h_shwr_direction_cut_zy", "h_shwr_direction_cut_zy", 100, -1, 1, 100, -1, 1);
	h_shwr_direction_y_vs_nearest_cosmic = new TH2D("h_shwr_direction_y_vs_nearest_cosmic", "h_shwr_direction_y_vs_nearest_cosmic", 30, 0, 120, 30, -1, 1);
	h_shwr_theta_phi = new TH2D("h_shwr_theta_phi", "h_shwr_theta_phi", 60, -100, 100, 60, -190, 190);
	h_shwr_cut_theta_phi = new TH2D("h_shwr_cut_theta_phi", "h_shwr_cut_theta_phi", 60, -100, 100, 60, -190, 190);

	h_shwr_length = new TH1D ("h_shwr_length", "h_shwr_length", 50, 0, 50);
	h_shwr_end_width = new TH1D ("h_shwr_end_width", "h_shwr_end_width", 35, 0, 35);
	h_shwr_open_angle = new TH1D ("h_shwr_open_angle", "h_shwr_open_angle", 60, 0, 180);

	c1 = new TCanvas();
	c1b = new TCanvas();
	c1c = new TCanvas();
	c1d = new TCanvas();
	c1e = new TCanvas();
	c2 = new TCanvas();
	c3 = new TCanvas();
	c3b = new TCanvas();
	c4 = new TCanvas();
	c5 = new TCanvas();
	c6 = new TCanvas();
	c7 = new TCanvas();
	c8 = new TCanvas();
	c9 = new TCanvas();
	c10 = new TCanvas();
	c11 = new TCanvas();
	c12 = new TCanvas();
	c13 = new TCanvas();
	c14 = new TCanvas();
	c15 = new TCanvas();
	c16 = new TCanvas();
	c17 = new TCanvas();
	c17b = new TCanvas();
	c18 = new TCanvas();
	c19 = new TCanvas();
	c19b = new TCanvas();
	c19c = new TCanvas();
	c19d = new TCanvas();
	c19e = new TCanvas();
	c20a = new TCanvas();
	c20b = new TCanvas();
	c21a = new TCanvas();
	c21b = new TCanvas();
	c21c = new TCanvas();
	c22 = new TCanvas();
	c23a = new TCanvas();
	c23b = new TCanvas();
	c24 = new TCanvas();
	c25a = new TCanvas();
	c25b = new TCanvas();
	c25c = new TCanvas();
	c25d = new TCanvas();
	c25e = new TCanvas();
	c26 = new TCanvas();
	c27a = new TCanvas();
	c27b = new TCanvas();
	c27c = new TCanvas();
	c27d = new TCanvas();
	c27e = new TCanvas();
	c27f = new TCanvas();
	c27g = new TCanvas();
	c28a = new TCanvas();
	c28b = new TCanvas();
	c28c = new TCanvas();

	h_numu_like_vtx_xy->Fill(1, 1);

	std::cout << "Finished Generating Histograms" << std::endl;

}

void h_manager::draw_save()
{

	std::cout << "Start Saving Histograms" << std::endl;

	c1->cd();
	h_nue_like_daughters->Draw("colz");
	h_nue_like_daughters->GetXaxis()->SetTitle("showers");
	h_nue_like_daughters->GetYaxis()->SetTitle("tracks");
	c1->Print("nue-like_daughters.pdf");
	c1b->cd();
	h_nue_like_trk_daughters->Draw();
	h_nue_like_trk_daughters->GetXaxis()->SetTitle("Tracks");
	h_nue_like_trk_daughters->GetYaxis()->SetTitle("Events");
	c1b->Print("nue-like_trk_daughters.pdf");
	c1c->cd();
	h_nue_like_daughters_cuts->Draw("colz");
	h_nue_like_daughters_cuts->GetXaxis()->SetTitle("showers");
	h_nue_like_daughters_cuts->GetYaxis()->SetTitle("tracks");
	c1c->Print("nue-like_daughters_cuts.pdf");
	c1d->cd();
	c1d->SetLogz();
	h_nue_like_daughters_logz->Draw("colz");
	h_nue_like_daughters_logz->GetXaxis()->SetTitle("showers");
	h_nue_like_daughters_logz->GetYaxis()->SetTitle("tracks");
	c1d->Print("nue-like_daughters_logz.pdf");
	c1e->cd();
	c1e->SetLogz();
	h_nue_like_daughters_cuts_logz->Draw("colz");
	h_nue_like_daughters_cuts_logz->GetXaxis()->SetTitle("showers");
	h_nue_like_daughters_cuts_logz->GetYaxis()->SetTitle("tracks");
	c1e->Print("nue-like_daughters_cuts_logz.pdf");

	c2->cd();
	h_numu_like_daughters->Draw("colz");
	h_numu_like_daughters->GetXaxis()->SetTitle("showers");
	h_numu_like_daughters->GetYaxis()->SetTitle("tracks");
	c2->Print("numu-like_daughters.pdf");

	c3->cd();
	h_nue_fv_cuts->Draw();
	h_nue_fv_cuts->GetXaxis()->SetTitle("Fiducial Volume Cut [cm]");
	h_nue_fv_cuts->GetYaxis()->SetTitle("Events in Volume");
	c3->Print("nue-like_fiducial_volume.pdf");
	c3b->cd();
	h_nue_fv_top_cuts->Draw();
	h_nue_fv_top_cuts->GetXaxis()->SetTitle("Fiducial Volume from Top [cm]");
	h_nue_fv_top_cuts->GetYaxis()->SetTitle("Events in Volume");
	c3b->Print("nue-like_fiducial_volume_y.pdf");
	c4->cd();
	h_numu_fv_cuts->Draw();
	h_numu_fv_cuts->GetXaxis()->SetTitle("Fiducial Volume Cut [cm]");
	h_numu_fv_cuts->GetYaxis()->SetTitle("Events in Volume");
	c4->Print("numu-like_fiducial_volume.pdf");

	c5->cd();
	h_nue_like_shwr_daughters_xy->Draw("colz");
	h_nue_like_shwr_daughters_xy->GetXaxis()->SetTitle("x [cm]");
	h_nue_like_shwr_daughters_xy->GetYaxis()->SetTitle("y [cm]");
	c5->Print("nue-like_shwr_daughters_xy.pdf");
	c6->cd();
	h_nue_like_shwr_daughters_yz->Draw("colz");
	h_nue_like_shwr_daughters_yz->GetXaxis()->SetTitle("z [cm]");
	h_nue_like_shwr_daughters_yz->GetYaxis()->SetTitle("y [cm]");
	c6->Print("nue-like_shwr_daughters_zy.pdf");
	c7->cd();
	h_nue_like_trk_daughters_xy->Draw("colz");
	h_nue_like_trk_daughters_xy->GetXaxis()->SetTitle("x [cm]");
	h_nue_like_trk_daughters_xy->GetYaxis()->SetTitle("y [cm]");
	c7->Print("nue-like_trk_daughters_xy.pdf");
	c8->cd();
	h_nue_like_trk_daughters_yz->Draw("colz");
	h_nue_like_trk_daughters_yz->GetXaxis()->SetTitle("z [cm]");
	h_nue_like_trk_daughters_yz->GetYaxis()->SetTitle("y [cm]");
	c8->Print("nue-like_trk_daughters_zy.pdf");
	c9->cd();
	h_nue_like_vtx_xy->Draw("colz");
	h_nue_like_vtx_xy->GetXaxis()->SetTitle("x [cm]");
	h_nue_like_vtx_xy->GetYaxis()->SetTitle("y [cm]");
	c9->Print("nue-like_vtx_xy.pdf");
	c10->cd();
	h_nue_like_vtx_yz->Draw("colz");
	h_nue_like_vtx_yz->GetXaxis()->SetTitle("z [cm]");
	h_nue_like_vtx_yz->GetYaxis()->SetTitle("y [cm]");
	c10->Print("nue-like_vtx_zy.pdf");

	c11->cd();
	h_numu_like_shwr_daughters_xy->Draw("colz");
	h_numu_like_shwr_daughters_xy->GetXaxis()->SetTitle("x [cm]");
	h_numu_like_shwr_daughters_xy->GetYaxis()->SetTitle("y [cm]");
	c11->Print("numu-like_shwr_daughters_xy.pdf");
	c12->cd();
	h_numu_like_shwr_daughters_yz->Draw("colz");
	h_numu_like_shwr_daughters_yz->GetXaxis()->SetTitle("z [cm]");
	h_numu_like_shwr_daughters_yz->GetYaxis()->SetTitle("y [cm]");
	c12->Print("numu-like_shwr_dauhters_zy.pdf");
	c13->cd();
	h_numu_like_trk_daughters_xy->Draw("colz");
	h_numu_like_trk_daughters_xy->GetXaxis()->SetTitle("x [cm]");
	h_numu_like_trk_daughters_xy->GetYaxis()->SetTitle("y [cm]");
	c13->Print("numu-like_trk_daughters_xy.pdf");
	c14->cd();
	h_numu_like_trk_daughters_yz->Draw("colz");
	h_numu_like_trk_daughters_yz->GetXaxis()->SetTitle("z [cm]");
	h_numu_like_trk_daughters_yz->GetYaxis()->SetTitle("y [cm]");
	c14->Print("numu-like_trk_daughters_zy.pdf");
	c15->cd();
	h_numu_like_vtx_xy->Draw("colz");
	h_numu_like_vtx_xy->GetXaxis()->SetTitle("x [cm]");
	h_numu_like_vtx_xy->GetYaxis()->SetTitle("y [cm]");
	c15->Print("numu-like_vtx_xy.pdf");
	c16->cd();
	h_numu_like_vtx_yz->Draw("colz");
	h_numu_like_vtx_yz->GetXaxis()->SetTitle("z [cm]");
	h_numu_like_vtx_yz->GetYaxis()->SetTitle("y [cm]");
	c16->Print("numu-like_vtx_zy.pdf");

	c17->cd();
	h_nue_cosmic_closest->Draw();
	h_nue_cosmic_closest->GetXaxis()->SetTitle("Distance [cm]");
	h_nue_cosmic_closest->GetYaxis()->SetTitle("Events");
	c17->Print("nue-like_cosmic_closest.pdf");
	c17b->cd();
	h_nue_shwr_cosmic_closest->Draw();
	h_nue_shwr_cosmic_closest->GetXaxis()->SetTitle("Distance [cm]");
	h_nue_shwr_cosmic_closest->GetYaxis()->SetTitle("Events");
	c17b->Print("nue-like_shwr_cosmic_closest.pdf");

	c18->cd();
	h_nue_shwr_vtx_dist->Draw();
	h_nue_shwr_vtx_dist->GetXaxis()->SetTitle("Distance [cm]");
	h_nue_shwr_vtx_dist->GetYaxis()->SetTitle("Events");
	c18->Print("nue-like_shwr_vtx_distance.pdf");

	c19->cd();
	h_nue_shwr_E->Draw();
	h_nue_shwr_E->GetXaxis()->SetTitle("Total Shower Energy [GeV]");
	h_nue_shwr_E->GetYaxis()->SetTitle("Events");
	c19->Print("nue-like_shwr_E.pdf");
	c19b->cd();
	h_nue_shwr_cosmic_closest_vs_E->Draw("colz");
	h_nue_shwr_cosmic_closest_vs_E->GetXaxis()->SetTitle("Distance to Nearest Cosmic Track [cm]");
	h_nue_shwr_cosmic_closest_vs_E->GetYaxis()->SetTitle("Total Shower Energy [GeV]");
	h_nue_shwr_cosmic_closest_vs_E->SetStats(kFALSE);
	c19b->Print("nue-like_shwr_vtx_distance_vs_E.pdf");
	c19c->cd();
	h_nue_shwr_cosmic_closest_vs_y->Draw("colz");
	h_nue_shwr_cosmic_closest_vs_y->GetXaxis()->SetTitle("Distance to nearest cosmic track [cm]");
	h_nue_shwr_cosmic_closest_vs_y->GetYaxis()->SetTitle("y [cm]");
	h_nue_shwr_cosmic_closest_vs_y->SetStats(kFALSE);
	c19c->Print("nue-like_shwr_cosmic_closest_vs_y.pdf");
	c19d->cd();
	h_nue_shwr_cosmic_closest_vs_E_zoom->Draw("colz");
	h_nue_shwr_cosmic_closest_vs_E_zoom->GetXaxis()->SetTitle("Distance to Nearest Cosmic Track [cm]");
	h_nue_shwr_cosmic_closest_vs_E_zoom->GetYaxis()->SetTitle("Total Shower Energy [GeV]");
	h_nue_shwr_cosmic_closest_vs_E_zoom->SetStats(kFALSE);
	c19d->Print("nue-like_shwr_vtx_distance_vs_E_zoom.pdf");
	c19e->cd();
	h_nue_shwr_cosmic_closest_vs_y_zoom->Draw("colz");
	h_nue_shwr_cosmic_closest_vs_y_zoom->GetXaxis()->SetTitle("Distance to nearest cosmic track [cm]");
	h_nue_shwr_cosmic_closest_vs_y_zoom->GetYaxis()->SetTitle("y [cm]");
	h_nue_shwr_cosmic_closest_vs_y_zoom->SetStats(kFALSE);
	c19e->Print("nue-like_shwr_cosmic_closest_vs_y_zoom.pdf");

	c20a->cd();
	h_cosmic_trk_length->Draw();
	h_cosmic_trk_length->GetXaxis()->SetTitle("Length [cm]");
	h_cosmic_trk_length->GetYaxis()->SetTitle("Events");
	c20a->Print("cosmic_trk_length.pdf");
	c20b->cd();
	h_nue_trk_length->Draw();
	h_nue_trk_length->GetXaxis()->SetTitle("Length [cm]");
	h_nue_trk_length->GetYaxis()->SetTitle("Events");
	c20b->Print("nue-like_trk_length.pdf");

	c21a->cd();
	h_nue_trk_closest->Draw();
	h_nue_trk_closest->GetXaxis()->SetTitle("Distance to nearest track [cm]");
	h_nue_trk_closest->GetYaxis()->SetTitle("Events");
	c21a->Print("nue-like_trk_closest.pdf");
	c21b->cd();
	h_nue_shwr_trk_closest->Draw();
	h_nue_shwr_trk_closest->GetXaxis()->SetTitle("Distance to nearest track [cm]");
	h_nue_shwr_trk_closest->GetYaxis()->SetTitle("Events");
	c21b->Print("nue-like_shwr_trk_closest.pdf");
	c21c->cd();
	h_nue_trk_closest_zoom->Draw();
	h_nue_trk_closest_zoom->GetXaxis()->SetTitle("Distance to nearest track [cm]");
	h_nue_trk_closest_zoom->GetYaxis()->SetTitle("Events");
	c21c->Print("nue-like_trk_closest_zoom.pdf");

	c22->cd();
	h_num_trks_nearby->Draw();
	h_num_trks_nearby->GetXaxis()->SetTitle("Number nearby tracks");
	h_num_trks_nearby->GetYaxis()->SetTitle("Events");
	c22->Print("nue-like_shwr_nearby_trks.pdf");

	c23a->cd();
	h_nue_shwr_cut_vtx_xy->Draw("colz");
	h_nue_shwr_cut_vtx_xy->GetXaxis()->SetTitle("x [cm]");
	h_nue_shwr_cut_vtx_xy->GetYaxis()->SetTitle("y [cm]");
	c23a->Print("nue-like_shwr_cut_vtx_xy.pdf");
	c23b->cd();
	h_nue_shwr_cut_vtx_zy->Draw("colz");
	h_nue_shwr_cut_vtx_zy->GetXaxis()->SetTitle("z [cm]");
	h_nue_shwr_cut_vtx_zy->GetYaxis()->SetTitle("y [cm]");
	c23b->Print("nue-like_shwr_cut_vtx_zy.pdf");

	c24->cd();
	h_num_nue_per_event->Draw();
	h_num_nue_per_event->GetXaxis()->SetTitle("Number of Nue-like per Event");
	h_num_nue_per_event->GetYaxis()->SetTitle("Events");
	c24->Print("num_nue_per_event.pdf");

	c25a->cd();
	h_nue_like_shwr_lrgDist_vtx_xy->Draw("colz");
	h_nue_like_shwr_lrgDist_vtx_xy->GetXaxis()->SetTitle("x [cm]");
	h_nue_like_shwr_lrgDist_vtx_xy->GetYaxis()->SetTitle("y [cm]");
	c25a->Print("nue-like_shwr_lrgDist_vtx_xy.pdf");
	c25b->cd();
	h_nue_like_shwr_lrgDist_vtx_zy->Draw("colz");
	h_nue_like_shwr_lrgDist_vtx_zy->GetXaxis()->SetTitle("z [cm]");
	h_nue_like_shwr_lrgDist_vtx_zy->GetYaxis()->SetTitle("y [cm]");
	c25b->Print("nue-like_shwr_lrgDist_vtx_zy.pdf");
	c25c->cd();
	h_nue_like_shwr_lrgDist_dist_to_track->Draw();
	h_nue_like_shwr_lrgDist_dist_to_track->GetXaxis()->SetTitle("Distance to Nearest PandoraNu Track [cm]");
	h_nue_like_shwr_lrgDist_dist_to_track->GetYaxis()->SetTitle("Events");
	c25c->Print("nue-like_shwr_lrgDist_dist_to_track.pdf");
	c25d->cd();
	h_nue_like_shwr_lrgDist_num_trks->Draw();
	h_nue_like_shwr_lrgDist_num_trks->GetXaxis()->SetTitle("Number of PandoraNu Tracks");
	h_nue_like_shwr_lrgDist_num_trks->GetYaxis()->SetTitle("Events");
	c25d->Print("nue-like_shwr_lrgDist_num_trks.pdf");
	c25e->cd();
	h_nue_like_shwr_lrgDist_dist_to_cosmic->Draw();
	h_nue_like_shwr_lrgDist_dist_to_cosmic->GetXaxis()->SetTitle("Distance to Nearest Cosmic Track [cm]");
	h_nue_like_shwr_lrgDist_dist_to_cosmic->GetYaxis()->SetTitle("Events");
	c25e->Print("nue-like_shwr_lrgDist_dist_to_cosmic.pdf");

	c26->cd();
	h_cylinder_vol->Draw();
	h_cylinder_vol->GetXaxis()->SetTitle("Fiducial Volume Loss [cm^3]");
	h_cylinder_vol->GetYaxis()->SetTitle("Events");
	c26->Print("nue-like_cylinder_vol_cut.pdf");

	c27a->cd();
	h_shwr_direction_xy->Draw("colz");
	h_shwr_direction_xy->GetXaxis()->SetTitle("x Dir Cos");
	h_shwr_direction_xy->GetYaxis()->SetTitle("y Dir Cos");
	c27a->Print("nue-like_shwr_dir_xy.pdf");
	c27b->cd();
	h_shwr_direction_zy->Draw("colz");
	h_shwr_direction_zy->GetXaxis()->SetTitle("z Dir Cos");
	h_shwr_direction_zy->GetYaxis()->SetTitle("y Dir Cos");
	c27b->Print("nue-like_shwr_dir_zy.pdf");
	c27c->cd();
	h_shwr_direction_cut_xy->Draw("colz");
	h_shwr_direction_cut_xy->GetXaxis()->SetTitle("x Dir Cos");
	h_shwr_direction_cut_xy->GetYaxis()->SetTitle("y Dir Cos");
	c27c->Print("nue-like_shwr_dir_cut_xy.pdf");
	c27d->cd();
	h_shwr_direction_cut_zy->Draw("colz");
	h_shwr_direction_cut_zy->GetXaxis()->SetTitle("z Dir Cos");
	h_shwr_direction_cut_zy->GetYaxis()->SetTitle("y Dir Cos");
	c27d->Print("nue-like_shwr_dir_cut_zy.pdf");
	c27e->cd();
	h_shwr_direction_y_vs_nearest_cosmic->Draw("colz");
	h_shwr_direction_y_vs_nearest_cosmic->GetYaxis()->SetTitle("y Dir Cos");
	h_shwr_direction_y_vs_nearest_cosmic->GetXaxis()->SetTitle("Distance to Nearest Cosmic Track [cm]");
	c27e->Print("nue-like_shwr_dir_y_vs_nearest_cosmic.pdf");
	c27f->cd();
	h_shwr_theta_phi->Draw("colz");
	h_shwr_theta_phi->GetYaxis()->SetTitle("Phi [Degrees]");
	h_shwr_theta_phi->GetXaxis()->SetTitle("Theta [Degrees]");
	c27f->Print("nue-like_shwr_theta_phi.pdf");
	c27g->cd();
	h_shwr_cut_theta_phi->Draw("colz");
	h_shwr_cut_theta_phi->GetYaxis()->SetTitle("Phi [Degrees]");
	h_shwr_cut_theta_phi->GetXaxis()->SetTitle("Theta [Degrees]");
	c27g->Print("nue-like_shwr_cut_theta_phi.pdf");

	c28a->cd();
	h_shwr_open_angle->Draw();
	h_shwr_open_angle->GetYaxis()->SetTitle("Events");
	h_shwr_open_angle->GetXaxis()->SetTitle("Opening Angle [Degrees]");
	c28a->Print("nue-like_shwr_open_angle.pdf");
	c28b->cd();
	h_shwr_length->Draw();
	h_shwr_length->GetYaxis()->SetTitle("Events");
	h_shwr_length->GetXaxis()->SetTitle("Shower Length [cm]");
	c28b->Print("nue-like_shwr_length.pdf");
	c28c->cd();
	h_shwr_end_width->Draw();
	h_shwr_end_width->GetYaxis()->SetTitle("Events");
	h_shwr_end_width->GetXaxis()->SetTitle("Shower End Width [cm]");
	c28c->Print("nue-like_shwr_end_width.pdf");


	std::cout << "Finished Saving Histograms" << std::endl;

}
