#ifndef LARLITE_EXAMPLE_ANA_CXX
#define LARLITE_EXAMPLE_ANA_CXX

#include "example_ana.h"

namespace galleryfmwk {

bool example_ana::inFV(double x_vtx, double y_vtx, double z_vtx,
                       double x1, double x2, double y1, double y2, double z1, double z2)
{
	//is vertex in given FV?
	const double x_boundary1 = 0;
	const double x_boundary2 = 256.35;
	const double y_boundary1 = -116.5;
	const double y_boundary2 = 116.5;
	const double z_boundary1 = 0;
	const double z_boundary2 = 1036.8;

	if(x_vtx > x_boundary1 + x2 &&
	   x_vtx < x_boundary2 - x1 &&
	   y_vtx > y_boundary1 + y2 &&
	   y_vtx < y_boundary2 - y1 &&
	   z_vtx > z_boundary1 + z1 &&
	   z_vtx < z_boundary2 - z2)
	{
		return true;
	}
	return false;

}



bool example_ana::initialize() {

	//
	// This function is called in the beggining of event loop
	// Do all variable initialization you wish to do here.
	// If you have a histogram to fill in the event loop, for example_ana,
	// here is a good place to create one on the heap (i.e. "new TH1D").
	//

	//gROOT->SetBatch();

	x_boundary1 = 0;
	x_boundary2 = 256.35;
	y_boundary1 = -116.5;
	y_boundary2 = 116.5;
	z_boundary1 = 0;
	z_boundary2 = 1036.8;

	ub_total_vol = (x_boundary2 - x_boundary1) * (y_boundary2 - y_boundary1) *(z_boundary2 - z_boundary1);

	//fiducial cut beyond TPC
	fromWall = 0;

	num_cosmic = 0;
	num_primary_pfp = 0;
	num_nue = 0;
	num_numu = 0;
	cosmic_vertex_cut_pass = 0;
	cosmic_vertex_shower_cut_pass = 0;

	//max fv cut for plotting (in cm)
	fv_cut_max = 50;
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
	h_shwr_cut_theta_phi = new TH2D("h_shwr_theta_phi", "h_shwr_theta_phi", 60, -100, 100, 60, -190, 190);

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

	return true;
}


bool example_ana::analyze(gallery::Event * ev) {

	num_cosmic++;

	// For each file, loop over all events.

	// Get all of the tracks from the event:
	art::InputTag tracks_tag(_track_producer);
	art::InputTag showers_tag(_shower_producer);
	double cut_distance_to_point(_cut);
	double fiducial_volume_x_right(_right);
	double fiducial_volume_x_left(_left);
	double fiducial_volume_y_up(_up);
	double fiducial_volume_y_down(_down);
	double fiducial_volume_z_back(_back);
	double fiducial_volume_z_front(_front);

	auto const & tracks
	        = ev->getValidHandle<std::vector <recob::Track> >(tracks_tag);
	auto const & showers
	        = ev->getValidHandle<std::vector <recob::Shower> >(showers_tag);
	auto const & pfp
	        = ev->getValidHandle<std::vector <recob::PFParticle> > ("pandoraNu");
	auto const & pfparticles(*pfp);
	auto const & cosmic_pfp
	        = ev->getValidHandle<std::vector < recob::PFParticle> > ("pandoraCosmic");
	auto const & cosmicpfps(*cosmic_pfp);

	art::FindMany<recob::Vertex> vertex_for_pfp(pfp, *ev, "pandoraNu");
	art::FindMany<recob::Track> track_for_pfp(pfp, *ev, "pandoraNu");
	art::FindMany<recob::Shower> shower_for_pfp(pfp, *ev, "pandoraNu");

	art::FindMany<recob::Track> cosmic_track_for_pfp(cosmic_pfp, *ev, "pandoraCosmic");
	art::FindMany<recob::Shower> cosmic_shower_for_pfp(cosmic_pfp, *ev, "pandoraCosmic");

	std::vector < geoalgo::Trajectory_t > track_trajectory_list;
	std::vector < geoalgo::Trajectory_t > cosmic_track_trajectory_list;
	std::vector < double > cosmic_track_length_list;
	std::vector < geoalgo::Point_t> nue_vertex_list;
	std::vector < geoalgo::Point_t> shwr_vertex_list;
	std::vector < double > shwr_energy_list;
	std::vector < geoalgo::Point_t> shwr_vertex_lrgDist_list;

	std::vector < double > shwr_dir_vector;
	std::vector < std::vector < double > > shwr_dir_list;

	const int num_pfps = pfparticles.size();
	const int num_cosmics = cosmicpfps.size();

	//vector < recob::Track> cosmic_track_obj_list;

	num_nue_per_event = 0;
	//pfp loop
	for(std::size_t this_pfp = 0; this_pfp < num_pfps; this_pfp++)
	{
		//******************************
		//check for reconstructed vertex
		//******************************
		std::vector<recob::Vertex const*> vertex;
		vertex_for_pfp.get(this_pfp,vertex);
		if(vertex.size() ==0 )
		{
			if(_verbose == true) {std::cout << "No vertex association found!" << std::endl; }
			return false;
		}
		//get vertex vector
		double xyz [3];
		vertex.at(0)->XYZ(xyz);
		if( inFV(xyz[0], xyz[1], xyz[2], _right, _left, _up, _down, _back, _front) == false)
		{
			if(_verbose == true)
			{
				std::cout << "Reco vertex outside fiducial volume!" << std::endl;
				std::cout << "Skipping event..." << std::endl;
			}
			continue;
		}

		//************************************
		//check if pfp is neutrino-like object
		//************************************
		auto const pfparts = pfparticles.at(this_pfp);
		int shwr_daughters = 0;
		int trk_daughters = 0;
		if(pfparts.IsPrimary() == true)
		{
			num_primary_pfp++;
			//nues!
			if(pfparts.PdgCode() == 12)
			{
				num_nue_per_event++;
				num_nue++;
				geoalgo::Point_t const nue_vtx (xyz[0], xyz[1], xyz[2]);
				nue_vertex_list.push_back(nue_vtx);

				h_nue_like_vtx_xy->Fill(xyz[0], xyz[1]);
				h_nue_like_vtx_yz->Fill(xyz[2], xyz[1]);

				auto const daughters = pfparts.Daughters();
				for(std::size_t const i : daughters)
				{
					auto const daughter = pfparticles.at(i);
					//let's get the vertex associations for the daughters
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					if(d_vertex.size() == 0 )
					{
						if(_verbose == true)
						{
							std::cout << "No vertex association found for daughter!" << std::endl;
						}
						return false;
					}
					//get vertex vector
					double d_xyz [3];
					d_vertex.at(0)->XYZ(d_xyz);

					//shwr daughters
					if(daughter.PdgCode() == 11)
					{
						shwr_daughters++;
						h_nue_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						h_nue_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);

						//let's get the shower associations
						std::vector<recob::Shower const*> shower;
						shower_for_pfp.get(i, shower);
						if(shower.size() == 0)
						{
							if(_verbose) {std::cout << "No shower for this pfp shower!" << std::endl; }
							continue;
						}

						//let's check the distance between the shwr vtx and the nue vtx
						const double dist_x = d_xyz[0] - xyz[0];
						const double dist_y = d_xyz[1] - xyz[1];
						const double dist_z = d_xyz[2] - xyz[2];
						const double dist =
						        sqrt((dist_x * dist_x)+
						             (dist_y * dist_y)+
						             (dist_z * dist_z));
						h_nue_shwr_vtx_dist->Fill(dist);

						//what does it mean when the distance between the nue and shwr vtx
						//is large?
						if(dist >= 5 )
						{
							h_nue_like_shwr_lrgDist_vtx_xy->Fill(d_xyz[0], d_xyz[1]);
							h_nue_like_shwr_lrgDist_vtx_zy->Fill(d_xyz[2], d_xyz[1]);
							geoalgo::Point_t const shwr_vtx_lrgDist (d_xyz[0], d_xyz[1], d_xyz[2]);
							shwr_vertex_lrgDist_list.push_back(shwr_vtx_lrgDist);
						}

						//let's get the energy! Energy() - GeV?
						double total_energy = 0;
						const std::vector < double > plane_energy = shower.at(0)->Energy();
						const int best_plane = shower.at(0)->best_plane();
						total_energy = plane_energy.at(best_plane);
						h_nue_shwr_E->Fill(total_energy);
						shwr_energy_list.push_back(total_energy);

						geoalgo::Point_t const shwr_vtx (d_xyz[0], d_xyz[1], d_xyz[2]);
						shwr_vertex_list.push_back(shwr_vtx);

						//let's look at the shower directions
						const double dir_x = shower.at(0)->Direction().X();
						const double dir_y = shower.at(0)->Direction().Y();
						const double dir_z = shower.at(0)->Direction().Z();
						const double shwr_theta = TMath::ASin(dir_y) * (180/3.1415);
						const double shwr_phi = TMath::ATan2(dir_x, dir_z) * (180/3.1415);
						h_shwr_direction_xy->Fill(dir_x, dir_y);
						h_shwr_direction_zy->Fill(dir_z, dir_y);
						h_shwr_theta_phi->Fill(shwr_theta, shwr_phi);
						shwr_dir_vector.push_back(dir_x);
						shwr_dir_vector.push_back(dir_y);
						shwr_dir_vector.push_back(dir_z);
						shwr_dir_list.push_back(shwr_dir_vector);
						if(!shwr_dir_vector.empty()) {shwr_dir_vector.clear(); }


					}//end shwr daughters
					 //trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						h_nue_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						h_nue_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);

						std::vector<recob::Track const*> track;
						track_for_pfp.get(i, track);
						if(track.size() == 0)
						{
							if(_verbose) {std::cout << "No track for this pfp track!" << std::endl; }
							continue;
						}
						//let's construct the path of the tracks
						std::vector<geoalgo::Point_t> track_path;
						for(int pts = 0; pts < track.at(0)->NPoints(); pts++)
						{
							geoalgo::Point_t const track_point (
							        track.at(0)->LocationAtPoint(pts).X(),
							        track.at(0)->LocationAtPoint(pts).Y(),
							        track.at(0)->LocationAtPoint(pts).Z());
							track_path.push_back(track_point);
						}
						const geoalgo::Trajectory_t trj = track_path;
						if(!track_path.empty()) {track_path.clear(); }
						track_trajectory_list.push_back(trj);

						//let's get the track length!
						const double track_length = track.at(0)->Length();
						h_nue_trk_length->Fill(track_length);

						//let's look at the track directions
						const double dir_x = track.at(0)->VertexDirection().X();
						const double dir_y = track.at(0)->VertexDirection().Y();
						const double dir_z = track.at(0)->VertexDirection().Z();
						//std::cout << dir_x << ", " << dir_y << ", " << dir_z << std::endl;

					}//end nue track daughters
				}//end nue daughters
				h_nue_like_trk_daughters->Fill(trk_daughters);
				h_nue_like_daughters->Fill(shwr_daughters, trk_daughters);
				h_nue_like_daughters_logz->Fill(shwr_daughters, trk_daughters);
				for(int fv_cut = 0; fv_cut < fv_cut_max; fv_cut++)
				{
					if(inFV(xyz[0], xyz[1], xyz[2], fv_cut, fv_cut, fv_cut, fv_cut, fv_cut, fv_cut) == true)
					{
						h_nue_fv_cuts->Fill(fv_cut);
					}
					//just fv cut from top
					if(inFV(xyz[0], xyz[1], xyz[2], 0, 0, fv_cut, 0, 0, 0) == true)
					{
						h_nue_fv_top_cuts->Fill(fv_cut);
					}
				}
			}
			//numus!
			if(pfparts.PdgCode() == 14)
			{
				num_numu++;
				h_numu_like_vtx_xy->Fill(xyz[0], xyz[1]);
				h_numu_like_vtx_yz->Fill(xyz[2], xyz[1]);

				auto const daughters = pfparts.Daughters();
				for(std::size_t const i : daughters)
				{
					auto const daughter = pfparticles.at(i);
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					if(d_vertex.size() ==0 )
					{
						if(_verbose == true) {std::cout << "No vertex association found for daughter!" << std::endl; }
						return false;
					}
					//get vertex vector
					double d_xyz [3];
					d_vertex.at(0)->XYZ(d_xyz);

					//shwr daughters
					if(daughter.PdgCode() == 11)
					{
						shwr_daughters++;
						h_numu_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						h_numu_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}
					//trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						h_numu_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						h_numu_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}

				}
				h_numu_like_daughters->Fill(shwr_daughters, trk_daughters);
				// for(int fv_cut = 0; fv_cut < fv_cut_max; fv_cut++)
				// {
				//      if(inFV(xyz[0], xyz[1], xyz[2], fv_cut, fv_cut, fv_cut, fv_cut, fv_cut, fv_cut) == true)
				//      {
				//              h_numu_fv_cuts->Fill(fv_cut);
				//      }
				//      //just fv cut from top
				//      if(inFV(xyz[0], xyz[1], xyz[2], 0, 0, fv_cut, 0, 0, 0) == true)
				//      {
				//              h_nue_fv_top_cuts->Fill(fv_cut);
				//      }
				// }
			}
		}//end if nu-like
	}//end loop pfps

	h_num_nue_per_event->Fill(num_nue_per_event);

	//**************************
	//loop over pandora cosmics
	//*************************
	for(std::size_t this_cosmic = 0; this_cosmic < num_cosmics; this_cosmic++)
	{
		auto const cosmic = cosmicpfps.at(this_cosmic);
		if(cosmic.PdgCode() == 13)
		{
			//let's get the cosmic to track associations
			std::vector<recob::Track const*> cosmic_track;
			cosmic_track_for_pfp.get(this_cosmic, cosmic_track);
			if(cosmic_track.size() == 0)
			{
				if(_verbose == true) {std::cout << "No track for pfp!" << std::endl; }
				continue;
			}
			std::vector<geoalgo::Point_t> cosmic_track_path;
			for(int pts = 0; pts < cosmic_track.at(0)->NPoints(); pts++)
			{
				geoalgo::Point_t const cosmic_track_point (
				        cosmic_track.at(0)->LocationAtPoint(pts).X(),
				        cosmic_track.at(0)->LocationAtPoint(pts).Y(),
				        cosmic_track.at(0)->LocationAtPoint(pts).Z());
				cosmic_track_path.push_back(cosmic_track_point);
			}
			const geoalgo::Trajectory_t trj = cosmic_track_path;
			if(!cosmic_track_path.empty()) {cosmic_track_path.clear(); }
			cosmic_track_trajectory_list.push_back(trj);

			//let's get the track length
			const double cosmic_length = cosmic_track.at(0)->Length();
			h_cosmic_trk_length->Fill(cosmic_length);
			cosmic_track_length_list.push_back(cosmic_length);

			//let's get the track energy
			const double cosmic_trk_energy = cosmic_track.at(0)->StartMomentum();

		} //end loop tracks
		if(cosmic.PdgCode() == 11)
		{
			//let's get the cosmic to shower associations
			std::vector<recob::Shower const*> cosmic_shower;
			cosmic_shower_for_pfp.get(this_cosmic, cosmic_shower);
			if(cosmic_shower.size() == 0)
			{
				if(_verbose) {std::cout << "No shower for pfp!" << std::endl; }
				continue;
			}

			//let's get the cosmic shower energy
			double cosmic_shwr_energy;
			const std::vector < double > plane_energy = cosmic_shower.at(0)->Energy();
			const int bestplane = cosmic_shower.at(0)->best_plane();
			const double total_energy = plane_energy.at(bestplane);
		}        //end loop showers
	}        //end loop cosmics

	//*********************************************************
	//******************* Geometry studies! *******************
	//*********************************************************
	std::vector < geoalgo::Point_t > cut_nue_vertex;
	std::vector < geoalgo::Point_t > cut_nue_shwr_vertex;

	//closest point between **nue vertex** and cosmic track
	for(int nNue = 0; nNue < nue_vertex_list.size(); nNue++)
	{
		double cosmic_closest_point = 0;
		if(!cosmic_track_trajectory_list.empty())
		{
			geoalgo::Point_t nue_vertex = nue_vertex_list.at(nNue);
			cosmic_closest_point = _geo_algo_instance.SqDist(nue_vertex, cosmic_track_trajectory_list);
			h_nue_cosmic_closest->Fill(cosmic_closest_point);
			//*****************************************************************
			//let's see what happens if we remove nue vtx close to tagged cosmic
			//*****************************************************************
			if(cosmic_closest_point >= cut_distance_to_point)
			{
				cut_nue_vertex.push_back(nue_vertex);
				cosmic_vertex_cut_pass++;
			}
		}
		//closest point between **nue vertex** and nue track
		if(!track_trajectory_list.empty())
		{
			geoalgo::Point_t nue_vertex = nue_vertex_list.at(nNue);
			double closest_point = _geo_algo_instance.SqDist(nue_vertex, track_trajectory_list);
			h_nue_trk_closest->Fill(closest_point);
			if(cosmic_closest_point >= cut_distance_to_point) {h_nue_trk_closest_zoom->Fill(closest_point); }
		}
	}
	//closest point between **nue shwr vertex** and cosmic track
	bool first = true;
	for (int nE = 0; nE < shwr_vertex_list.size(); nE++)
	{
		if(!cosmic_track_trajectory_list.empty())
		{
			geoalgo::Point_t shwr_vertex = shwr_vertex_list.at(nE);
			double closest_point = _geo_algo_instance.SqDist(shwr_vertex, cosmic_track_trajectory_list);
			h_nue_shwr_cosmic_closest->Fill(closest_point);
			h_nue_shwr_cosmic_closest_vs_E->Fill(closest_point, shwr_energy_list.at(nE));
			h_nue_shwr_cosmic_closest_vs_y->Fill(closest_point, shwr_vertex[1]);
			h_nue_shwr_cosmic_closest_vs_E_zoom->Fill(closest_point, shwr_energy_list.at(nE));
			h_nue_shwr_cosmic_closest_vs_y_zoom->Fill(closest_point, shwr_vertex[1]);
			h_shwr_direction_y_vs_nearest_cosmic->Fill(closest_point, shwr_dir_list.at(nE).at(1));

			//***********************************************************************
			//let's see what happens if we remove nue shwr vtx close to tagged cosmic
			//***********************************************************************
			if(closest_point >= cut_distance_to_point)
			{
				cut_nue_shwr_vertex.push_back(shwr_vertex);
				cosmic_vertex_shower_cut_pass++;
				if(first == true)
				{
					first = false;
					h_nue_like_daughters_cuts->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
					h_nue_like_daughters_cuts_logz->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
				}
				h_shwr_direction_cut_xy->Fill(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(1));
				h_shwr_direction_cut_zy->Fill(shwr_dir_list.at(nE).at(2), shwr_dir_list.at(nE).at(1));
				const double shwr_cut_theta = TMath::ASin(shwr_dir_list.at(nE).at(1)) * (180/3.1415);
				const double shwr_cut_phi   = TMath::ATan2(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(2)) * (180/3.1415);
				h_shwr_cut_theta_phi->Fill(shwr_cut_theta, shwr_cut_phi);
			}
			//let's loop over all of the nue shwrs and the cosmic tracks
			//this will take the length to compute a cut out "fiducial" volume
			double closest_track = 100;
			int track_counter = 0;
			int this_track_counter = 0;
			for( auto cosmic_track_traj : cosmic_track_trajectory_list)
			{
				//which track is closest?
				track_counter++;
				double close_track = _geo_algo_instance.SqDist(shwr_vertex, cosmic_track_traj);
				if(close_track < closest_track)
				{
					closest_track = close_track;
					this_track_counter = track_counter;
				}
			}
			double closest_track_length = cosmic_track_length_list.at(this_track_counter);
			//calculate the volume of the cylinder
			double cylinder_vol = closest_track_length * 3.1415 * cut_distance_to_point * cut_distance_to_point;
			h_cylinder_vol->Fill(cylinder_vol/ub_total_vol);
		}
		//closest point between **nue shwr vertex** and nue track
		if(!track_trajectory_list.empty())
		{
			geoalgo::Point_t shwr_vertex = shwr_vertex_list.at(nE);
			double closest_point = _geo_algo_instance.SqDist(shwr_vertex, track_trajectory_list);
			h_nue_shwr_trk_closest->Fill(closest_point);
			int num_trks_nearby = 0;
			for(auto this_track : track_trajectory_list)
			{
				double this_closest_point = _geo_algo_instance.SqDist(shwr_vertex, this_track);
				if(this_closest_point <= 5) {num_trks_nearby++; }
			}
			h_num_trks_nearby->Fill(num_trks_nearby);
		}//end if event has track
	}//end loop nue shwr vtx
	 //***************************************************************
	 //loop for showers with a large distance from the reco nue vertex
	for(int nE = 0; nE < shwr_vertex_lrgDist_list.size(); nE++)
	{
		const int num_trks_lrgDist = track_trajectory_list.size();
		if(!track_trajectory_list.empty())
		{
			geoalgo::Point_t shwr_vertex_lrgDist = shwr_vertex_lrgDist_list.at(nE);

			double this_closest_point_lrgDist = _geo_algo_instance.SqDist(shwr_vertex_lrgDist, track_trajectory_list);
			h_nue_like_shwr_lrgDist_dist_to_track->Fill(this_closest_point_lrgDist);

			double this_closest_point_cosmic_lrgDist = _geo_algo_instance.SqDist(shwr_vertex_lrgDist, cosmic_track_trajectory_list);
			h_nue_like_shwr_lrgDist_dist_to_cosmic->Fill(this_closest_point_cosmic_lrgDist);
		}
		h_nue_like_shwr_lrgDist_num_trks->Fill(num_trks_lrgDist);
	}

	//cut on sqdist to remove small distances - does this change vertex ave. position
	for(auto vtx : cut_nue_shwr_vertex)
	{
		h_nue_shwr_cut_vtx_xy->Fill(vtx[0], vtx[1]);
		h_nue_shwr_cut_vtx_zy->Fill(vtx[2], vtx[1]);
	}

//*****************************
// *End Geometry Calculations*
//*****************************

	return true;
}



bool example_ana::finalize() {

	/*********************************
	** Histogram Saving and Editing **
	*///******************************
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

	std::cout << "Number of events: " << num_cosmic << std::endl;
	std::cout << "Number of primary pfps: " << num_primary_pfp << std::endl;
	std::cout << "Number of nue-like: " << num_nue << std::endl;
	std::cout << "Number of numu-like: " << num_numu << std::endl;
	std::cout << "Nue-like Showers Remaining after a " << _cut << " cm cut: " << cosmic_vertex_shower_cut_pass << std::endl;
	std::cout << "Nue-like Events  Remaining after a " << _cut << " cm cut: " << cosmic_vertex_cut_pass << std::endl;

	return true;
}


}
#endif
