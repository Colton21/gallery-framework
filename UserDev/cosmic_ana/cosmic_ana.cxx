#ifndef GALLERY_FMWK_COSMIC_ANA_CXX
#define GALLERY_FMWK_COSMIC_ANA_CXX

#include "cosmic_ana.h"

namespace galleryfmwk {

bool cosmic_ana::initialize() {

	//
	// This function is called in the beginning of event loop
	// Do all variable initialization you wish to do here.
	// If you have a histogram to fill in the event loop, for example_ana,
	// here is a good place to create one on the heap (i.e. "new TH1D").
	//

	_h_manager_instance.gen_histograms(gDirectory);

	double fiducial_volume_x_right(_right);
	double fiducial_volume_x_left(_left);
	double fiducial_volume_y_up(_up);
	double fiducial_volume_y_down(_down);
	double fiducial_volume_z_back(_back);
	double fiducial_volume_z_front(_front);

	x_boundary1 = 0      + fiducial_volume_x_left;
	x_boundary2 = 256.35 - fiducial_volume_x_right;
	y_boundary1 = -116.5 + fiducial_volume_y_down;
	y_boundary2 = 116.5  - fiducial_volume_y_up;
	z_boundary1 = 0      + fiducial_volume_z_back;
	z_boundary2 = 1036.8 - fiducial_volume_z_front;

	ub_total_vol = (x_boundary2 - x_boundary1) * (y_boundary2 - y_boundary1) *(z_boundary2 - z_boundary1);

	//fiducial cut beyond TPC
	fromWall = 0;

	num_cosmic = 0;
	num_primary_pfp = 0;
	num_nue = 0;
	num_numu = 0;
	cosmic_vertex_cut_pass = 0;
	cosmic_vertex_shower_cut_pass = 0;

	fv_cut_max = 50;


	return true;
}


bool cosmic_ana::analyze(gallery::Event * ev) {

	num_cosmic++;

	// For each file, loop over all events.
	// Get all of the tracks from the event:
	art::InputTag pfp_tag(_pfp_tag);
	art::InputTag pfp_cosmic_tag(_pfp_cosmic_tag);
	double cut_distance_to_point(_cut);

	auto const & pfp
	        = ev->getValidHandle<std::vector <recob::PFParticle> > (pfp_tag);
	auto const & pfparticles(*pfp);
	auto const & cosmic_pfp
	        = ev->getValidHandle<std::vector < recob::PFParticle> > (pfp_cosmic_tag);
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
	std::vector < double > trk_energy_list;
	std::vector < geoalgo::Point_t> shwr_vertex_lrgDist_list;

	std::vector < double > shwr_dir_vector;
	std::vector < std::vector < double > > shwr_dir_list;

	num_pfps = 0;
	num_cosmics = 0;
	num_nue_per_event = 0;

	//pfp loop
	for(auto pfparts : pfparticles)
	{
		//******************************
		//check for reconstructed vertex
		//******************************
		std::vector<recob::Vertex const*> vertex;

		vertex_for_pfp.get(num_pfps,vertex);
		//get vertex vector
		double xyz [3];
		vertex.at(0)->XYZ(xyz);

		//************************************
		//check if pfp is neutrino-like object
		//************************************
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

				_h_manager_instance.h_nue_like_vtx_xy->Fill(xyz[0], xyz[1]);
				_h_manager_instance.h_nue_like_vtx_yz->Fill(xyz[2], xyz[1]);

				for(std::size_t const i : pfparts.Daughters())
				{

					auto const daughter = pfparticles.at(i);
					//let's get the vertex associations for the daughters
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					//get vertex vector
					double d_xyz [3];
					d_vertex.at(0)->XYZ(d_xyz);

					//shwr daughters
					if(daughter.PdgCode() == 11)
					{
						shwr_daughters++;
						//let's get the shower associations
						std::vector<recob::Shower const*> shower;
						shower_for_pfp.get(i, shower);

						//shower daughter vertices
						_h_manager_instance.h_nue_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						_h_manager_instance.h_nue_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
						geoalgo::Point_t const shwr_vtx (d_xyz[0], d_xyz[1], d_xyz[2]);
						shwr_vertex_list.push_back(shwr_vtx);

						//let's check the distance between the shwr vtx and the nue vtx
						const double dist = _utility_instance.geo_distance(d_xyz[0], xyz[0], d_xyz[1], xyz[1], d_xyz[2], xyz[2]);
						_h_manager_instance.h_nue_shwr_vtx_dist->Fill(dist);

						//what does it mean when the distance between the nue and shwr vtx
						//is large?
						if(dist >= 5 )
						{
							_h_manager_instance.h_nue_like_shwr_lrgDist_vtx_xy->Fill(d_xyz[0], d_xyz[1]);
							_h_manager_instance.h_nue_like_shwr_lrgDist_vtx_zy->Fill(d_xyz[2], d_xyz[1]);
							geoalgo::Point_t const shwr_vtx_lrgDist (d_xyz[0], d_xyz[1], d_xyz[2]);
							shwr_vertex_lrgDist_list.push_back(shwr_vtx_lrgDist);
						}

						//let's get the energy! Energy() - GeV?
						const std::vector < double > plane_energy = shower.at(0)->Energy();
						const int best_plane = shower.at(0)->best_plane();
						const double total_energy = plane_energy.at(best_plane);
						_h_manager_instance.h_nue_shwr_E->Fill(total_energy);
						shwr_energy_list.push_back(total_energy);

						//let's look at the shower directions
						const double dir_x = shower.at(0)->Direction().X();
						const double dir_y = shower.at(0)->Direction().Y();
						const double dir_z = shower.at(0)->Direction().Z();
						const double shwr_theta = TMath::ASin(dir_y) * (180/3.1415);
						const double shwr_phi = TMath::ATan2(dir_x, dir_z) * (180/3.1415);
						_h_manager_instance.h_shwr_direction_xy->Fill(dir_x, dir_y);
						_h_manager_instance.h_shwr_direction_zy->Fill(dir_z, dir_y);
						_h_manager_instance.h_shwr_theta_phi->Fill(shwr_theta, shwr_phi);
						shwr_dir_vector.push_back(dir_x);
						shwr_dir_vector.push_back(dir_y);
						shwr_dir_vector.push_back(dir_z);
						shwr_dir_list.push_back(shwr_dir_vector);
						if(!shwr_dir_vector.empty()) {shwr_dir_vector.clear(); }

						const double open_angle = shower.at(0)->OpenAngle() * (180 / 3.1415);
						_h_manager_instance.h_shwr_open_angle->Fill(open_angle);
						const double shwr_length = shower.at(0)->Length();
						_h_manager_instance.h_shwr_length->Fill(shwr_length);
						const double end_width = _utility_instance.calc_end_width(shwr_length, open_angle);
						_h_manager_instance.h_shwr_end_width->Fill(end_width);

					}//end shwr daughters
					 //start trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						std::vector<recob::Track const*> track;
						track_for_pfp.get(i, track);

						auto const this_track = track.at(0);

						//track vertices
						_h_manager_instance.h_nue_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						_h_manager_instance.h_nue_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);

						//let's construct the path of the tracks
						std::vector<geoalgo::Point_t> track_path;
						const int track_points = this_track->NPoints();
						for(int pts = 0; pts < track_points; pts++)
						{
							geoalgo::Point_t const track_point (
							        this_track->LocationAtPoint(pts).X(),
							        this_track->LocationAtPoint(pts).Y(),
							        this_track->LocationAtPoint(pts).Z());
							track_path.push_back(track_point);
						}
						const geoalgo::Trajectory_t trj = track_path;
						if(!track_path.empty()) {track_path.clear(); }
						track_trajectory_list.push_back(trj);

						//let's get the track length!
						const double track_length = this_track->Length();
						_h_manager_instance.h_nue_trk_length->Fill(track_length);

						//let's look at the track directions
						//const double dir_x = track.at(0)->VertexDirection().X();
						//const double dir_y = track.at(0)->VertexDirection().Y();
						//const double dir_z = track.at(0)->VertexDirection().Z();

						//let's get the momentum at the start of the track
						double trk_start_m = this_track->VertexMomentum();
						//std::cout << this_track->StartMomentumVector().X()+this_track->StartMomentumVector().Y()+this_track->StartMomentumVector().Z() << std::endl;
						//std::cout << trk_start_m << std::endl;
						trk_energy_list.push_back(trk_start_m);

					}//end nue track daughters
				}//end nue daughters
				_h_manager_instance.h_nue_like_trk_daughters->Fill(trk_daughters);
				_h_manager_instance.h_nue_like_daughters->Fill(shwr_daughters, trk_daughters);
				_h_manager_instance.h_nue_like_daughters_logz->Fill(shwr_daughters, trk_daughters);

				//sum the energy of the nue event
				double shwr_E = 0;
				double trk_E  = 0;
				for(auto const shwr_energy : shwr_energy_list) {shwr_E += shwr_energy; }
				for(auto const trk_energy  : trk_energy_list ) {trk_E += trk_energy; }
				const double total_nue_momentum = shwr_E + trk_E;

			}//end nues
			 //numus!
			if(pfparts.PdgCode() == 14)
			{
				num_numu++;
				_h_manager_instance.h_numu_like_vtx_xy->Fill(xyz[0], xyz[1]);
				_h_manager_instance.h_numu_like_vtx_yz->Fill(xyz[2], xyz[1]);


				for(std::size_t const i : pfparts.Daughters())
				{
					auto const daughter = pfparticles.at(i);
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					//get vertex vector
					double d_xyz [3];
					d_vertex.at(0)->XYZ(d_xyz);

					//shwr daughters
					if(daughter.PdgCode() == 11)
					{
						shwr_daughters++;
						_h_manager_instance.h_numu_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						_h_manager_instance.h_numu_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}
					//trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						_h_manager_instance.h_numu_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						_h_manager_instance.h_numu_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}

				}
				_h_manager_instance.h_numu_like_daughters->Fill(shwr_daughters, trk_daughters);

			}//end if numu-like
		}//end if nu-like
		num_pfps++;

	}//end loop pfps

	_h_manager_instance.h_num_nue_per_event->Fill(num_nue_per_event);

	//**************************
	//loop over pandora cosmics
	//*************************
	for(auto const cosmic : cosmicpfps)
	{
		if(cosmic.PdgCode() == 13)
		{
			//let's get the cosmic to track associations
			std::vector<recob::Track const*> cosmic_track;
			cosmic_track_for_pfp.get(num_cosmics, cosmic_track);
			if(cosmic_track.size() == 0)
			{
				if(_verbose) {std::cout << "No track for pfp!" << std::endl; }
				continue;
			}
			std::vector<geoalgo::Point_t> cosmic_track_path;
			int pts = 0;
			for(auto this_point : cosmic_track)
			{
				geoalgo::Point_t const cosmic_track_point (
				        this_point->LocationAtPoint(pts).X(),
				        this_point->LocationAtPoint(pts).Y(),
				        this_point->LocationAtPoint(pts).Z());
				cosmic_track_path.push_back(cosmic_track_point);
				pts++;
			}
			const geoalgo::Trajectory_t trj = cosmic_track_path;
			if(!cosmic_track_path.empty()) {cosmic_track_path.clear(); }
			cosmic_track_trajectory_list.push_back(trj);

			//let's get the track length
			const double cosmic_length = cosmic_track.at(0)->Length();
			_h_manager_instance.h_cosmic_trk_length->Fill(cosmic_length);
			cosmic_track_length_list.push_back(cosmic_length);

			//let's get the track energy
			const double cosmic_trk_energy = cosmic_track.at(0)->StartMomentum();

		} //end loop tracks
		if(cosmic.PdgCode() == 11)
		{
			//let's get the cosmic to shower associations
			std::vector<recob::Shower const*> cosmic_shower;
			cosmic_shower_for_pfp.get(num_cosmics, cosmic_shower);
			if(cosmic_shower.size() == 0)
			{
				if(_verbose) {std::cout << "No shower for pfp!" << std::endl; }
				continue;
			}

			//let's get the cosmic shower energy
			const std::vector < double > plane_energy = cosmic_shower.at(0)->Energy();
			const int bestplane = cosmic_shower.at(0)->best_plane();
			const double cosmic_shwr_energy = plane_energy.at(bestplane);
		} //end loop showers
		num_cosmics++;
	} //end loop cosmics

	//*********************************************************
	//******************* Geometry studies! *******************
	//*********************************************************
	std::vector < geoalgo::Point_t > cut_nue_vertex;
	std::vector < geoalgo::Point_t > cut_nue_shwr_vertex;

	//closest point between **nue vertex** and cosmic track
	for(geoalgo::Point_t const nue_vertex : nue_vertex_list)
	{
		double cosmic_closest_point = 0;
		if(!cosmic_track_trajectory_list.empty())
		{
			cosmic_closest_point = _geo_algo_instance.SqDist(nue_vertex, cosmic_track_trajectory_list);
			_h_manager_instance.h_nue_cosmic_closest->Fill(cosmic_closest_point);
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
			double closest_point = _geo_algo_instance.SqDist(nue_vertex, track_trajectory_list);
			_h_manager_instance.h_nue_trk_closest->Fill(closest_point);
			if(cosmic_closest_point >= cut_distance_to_point) {_h_manager_instance.h_nue_trk_closest_zoom->Fill(closest_point); }
		}
	}
	//closest point between **nue shwr vertex** and cosmic track
	bool first = true;
	int nE = 0;
	for (geoalgo::Point_t shwr_vertex : shwr_vertex_list)
	{
		if(!cosmic_track_trajectory_list.empty())
		{
			const double closest_point = _geo_algo_instance.SqDist(shwr_vertex, cosmic_track_trajectory_list);
			_h_manager_instance.h_nue_shwr_cosmic_closest->Fill(closest_point);
			_h_manager_instance.h_nue_shwr_cosmic_closest_vs_E->Fill(closest_point, shwr_energy_list.at(nE));
			_h_manager_instance.h_nue_shwr_cosmic_closest_vs_y->Fill(closest_point, shwr_vertex[1]);
			_h_manager_instance.h_nue_shwr_cosmic_closest_vs_E_zoom->Fill(closest_point, shwr_energy_list.at(nE));
			_h_manager_instance.h_nue_shwr_cosmic_closest_vs_y_zoom->Fill(closest_point, shwr_vertex[1]);
			_h_manager_instance.h_shwr_direction_y_vs_nearest_cosmic->Fill(closest_point, shwr_dir_list.at(nE).at(1));

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
					_h_manager_instance.h_nue_like_daughters_cuts->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
					_h_manager_instance.h_nue_like_daughters_cuts_logz->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
				}
				_h_manager_instance.h_shwr_direction_cut_xy->Fill(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(1));
				_h_manager_instance.h_shwr_direction_cut_zy->Fill(shwr_dir_list.at(nE).at(2), shwr_dir_list.at(nE).at(1));
				const double shwr_cut_theta = TMath::ASin(shwr_dir_list.at(nE).at(1)) * (180/3.1415);
				const double shwr_cut_phi   = TMath::ATan2(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(2)) * (180/3.1415);
				_h_manager_instance.h_shwr_cut_theta_phi->Fill(shwr_cut_theta, shwr_cut_phi);
			}
			//let's loop over all of the nue shwrs and the cosmic tracks
			//this will take the length to compute a cut out "fiducial" volume
			double cylinder_vol = _utility_instance.cylinder_fid_vol(
			        cosmic_track_trajectory_list,
			        cosmic_track_length_list,
			        _geo_algo_instance,
			        shwr_vertex,
			        cut_distance_to_point
			        );
			_h_manager_instance.h_cylinder_vol->Fill(cylinder_vol/ub_total_vol);
		}
		//closest point between **nue shwr vertex** and nue track
		if(!track_trajectory_list.empty())
		{
			const double closest_point = _geo_algo_instance.SqDist(shwr_vertex, track_trajectory_list);
			_h_manager_instance.h_nue_shwr_trk_closest->Fill(closest_point);
			int num_trks_nearby = 0;
			for(auto this_track : track_trajectory_list)
			{
				const double this_closest_point = _geo_algo_instance.SqDist(shwr_vertex, this_track);
				if(this_closest_point <= 5) {num_trks_nearby++; }
			}
			_h_manager_instance.h_num_trks_nearby->Fill(num_trks_nearby);
		}//end if event has track
		nE++;
	}//end loop nue shwr vtx
	 //***************************************************************
	 //loop for showers with a large distance from the reco nue vertex
	for(geoalgo:: Point_t const shwr_vertex_lrgDist : shwr_vertex_lrgDist_list)
	{
		if(!track_trajectory_list.empty())
		{
			const double this_closest_point_lrgDist = _geo_algo_instance.SqDist(shwr_vertex_lrgDist, track_trajectory_list);
			_h_manager_instance.h_nue_like_shwr_lrgDist_dist_to_track->Fill(this_closest_point_lrgDist);

			const double this_closest_point_cosmic_lrgDist = _geo_algo_instance.SqDist(shwr_vertex_lrgDist, cosmic_track_trajectory_list);
			_h_manager_instance.h_nue_like_shwr_lrgDist_dist_to_cosmic->Fill(this_closest_point_cosmic_lrgDist);
		}
		const int num_trks_lrgDist = track_trajectory_list.size();
		_h_manager_instance.h_nue_like_shwr_lrgDist_num_trks->Fill(num_trks_lrgDist);
	}

	//cut on sqdist to remove small distances - does this change vertex ave. position
	for(auto vtx : cut_nue_shwr_vertex)
	{
		_h_manager_instance.h_nue_shwr_cut_vtx_xy->Fill(vtx[0], vtx[1]);
		_h_manager_instance.h_nue_shwr_cut_vtx_zy->Fill(vtx[2], vtx[1]);
	}

//*****************************
// *End Geometry Calculations*
//*****************************

	return true;
}



bool cosmic_ana::finalize() {

	/*********************************
	** Histogram Saving and Editing **
	*///******************************
	_h_manager_instance.draw_save();


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
