#ifndef GALLERY_FMWK_COSMIC_ANA_CXX
#define GALLERY_FMWK_COSMIC_ANA_CXX

#include "cosmic_ana.h"
#include "histo_manager.h"

namespace galleryfmwk {

bool cosmic_ana::inFV(double x_vtx, double y_vtx, double z_vtx,
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

//calculates lost fiducial volume as a result of the closest cosmic track to a shower vertex
double cosmic_ana::cylinder_fid_vol(
        std::vector < geoalgo::Trajectory_t > cosmic_track_trajectory_list,
        std::vector < double > cosmic_track_length_list,
        geoalgo::GeoAlgo const _geo_algo_instance,
        geoalgo::Point_t shwr_vertex,
        double cut_distance_to_point
        )
{
	double closest_track = 100;
	int track_counter = 0;
	int this_track_counter = 0;
	for( auto cosmic_track_traj : cosmic_track_trajectory_list)
	{
		//which track is closest?
		track_counter++;
		const double close_track = _geo_algo_instance.SqDist(shwr_vertex, cosmic_track_traj);
		if(close_track < closest_track)
		{
			closest_track = close_track;
			this_track_counter = track_counter;
		}
	}
	const double closest_track_length = cosmic_track_length_list.at(this_track_counter);
	//calculate the volume of the cylinder
	return closest_track_length * 3.1415 * cut_distance_to_point * cut_distance_to_point;
}

double cosmic_ana::geo_distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2) const
{
	const double dist_x = x2 - x1;
	const double dist_y = y2 - y1;
	const double dist_z = z2 - z1;
	const double dist =
	        sqrt((dist_x * dist_x)+
	             (dist_y * dist_y)+
	             (dist_z * dist_z));
	return dist;
}

double shower_ana::calc_end_width(const double length, const double open_angle) const
{
	const double width = 2 * length * TMath::Tan(open_angle / 2);
	return width;
}


bool cosmic_ana::initialize() {

	//
	// This function is called in the beginning of event loop
	// Do all variable initialization you wish to do here.
	// If you have a histogram to fill in the event loop, for example_ana,
	// here is a good place to create one on the heap (i.e. "new TH1D").
	//

	histogram::gen_histograms();

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
	// art::InputTag tracks_tag(_track_producer);
	// art::InputTag showers_tag(_shower_producer);
	art::InputTag pfp_tag(_pfp_tag);
	art::InputTag pfp_cosmic_tag(_pfp_cosmic_tag);
	double cut_distance_to_point(_cut);


	// auto const & tracks
	//         = ev->getValidHandle<std::vector <recob::Track> >(tracks_tag);
	// auto const & showers
	//         = ev->getValidHandle<std::vector <recob::Shower> >(showers_tag);
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
	std::vector < geoalgo::Point_t> shwr_vertex_lrgDist_list;

	std::vector < double > shwr_dir_vector;
	std::vector < std::vector < double > > shwr_dir_list;

	int num_pfps = 0;
	int num_cosmics = 0;

	//vector < recob::Track> cosmic_track_obj_list;

	num_nue_per_event = 0;
	//pfp loop
	for(auto pfparts : pfparticles)
	{

		//******************************
		//check for reconstructed vertex
		//******************************
		std::vector<recob::Vertex const*> vertex;
		vertex_for_pfp.get(num_pfps,vertex);
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
			if(_verbose) {std::cout << "Reco vertex outside fiducial volume!" << std::endl; }
			continue;
		}

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

				histogram::h_nue_like_vtx_xy->Fill(xyz[0], xyz[1]);
				histogram::h_nue_like_vtx_yz->Fill(xyz[2], xyz[1]);

				for(std::size_t const i : pfparts.Daughters())
				{
					auto const daughter = pfparticles.at(i);
					//let's get the vertex associations for the daughters
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					//if no daughter vertex
					if(d_vertex.size() == 0 )
					{
						if(_verbose) {std::cout << "No vertex association found for daughter!" << std::endl; }
						return false;
					}
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
						if(shower.size() == 0)
						{
							if(_verbose) {std::cout << "No shower for this pfp shower!" << std::endl; }
							continue;
						}

						//shower daughter vertices
						histogram::h_nue_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						histogram::h_nue_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
						geoalgo::Point_t const shwr_vtx (d_xyz[0], d_xyz[1], d_xyz[2]);
						shwr_vertex_list.push_back(shwr_vtx);

						//let's check the distance between the shwr vtx and the nue vtx
						const double dist = cosmic_ana::geo_distance(d_xyz[0], xyz[0], d_xyz[1], xyz[1], d_xyz[2], xyz[2]);
						histogram::h_nue_shwr_vtx_dist->Fill(dist);

						//what does it mean when the distance between the nue and shwr vtx
						//is large?
						if(dist >= 5 )
						{
							histogram::h_nue_like_shwr_lrgDist_vtx_xy->Fill(d_xyz[0], d_xyz[1]);
							histogram::h_nue_like_shwr_lrgDist_vtx_zy->Fill(d_xyz[2], d_xyz[1]);
							geoalgo::Point_t const shwr_vtx_lrgDist (d_xyz[0], d_xyz[1], d_xyz[2]);
							shwr_vertex_lrgDist_list.push_back(shwr_vtx_lrgDist);
						}

						//let's get the energy! Energy() - GeV?
						const std::vector < double > plane_energy = shower.at(0)->Energy();
						const int best_plane = shower.at(0)->best_plane();
						const double total_energy = plane_energy.at(best_plane);
						histogram::h_nue_shwr_E->Fill(total_energy);
						shwr_energy_list.push_back(total_energy);

						//let's look at the shower directions
						const double dir_x = shower.at(0)->Direction().X();
						const double dir_y = shower.at(0)->Direction().Y();
						const double dir_z = shower.at(0)->Direction().Z();
						const double shwr_theta = TMath::ASin(dir_y) * (180/3.1415);
						const double shwr_phi = TMath::ATan2(dir_x, dir_z) * (180/3.1415);
						histogram::h_shwr_direction_xy->Fill(dir_x, dir_y);
						histogram::h_shwr_direction_zy->Fill(dir_z, dir_y);
						histogram::h_shwr_theta_phi->Fill(shwr_theta, shwr_phi);
						shwr_dir_vector.push_back(dir_x);
						shwr_dir_vector.push_back(dir_y);
						shwr_dir_vector.push_back(dir_z);
						shwr_dir_list.push_back(shwr_dir_vector);
						if(!shwr_dir_vector.empty()) {shwr_dir_vector.clear(); }

						const double open_angle = shower.at(0)->OpenAngle();
						histogram::h_shwr_open_angle->Fill(open_angle);
						const double shwr_length = shower.at(0)->Length();
						histogram::h_shwr_length->Fill(shwr_length);
						const double end_width = shower_ana::calc_end_width(shwr_length, open_angle);
						histogram::h_shwr_end_width->Fill(end_width);

					}//end shwr daughters
					 //start trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						std::vector<recob::Track const*> track;
						track_for_pfp.get(i, track);
						if(track.size() == 0)
						{
							if(_verbose) {std::cout << "No track for this pfp track!" << std::endl; }
							continue;
						}

						//track vertices
						histogram::h_nue_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						histogram::h_nue_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);

						//let's construct the path of the tracks
						std::vector<geoalgo::Point_t> track_path;
						const int track_points = track.at(0)->NPoints();
						for(int pts = 0; pts < track_points; pts++)
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
						histogram::h_nue_trk_length->Fill(track_length);

						//let's look at the track directions
						//const double dir_x = track.at(0)->VertexDirection().X();
						//const double dir_y = track.at(0)->VertexDirection().Y();
						//const double dir_z = track.at(0)->VertexDirection().Z();

					}//end nue track daughters
				}//end nue daughters
				histogram::h_nue_like_trk_daughters->Fill(trk_daughters);
				histogram::h_nue_like_daughters->Fill(shwr_daughters, trk_daughters);
				histogram::h_nue_like_daughters_logz->Fill(shwr_daughters, trk_daughters);

				for(int fv_cut = 0; fv_cut < fv_cut_max; fv_cut++)
				{
					if(inFV(xyz[0], xyz[1], xyz[2], fv_cut, fv_cut, fv_cut, fv_cut, fv_cut, fv_cut) == true)
					{histogram::h_nue_fv_cuts->Fill(fv_cut); }
					//just fv cut from top
					if(inFV(xyz[0], xyz[1], xyz[2], 0, 0, fv_cut, 0, 0, 0) == true)
					{histogram::h_nue_fv_top_cuts->Fill(fv_cut); }
				}
			}
			//numus!
			if(pfparts.PdgCode() == 14)
			{
				num_numu++;
				histogram::h_numu_like_vtx_xy->Fill(xyz[0], xyz[1]);
				histogram::h_numu_like_vtx_yz->Fill(xyz[2], xyz[1]);

				for(std::size_t const i : pfparts.Daughters())
				{
					auto const daughter = pfparticles.at(i);
					std::vector<recob::Vertex const*> d_vertex;
					vertex_for_pfp.get(i, d_vertex);
					if(d_vertex.size() ==0 )
					{
						if(_verbose) {std::cout << "No vertex association found for daughter!" << std::endl; }
						return false;
					}
					//get vertex vector
					double d_xyz [3];
					d_vertex.at(0)->XYZ(d_xyz);

					//shwr daughters
					if(daughter.PdgCode() == 11)
					{
						shwr_daughters++;
						histogram::h_numu_like_shwr_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						histogram::h_numu_like_shwr_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}
					//trk daughters
					if(daughter.PdgCode() == 13)
					{
						trk_daughters++;
						histogram::h_numu_like_trk_daughters_xy->Fill(d_xyz[0], d_xyz[1]);
						histogram::h_numu_like_trk_daughters_yz->Fill(d_xyz[2], d_xyz[1]);
					}

				}
				histogram::h_numu_like_daughters->Fill(shwr_daughters, trk_daughters);
			}//end if numu-like
		}//end if nu-like
		num_pfps++;
	}//end loop pfps

	histogram::h_num_nue_per_event->Fill(num_nue_per_event);

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
			histogram::h_cosmic_trk_length->Fill(cosmic_length);
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
		}        //end loop showers
		num_cosmics++;
	}        //end loop cosmics

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
			histogram::h_nue_cosmic_closest->Fill(cosmic_closest_point);
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
			histogram::h_nue_trk_closest->Fill(closest_point);
			if(cosmic_closest_point >= cut_distance_to_point) {histogram::h_nue_trk_closest_zoom->Fill(closest_point); }
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
			histogram::h_nue_shwr_cosmic_closest->Fill(closest_point);
			histogram::h_nue_shwr_cosmic_closest_vs_E->Fill(closest_point, shwr_energy_list.at(nE));
			histogram::h_nue_shwr_cosmic_closest_vs_y->Fill(closest_point, shwr_vertex[1]);
			histogram::h_nue_shwr_cosmic_closest_vs_E_zoom->Fill(closest_point, shwr_energy_list.at(nE));
			histogram::h_nue_shwr_cosmic_closest_vs_y_zoom->Fill(closest_point, shwr_vertex[1]);
			histogram::h_shwr_direction_y_vs_nearest_cosmic->Fill(closest_point, shwr_dir_list.at(nE).at(1));

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
					histogram::h_nue_like_daughters_cuts->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
					histogram::h_nue_like_daughters_cuts_logz->Fill(shwr_vertex_list.size(), track_trajectory_list.size());
				}
				histogram::h_shwr_direction_cut_xy->Fill(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(1));
				histogram::h_shwr_direction_cut_zy->Fill(shwr_dir_list.at(nE).at(2), shwr_dir_list.at(nE).at(1));
				const double shwr_cut_theta = TMath::ASin(shwr_dir_list.at(nE).at(1)) * (180/3.1415);
				const double shwr_cut_phi   = TMath::ATan2(shwr_dir_list.at(nE).at(0), shwr_dir_list.at(nE).at(2)) * (180/3.1415);
				histogram::h_shwr_cut_theta_phi->Fill(shwr_cut_theta, shwr_cut_phi);
			}
			//let's loop over all of the nue shwrs and the cosmic tracks
			//this will take the length to compute a cut out "fiducial" volume
			double cylinder_vol = cosmic_ana::cylinder_fid_vol(
			        cosmic_track_trajectory_list,
			        cosmic_track_length_list,
			        _geo_algo_instance,
			        shwr_vertex,
			        cut_distance_to_point
			        );
			histogram::h_cylinder_vol->Fill(cylinder_vol/ub_total_vol);
		}
		//closest point between **nue shwr vertex** and nue track
		if(!track_trajectory_list.empty())
		{
			const double closest_point = _geo_algo_instance.SqDist(shwr_vertex, track_trajectory_list);
			histogram::h_nue_shwr_trk_closest->Fill(closest_point);
			int num_trks_nearby = 0;
			for(auto this_track : track_trajectory_list)
			{
				const double this_closest_point = _geo_algo_instance.SqDist(shwr_vertex, this_track);
				if(this_closest_point <= 5) {num_trks_nearby++; }
			}
			histogram::h_num_trks_nearby->Fill(num_trks_nearby);
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
			histogram::h_nue_like_shwr_lrgDist_dist_to_track->Fill(this_closest_point_lrgDist);

			const double this_closest_point_cosmic_lrgDist = _geo_algo_instance.SqDist(shwr_vertex_lrgDist, cosmic_track_trajectory_list);
			histogram::h_nue_like_shwr_lrgDist_dist_to_cosmic->Fill(this_closest_point_cosmic_lrgDist);
		}
		const int num_trks_lrgDist = track_trajectory_list.size();
		histogram::h_nue_like_shwr_lrgDist_num_trks->Fill(num_trks_lrgDist);
	}

	//cut on sqdist to remove small distances - does this change vertex ave. position
	for(auto vtx : cut_nue_shwr_vertex)
	{
		histogram::h_nue_shwr_cut_vtx_xy->Fill(vtx[0], vtx[1]);
		histogram::h_nue_shwr_cut_vtx_zy->Fill(vtx[2], vtx[1]);
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
	histogram::draw_save();


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
