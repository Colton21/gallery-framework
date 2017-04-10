#ifndef GALLERY_FMWK_INFV_FILTER_CXX
#define GALLERY_FMWK_INFV_FILTER_CXX

#include "inFV_filter.h"

//also add containment of the showers and tracks!
//test that this works for nues
//also once tested, can remove excess comment lines in cosmic_ana

namespace galleryfmwk {

bool inFV_filter::initialize() {

	_h_manager_instance.gen_histograms_fv(gDirectory);

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


bool inFV_filter::analyze(gallery::Event * ev) {

	num_cosmic++;

	// For each file, loop over all events.

	bool _debug = false;


	// Get all of the tracks from the event:
	art::InputTag pfp_tag(_pfp_tag);
	art::InputTag pfp_cosmic_tag(_pfp_cosmic_tag);
	//double cut_distance_to_point(_cut);

	auto const & pfp
	        = ev->getValidHandle<std::vector <recob::PFParticle> > (pfp_tag);
	auto const & pfparticles(*pfp);

	art::FindMany<recob::Vertex> vertex_for_pfp(pfp, *ev, pfp_tag);
	art::FindMany<recob::Track> track_for_pfp(pfp, *ev, "pandoraNu");
	art::FindMany<recob::Shower> shower_for_pfp(pfp, *ev, "pandoraNu");

	if(_verbose) {std::cout << "FV Filter Event" << std::endl; }

	num_pfps = 0;

	//pfp loop
	for(auto pfparts : pfparticles)
	{

		if(_debug) {std::cout << "a" << std::endl; }

		//******************************
		//check for reconstructed vertex
		//******************************
		std::vector<recob::Vertex const*> vertex;
		vertex_for_pfp.get(num_pfps,vertex);

		if(vertex.size() == 0 )
		{
			if(_verbose == true) {std::cout << "No vertex association found!" << std::endl; }
			return false;
		}

		//get vertex vector
		double xyz [3];
		vertex.at(0)->XYZ(xyz);
		if(_utility_instance.inFV(xyz[0], xyz[1], xyz[2], _right, _left, _up, _down, _back, _front) == false)
		{
			if(_verbose) {std::cout << "Reco vertex outside fiducial volume!" << std::endl; }
			return false;
		}

		if(_debug) {std::cout << "b" << std::endl; }

		//************************************
		//check if pfp is neutrino-like object
		//************************************
		int shwr_daughters = 0;
		int trk_daughters = 0;
		if(pfparts.IsPrimary() == true)
		{
			num_primary_pfp++;
			//nues!

			if(_debug) {std::cout << "c" << std::endl; }

			if(pfparts.PdgCode() == 12)
			{
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
					//check if daughter vertex is inFV
					if(_utility_instance.inFV(d_xyz[0], d_xyz[1], d_xyz[2], _right, _left, _up, _down, _back, _front) == false) {return false; }

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
							return false;
						}

						// const double open_angle = shower.at(0)->OpenAngle() * (180 / 3.1415);
						// _h_manager_instance.h_shwr_open_angle->Fill(open_angle);
						// const double shwr_length = shower.at(0)->Length();
						// _h_manager_instance.h_shwr_length->Fill(shwr_length);
						// const double end_width = _utility_instance.calc_end_width(shwr_length, open_angle);
						// _h_manager_instance.h_shwr_end_width->Fill(end_width);

					}//end shwr daughters
					 //start trk daughters
					if(daughter.PdgCode() == 13)
					{

						if(_debug) {std::cout << "d" << std::endl; }

						trk_daughters++;
						std::vector<recob::Track const*> track;
						track_for_pfp.get(i, track);
						if(track.size() == 0)
						{
							if(_verbose) {std::cout << "No track for this pfp track!" << std::endl; }
							return false;
						}

						//auto const this_track = track.at(0);

						//let's construct the path of the tracks

						//let's get the track length!
						//const double track_length = this_track->Length();

						//let's look at the track directions
						//const double dir_x = this_track->VertexDirection().X();
						//const double dir_y = this_track->VertexDirection().Y();
						//const double dir_z = this_track->VertexDirection().Z();

					}//end nue track daughters
				}//end nue daughters

				if(_debug) {std::cout << "e" << std::endl; }

				for(int fv_cut = 0; fv_cut < fv_cut_max; fv_cut++)
				{
					if(_utility_instance.inFV(xyz[0], xyz[1], xyz[2], fv_cut, fv_cut, fv_cut, fv_cut, fv_cut, fv_cut) == true)
					{_h_manager_instance.h_nue_fv_cuts->Fill(fv_cut); }
					//just fv cut from top
					if(_utility_instance.inFV(xyz[0], xyz[1], xyz[2], 0, 0, fv_cut, 0, 0, 0) == true)
					{_h_manager_instance.h_nue_fv_top_cuts->Fill(fv_cut); }
				}
			}//end nues

			if(_debug) {std::cout << "f" << std::endl; }

			//numus!
			if(pfparts.PdgCode() == 14)
			{
				num_numu++;

				if(_debug) {std::cout << "g" << std::endl; }

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

					if(_utility_instance.inFV(d_xyz[0], d_xyz[1], d_xyz[2], _right, _left, _up, _down, _back, _front) == false) {return false; }
				}
			}//end if numu-like
		}//end if nu-like
		num_pfps++;

		if(_debug) {std::cout << "h" << std::endl; }

	}//end loop pfps

	if(_debug) {std::cout << "i" << std::endl; }

	return true;
}



bool inFV_filter::finalize() {

	std::cout << "Finished inFV Filter" << std::endl;

	_h_manager_instance.draw_save_fv();

	return true;
}


}
#endif
