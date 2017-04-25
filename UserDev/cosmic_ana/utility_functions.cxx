#include "utility_functions.h"


bool utility::inFV(double x_vtx, double y_vtx, double z_vtx,
                   double x1, double x2, double y1, double y2, double z1, double z2)
{
	//is vertex in given FV?

	if(x_vtx < x2 ||
	   x_vtx > x1 ||
	   y_vtx < y2 ||
	   y_vtx > y1 ||
	   z_vtx < z1 ||
	   z_vtx > z2)
	{
		return false;
	}
	return true;
}

//calculates lost fiducial volume as a result of the closest cosmic track to a shower vertex
double utility::cylinder_fid_vol(
        std::vector < geoalgo::Trajectory_t > cosmic_track_trajectory_list,
        std::vector < double > cosmic_track_length_list,
        geoalgo::GeoAlgo const _geo_algo_instance,
        geoalgo::Point_t shwr_vertex,
        double cut_distance_to_point
        )
{
	bool _verbose = false;
	double closest_track = 100;
	int track_counter = 0;
	int this_track_counter = 0;
	for( auto cosmic_track_traj : cosmic_track_trajectory_list)
	{
		//which track is closest?
		track_counter++;
		const double close_track = _geo_algo_instance.SqDist(shwr_vertex, cosmic_track_traj);
		if(close_track <= closest_track)
		{
			closest_track = close_track;
			this_track_counter = track_counter;
		}
	}
	double closest_track_length = 0;
	try
	{
		closest_track_length = cosmic_track_length_list.at(this_track_counter);
	}
	catch(...)
	{
		if(_verbose) {std::cout << "Track counter not same length as track length vector" << std::endl; }
		//return 0;
	}
	//calculate the volume of the cylinder
	const double cylinder_volume = (closest_track_length * 3.1415 * cut_distance_to_point * cut_distance_to_point);
	return cylinder_volume;
}

double utility::geo_distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2) const
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

double utility::calc_end_width(const double length, const double open_angle) const
{
	const double width = 2 * length * TMath::Tan(open_angle / 2);
	return width;
}
