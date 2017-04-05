#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H


#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoTrajectory.h"

class utility
{

public:

bool inFV(double x_vtx, double y_vtx, double z_vtx,
          double x1, double x2, double y1, double y2, double z1, double z2);
double cylinder_fid_vol(
        std::vector < geoalgo::Trajectory_t > cosmic_track_trajectory_list,
        std::vector < double > cosmic_track_length_list,
        geoalgo::GeoAlgo const _geo_algo_instance,
        geoalgo::Point_t shwr_vertex,
        double cut_distance_to_point
        );
double geo_distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2) const;
double calc_end_width(const double length, const double open_angle) const;

};


#endif
