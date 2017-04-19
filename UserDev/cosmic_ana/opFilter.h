/**
 * \file opFilter.h
 *
 * \ingroup nuexsec_analysis
 *
 * \brief Class def header for a class opFilter
 *
 * @author cadams
 */

/** \addtogroup nuexsec_analysis

    @{*/

#ifndef GALLERY_FMWK_OPTFILTER_H
#define GALLERY_FMWK_OPTFILTER_H

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "Analysis/ana_base.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoTrajectory.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"

namespace galleryfmwk {

/**
   \class example_ana
   User custom analysis class made by SHELL_USER_NAME
 */
class optFilter : galleryfmwk::ana_base {

public:

/// Default constructor
optFilter() {
	_verbose = false;
}

/// Default destructor
// ~example_ana() {}


bool initialize();


bool analyze(gallery::Event * ev);


bool finalize();

void setPfpProducer(std::string s) {
	_pfp_tag = s;
}

void setTrackProducer(std::string s) {
	_track_producer = s;
}
void setShowerProducer(std::string s) {
	_shower_producer = s;
}
void setFlashProducer(std::string s) {
	_flash_producer = s;
}
void setPEThreshold(double threshold){
	_pe_threshold = threshold;
}
void togglePlotting(bool yn){
	_threshold_plotting = yn;
}

void setMinTime(double d){
	_min_time = d;
}
void setMaxTime(double d){
	_max_time = d;
}

void setVerbose(bool b){
	_verbose = b;
}

protected:

bool _threshold_plotting;
std::string _pfp_tag;
std::string _track_producer;
std::string _shower_producer;
std::string _flash_producer;
double _pe_threshold;
double _min_time;
double _max_time;

double x_boundary1;
double x_boundary2;
double y_boundary1;
double y_boundary2;
double z_boundary1;
double z_boundary2;
double fromWall;

double fv_cut_max;

bool _verbose;

int flash_pass_counter;
int total_flash_counter;

TH2D * h_flash_zywidth;
TH2D * h_flash_zycenter;
TH1D * h_flash_threshold;

TH1D * h_largest_flash;
TH2D * h_largest_flash_y;
TH1D * h_nue_shwr_vtx_flash_dist_zy;
TH2D * h_nue_shwr_vtx_flash_width_zy;

};

}

#endif

//**************************************************************************
//
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group
