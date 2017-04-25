/**
 * \file example_ana.h
 *
 * \ingroup nuexsec_analysis
 *
 * \brief Class def header for a class example_ana
 *
 * @author cadams
 */

/** \addtogroup nuexsec_analysis

    @{*/

#ifndef GALLERY_FMWK_COSMIC_FILTER_H
#define GALLERY_FMWK_COSMIC_FILTER_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "Analysis/ana_base.h"
#include "histo_manager.h"
#include "utility_functions.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"


#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"

#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoTrajectory.h"

namespace galleryfmwk {

/**
   \class example_ana
   User custom analysis class made by SHELL_USER_NAME
 */
class cosmic_filter : galleryfmwk::ana_base {

//class instances
geoalgo::GeoAlgo const _geo_algo_instance;
h_manager _h_manager_instance;
utility _utility_instance;

public:

/// Default constructor
cosmic_filter() {
	_verbose = false;
}

/// Default destructor
// ~example_ana() {}


bool initialize();


bool analyze(gallery::Event * ev);


bool finalize();

void setWantCC(bool b) {
	wantCC = b;
}

void setTrackProducer(std::string s) {
	_track_producer = s;
}
void setShowerProducer(std::string s) {
	_shower_producer = s;
}
void setVerbose(bool b){
	_verbose = b;
}
void setNearestCutDist(double cut){
	_cut = cut;
}
void fiducial_volume_x_right(double right){
	_right = right;
}
void fiducial_volume_x_left(double left){
	_left = left;
}
void fiducial_volume_y_up(double up){
	_up = up;
}
void fiducial_volume_y_down(double down){
	_down = down;
}
void fiducial_volume_z_back(double back){
	_back = back;
}
void fiducial_volume_z_front(double front){
	_front = front;
}
void setMCProducer(std::string s){
	_mc_part_tag = s;
}
void setPfpCosmicProducer(std::string s){
	_pfp_cosmic_tag = s;
}

protected:

bool _verbose;
std::string _track_producer;
std::string _shower_producer;
std::string _mc_part_tag;
std::string _pfp_cosmic_tag;
double _cut;
double _right;
double _left;
double _up;
double _down;
double _back;
double _front;
bool wantCC;

double x_boundary1;
double x_boundary2;
double y_boundary1;
double y_boundary2;
double z_boundary1;
double z_boundary2;

int num_events;
int num_events_remaining;

TH1D * h_nu_vtx_x;
TH1D * h_nu_vtx_y;
TH1D * h_nu_vtx_z;

TH2D * h_nu_vtx_xy_outside;
TH2D * h_nu_vtx_zy_outside;
TH2D * h_nu_vtx_xy_inside;
TH2D * h_nu_vtx_zy_inside;

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
