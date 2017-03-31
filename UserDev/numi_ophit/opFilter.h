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

#ifndef GALLERY_FMWK_OPFILTER_H
#define GALLERY_FMWK_OPFILTER_H

//some root includes
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "Analysis/ana_base.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "GeoAlgo/GeoAlgo.h"
#include "GeoAlgo/GeoVector.h"
#include "GeoAlgo/GeoTrajectory.h"

namespace galleryfmwk {

/**
   \class example_ana
   User custom analysis class made by SHELL_USER_NAME
 */
class opFilter : galleryfmwk::ana_base {

const raw::ubdaqSoftwareTriggerData * ubTrigData;

public:

/// Default constructor
opFilter() {
	_verbose = false;
}

/// Default destructor
// ~example_ana() {}


bool initialize();


bool analyze(gallery::Event * ev);


bool finalize();

void setFlashProducer(std::string s) {
	_flash_producer = s;
}
void setVerbose(bool b){
	_verbose = b;
}
void setSoftWareProducer(std::string s) {
	_sftwr_producer = s;
}
void setTrigAlgProducer(std::string s) {
	_trigAlg_producer = s;	
}

protected:

TFile * datafile = new TFile("ext_numi_ana.root", "RECREATE", "timing");
TTree * datatree = new TTree("datatree", "data tree");
double numi_time;

std::string _flash_producer;
std::string _sftwr_producer;
std::string _trigAlg_producer;


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

int sftwTrig_counter;
int noSftwTrig_counter;

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
