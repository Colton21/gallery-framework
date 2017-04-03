#ifndef LARLITE_OPFILTER_CXX
#define LARLITE_OPFILTER_CXX

#include "opFilter.h"

namespace galleryfmwk {

bool optFilter::initialize() {

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

	//fiducial cut beyond TPC
	fromWall = 0;

	flash_pass_counter = 0;
	total_flash_counter = 0;

	return true;
}


bool optFilter::analyze(gallery::Event * ev) {

	//num_cosmic++;

	// For each file, loop over all events.
	//
	// Determine criteria for rejecting muon like events
	// Save metadata for each event (neutrino pdg, energy, vertex)
	// as well as filter results.


	// Get all of the tracks from the event:
	art::InputTag tracks_tag(_track_producer);
	art::InputTag showers_tag(_shower_producer);
	art::InputTag flash_tag(_flash_producer);

	auto const & opf = ev->getValidHandle<std::vector < recob::OpFlash> >(flash_tag);
	auto const & opflashes(*opf);

	const int num_flashes = opflashes.size();
	total_flash_counter++;

	for(int this_flash = 0; this_flash < num_flashes; this_flash++)
	{
		auto const opflsh = opflashes.at(this_flash);
		if(opflsh.Time() >= 3 && opflsh.Time() <= 5)
		{
			if(opflsh.TotalPE() >= 50)
			{
				flash_pass_counter++;
				return true;
			}
		}
	}

	//default filters
	return false;
}

//for timing: marco does 3-5 us window, and total PE over all opdets >= 50 PE

bool optFilter::finalize() {

	// If you need, you can store your ROOT class instance in the output
	// file. You have an access to the output file through "_fout" pointer.
	//
	// Say you made a histogram pointer h1 to store. You can do this:
	//
	// if (_fout) { _fout->cd(); _tree->Write(); }

	std::cout << "Total Events: " << total_flash_counter << std::endl;
	std::cout << "Remaining Events: " << flash_pass_counter << std::endl;

	return true;
}


}
#endif
