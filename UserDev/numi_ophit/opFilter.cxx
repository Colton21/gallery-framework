#ifndef LARLITE_OPFILTER_CXX
#define LARLITE_OPFILTER_CXX

#include "opFilter.h"

namespace galleryfmwk {

bool opFilter::initialize() {

	//
	// This function is called in the beggining of event loop
	// Do all variable initialization you wish to do here.
	// If you have a histogram to fill in the event loop, for example_ana,
	// here is a good place to create one on the heap (i.e. "new TH1D").
	//

	//gROOT->SetBatch();

	datatree->Branch("numi_time", &numi_time, "numi_time/D");

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

	sftwTrig_counter = 0;
	noSftwTrig_counter = 0;

	return true;
}


bool opFilter::analyze(gallery::Event * ev) {

	//num_cosmic++;

	// For each file, loop over all events.
	//
	// Determine criteria for rejecting muon like events
	// Save metadata for each event (neutrino pdg, energy, vertex)
	// as well as filter results.


	// Get all of the tracks from the event:
	art::InputTag flash_tag(_flash_producer);
	art::InputTag sftwr_tag(_sftwr_producer);
	std::string trigAlg_tag(_trigAlg_producer);

	auto const & opf = ev->getValidHandle<std::vector < recob::OpFlash> >(flash_tag);
	auto const & opflashes(*opf);

	total_flash_counter++;

	auto const & sftwrTrig = ev->getValidHandle<raw::ubdaqSoftwareTriggerData> (sftwr_tag);
	ubTrigData = sftwrTrig.product();

	//check if event passes software trigger
	sftwTrig_counter ++;
	
	
        const int numTriggers = ubTrigData->getNumberOfAlgorithms();
	std::vector < std::string> listOfTriggers = ubTrigData->getListOfAlgorithms();
	if(_verbose){std::cout << "Number of Triggers: " << numTriggers << std::endl;}

	for(auto triggers : listOfTriggers)
	{
	  if(_verbose){std::cout << triggers << std::endl;}
	}
	
	if(ubTrigData->passedAlgo(trigAlg_tag) == false)
	{
	    if(_verbose == true){std::cout << "Incorrect Algorithm: " << trigAlg_tag << std::endl;}
	    noSftwTrig_counter ++;
	    return false;
	}

	bool event_pass = false;
	for(auto this_flash : opflashes)
	{
		if(this_flash.TotalPE() >= 50)
		{
			if(event_pass == false){flash_pass_counter++; event_pass = true;}
			if(_verbose){std::cout << this_flash.Time() << std::endl;}
			if(this_flash.Time() >= -10 && this_flash.Time() <= 25)
			{
			  numi_time = this_flash.Time();
			  datatree->Fill();  
			}
			//return true;
		}
	}

	//default filters
	return false;
}


bool opFilter::finalize() {

	// If you need, you can store your ROOT class instance in the output
	// file. You have an access to the output file through "_fout" pointer.
	//
	// Say you made a histogram pointer h1 to store. You can do this:
	//
	// if (_fout) { _fout->cd(); _tree->Write(); }

	std::cout << "Total Events: " << total_flash_counter << std::endl;
	std::cout << "Remaining Events: " << flash_pass_counter << std::endl;
	std::cout << "Software Triggers: " << sftwTrig_counter << std::endl;
	std::cout << "Failed Software Triggers: " << noSftwTrig_counter << std::endl;

	datafile->Write();
	return true;
}


}
#endif
