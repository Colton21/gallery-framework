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

	h_flash_zywidth = new TH2D("h_flash_zywidth", "h_flash_zywidth", 50, 0, 1050, 60, -120, 120);
	h_flash_zycenter = new TH2D("h_flash_zycenter", "h_flash_zycenter", 50, 0, 50, 50, 0, 50);

	bool threshold_plotting(_threshold_plotting);
	if(threshold_plotting == true)
	{
		h_flash_threshold = new TH1D ("h_flash_threshold", "h_flash_threshold", 50, 0, 100);
	}

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
	double pe_threshold(_pe_threshold);
	bool threshold_plotting(_threshold_plotting);

	auto const & opf = ev->getValidHandle<std::vector < recob::OpFlash> >(flash_tag);
	auto const & opflashes(*opf);

	total_flash_counter++;

	if(threshold_plotting == true)
	{
		for(auto opflsh : opflashes)
		{
			if(opflsh.Time() >= 3 && opflsh.Time() <= 5)
			{
				for(int i = 0; i < 100; i+=2)
				{
					if(opflsh.TotalPE() >= i)
					{
						h_flash_threshold->Fill(i);
					}
				}
			}
		}
	}


	for(auto opflsh : opflashes)
	{
		if(opflsh.Time() >= 3 && opflsh.Time() <= 5)
		{
			if(opflsh.TotalPE() >= pe_threshold)
			{
				flash_pass_counter++;

				h_flash_zycenter->Fill(opflsh.ZCenter(), opflsh.YCenter());
				h_flash_zywidth->Fill(opflsh.ZWidth(), opflsh.YWidth());


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

	TCanvas * c1a = new TCanvas();
	c1a->cd();
	h_flash_zywidth->Draw("colz");
	h_flash_zywidth->GetXaxis()->SetTitle("z [cm]");
	h_flash_zywidth->GetYaxis()->SetTitle("y [cm]");
	c1a->Print("flahsh_zy_width.pdf");

	TCanvas * c1b = new TCanvas();
	c1b->cd();
	h_flash_zycenter->Draw("colz");
	h_flash_zycenter->GetXaxis()->SetTitle("z [cm]");
	h_flash_zycenter->GetYaxis()->SetTitle("y [cm]");
	c1b->Print("flash_zy_center.pdf");

	bool threshold_plotting(_threshold_plotting);
	if(threshold_plotting == true)
	{
		TCanvas * c2 = new TCanvas();
		h_flash_threshold->Draw();
		h_flash_threshold->GetXaxis()->SetTitle("PE Threshold");
		h_flash_threshold->GetYaxis()->SetTitle("Events");
		c2->Print("flash_threshold.pdf");
	}


	std::cout << "Total Events: " << total_flash_counter << std::endl;
	std::cout << "Remaining Events: " << flash_pass_counter << std::endl;

	return true;
}


}
#endif
