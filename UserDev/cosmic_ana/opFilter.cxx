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
	flash_center_cut_counter = 0;

	TH1D::AddDirectory(kFALSE);

	h_flash_zycenter = new TH2D("h_flash_zycenter", "h_flash_zycenter", 50, 0, 1050, 60, -120, 120);
	h_flash_zywidth = new TH2D("h_flash_zywidth", "h_flash_zywidth", 50, 0, 200, 50, 0, 100);

	h_largest_flash = new TH1D("h_largest_flash", "h_largest_flash", 50, 0, 1500);
	h_largest_flash_y = new TH2D("h_largest_flash_y", "h_largest_flash_y", 50, 0, 1500, 50, -120, 120);
	h_nue_shwr_vtx_flash_dist_zy = new TH1D("h_nue_shwr_vtx_flash_dist_zy", "h_nue_shwr_vtx_flash_dist_zy", 50, 0, 200);
	h_nue_shwr_vtx_flash_width_zy = new TH2D("h_nue_shwr_vtx_flash_width_zy", "h_nue_shwr_vtx_flash_width_zy", 50, -150, 200, 50, -150, 300);

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
	art::InputTag pfp_tag(_pfp_tag);
	art::InputTag tracks_tag(_track_producer);
	art::InputTag showers_tag(_shower_producer);
	art::InputTag flash_tag(_flash_producer);
	double pe_threshold(_pe_threshold);
	double min_time(_min_time);
	double max_time(_max_time);
	bool threshold_plotting(_threshold_plotting);
	bool flash_center_cut(_flash_center_cut);
	double flash_center_cut_distance(_flash_center_cut_distance);

	auto const & pfp
	        = ev->getValidHandle<std::vector <recob::PFParticle> > (pfp_tag);
	auto const & pfparticles(*pfp);

	auto const & opf = ev->getValidHandle<std::vector < recob::OpFlash> >(flash_tag);
	auto const & opflashes(*opf);

	art::FindMany<recob::Vertex> vertex_for_pfp(pfp, *ev, "pandoraNu");


	bool sufficient_flash = false;

	total_flash_counter++;

	if(threshold_plotting == true)
	{
		for(auto opflsh : opflashes)
		{
			if(opflsh.Time() >= min_time && opflsh.Time() <= max_time)
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

	double largestPE = 0;
	bool this_event_flash = false;
	recob::OpFlash this_flash;
	for(auto opflsh : opflashes)
	{
		if(opflsh.Time() >= min_time && opflsh.Time() <= max_time)
		{
			if(opflsh.TotalPE() >= pe_threshold)
			{
				if(this_event_flash == false) {this_event_flash = true; }

				h_flash_zycenter->Fill(opflsh.ZCenter(), opflsh.YCenter());
				h_flash_zywidth->Fill(opflsh.ZWidth(), opflsh.YWidth());

				sufficient_flash = true;

				//let's get the largest flash
				if(opflsh.TotalPE() > largestPE)
				{
					largestPE = opflsh.TotalPE();
					this_flash = opflsh;
				}
			}//end if over PE threshold
		}//end if in time window
	}//end looping flashes



	//default filters
	if(sufficient_flash == true)
	{
		//marco often considers the largest flash the primary
		h_largest_flash->Fill(largestPE);
		h_largest_flash_y->Fill(largestPE, this_flash.YCenter());

		for(auto pfparts : pfparticles)
		{
			if(pfparts.IsPrimary() == true)
			{
				if(pfparts.PdgCode() == 12)
				{
					for(std::size_t const i : pfparts.Daughters())
					{
						auto const daughter = pfparticles.at(i);
						//let's get the vertex associations for the daughters
						std::vector<recob::Vertex const*> d_vertex;
						vertex_for_pfp.get(i, d_vertex);
						if(d_vertex.size() == 0 )
						{
							if(_verbose) {std::cout << "No vertex association found for daughter!" << std::endl; }
							return false;
						}
						//get vertex vector
						double d_xyz [3];
						d_vertex.at(0)->XYZ(d_xyz);
						if(daughter.PdgCode() == 11)
						{
							//let's get the distance between the nue-like shower vertex and the largest flash
							const double zy_dist = sqrt((this_flash.ZCenter() - d_xyz[2]) * (this_flash.ZCenter() - d_xyz[2])
							                            + (this_flash.YCenter() - d_xyz[1]) * (this_flash.YCenter() - d_xyz[1]));
							h_nue_shwr_vtx_flash_dist_zy->Fill(zy_dist);
							h_nue_shwr_vtx_flash_width_zy->Fill(this_flash.ZWidth() - (d_xyz[2] - this_flash.ZCenter()), this_flash.YWidth() - (d_xyz[1] - this_flash.YCenter()));

							//cut on the distance from the center of the largest flash
							//to the reconstructed shower vertex
							//a guess from the plots is that ~100 cm should be good
							if(flash_center_cut == true)
							{
								if(zy_dist >= flash_center_cut_distance)
								{
									if(_verbose)
									{
										std::cout << "Shower vertex too far from flash center" << std::endl;
									}
									flash_center_cut_counter++;
									return false;
								}
							}


							// //let's get the shower associations
							// std::vector<recob::Shower const*> shower;
							// shower_for_pfp.get(i, shower);
							// if(shower.size() == 0)
							// {
							//      if(_verbose) {std::cout << "No shower for this pfp shower!" << std::endl; }
							//      return false;
							// }
						}
					}//end daughters
				}//end if nue-like
			}//end if primary
		}//end looping pfp

		flash_pass_counter++;
		return true;
	}
	else{return false; }
}

//for timing: marco does 3-5 us window, and total PE over all opdets >= 50 PE

bool optFilter::finalize() {

	// If you need, you can store your ROOT class instance in the output
	// file. You have an access to the output file through "_fout" pointer.
	//
	// Say you made a histogram pointer h1 to store. You can do this:
	//
	// if (_fout) { _fout->cd(); _tree->Write(); }

	bool flash_center_cut(_flash_center_cut);

	TCanvas * c1a = new TCanvas();
	c1a->cd();
	h_flash_zywidth->Draw("colz");
	h_flash_zywidth->GetXaxis()->SetTitle("z [cm]");
	h_flash_zywidth->GetYaxis()->SetTitle("y [cm]");
	c1a->Print("flash_zy_width.pdf");

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
		c2->cd();
		h_flash_threshold->Draw();
		h_flash_threshold->GetXaxis()->SetTitle("PE Threshold");
		h_flash_threshold->GetYaxis()->SetTitle("Events");
		c2->Print("flash_threshold.pdf");
	}

	TCanvas * c3a = new TCanvas();
	c3a->cd();
	h_largest_flash->Draw();
	h_largest_flash->GetXaxis()->SetTitle("Largest Flash (always >= 50)");
	h_largest_flash->GetYaxis()->SetTitle("Counts");
	c3a->Print("largest_flash.pdf");
	TCanvas * c3b = new TCanvas();
	c3b->cd();
	h_largest_flash_y->Draw("colz");
	h_largest_flash_y->GetXaxis()->SetTitle("Largest Flash (always >= 50)");
	h_largest_flash_y->GetYaxis()->SetTitle("Flash y [cm]");
	c3b->Print("largest_flash_y.pdf");
	TCanvas * c3c = new TCanvas();
	h_nue_shwr_vtx_flash_dist_zy->Draw();
	h_nue_shwr_vtx_flash_dist_zy->GetXaxis()->SetTitle("Nue shwr vtx to Flash Center [cm]");
	h_nue_shwr_vtx_flash_dist_zy->GetYaxis()->SetTitle("Counts");
	c3c->Print("nue-like_shwr_vtx_flashCenter.pdf");
	TCanvas * c3d = new TCanvas();
	h_nue_shwr_vtx_flash_width_zy->Draw("colz");
	h_nue_shwr_vtx_flash_width_zy->GetXaxis()->SetTitle("Z Flash Width - Vertex Dist [cm]");
	h_nue_shwr_vtx_flash_width_zy->GetYaxis()->SetTitle("Y Flash Width - Vertex Dist [cm]");
	c3d->Print("nue-like_shwr_vtx_flashCenter_width.pdf");

	std::cout << "Total Events     (OpFilter):     " << total_flash_counter << std::endl;
	if(flash_center_cut == false) {std::cout << "No cut on flash center to vtx" << std::endl; }
	if(flash_center_cut == true) {std::cout << "Events Cut: flash center to vtx: " << flash_center_cut_counter << std::endl; }
	std::cout << "Remaining Events (OpFilter):     " << flash_pass_counter << std::endl;

	return true;
}


}
#endif
