#ifndef GALLERY_FMWK_COSMIC_FILTER_CXX
#define GALLERY_FMWK_COSMIC_FILTER_CXX

#include "cosmic_filter.h"

namespace galleryfmwk {

bool cosmic_filter::initialize() {

	TH1D::AddDirectory(kFALSE);

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

	h_nu_vtx_x = new TH1D ("h_nu_vtx_x", "True Nu Vtx X", 50, -50, 300);
	h_nu_vtx_y = new TH1D ("h_nu_vtx_y", "True Nu Vtx Y", 50, -150, 150);
	h_nu_vtx_z = new TH1D ("h_nu_vtx_z", "True Nu Vtx Z", 50, -50, 1100);

	h_nu_vtx_xy_outside = new TH2D ("h_nu_vtx_xy_outside", "True Nu Vtx XY - Outside TPC", 50, -50, 300, 50, -150, 150);
	h_nu_vtx_zy_outside = new TH2D ("h_nu_vtx_zy_outside", "True Nu Vtx ZY - Outside TPC", 50, -50, 1100, 50, -150, 150);

	h_nu_vtx_xy_inside = new TH2D ("h_nu_vtx_xy_inside", "True Nu Vtx XY - Inside TPC", 50, -50, 300, 50, -150, 150);
	h_nu_vtx_zy_inside = new TH2D ("h_nu_vtx_zy_inside", "True Nu Vtx ZY - Inside TPC", 50, -50, 1100, 50, -150, 150);

	num_events = 0;
	num_events_remaining = 0;

	return true;

}

bool cosmic_filter::analyze(gallery::Event * ev) {

	bool setWantCC(wantCC);

	// auto const & mcparticle
	//         = ev->getValidHandle<std::vector <simb::MCParticle> > ("largeant");
	// auto const & mcparts(*mcparticle);

	auto const & mctruth = ev->getValidHandle< std::vector < simb::MCTruth> > ("generator");
	auto const & mctrue(*mctruth);

	//auto const & mcneutrino = ev->getValidHandle<std::vector<simb::MCNeutrino> > (_mc_part_tag);
	//auto const & mcnus(*mcneutrino);

	for(auto mct : mctrue)
	{
		auto const mcnu = mct.GetNeutrino();
		auto const mcpart = mcnu.Nu();

		const double nu_vtx_x = mcpart.Vx();
		const double nu_vtx_y = mcpart.Vy();
		const double nu_vtx_z = mcpart.Vz();

		h_nu_vtx_x->Fill(nu_vtx_x);
		h_nu_vtx_y->Fill(nu_vtx_y);
		h_nu_vtx_z->Fill(nu_vtx_z);

		//if true vertex is outside the volume
		if(_utility_instance.inFV(nu_vtx_x, nu_vtx_y, nu_vtx_z, x_boundary2, x_boundary1, y_boundary2, y_boundary1, z_boundary1, z_boundary2) == false)
		{
			h_nu_vtx_xy_outside->Fill(nu_vtx_x, nu_vtx_y);
			h_nu_vtx_zy_outside->Fill(nu_vtx_z, nu_vtx_y);
			//return false;
		}
		//if true vertex is inside the volume
		if(_utility_instance.inFV(nu_vtx_x, nu_vtx_y, nu_vtx_z, x_boundary2, x_boundary1, y_boundary2, y_boundary1, z_boundary1, z_boundary2) == true)
		{
			h_nu_vtx_xy_inside->Fill(nu_vtx_x, nu_vtx_y);
			h_nu_vtx_zy_inside->Fill(nu_vtx_z, nu_vtx_y);
			//return false;
		}

		num_events++;
		//if CC
		if(mcnu.CCNC() == false)
		{
			if(setWantCC == true)
			{
				num_events_remaining++;
				//if I only want events inside the volume
				if(_utility_instance.inFV(nu_vtx_x, nu_vtx_y, nu_vtx_z, x_boundary2, x_boundary1, y_boundary2, y_boundary1, z_boundary1, z_boundary2) == false)
				{
					return false;
				}
				//if I only want events outside the volume
				// if(_utility_instance.inFV(nu_vtx_x, nu_vtx_y, nu_vtx_z, x_boundary2, x_boundary1, y_boundary2, y_boundary1, z_boundary1, z_boundary2) == true)
				// {
				//      return false;
				// }
				return true;
			}
			if(setWantCC == false)
			{
				return false;
			}
		}
		//if NC
		if(mcnu.CCNC() == true)
		{
			if(setWantCC == true)
			{
				return false;
			}
			if(setWantCC == false)
			{
				return true;
			}
		}
	}//end loop mcnus

	return false;

}

bool cosmic_filter::finalize() {

	TCanvas * c1 = new TCanvas();
	c1->cd();
	h_nu_vtx_x->Draw();
	h_nu_vtx_x->GetXaxis()->SetTitle("X Vertex [cm]");
	h_nu_vtx_x->GetXaxis()->SetTitle("Events");
	c1->Print("true_nu_vtx_x.pdf");
	TCanvas * c2 = new TCanvas();
	c2->cd();
	h_nu_vtx_y->Draw();
	h_nu_vtx_y->GetXaxis()->SetTitle("Y Vertex [cm]");
	h_nu_vtx_y->GetXaxis()->SetTitle("Events");
	c2->Print("true_nu_vtx_y.pdf");
	TCanvas * c3 = new TCanvas();
	c3->cd();
	h_nu_vtx_z->Draw();
	h_nu_vtx_z->GetXaxis()->SetTitle("Z Vertex [cm]");
	h_nu_vtx_z->GetXaxis()->SetTitle("Events");
	c3->Print("true_nu_vtx_z.pdf");

	TCanvas * c4a = new TCanvas();
	c4a->cd();
	h_nu_vtx_xy_outside->Draw("colz");
	h_nu_vtx_xy_outside->GetXaxis()->SetTitle("X Vertex [cm]");
	h_nu_vtx_xy_outside->GetYaxis()->SetTitle("Y Vertex [cm]");
	c4a->Print("true_nu_vtx_xy_outside.pdf");
	TCanvas * c4b = new TCanvas();
	c4b->cd();
	h_nu_vtx_zy_outside->Draw("colz");
	h_nu_vtx_zy_outside->GetXaxis()->SetTitle("Z Vertex [cm]");
	h_nu_vtx_zy_outside->GetYaxis()->SetTitle("Y Vertex [cm]");
	c4b->Print("true_nu_vtx_zy_outside.pdf");
	TCanvas * c4c = new TCanvas();
	c4c->cd();
	h_nu_vtx_xy_inside->Draw("colz");
	h_nu_vtx_xy_inside->GetXaxis()->SetTitle("X Vertex [cm]");
	h_nu_vtx_xy_inside->GetYaxis()->SetTitle("Y Vertex [cm]");
	c4c->Print("true_nu_vtx_xy_inside.pdf");
	TCanvas * c4d = new TCanvas();
	c4d->cd();
	h_nu_vtx_zy_inside->Draw("colz");
	h_nu_vtx_zy_inside->GetXaxis()->SetTitle("Z Vertex [cm]");
	h_nu_vtx_zy_inside->GetYaxis()->SetTitle("Y Vertex [cm]");
	c4d->Print("true_nu_vtx_zy_inside.pdf");

	std::cout << "Total Events CCNC Filter:     " << num_events << std::endl;
	std::cout << "Remaining Events CCNC Filter: " << num_events_remaining << std::endl;

	return true;

}
}

#endif
