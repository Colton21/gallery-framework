#ifndef GALLERY_FMWK_TRUTH_FIG_CXX
#define GALLERY_FMWK_TRUTH_FIG_CXX

#include "truth_fig.h"

namespace galleryfmwk {

bool truth_fig::initialize() {

	//want to save a root file which has the histograms I need,
	//so I can create a TStack later with different files

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

	num_events = 0;
	num_events_remaining = 0;


	t1.Branch("mode0Ecc", &mode0Ecc, "mode0Ecc/D");
	t1.Branch("mode1Ecc", &mode1Ecc, "mode1Ecc/D");
	t1.Branch("mode2Ecc", &mode2Ecc, "mode2Ecc/D");
	t1.Branch("mode3Ecc", &mode3Ecc, "mode3Ecc/D");
	t1.Branch("mode10Ecc", &mode10Ecc, "mode10Ecc/D");


	t1.Branch("mode0Enc", &mode0Enc, "mode0Enc/D");
	t1.Branch("mode1Enc", &mode1Enc, "mode1Enc/D");
	t1.Branch("mode2Enc", &mode2Enc, "mode2Enc/D");
	t1.Branch("mode3Enc", &mode3Enc, "mode3Enc/D");
	t1.Branch("mode10Enc", &mode10Enc, "mode10Enc/D");
	return true;

}

bool truth_fig::analyze(gallery::Event * ev) {

	// auto const & mcparticle
	//         = ev->getValidHandle<std::vector <simb::MCParticle> > ("largeant");
	// auto const & mcparts(*mcparticle);

	auto const & mctruth = ev->getValidHandle< std::vector < simb::MCTruth> > ("generator");
	auto const & mctrue(*mctruth);

	for(auto mct : mctrue)
	{
		auto const mcnu = mct.GetNeutrino();
		auto const mcpart = mcnu.Nu();

		const double nu_vtx_x = mcpart.Vx();
		const double nu_vtx_y = mcpart.Vy();
		const double nu_vtx_z = mcpart.Vz();


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
			//"qe-like"
			if(mcnu.Mode() == 0)
			{

			}
			//resonant
			if(mcnu.Mode() == 1)
			{

			}
			//DIS
			if(mcnu.Mode() == 2)
			{

			}
			//Coherent
			if(mcnu.Mode() == 3)
			{

			}
			//MEC - meson exchange current
			if(mcnu.Mode() == 10)
			{

			}

		}
		//if NC
		if(mcnu.CCNC() == true)
		{

			//"qe-like"
			if(mcnu.Mode() == 0)
			{

			}
			if(mcnu.Mode() == 1)
			{

			}
			if(mcnu.Mode() == 2)
			{

			}
			if(mcnu.Mode() == 3)
			{

			}
			if(mcnu.Mode() == 10)
			{

			}

		}
	}//end loop mcnus

	return false;

}

bool truth_fig::finalize() {


	std::cout << "Total Events CCNC Filter:     " << num_events << std::endl;
	std::cout << "Remaining Events CCNC Filter: " << num_events_remaining << std::endl;

	return true;

}
}

#endif
