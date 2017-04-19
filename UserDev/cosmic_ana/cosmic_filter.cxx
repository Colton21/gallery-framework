#ifndef GALLERY_FMWK_COSMIC_FILTER_CXX
#define GALLERY_FMWK_COSMIC_FILTER_CXX

#include "cosmic_filter.h"

namespace galleryfmwk {

bool cosmic_filter::initialize() {

	num_events = 0;
	num_events_remaining = 0;

	return true;

}

bool cosmic_filter::analyze(gallery::Event * ev) {

	bool setWantCC(wantCC);

	auto const & mcparticle
	        = ev->getValidHandle<std::vector <simb::MCParticle> > ("largeant");
	auto const & mcparts(*mcparticle);

	auto const & mctruth = ev->getValidHandle< std::vector < simb::MCTruth> > ("generator");
	auto const & mctrue(*mctruth);

	//auto const & mcneutrino = ev->getValidHandle<std::vector<simb::MCNeutrino> > (_mc_part_tag);
	//auto const & mcnus(*mcneutrino);

	for(auto mct : mctrue)
	{
		auto const mcnu = mct.GetNeutrino();
		num_events++;
		//if CC
		if(mcnu.CCNC() == false)
		{
			if(setWantCC == true)
			{
				num_events_remaining++;
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

	std::cout << "Total Events: " << num_events << std::endl;
	std::cout << "Remaining Events: " << num_events_remaining << std::endl;

	return true;

}
}

#endif
