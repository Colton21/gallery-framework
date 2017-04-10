import ROOT
from ROOT import gallery, galleryfmwk, larutil

import sys


def process_file(_file):

    # Create ana_processor instance
    my_proc = galleryfmwk.ana_processor()

    myfile = open(_file, 'r')

    # Set input root file
    for _f in myfile:
        # print _f
        my_proc.add_input_file(_f[:-1])

    # make root run in batch mode
    ROOT.gROOT.SetBatch(1)

    # Specify output root file name
    #my_proc.set_ana_output_file(_file.replace('.root', '') + "_ana.root")
    # my_proc.set_output_file("")

    use_opFilter = False
    use_inFVFilter = True
    use_cosmic_ana = True

    if(use_opFilter == False and use_cosmic_ana == False and use_inFVFilter == False):
        print 'No module was selected ... exiting'
        exit(1)

    opFilterModule = galleryfmwk.optFilter()
    opFilterModule.setTrackProducer("pandoraNuKHit")
    opFilterModule.setShowerProducer("showerrecopandora")
    opFilterModule.setFlashProducer("simpleFlashBeam")
    opFilterModule.togglePlotting(True)
    opFilterModule.setPEThreshold(50)
    opFilterModule.setVerbose(False)

    inFVModule = galleryfmwk.inFV_filter()
    inFVModule.setPfpProducer("pandoraNu")
    inFVModule.setPfpCosmicProducer("pandoraCosmic")
    inFVModule.setNearestCutDist(5)
    inFVModule.fiducial_volume_x_right(0)
    inFVModule.fiducial_volume_x_left(0)
    inFVModule.fiducial_volume_y_up(0)
    inFVModule.fiducial_volume_y_down(0)
    inFVModule.fiducial_volume_z_back(0)
    inFVModule.fiducial_volume_z_front(0)
    inFVModule.setVerbose(False)

    cosmicanaModule = galleryfmwk.cosmic_ana()
    cosmicanaModule.setPfpProducer("pandoraNu")
    cosmicanaModule.setPfpCosmicProducer("pandoraCosmic")
    cosmicanaModule.setNearestCutDist(5)
    cosmicanaModule.fiducial_volume_x_right(0)
    cosmicanaModule.fiducial_volume_x_left(0)
    cosmicanaModule.fiducial_volume_y_up(0)
    cosmicanaModule.fiducial_volume_y_down(0)
    cosmicanaModule.fiducial_volume_z_back(0)
    cosmicanaModule.fiducial_volume_z_front(0)
    cosmicanaModule.setVerbose(False)

    # Attach an analysis unit ... here we use a base class which do
    if(use_opFilter == True):
        my_proc.add_process(opFilterModule)
    if(use_inFVFilter == True):
        my_proc.add_process(inFVModule)
    if(use_cosmic_ana == True):
        my_proc.add_process(cosmicanaModule)

    my_proc.run()


def main():

    if len(sys.argv) < 2:
        print "Error: must include an input file."
        exit()

    _file = sys.argv[-1]
    process_file(_file)


if __name__ == '__main__':
    main()
