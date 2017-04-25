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

    use_opFilter = True
    use_flashCenterCut = True
    use_ccnc = True
    use_inFVFilter = True
    use_cosmic_ana = True

    # fiducial volume settings
    right = 10
    left = 10
    up = 30
    down = 10
    back = 10
    front = 10

    if(use_opFilter == False and use_cosmic_ana == False and use_inFVFilter == False and use_ccnc == False):
        print 'No module was selected ... exiting'
        exit(1)

    opFilterModule = galleryfmwk.optFilter()
    opFilterModule.setTrackProducer("pandoraNu")
    opFilterModule.setPfpProducer("pandoraNu")
    opFilterModule.setShowerProducer("pandoraNu")
    opFilterModule.setFlashProducer("simpleFlashBeam")
    opFilterModule.togglePlotting(True)
    opFilterModule.setPEThreshold(50)
    # For BNB use 3 and 5
    # For NuMI use 5 16 (5.03-16.75)
    min_time = 5
    max_time = 16
    opFilterModule.setMinTime(min_time)
    opFilterModule.setMaxTime(max_time)
    opFilterModule.setVerbose(False)
    opFilterModule.flashCenterCut(use_flashCenterCut)
    opFilterModule.flashCenterCutDistance(100)

    cosmicfilterModule = galleryfmwk.cosmic_filter()
    cosmicfilterModule.setMCProducer("generator")
    cosmicfilterModule.fiducial_volume_x_right(right)
    cosmicfilterModule.fiducial_volume_x_left(left)
    cosmicfilterModule.fiducial_volume_y_up(up)
    cosmicfilterModule.fiducial_volume_y_down(down)
    cosmicfilterModule.fiducial_volume_z_back(back)
    cosmicfilterModule.fiducial_volume_z_front(front)
    cosmicfilterModule.setWantCC(True)

    if(use_opFilter == True):
        print "Time Window: ", min_time, max_time

    inFVModule = galleryfmwk.inFV_filter()
    inFVModule.setPfpProducer("pandoraNu")
    inFVModule.setPfpCosmicProducer("pandoraCosmic")
    inFVModule.setNearestCutDist(5)
    inFVModule.fiducial_volume_x_right(right)
    inFVModule.fiducial_volume_x_left(left)
    inFVModule.fiducial_volume_y_up(up)
    inFVModule.fiducial_volume_y_down(down)
    inFVModule.fiducial_volume_z_back(back)
    inFVModule.fiducial_volume_z_front(front)
    inFVModule.setVerbose(False)

    cosmicanaModule = galleryfmwk.cosmic_ana()
    cosmicanaModule.setPfpProducer("pandoraNu")
    cosmicanaModule.setPfpCosmicProducer("pandoraCosmic")
    cosmicanaModule.setNearestCutDist(5)
    cosmicanaModule.fiducial_volume_x_right(right)
    cosmicanaModule.fiducial_volume_x_left(left)
    cosmicanaModule.fiducial_volume_y_up(up)
    cosmicanaModule.fiducial_volume_y_down(down)
    cosmicanaModule.fiducial_volume_z_back(back)
    cosmicanaModule.fiducial_volume_z_front(front)
    cosmicanaModule.setVerbose(False)

    # Attach an analysis unit ... here we use a base class which do
    if(use_ccnc == True):
        my_proc.add_process(cosmicfilterModule)
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
