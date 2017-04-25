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

    # print my_proc._input_files
    # exit()

    # my_proc.add_input_file(_file)

    # make root run in batch mode
    ROOT.gROOT.SetBatch(1)

    # Specify output root file name
    #my_proc.set_ana_output_file(_file.replace('.root', '') + "_ana.root")
    # my_proc.set_output_file("")

    opFilterModule = galleryfmwk.opFilter()
    #opFilterModule.setTrackProducer("pandoraNuKHit")
    #opFilterModule.setShowerProducer("showerrecopandora")
    #opFilterModule.setFlashProducer("opflashSat")
    opFilterModule.setFlashProducer("simpleFlashBeam")
    opFilterModule.setSoftWareProducer("daq")
    opFilterModule.setTrigAlgProducer("EXT_NUMIwin_FEMBeamTriggerAlgo")
    opFilterModule.setVerbose(False)

    # Attach an analysis unit ... here we use a base class which do
    my_proc.add_process(opFilterModule)

    my_proc.run()


def main():

    if len(sys.argv) < 2:
        print "Error: must include an input file."
        exit()

    _file = sys.argv[-1]
    process_file(_file)


if __name__ == '__main__':
    main()