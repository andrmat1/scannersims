#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description=
"""Process HPGe detector simulations from one or more G4Simple ROOT files
with a single detector geometry from the LEGEND200 metadata repo.""")

parser.add_argument('file', nargs = '+',
                    help="Input MaGe file(s); can use wildcards")
parser.add_argument('-o', '--output', default='', type=str,
                    help="Name of output ROOT file. Default to 'mpp_[input file name].root'")
parser.add_argument('-v', '--verbose', default=False, action='store_true',
                    help="Verbose output")
parser.add_argument('-i', '--mageid', default=0, type=int,
                    help="MaGe detector ID to process. By default, process all detectors with ID starting with 1; this will use the same geometry for all!")
parser.add_argument('-g', '--geometry', default=None, type=str,
                    help="Detector geometry style, using the L200 ICPC file specification in the L200 metadata repo. By default, use a generic ICPC.")
parser.add_argument('-d', '--dlparams', default=None, type=float, nargs=3,
                    help="Parameters for the ExLinT transition dead layer model. Takes three arguments, each in mm. Default is (1, 0.5, 0.3).")
parser.add_argument('-E', '--energyres', default=None, type=float, nargs=6,
                    help="Energy resolution parameters. Takes 6 args: (sig0, sig1, sig2, h_tail, tau0, tau1) in units of some power of keV.")
parser.add_argument('-t', '--threshold', default=5, type=float,
                    help="Energy threshold in keV. Default is 5 keV.")
parser.add_argument('-z', '--zcenter', default=False, action='store_true',
                    help="Put origin at center of crystal.")
parser.add_argument('-n', '--nevents', default=-1, type=int,
                    help="Number of events to process. By default process all")
parser.add_argument('-T', '--threads', default=1, type=int,
                    help="Number of threads to use for processing")

args = parser.parse_args()
output = args.output if args.output != '' else 'mpp_'+args.file[0].split('/')[-1]

import time
import magepostproc as mpp
import ROOT

mm = mpp.CLHEP.mm
um = mpp.CLHEP.micrometer
us = mpp.CLHEP.us
keV = mpp.CLHEP.keV

########################### Load detector geometry ############################
crystal = mpp.MGTL200ICPCData()
if args.geometry:
    crystal.ReadDetectorData(args.geometry)
else:
    meas = mpp.MGTL200ICPCData.Measurements()
    meas.diameter           = 74*mm
    meas.height             = 80*mm
    meas.contact_radius     = 0
    meas.contact_depth      = 0
    meas.groove_width       = 4.20*mm
    meas.groove_radius      = 9.65*mm
    meas.groove_depth       = 2.00*mm
    meas.well_radius        = 10*mm
    meas.well_depth         = 52*mm
    meas.taper_outer_radius = 0
    meas.taper_outer_height = 0
    meas.taper_inner_radius = 0
    meas.taper_inner_height = 0
    crystal.SetMeasurements(meas)

if args.zcenter:
    crystal.SetZOffset(crystal.GetCrystalHeight()/2)

############################### Set Dead Layer ################################
dead_layer = mpp.MGTDLExLinT(1*mm, 0.5*mm, 0.3*mm)
if args.dlparams:
    dead_layer.SetDLPars(*args.dlparams)
crystal.SetNPlusOuterDLModel(dead_layer)
crystal.SetNPlusInnerDLModel(dead_layer)

############################ Set Energy Resolution ############################
EParams = mpp.MPPEnergyCalculator.EnergyParams() # default values in .hh file
if args.energyres:
    EParams.fSigma0 = args.energyres[0]
    EParams.fSigma1 = args.energyres[1]
    EParams.fSigma2 = args.energyres[2]
    EParams.fHTail0 = args.energyres[3]
    EParams.fTau0 = args.energyres[4]
    EParams.fTau1 = args.energyres[5]

#--------------------------------------------------------------------------
#---------------------Set up post-processing sequence----------------------
#--------------------------------------------------------------------------

selector = mpp.TAMSelector()

# Read input leafs
inputDataPoster = mpp.MPPInputDataPoster()
inputDataPoster.PostInputLeaf['int']("event", "eventid")
selector.AddInput(inputDataPoster)

# Add windower
windowingTime = 10.*us
windower = mpp.MPPStepsWindower(windowingTime, mpp.MPPStepsWindower.kG4SimpleStepsReader)
if args.mageid:
    windower.AddSensVolGroup("ge_id", "ge_steps", [args.mageid])
    windower.SetDefaultSensVolIDs("unused_id")
    windower.SetDefaultStepsList("unused_steps")
else:
    windower.AddSensVolGroup("ge_id", "ge_steps", 1)
    windower.SetDefaultSensVolIDs("unused_id")
    windower.SetDefaultStepsList("unused_steps")
windower.SetVerbosity(args.verbose)
selector.AddInput(windower)

# Dead layer processor
deadLayerProcessor = mpp.MPPGeDeadLayerProcessor("ge_id", "ge_steps", crystal)
selector.AddInput(deadLayerProcessor)

# Energy smearer
energyCalc = mpp.MPPEnergyCalculator("ge_id", "ge_steps")
energyCalc.SetDefaultSensVolParams(EParams)
energyCalc.activeness = deadLayerProcessor.activeness
energyCalc.SetVerbosity(args.verbose)
selector.AddInput(energyCalc)


# Step clusterer (1 mm clusters)
stepsClusterer = mpp.MPPStepsClusterer("ge_id", "ge_steps", 1.*mm)
stepsClusterer.activeness = energyCalc.activeness
stepsClusterer.windowT0  = windower.windowT0
stepsClusterer.SetVerbosity(args.verbose)
selector.AddInput(stepsClusterer)


# Tree writer
treeWriter = mpp.MPPWindowedDataTreeWriter(output, "ge_id")
treeWriter.AddPostedObject("eventid")
treeWriter.AddPostedObject(windower.windowT0)
treeWriter.AddPostedObject(energyCalc.sensVolE)
treeWriter.AddPostedObject(energyCalc.sensVolEdep)
treeWriter.AddPostedObject(energyCalc.sumE)
treeWriter.AddPostedObject(stepsClusterer.nClusters)
treeWriter.AddPostedObject(stepsClusterer.clusters)
selector.AddInput(treeWriter)


#--------------------------------------------------------------------------
#----------------------Load input tree and process-------------------------
#--------------------------------------------------------------------------

print('-'*65)
ch = ROOT.TChain("g4sntuple")
for f in args.file:
    print('processing', f)
    if not ch.Add(f):
        print("File", f, "does not contain a TTree fTree")
print("File(s) have", ch.GetEntries(), "entries")
print('-'*65)

start_time = time.time()
if args.threads > 1: # multiprocessing
    pool = ROOT.TTreeProcessorMP(args.threads)
    if args.nevents<0:
        pool.Process(ch, selector)
    else:
        pool.Process(ch, selector, "", args.nevents)
else: # single thread
    if args.nevents<0:
        ch.Process(selector)
    else:
        ch.Process(selector, "", args.nevents)
print("Took", time.time()-start_time, "s")

