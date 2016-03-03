import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("razorTuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/70001/60962869-5815-E511-BA9E-02163E014297.root'
        'file:60962869-5815-E511-BA9E-02163E014297.root'
        )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file:test.root")
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#global tag for PHYS14 asymptotic 25ns scenario
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

from RecoMET.METProducers.PFMET_cfi import pfMet
process.packedPFCandidates30 = cms.EDFilter("CandPtrSelector",
                                            src = cms.InputTag("particleFlow"),
                                            cut = cms.string("abs(eta) < 3.0"))
process.pfMet30                      = pfMet.clone();
process.pfMet30.src                  = cms.InputTag('packedPFCandidates30')

process.ntuples = cms.EDAnalyzer('RazorTuplizer',
                                 isData = cms.bool(False),
                                 useGen = cms.bool(True),
                                 enableTriggerInfo = cms.bool(False),

                                 triggerPathNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorHLTPathnames.dat"),
                                 eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
                                 muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
                                 photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

                                 triggerBits = cms.InputTag("TriggerResults","","HLT"),
                                 triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"),

                                 vertices = cms.InputTag("offlinePrimaryVertices"),
                                 muons = cms.InputTag("muons"),
                                 electrons = cms.InputTag("gedGsfElectrons"),
                                 taus = cms.InputTag("hpsPFTauProducer"),
                                 photons = cms.InputTag("gedPhotons"),
                                 jets = cms.InputTag("ak4PFJetsCHS"), # or not CHS?
                                 jetsPuppi = cms.InputTag("ak4PFJets"), # not PUPPI
                                 jetsAK8 = cms.InputTag("ak8PFJetsCHS"),# or not CHS?

                                 mets = cms.InputTag("pfMet"),
                                 metsNoHF = cms.InputTag("pfMet30"),
                                 metsPuppi = cms.InputTag("pfMet"),
                                 pfCands = cms.InputTag("PFCandidates"),

                                 genParticles = cms.InputTag("genParticles"),
                                 genJets = cms.InputTag("ak4GenJets"),

                                 lheInfo = cms.InputTag("externalLHEProducer", "", "LHE"),
                                 genInfo = cms.InputTag("generator", "", "SIM"),
                                 puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
                                 #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup

                                 rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
                                 rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
                                 rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
                                 rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
                                 rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
                                 rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

                                 )

process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
    process.packedPFCandidates30*
    process.pfMet30*
    process.ntuples)
