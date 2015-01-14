import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("CheckFlattening")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load("RecoHI.HiCentralityAlgos/HiCentrality_cfi")
process.load("HeavyIonsAnalysis.HiEvtPlaneCalib/checkflattening_cfi")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

process.GlobalTag.globaltag = 'GR_R_53_LV6::All'
process.MessageLogger.cerr.FwkReport.reportEvery=1000

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    centralitySrc = cms.InputTag("hiCentrality")
    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd.unl.edu//store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/C8D19E59-F1AF-E311-A910-FA163EC7F60B.root'
  )
)

##
## Uncomment the following and comment the subsequent toGet for local db file
##

#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.CondDBCommon.connect = "sqlite_file:flatparms.db"
#process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
#                                      process.CondDBCommon,
#                                      toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
#                                                                 tag = cms.string('EPFlattening_HIRun2011_v5320devx02_offline')
#                                                                 )
#                                                        )
#                                      )

process.GlobalTag.toGet.extend([
        cms.PSet(record = cms.string("HeavyIonRPRcd"),
                 tag = cms.string('HeavyIonRPRcd_PbPb2011_5320_v01_offline'),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_PAT_000")
                 )
        ])

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("rpflat.root")
)

# Minimum bias trigger selection (later runs)
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasHFOrBSC = process.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v1"]

process.hiEvtPlane.trackCollection_ = cms.InputTag("hiGeneralTracks");
process.hiEvtPlane.minvtx_ = cms.untracked.double(-25.)
process.hiEvtPlane.maxvtx_ = cms.untracked.double(25.)


process.p = cms.Path(process.collisionEventSelection*process.hltMinBiasHFOrBSC*process.hiEvtPlane*process.hiEvtPlaneFlat*process.checkflattening)

