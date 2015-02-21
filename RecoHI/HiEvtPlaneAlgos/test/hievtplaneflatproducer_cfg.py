import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("FlatCalib")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.GlobalTag.globaltag='GR_R_74_V0A::All'
process.MessageLogger.cerr.FwkReport.reportEvery=1000
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.GlobalTag.toGet.extend([

 cms.PSet(record = cms.string("HeavyIonRcd"),
  tag = cms.string("CentralityTable_HFtowers200_Glauber2010A_v5315x02_offline"),
  connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
  label = cms.untracked.string("HFtowers")
 )
])

#readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/00726191-ACA8-E411-AA70-02163E00E994.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/007E7AE0-ABA8-E411-9F36-02163E00E91E.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/02044B3B-B2A8-E411-B5B5-02163E00CAB8.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0218BFBE-ABA8-E411-AF64-0025904B2018.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/022EDFE2-ABA8-E411-98C2-002590494D1A.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/044E159A-ABA8-E411-AFFF-02163E00B70E.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/04A21CD0-A6A8-E411-9741-02163E00E880.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/04AB20AA-ABA8-E411-9D04-02163E00E969.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/04DFC312-A8A8-E411-862C-02163E00E713.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0666C395-ABA8-E411-82DA-02163E00F2CF.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/066F0B4D-A8A8-E411-A743-02163E00C06B.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/08522997-ABA8-E411-A4FF-02163E00E806.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/08777C73-ACA8-E411-9FA4-02163E00BFD0.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0CE8F940-ABA8-E411-A355-02163E00F524.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0CF13D7B-ABA8-E411-8F5C-02163E00E696.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0E0F90E0-ABA8-E411-B752-02163E00E5B8.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/0EEDC179-ADA8-E411-9588-02163E00EA03.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/10189D91-ABA8-E411-AE2C-02163E00CD63.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/1031BB03-33AA-E411-A905-02163E00CE07.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/104DF356-ABA8-E411-B9FE-02163E00E9F1.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/126925A1-ACA8-E411-8110-02163E00E982.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/12BB9E7C-ABA8-E411-A604-02163E00E9F6.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/14CD6A14-ACA8-E411-A71A-02163E00BBA9.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/16A4EA51-A7A8-E411-8739-02163E00F472.root',
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre6/HIMinBiasUPC/RECO/GR_R_74_V0A_RelVal_hi2011-v1/00000/16E0EFF0-A8A8-E411-8A19-02163E00BA2F.root',
),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            inputCommands=cms.untracked.vstring(
        'keep *',
        'drop *_hiEvtPlane_*_*'
        ),
                            dropDescendantsOfDroppedBranches=cms.untracked.bool(False)
                            )



#process.CondDBCommon.connect = "sqlite_file:HeavyIonRPRcd_PbPb2011_74X_v01_offline.db"
#process.PoolDBESSource = cms.ESSource("PoolDBESSource",
#                                       process.CondDBCommon,
#                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
#                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2011_74X_v01_offline')
#                                                                  )
#                                                         )
#                                      )

process.GlobalTag.toGet.extend([
        cms.PSet(record = cms.string("HeavyIonRPRcd"),
                 tag = cms.string('HeavyIonRPRcd_PbPb2011_74X_v01_offline'),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_PAT_000")
                 )
        ])


process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.hiEvtPlane.loadDB_ = cms.untracked.bool(True)
process.hiEvtPlaneFlat.genFlatPsi_ = cms.untracked.bool(True)
process.p = cms.Path(process.centralityBin*process.hiEvtPlane*process.hiEvtPlaneFlat)


                        

