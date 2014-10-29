import FWCore.ParameterSet.Config as cms

hiEvtPlaneFlatCalib = cms.EDAnalyzer('HiEvtPlaneFlatCalib',
                                     vtxCollection_=cms.InputTag("hiSelectedVertex"),
                                     inputPlanes_ = cms.InputTag("hiEvtPlane","recoLevel"),
                                     genFlatPsi_ = cms.untracked.bool(True),
                                     FlatOrder_ = cms.untracked.int32(9),
                                     NumFlatCentBins_ = cms.untracked.int32(50),
                                     CentBinCompression_ = cms.untracked.int32(4)
                                     )
