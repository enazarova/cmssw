import FWCore.ParameterSet.Config as cms

hiEvtPlaneFlat = cms.EDProducer('HiEvtPlaneFlatProducer',
                                vtxCollection_=cms.InputTag("hiSelectedVertex"),
                                inputPlanes_ = cms.InputTag("hiEvtPlane","recoLevel"),
                                FlatOrder_ = cms.untracked.int32(9),
                                NumFlatCentBins_ = cms.untracked.int32(50),
                                CentBinCompression_ = cms.untracked.int32(4)
                                )
