import FWCore.ParameterSet.Config as cms

hiEvtPlaneFlat = cms.EDProducer('HiEvtPlaneFlatProducer',
                                vertexTag_=cms.InputTag("hiSelectedVertex"),
                                centralityTag_=cms.InputTag("hiCentrality"),
                                inputPlanesTag_ = cms.InputTag("hiEvtPlane","recoLevel"),
                                FlatOrder_ = cms.untracked.int32(9),
                                NumFlatBins_ = cms.untracked.int32(40),
                                CentBinCompression_ = cms.untracked.int32(5),
                                caloCentRef_ = cms.untracked.double(80.),
                                caloCentRefWidth_ = cms.untracked.double(5.0),
                                HFEtScale_ = cms.untracked.int32(3800),
                                useOffsetPsi_ = cms.untracked.bool(True)
                                )
