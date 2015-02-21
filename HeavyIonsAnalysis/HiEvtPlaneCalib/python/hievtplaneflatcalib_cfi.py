import FWCore.ParameterSet.Config as cms

hiEvtPlaneFlatCalib = cms.EDAnalyzer('HiEvtPlaneFlatCalib',
                                     vertexTag_=cms.InputTag("hiSelectedVertex"),
                                     centralityTag_=cms.InputTag("hiCentrality"),
                                     centralityBinTag_ = cms.InputTag("centralityBin","HFtowers"),
                                     centralityVariable = cms.string("HFtowers"),
                                     nonDefaultGlauberModel = cms.string(""),
                                     caloTag_=cms.InputTag("towerMaker"),
                                     castorTag_=cms.InputTag("CastorTowerReco"),
                                     trackTag_=cms.InputTag("hiGeneralTracks"),                           
                                     inputPlanesTag_ = cms.InputTag("hiEvtPlane","recoLevel"),
                                     genFlatPsi_ = cms.untracked.bool(True),
                                     FlatOrder_ = cms.untracked.int32(9),
                                     NumFlatBins_ = cms.untracked.int32(40),
                                     CentBinCompression_ = cms.untracked.int32(5),
                                     caloCentRef_ = cms.untracked.double(80.),
                                     caloCentRefWidth_ = cms.untracked.double(5.0),
                                     HFEtScale_ = cms.untracked.int32(3800),
                                     useOffsetPsi_ = cms.untracked.bool(True),
                                     useECAL_ = cms.untracked.bool(True),
                                     useHCAL_ = cms.untracked.bool(True),
                                     useTrackPtWeight_ = cms.untracked.bool(True),
                                     useMomentumCorrV1_ = cms.untracked.bool(True),
                                     minet_ = cms.untracked.double(-1.),
                                     maxet_ = cms.untracked.double(-1.),
                                     effm_ = cms.untracked.double(0.0),
                                     minpt_ = cms.untracked.double(0.5),
                                     maxpt_ = cms.untracked.double(3.0),
                                     minvtx_ = cms.untracked.double(-25.),
                                     maxvtx_ = cms.untracked.double(25.),
                                     dzerr_ = cms.untracked.double(10.),
                                     chi2_ = cms.untracked.double(40.),
                                     Noffmin_ = cms.untracked.int32 (-1),
                                     Noffmax_ = cms.untracked.int32 (10000),
                                     minrun_ = cms.untracked.int32 (210490),
                                     maxrun_ = cms.untracked.int32 (211631)
                                     )
