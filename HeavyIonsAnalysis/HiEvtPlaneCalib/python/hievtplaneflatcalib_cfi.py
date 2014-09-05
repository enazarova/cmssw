import FWCore.ParameterSet.Config as cms

hiEvtPlaneFlatCalib = cms.EDAnalyzer('HiEvtPlaneFlatCalib',
                                     vtxCollection_=cms.InputTag("hiSelectedVertex"),
                                     inputPlanes_ = cms.InputTag("hiEvtPlane","recoLevel"),
                                     centrality_  = cms.InputTag("hiCentrality"),
                                     genFlatPsi_ = cms.untracked.bool(True),
                                     useOffsetPsi_ = cms.untracked.bool(True),
                                     caloCollection_=cms.InputTag("towerMaker"),
                                     trackCollection_=cms.InputTag("hiGeneralTracks"),                           
                                     useECAL_ = cms.untracked.bool(True),
                                     useHCAL_ = cms.untracked.bool(True),
                                     useTrackPtWeight_ = cms.untracked.bool(True),
                                     useMomentumCorrV1_ = cms.untracked.bool(True),
                                     minet_ = cms.untracked.double(0.5),
                                     maxet_ = cms.untracked.double(80.0),
                                     effm_ = cms.untracked.double(0.0),
                                     minpt_ = cms.untracked.double(0.7),
                                     maxpt_ = cms.untracked.double(2.6),
                                     minvtx_ = cms.untracked.double(-25.),
                                     maxvtx_ = cms.untracked.double(25.),
                                     dzerr_ = cms.untracked.double(10.),
                                     chi2_ = cms.untracked.double(40.)
                                     )
