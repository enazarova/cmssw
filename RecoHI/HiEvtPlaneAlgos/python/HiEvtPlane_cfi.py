import FWCore.ParameterSet.Config as cms

hiEvtPlane = cms.EDProducer("EvtPlaneProducer",
                            vtxCollection_=cms.InputTag("hiSelectedVertex"),
                            caloCollection_=cms.InputTag("towerMaker"),
                            trackCollection_=cms.InputTag("hiGeneralTracks"),                           
                            useECAL_ = cms.untracked.bool(True),
                            useHCAL_ = cms.untracked.bool(True),
                            useTrackPtWeight_ = cms.untracked.bool(True),
                            minet_ = cms.untracked.double(0.3),
                            maxet_ = cms.untracked.double(80.0),
                            minpt_ = cms.untracked.double(0.3),
                            maxpt_ = cms.untracked.double(2.6),
                            minvtx_ = cms.untracked.double(-25.),
                            maxvtx_ = cms.untracked.double(25.),
                            dzerr_ = cms.untracked.double(10.),
                            chi2_ = cms.untracked.double(40.)
                            )
                            




    
