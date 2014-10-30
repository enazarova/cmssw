import FWCore.ParameterSet.Config as cms

checkflattening = cms.EDAnalyzer("CheckFlattening",
                            vtxCollection_=cms.InputTag("hiSelectedVertex"),
                            caloCollection_=cms.InputTag("towerMaker"),
                            inputPlanes_ = cms.InputTag("hiEvtPlane","recoLevel"),
                            trackCollection_=cms.InputTag("hiGeneralTracks"),                           
                            centrality_ = cms.InputTag("centralityBin")
 )
