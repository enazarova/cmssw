import FWCore.ParameterSet.Config as cms

checkflattening = cms.EDAnalyzer("CheckFlattening",
                            vtxCollection_=cms.InputTag("hiSelectedVertex"),
                            caloCollection_=cms.InputTag("towerMaker"),
                            hiCentrality_ = cms.InputTag("hiCentrality"),  
                            inputPlanes_ = cms.InputTag("hiEvtPlane","recoLevel"),
                            inputPlanesFlat_ = cms.InputTag("hiEvtPlaneFlat","recoLevel"),
                            trackCollection_=cms.InputTag("hiGeneralTracks"),                           
                            centrality_ = cms.InputTag("centralityBin")
 )
