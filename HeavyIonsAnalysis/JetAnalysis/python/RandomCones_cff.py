from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoHI.HiJetAlgos.HiPFJetParameters_cff import *
from RecoHI.HiJetAlgos.HiCaloJetParameters_cff import *

akPu1PFCones = cms.EDProducer(
    "JetAlgorithmAnalyzer",
    HiPFJetParameters,
    AnomalousCellParameters,
    MultipleAlgoIteratorBlock,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.1),
    useInputJets = cms.untracked.bool(True),
    useTowersForBkg = cms.untracked.bool(True),
    centralityTag = cms.InputTag("hiCentrality"),
    evtPlaneTag = cms.InputTag("hiEvtPlane","recoLevel"),
    avoidNegative = cms.bool(False),
    patJetSrc = cms.untracked.InputTag("akPu1PFpatJets"),
    evtPlaneIndex = cms.untracked.int32(21),
    doBackToBack  = cms.untracked.bool(True),
    centrality  = cms.untracked.int32(-1),
    doRecoEvtPlane  = cms.untracked.bool(False),        
    doAnalysis  = cms.untracked.bool(False),
    puPtMin = cms.double(15.0)    
    )

akPu1PFCones.doPUOffsetCorr = True
akPu1PFCones.jetType = 'BasicJet'
akPu1PFCones.src = cms.InputTag("PFTowers")

akPu1PFCones = cms.PSet(
    initialSeed = cms.untracked.uint32(6),
    engineName = cms.untracked.string('TRandom3')
  )

akPu2PFCones = akPu1PFCones.clone(rParam = cms.double(0.2), patJetSrc = cms.untracked.InputTag("akPu2PFpatJets"))
akPu3PFCones = akPu1PFCones.clone(rParam = cms.double(0.3), patJetSrc = cms.untracked.InputTag("akPu3PFpatJets"))
akPu4PFCones = akPu1PFCones.clone(rParam = cms.double(0.4), patJetSrc = cms.untracked.InputTag("akPu4PFpatJets"))
akPu5PFCones = akPu1PFCones.clone(rParam = cms.double(0.5), patJetSrc = cms.untracked.InputTag("akPu5PFpatJets"))
akPu6PFCones = akPu1PFCones.clone(rParam = cms.double(0.6), patJetSrc = cms.untracked.InputTag("akPu6PFpatJets"))
#akPu7PFCones = akPu1PFCones.clone(rParam = cms.double(0.7), patJetSrc = cms.untracked.InputTag("akPu7PFpatJets"))

akPu1CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.1), patJetSrc = cms.untracked.InputTag("akPu1CalopatJets"))
akPu2CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.2), patJetSrc = cms.untracked.InputTag("akPu2CalopatJets"))
akPu3CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.3), patJetSrc = cms.untracked.InputTag("akPu3CalopatJets"))
akPu4CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.4), patJetSrc = cms.untracked.InputTag("akPu4CalopatJets"))
akPu5CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.5), patJetSrc = cms.untracked.InputTag("akPu5CalopatJets"))
akPu6CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.6), patJetSrc = cms.untracked.InputTag("akPu6CalopatJets"))
#akPu7CaloCones = akPu1PFCones.clone(src = cms.InputTag("towerMaker"), rParam = cms.double(0.7), patJetSrc = cms.untracked.InputTag("akPu7CalopatJets"))


akPu1PFCones.puPtMin = cms.double(5.0)
akPu2PFCones.puPtMin = cms.double(10.0)
akPu3PFCones.puPtMin = cms.double(15.0)
akPu4PFCones.puPtMin = cms.double(20.0)
akPu5PFCones.puPtMin = cms.double(25.0)
akPu6PFCones.puPtMin = cms.double(30.0)
#akPu7PFCones.puPtMin = cms.double(35.0)

akPu1CaloCones.puPtMin = cms.double(2.0)
akPu2CaloCones.puPtMin = cms.double(4.0)
akPu3CaloCones.puPtMin = cms.double(6.0)
akPu4CaloCones.puPtMin = cms.double(8.0)
akPu5CaloCones.puPtMin = cms.double(10.0)
akPu6CaloCones.puPtMin = cms.double(12.0)
#akPu7CaloCones.puPtMin = cms.double(14.0)

randomCones = cms.Sequence(
    akPu2PFCones+
    akPu3PFCones+
    akPu4PFCones+
    akPu5PFCones+
    akPu6PFCones+
    #akPu7PFCones+
    akPu2CaloCones+
    akPu3CaloCones+
    akPu4CaloCones+
    akPu5CaloCones+
    akPu6CaloCones
    #akPu7CaloCones
    )

randomCones2to5 = cms.Sequence(
    akPu2PFCones+
    akPu3PFCones+
    akPu4PFCones+
    akPu5PFCones+
    akPu2CaloCones+
    akPu3CaloCones+
    akPu4CaloCones+
    akPu5CaloCones
    )
