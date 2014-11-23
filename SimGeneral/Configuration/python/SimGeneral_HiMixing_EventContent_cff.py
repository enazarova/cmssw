import FWCore.ParameterSet.Config as cms

HiMixRAW = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_mixGen_*_*',
        'keep *_hiGenParticles_*_*'
    )
)

HiMixRECO = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_mixGen_*_*',
	'keep *_hiGenParticles_*_*'
    )
)

HiMixAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
	'keep *_hiGenParticles_*_*'
    )
)

