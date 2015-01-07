
#!/bin/bash
#rm -rf data/*.root
#rm -rf RescorTables

#hadd -f data/runs_181611_181670.root ../HiEvtPlaneFlatCalib/runs_181611_181670/res/*.root
#hadd -f data/runs_181671_181700.root ../HiEvtPlaneFlatCalib/runs_181671_181700/res/*.root
#hadd -f data/runs_181701_181800.root ../HiEvtPlaneFlatCalib/runs_181701_181800/res/*.root
#hadd -f data/runs_181801_182000.root ../HiEvtPlaneFlatCalib/runs_181801_182000/res/*.root
#hadd -f temp1 ../HiEvtPlaneFlatCalib/runs_182001_182200/res/*.root
#hadd -f temp2 ../HiEvtPlaneFlatCalib/runs_182201_182400/res/*.root
#hadd -f temp3 ../HiEvtPlaneFlatCalib/runs_182401_182600/res/*.root
#hadd -f temp4 ../HiEvtPlaneFlatCalib/runs_182601_182800/res/*.root
#hadd -f temp5 ../HiEvtPlaneFlatCalib/runs_182801_182950/res/*.root
#hadd -f temp6 ../HiEvtPlaneFlatCalib/runs_182951_183100/res/*.root
#hadd -f data/runs_182001_183100.root temp1 temp2 temp3 temp4 temp5 temp6
#rm -rf temp1
#rm -rf temp2
#rm -rf temp3
#rm -rf temp4
#rm -rf temp5
#rm -rf temp6
#hadd -f data/rpflat_combined.root data/runs*.root
#mkdir RescorTables
root -l -b -q rescor.C+
mv data/rpflat_combined.root data/rpflat_all.root

cp data/runs_181611_181670.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py
cp data/rpflat_combined.root ../HiEvtPlaneFlatCalib/runs_181611_181670.root
cp EPFlattening_HIRun2011_v5320devx01_offline.db ../HiEvtPlaneFlatCalib/EPFlattening_HIRun2011_v5320devx01_offline_1_181670.db
mv EPFlattening_HIRun2011_v5320devx01_offline.db  DBFiles/EPFlattening_HIRun2011_v5320devx01_offline_1_181670.db

cp data/runs_181671_181700.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py
cp data/rpflat_combined.root ../HiEvtPlaneFlatCalib/runs_181671_181700.root
cp EPFlattening_HIRun2011_v5320devx01_offline.db ../HiEvtPlaneFlatCalib/EPFlattening_HIRun2011_v5320devx01_offline_181671_181700.db
mv EPFlattening_HIRun2011_v5320devx01_offline.db  DBFiles/EPFlattening_HIRun2011_v5320devx01_offline_181671_181700.db

cp data/runs_181701_181800.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py
cp data/rpflat_combined.root ../HiEvtPlaneFlatCalib/runs_181701_181800.root
cp EPFlattening_HIRun2011_v5320devx01_offline.db ../HiEvtPlaneFlatCalib/EPFlattening_HIRun2011_v5320devx01_offline_181701_181800.db
mv EPFlattening_HIRun2011_v5320devx01_offline.db  DBFiles/EPFlattening_HIRun2011_v5320devx01_offline_181701_181800.db

cp data/runs_181801_182000.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py
cp data/rpflat_combined.root ../HiEvtPlaneFlatCalib/runs_181801_182000.root
cp EPFlattening_HIRun2011_v5320devx01_offline.db ../HiEvtPlaneFlatCalib/EPFlattening_HIRun2011_v5320devx01_offline_181801_182000.db
mv EPFlattening_HIRun2011_v5320devx01_offline.db  DBFiles/EPFlattening_HIRun2011_v5320devx01_offline_181801_182000.db

cp data/runs_182001_183100.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py
cp data/rpflat_combined.root ../HiEvtPlaneFlatCalib/runs_182001_183100.root
cp EPFlattening_HIRun2011_v5320devx01_offline.db ../HiEvtPlaneFlatCalib/EPFlattening_HIRun2011_v5320devx01_offline_182001_183100.db
mv EPFlattening_HIRun2011_v5320devx01_offline.db  DBFiles/EPFlattening_HIRun2011_v5320devx01_offline_182001_183100.db

cp data/rpflat_all.root data/rpflat_combined.root
