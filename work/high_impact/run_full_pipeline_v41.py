#!/usr/bin/env python3
import subprocess,sys,time,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
steps=[
 [sys.executable,'work/upgrade/run_meta_analysis_v2.py'],
 [sys.executable,'work/upgrade/run_pathway_robustness_v2.py'],
 [sys.executable,'work/high_impact/run_human_bone_localization.py'],
 [sys.executable,'work/high_impact/run_localization_study_unit_audit.py'],
 [sys.executable,'work/high_impact/run_localization_overlap_sensitivity.py'],
 [sys.executable,'work/high_impact/run_localization_covariate_robustness.py'],
 [sys.executable,'work/high_impact/run_independent_osteoblast_validation.py'],
 [sys.executable,'work/high_impact/run_itga11_integration.py'],
 [sys.executable,'work/high_impact/run_atac_validation.py'],
 [sys.executable,'work/high_impact/run_atac_bidirectional.py'],
 [sys.executable,'work/high_impact/run_atac_tpm_sensitivity.py'],
 [sys.executable,'work/high_impact/run_genetics_validation.py'],
 [sys.executable,'work/high_impact/run_representative_contrast_sensitivity.py'],
 ['Rscript','work/high_impact/preprocess_gse156508.R'],
 [sys.executable,'work/high_impact/run_osteoporosis_osteoblast_validation.py'],
 ['Rscript','work/high_impact/make_v4_figures.R'],
 [sys.executable,'work/high_impact/make_figure4.py']]
start_index=int(sys.argv[1]) if len(sys.argv)>1 else 0
log=[];start=time.time()
for cmd in steps[start_index:]:
 t=time.time();print('RUN',' '.join(cmd),flush=True)
 q=subprocess.run(cmd,cwd=ROOT)
 log.append({'command':cmd,'exit_code':q.returncode,'seconds':round(time.time()-t,2)})
 if q.returncode: raise SystemExit(q.returncode)
out=ROOT/'outputs/\u7814\u7a76\u5347\u7ea7/high_impact_v4/full_pipeline_run_log.json'
out.write_text(json.dumps({'completed':True,'start_index':start_index,'total_seconds':round(time.time()-start,2),'steps':log},indent=2))
print('WROTE',out)
