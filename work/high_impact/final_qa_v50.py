#!/usr/bin/env python3
from pathlib import Path
from docx import Document
from PIL import Image
import re,json,zipfile,hashlib,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs/研究升级/Communications_Biology_v4';errs=[]
def ck(x,msg):
 if not x:errs.append(msg)
d=Document(O/'Manuscript_Communications_Biology.docx');ps=[p.text for p in d.paragraphs]
title=ps[0];abstract=ps[4];ck(len(re.findall(r"\b[\w-]+\b",title))<=15,'title >15 words');ck(len(re.findall(r"\b[\w-]+\b",abstract))<150,'abstract >=150 words')
ck('Monte Carlo P=0.259' in '\n'.join(ps),'combined magnitude missing');ck('12,757' in '\n'.join(ps),'human gene count missing');ck('rho=0.810' in '\n'.join(ps),'human concordance missing')
for h in ['Results','Discussion','Methods','Data availability','Code availability','References']:ck(ps.count(h)==1,f'heading {h} count {ps.count(h)}')
# placeholders are intentionally limited to author-supplied information
text='\n'.join(ps);found=re.findall(r'\[[^\]]*(?:AUTHOR|VERIFY)[^\]]*\]',text)
ck(len(found)==0,f'unresolved manuscript placeholder count {len(found)}')
ck('10.5281/zenodo.22726636' in text,'archived release DOI missing')
ck('[AUTHOR TO PROVIDE REPOSITORY URL AND DOI]' not in text,'repository DOI placeholder remains')
for expected in ['Tingting Tan','Xia Chen','Zhengxiao Ouyang','These authors contributed equally','chenxiacxv99@csu.edu.cn','ouyangzhengxiao@csu.edu.cn','0000-0002-8997-0446','82371600','82571826','2026JJ20072','2025JJ50614']:
 ck(expected in text,f'manuscript metadata missing: {expected}')
all_docx_text=''
for n in ['Manuscript_Communications_Biology','Supplementary_Information_Communications_Biology','Cover_Letter_Communications_Biology','PRISMA_2020_Checklist','Evidence_Ceiling_Audit_v5']:
 x=Document(O/f'{n}.docx')
 all_docx_text+='\n'.join(p.text for p in x.paragraphs)+'\n'+'\n'.join(p.text for t in x.tables for row in t.rows for c in row.cells for p in c.paragraphs)
ck(not re.search(r'\[(?:AUTHOR|CORRESPONDING|VERIFY|COMPLETE|THE AUTHORS|ACKNOWLEDGEMENTS|DATE|CONFIRM)',all_docx_text),'unresolved submission placeholder remains')
mc=json.loads((R/'outputs/研究升级/high_impact_v4/multicohort_osteoporosis_validation/multicohort_osteoporosis_validation_summary.json').read_text())
hu=json.loads((R/'outputs/研究升级/high_impact_v4/human_only_sensitivity/human_only_sensitivity_summary.json').read_text())
ck(abs(mc['combined']['abs']['monte_carlo_p']-.25886948226103546)<1e-12,'MC mismatch');ck(hu['genes_fdr']==0 and hu['genes_k_ge_3']==12757,'human summary mismatch')
for i in range(1,5):
 for ext in ['png','pdf','tif']:ck((O/f'Figure_{i}.{ext}').exists(),f'missing Figure_{i}.{ext}')
 im=Image.open(O/f'Figure_{i}.png');ck(im.width>=1800,f'Figure_{i}.png low width {im.width}')
for n in ['Manuscript_Communications_Biology','Supplementary_Information_Communications_Biology','Cover_Letter_Communications_Biology','PRISMA_2020_Checklist','Evidence_Ceiling_Audit_v5']:
 ck((O/f'{n}.docx').exists(),f'missing {n}.docx');ck((O/f'{n}.pdf').exists(),f'missing {n}.pdf')
for zname in ['S1_Derived_Data_and_Code.zip','Communications_Biology_v4_投稿包.zip']:
 with zipfile.ZipFile(O/zname) as z:bad=z.testzip();ck(bad is None,f'{zname} corrupt {bad}')
with zipfile.ZipFile(O/'S1_Derived_Data_and_Code.zip') as z:
 names=z.namelist();ck(any(x.endswith('human_only_sensitivity_summary.json') for x in names),'human result absent S1');ck(any(x.endswith('multicohort_osteoporosis_validation_protocol_v1.md') for x in names),'protocol absent S1')
for line in (O/'SHA256SUMS.txt').read_text().splitlines():
 h,n=line.split('  ',1);ck(hashlib.sha256((O/n).read_bytes()).hexdigest()==h,f'checksum mismatch {n}')
report={'status':'PASS' if not errs else 'FAIL','title_words':len(title.split()),'abstract_words':len(re.findall(r"\b[\w-]+\b",abstract)),'manuscript_paragraphs':len(ps),'intentional_author_placeholders':found,'package_files':len(zipfile.ZipFile(O/'Communications_Biology_v4_投稿包.zip').namelist()),'S1_entries':len(zipfile.ZipFile(O/'S1_Derived_Data_and_Code.zip').namelist()),'errors':errs}
(O/'final_QA_v5.0.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if errs:raise SystemExit(1)
