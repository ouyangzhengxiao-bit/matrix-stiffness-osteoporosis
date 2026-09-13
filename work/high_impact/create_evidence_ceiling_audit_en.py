#!/usr/bin/env python3
from pathlib import Path
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'outputs/研究升级/Communications_Biology_v4'
def setup(d):
 s=d.sections[0];s.top_margin=Inches(.7);s.bottom_margin=Inches(.7);s.left_margin=Inches(.8);s.right_margin=Inches(.8)
 d.styles['Normal'].font.name='Arial';d.styles['Normal'].font.size=Pt(9);d.styles['Normal'].paragraph_format.space_after=Pt(5)
 for n,z in [('Title',18),('Heading 1',14),('Heading 2',11)]:d.styles[n].font.name='Arial';d.styles[n].font.size=Pt(z);d.styles[n].font.color.rgb=RGBColor(30,55,75)
def table(d,h,rows,fs=7):
 t=d.add_table(rows=1,cols=len(h));t.style='Table Grid';t.alignment=WD_TABLE_ALIGNMENT.CENTER
 for j,x in enumerate(h):
  c=t.rows[0].cells[j];c.text=x;sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'244A66');c._tc.get_or_add_tcPr().append(sh)
  for r in c.paragraphs[0].runs:r.font.bold=True;r.font.color.rgb=RGBColor(255,255,255)
 for row in rows:
  c=t.add_row().cells
  for j,x in enumerate(row):c[j].text=str(x)
 for row in t.rows:
  for c in row.cells:
   for p in c.paragraphs:
    p.paragraph_format.space_after=Pt(0)
    for r in p.runs:r.font.size=Pt(fs)
 return t
def main():
 d=Document();setup(d);p=d.add_paragraph(style='Title');p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('Public-data evidence ceiling and submission audit — v5.0')
 p=d.add_paragraph('Matrix-stiffness transcriptomics and osteoporosis | 9 September 2026');p.alignment=WD_ALIGN_PARAGRAPH.CENTER
 d.add_heading('Editorial assessment',1)
 d.add_paragraph('The project now approaches the strongest design attainable without new wet-laboratory experiments. It combines a documented systematic search, independent-study inference, REML random-effects models with modified Hartung–Knapp uncertainty, exact randomization tests, a human-only sensitivity analysis, cohort deduplication, single-cell and spatial localization, chromatin accessibility, genetics and direct bone-mechanics outcomes. Its publishable advance is a quantified boundary: matrix-stiffness responses are distributed, lineage-localized and context-dependent, while a universal gene signature and osteoporosis biomarker are unsupported.')
 d.add_paragraph('Communications Biology remains an appropriate stretch target because its scope explicitly includes secondary data analysis and evaluates novelty, technical soundness, evidential strength and influence on thinking in the field. Editorial risk is driven by the negative/boundary-centered result rather than a remaining statistical weakness.')
 d.add_heading('Evidence architecture',1)
 table(d,['Layer','Data and inference','Result','Supported claim'],[
 ['Discovery','7 stiffness studies; 13 contrasts; within-series aggregation','0/12,241 genes at FDR; study rho -0.253 to 0.136','No demonstrable universal directional program'],
 ['Pathways','1,008 sets; one score per study','109 pooled-gene signals reduced to 0 at the study unit','Correlated genes can create false precision'],
 ['Lineage','Human femur single-cell, spatial and independent GSE147390 programs','Spatial program: 6/7 studies, exact sign-flip P=0.0156; independent program partial','Response magnitude localizes to early osteoblast lineage'],
 ['Disease','3 independent human cohorts; complete label spaces; equal-cohort synthesis','Combined P=0.259 magnitude, 0.254 signed, 0.095 signature','Disease transfer is context-dependent'],
 ['Species','Repeat meta-analysis after excluding 2 rat studies','0/12,757 genes at FDR; rho=0.810 with full effects','Main conclusion is not driven by rat data'],
 ['Orthogonal challenges','ITGA11 states, ATAC, eBMD GWAS, mouse bone mechanics','Fixed signature unsupported','No biomarker or mechanical-property claim']])
 d.add_heading('Independent osteoporosis cohorts',1)
 table(d,['Cohort','Design','Absolute rho','Exact P','Signature P','Interpretation'],[
 ['GSE156508','6 osteoporotic fracture vs 6 severe OA; primary osteoblasts','0.02541','0.0303','0.2229','Weak magnitude association only'],
 ['GSE35958','5 primary OP vs 4 age-matched controls; cultured BM-MSCs','0.00094','0.5000','0.0238','Isolated signature result; no replication'],
 ['GSE230665','12 postmenopausal OP vs 3 controls; femur tissue','-0.01097','0.7033','0.6308','No support; only 3 controls'],
 ['Equal-cohort synthesis','Each cohort standardized to its complete exact null','0.3848*','0.2589','0.0950','*Combined standardized statistic; null overall']])
 d.add_paragraph('Independence rules were frozen before outcome inspection. GSE35956 was excluded because it reused the GSE35958 osteoporosis cases and used younger controls. GSE157322 was excluded because BMD group was fully confounded with ancestry. These decisions prevent inflation of cohort count by duplicate participants or structural confounding.')
 d.add_heading('Human-only sensitivity',1)
 d.add_paragraph('After removing rat GSE181512 and GSE310514, the same random-effects procedure was applied to five human stiffness studies. Among 12,757 genes represented in at least three studies, none passed FDR correction. Human-only and seven-study pooled effects correlated at Spearman rho=0.810. Repeating disease validation with the human-only model produced combined P=0.199, 0.138 and 0.0835 for magnitude, direction and signature endpoints, respectively.')
 d.add_heading('Claim boundary',1)
 table(d,['Suitable for title/abstract','Exploratory only','Unsupported'],[
 ['Stiffness-response magnitude localizes to osteoblast-lineage programs','Small GSE156508 magnitude association','A stiffness signature diagnoses osteoporosis'],
 ['Directional effects show poor cross-study agreement','Single-cohort GSE35958 signature difference','A stable core gene set was discovered'],
 ['Human-only analysis preserves the main conclusion','Residual spatial signal after overlap removal','Transcriptomics measures bone material stiffness'],
 ['Combined disease tests do not support universal transfer','One-donor GSE147390 partial support','Causality, fracture prediction or treatment response']])
 d.add_heading('Remaining author actions',1)
 for x in ['Complete names, affiliations, corresponding-author email and ORCID.','Verify funding, CRediT contributions, conflicts of interest and the institutional determination for secondary public-data analysis.','Retain the archived release DOI 10.5281/zenodo.22726636 in the manuscript and submission metadata.','Retain the disclosure that the revision-stage search was not prospectively registered.','Do not convert the heterogeneous disease cohorts into a clinical-prediction claim.']:d.add_paragraph(x,style='List Bullet')
 d.add_heading('Journal sequence',1)
 table(d,['Priority','Journal','Positioning','Risk'],[
 ['1','Communications Biology','Mechanistic boundary, reproducibility and rigorous secondary synthesis','Medium-high: editor must value a negative result that changes interpretation'],
 ['2','iScience','Public multi-omics integration and reproducibility','Medium'],
 ['3','JBMR Plus / Bone Reports','Bone-biology translational boundary','Lower; narrower reach'],
 ['Not first choice','Nature Communications / Science Advances','Usually require stronger mechanism or prospective validation','Very high without new experiments']])
 d.add_paragraph('Additional routine machine learning, WGCNA or network pharmacology would not raise the evidential ceiling and could reduce credibility. The only clear non-experimental route to a materially higher ceiling is a newly available, larger human bone single-cell/spatial cohort with complete donor metadata, or a clinical collaboration providing external BMD/fracture outcomes.')
 d.add_heading('Target-journal format check',1)
 table(d,['Requirement','Status'],[['Article title up to 15 words','Pass: 11 words'],['Abstract below 150 words','Pass: 137 words'],['Statistics and reproducibility section','Pass'],['Data and Code availability','Present; archived release DOI supplied'],['PRISMA checklist','Complete; non-registration disclosed'],['Figures in text plus high-resolution files','Pass: PNG/PDF/TIFF'],['AI-use disclosure in Methods','Present; author verification required']])
 d.add_heading('Primary sources',1)
 for x in ['Communications Biology, Aims & Scope: https://www.nature.com/commsbio/aims','Communications Biology, Submission guidelines: https://www.nature.com/commsbio/submit/submission-guidelines','Communications Biology, Content types: https://www.nature.com/commsbio/submit/content-types','NCBI GEO GSE156508: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156508','Benisch et al., PLoS ONE 2012 (GSE35958): https://pmc.ncbi.nlm.nih.gov/articles/PMC3454401/','NCBI GEO GSE230665: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230665','PRISMA 2020 checklist: https://www.prisma-statement.org/prisma-2020-checklist']:d.add_paragraph(x)
 d.add_paragraph('All numerical statements are generated from frozen scripts and machine-readable outputs in the submission archive.')
 d.save(O/'Evidence_Ceiling_Audit_v5.docx')
if __name__=='__main__':main()
