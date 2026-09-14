#!/usr/bin/env python3
from pathlib import Path
import json, re
import pandas as pd
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'outputs/\u7814\u7a76\u5347\u7ea7'
V4=BASE/'high_impact_v4'
OUT=BASE/'Communications_Biology_v4'
BLACK=RGBColor(0,0,0); RED=RGBColor(170,0,0); BLUE='244A66'

def style_doc(d):
    sec=d.sections[0]; sec.top_margin=Inches(.72); sec.bottom_margin=Inches(.72); sec.left_margin=Inches(.82); sec.right_margin=Inches(.82)
    st=d.styles['Normal']; st.font.name='Arial'; st.font.size=Pt(9.5); st.font.color.rgb=BLACK
    st.paragraph_format.space_after=Pt(4); st.paragraph_format.line_spacing=1.08
    for name,size in [('Title',17),('Heading 1',13),('Heading 2',10.5),('Heading 3',9.7)]:
        s=d.styles[name]; s.font.name='Arial'; s.font.size=Pt(size); s.font.bold=True; s.font.color.rgb=BLACK
        s.paragraph_format.space_before=Pt(8); s.paragraph_format.space_after=Pt(4); s.paragraph_format.keep_with_next=True

def replace_start(d,start,text):
    hits=[p for p in d.paragraphs if p.text.startswith(start)]
    if len(hits)!=1: raise RuntimeError(f'replace {start!r}: {len(hits)} hits')
    hits[0].text=text
    return hits[0]

def insert_before(d,start,blocks):
    hits=[p for p in d.paragraphs if p.text==start]
    if not hits: hits=[p for p in d.paragraphs if p.text.startswith(start)]
    if len(hits)!=1: raise RuntimeError(f'insert {start!r}: {len(hits)} hits')
    target=hits[0]
    for style,text in blocks:
        p=d.add_paragraph(text,style=style)
        target._p.addprevious(p._p)

def add_table(d,headers,rows,font=7.0):
    t=d.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.style='Table Grid'
    for j,h in enumerate(headers):
        c=t.rows[0].cells[j];c.text=str(h)
        shd=OxmlElement('w:shd');shd.set(qn('w:fill'),BLUE);c._tc.get_or_add_tcPr().append(shd)
        for r in c.paragraphs[0].runs:r.font.bold=True;r.font.color.rgb=RGBColor(255,255,255);r.font.size=Pt(font)
    for row in rows:
        cells=t.add_row().cells
        for j,v in enumerate(row):
            cells[j].text='' if pd.isna(v) else str(v)
            for p in cells[j].paragraphs:
                p.paragraph_format.space_after=Pt(0)
                for r in p.runs:r.font.size=Pt(font)
    return t

def flow_figure():
    fig,ax=plt.subplots(figsize=(6.7,6.8)); ax.set_xlim(0,10);ax.set_ylim(0,12);ax.axis('off')
    boxes=[
      (5,11,'Records identified through database searches\\nGEO n = 42; PubMed n = 3'),
      (5,9.3,'Additional record from citation/accession checking\\nn = 1'),
      (5,7.6,'Records screened\\nn = 46'),
      (8.2,7.6,'Excluded at title/summary\\nn = 34'),
      (5,5.7,'Full records assessed for eligibility\\nn = 12'),
      (8.2,5.7,'Excluded after full assessment\\nn = 5\\nOne donor (1); split series (2);\\nglass only (1); strain/vibration (1)'),
      (5,3.6,'Studies included in quantitative synthesis\\nn = 7'),
      (5,1.7,'Independent validation source\\nGSE147390; one donor; n = 5,329 osteoblasts')
    ]
    for x,y,txt in boxes:
        ax.text(x,y,txt,ha='center',va='center',fontsize=8.4,
                bbox=dict(boxstyle='round,pad=.45',fc='white',ec='#244A66',lw=1.2))
    def arrow(x1,y1,x2,y2): ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='->',lw=1,color='#555'))
    arrow(5,10.55,5,9.75);arrow(5,8.85,5,8.05);arrow(5,7.15,5,6.15);arrow(5,5.25,5,4.05);arrow(5,3.15,5,2.15)
    arrow(6.2,7.6,7.1,7.6);arrow(6.2,5.7,7.1,5.7)
    ax.text(.2,11.85,'Supplementary Fig. 1 | Public-data identification and eligibility flow',fontsize=10,fontweight='bold',va='top')
    fig.tight_layout();fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.png',dpi=320,bbox_inches='tight');fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.tif',dpi=320,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});fig.savefig(OUT/'Supplementary_Figure_1_PRISMA.pdf',bbox_inches='tight');plt.close(fig)

def manuscript():
    p=OUT/'Manuscript_Communications_Biology.docx';d=Document(p);style_doc(d)
    abstract=('Matrix stiffness regulates mesenchymal stromal cell behavior, but the extent to which transcriptomic effects generalize across experiments and connect to osteoporosis is unclear. We systematically identified seven eligible public stiffness studies comprising 13 contrasts and 12,241 genes, while preserving each GEO series as the unit of inference. No gene passed false-discovery-rate correction; study-level gene-effect correlations ranged from -0.253 to 0.136; and 109 pooled-gene pathway signals fell to zero under study-unit inference. Response magnitude nevertheless localized to osteoblast-lineage programs from human femoral single-cell and spatial data. The spatial program remained supported by exact study-level sign-flip inference after stricter matching for coverage, uncertainty and gene length (6/7 studies; P=0.0156). A truly independent human osteoblast program from GSE147390 showed partial support (gene-level matched P=0.0392; 6/7 studies; exact sign-flip P=0.0391, exact sign-count P=0.0625), with the important limitation that its source comprised one donor. A fixed 200-gene signature did not enrich bone-mineral-density associations, identify an ITGA11-high osteogenic state, or show specific accessibility reversal after mechanical interventions. Public data therefore support distributed, lineage-localized mechanosensitivity, while arguing against a universal stiffness gene signature or an osteoporosis biomarker.')
    replace_start(d,'Matrix stiffness regulates mesenchymal stromal cell behavior, yet',abstract)
    replace_start(d,'We tested these alternatives','We tested these alternatives in a systematic public-data synthesis of seven MSC stiffness datasets. We asked whether genes or pathways transferred across independent studies, localized response magnitude in human femoral cell programs, and then tested a fully independent osteoblast program from GSE147390. Orthogonal donor-resolved osteogenic states, mechanical-intervention ATAC-seq and large-scale bone-density genetics were used to challenge a fixed stiffness signature. The analysis separates lineage localization from gene-level transfer and retains negative and borderline results as constraints.')
    replace_start(d,'A prespecified replication used','A prespecified second tissue representation used human bone spatial programs from the same multiscale study as the single-cell atlas [9]. The union of MSC, Osteo-MSC and osteoblast programs was enriched (184 mapped genes; standardized matched effect 2.91; permutation P=0.00260; Fig. 2b). The three spatial states did not survive correction separately. The single-cell and spatial unions shared 125 genes; after removal, 64 mapped spatial genes retained a positive, non-significant effect (P=0.107). This is cross-modality corroboration within one study, rather than independent cohort replication.')
    replace_start(d,'To prevent correlated genes','To prevent correlated genes from serving as replicates, each GEO series contributed one matched enrichment effect. Under the original coverage- and variance-matched analysis, the single-cell program was positive in 5/7 studies (mean 0.1319; 95% t interval -0.0225 to 0.2864; exact one-sided sign-flip P=0.0469; sign-count P=0.2266; t-test sensitivity P=0.0408). The spatial program was positive in 7/7 studies (mean 0.2552; 95% t interval 0.0653 to 0.4450; exact sign-flip P=0.00781; sign-count P=0.00781; t-test sensitivity P=0.00832; Fig. 2c). All leave-one-study-out mean effects remained positive.')
    insert_before(d,'A fixed gene signature does not transfer',[
      ('Heading 2','Independent human osteoblast markers provide partial external support'),
      ('Normal','We froze eligibility and matching rules before merging a human osteoblast marker program from GSE147390, a dataset outside the multiscale atlas [17]. Of 2,235 source markers, 2,074 mapped to the stiffness meta-analysis. Their mean absolute meta z score exceeded sets matched for study coverage, standard error and gene length (standardized effect 1.76; 20,000-permutation P=0.0392). Study-level effects were positive in 6/7 series (mean 0.0511; 95% t interval -0.00520 to 0.1073; exact sign-flip P=0.0391; exact sign-count P=0.0625; t-test sensitivity P=0.0341), and every leave-one-study-out mean remained positive. Because strict sign-count inference did not cross 0.05 and the source atlas contained one donor, we classify this result as partial independent support rather than definitive replication.'),
      ('Heading 2','Post-review sensitivity analyses bound the localization claim'),
      ('Normal','A dated post hoc analysis added GENCODE v49 gene-length deciles to the original matching variables and corrected the two union tests by Bonferroni. Gene-level enrichment remained for the single-cell union (1,004 genes; adjusted P=0.00210) and spatial union (184 genes; adjusted P=0.00250). At the study level, the spatial union remained supported (6/7 positive; exact sign-flip P=0.0156), whereas the single-cell union became borderline (5/7; P=0.0547). Separately, retaining only one representative contrast per GEO series produced no FDR-significant gene and preserved the ordering of pooled effects (Spearman rho=0.814 versus the primary analysis).')
    ])
    replace_start(d,'Figure 3  Orthogonal datasets','Figure 3  Independent and orthogonal tests constrain transferability. a, One gene-length- and variance-matched effect per stiffness study for an independent human osteoblast marker program from GSE147390. b, Author-normalized ATAC signals reproduce global promoter shifts without separating the fixed positive and negative stiffness gene sets. c, Observed and expected numbers of eBMD-significant genes in the fixed signature.')
    replace_start(d,'The localization was strongest','The localization was strongest in two MSC states and pre-osteoblasts, whereas mature osteoblasts were not enriched. The spatial program was positive across studies, but it came from the same multiscale publication as the single-cell atlas and overlap removal weakened gene-level evidence. The independent GSE147390 program added a second source and retained a positive exact sign-flip result, while its sign-count P value and one-donor design prevent a claim of definitive replication. Gene-length matching further reduced the single-cell study-unit result to borderline evidence. Together, these results support a shared lineage signal with measured uncertainty.')
    replace_start(d,'This study is limited','This study is limited by small source studies, heterogeneous processed data, incomplete reporting of matrix chemistry and biological replication, and only seven independent series. Moderator analysis was therefore underpowered. Multi-contrast series required within-series aggregation, although a one-contrast-per-series sensitivity gave the same gene-level conclusion. The GSE147390 marker source contained one donor. The ATAC perturbation contained two biological samples per group. Single-cell and spatial programs from Chai et al. represent distinct datasets within one publication rather than independent replication. Localization, TPM and gene-length analyses added after the initial synthesis are labelled by timing. None of the integrated analyses establish causality or clinical prediction.')
    replace_start(d,'A protocol was frozen','The gene-level synthesis protocol was frozen before analysis. A revision-stage systematic search was run on 9 September 2026 and was not prospectively registered. Eligible studies imposed at least two extracellular-matrix stiffness levels on MSCs or mesenchymal progenitors, included at least two biological samples per condition and provided genome-wide expression data. One GEO series was one independent study. Times, coatings, passages and library preparations within a series were dependent contrasts. Later human-bone, ATAC and genetics protocols were dated before their corresponding dataset merges; subsequent sensitivities are explicitly labelled post hoc.')
    insert_before(d,'Expression processing and within-study effects',[
      ('Heading 2','Systematic search and screening'),
      ('Normal','NCBI GEO DataSets was searched on 9 September 2026 with: gse[ETYP] AND (stiffness[All Fields] OR "matrix elasticity"[All Fields] OR "substrate stiffness"[All Fields] OR "mechanical memory"[All Fields]) AND (mesenchymal[All Fields] OR BMSC[All Fields]). PubMed was searched with: (("substrate stiffness"[Title/Abstract] OR "matrix stiffness"[Title/Abstract] OR "matrix elasticity"[Title/Abstract]) AND ("mesenchymal stromal cell"[Title/Abstract] OR "mesenchymal stem cell"[Title/Abstract] OR BMSC[Title/Abstract]) AND (transcriptom*[Title/Abstract] OR RNA-seq[Title/Abstract] OR RNAseq[Title/Abstract] OR microarray[Title/Abstract])). Citation and known-accession checking was added. Forty-six records were screened, 12 received full assessment, five were excluded after full assessment and seven entered synthesis (Supplementary Fig. 1; screening log in the archive).')
    ])
    replace_start(d,'Genes measured in at least four','Genes measured in at least four studies entered random-effects meta-analysis. Between-study variance used restricted maximum likelihood. Two-sided confidence intervals and P values used modified Hartung-Knapp uncertainty with scale not below one and a t distribution with k-1 degrees of freedom. Benjamini-Hochberg correction covered 12,241 genes. Hallmark and Reactome sets from MSigDB 2025.1 were restricted to 15-500 genes. Pooled-gene rank tests were followed by one score per GEO series. Study-level pathway t tests were two-sided sensitivities, with FDR controlled across 1,008 pathways. A post hoc aggregation sensitivity retained the first ordered eligible contrast for each series.')
    replace_start(d,'For study-unit audit','For the study-unit audit, each GEO series contributed one difference between the mean absolute study effect in a fixed program and its matched null mean. The primary post-review inference was the exact one-sided sign-flip test over all 27 study-level sign assignments. Exact one-sided binomial sign counts and one-sided one-sample t tests were robustness summaries; 95% t intervals were descriptive. No normality assumption was required for permutation, sign-flip or sign-count tests. A dated post hoc sensitivity added gene-length deciles to the matching variables and used Bonferroni correction across the two lineage unions.')
    insert_before(d,'Mechanical-intervention ATAC-seq',[
      ('Heading 2','Independent osteoblast validation'),
      ('Normal','GSE147390 markers were taken from Gong et al. [17]. Source-table genes with adjusted P<0.05 and average fold change >1.2 were fixed before merger. Mean absolute meta z score was compared with 20,000 random sets matched by number of contributing studies, Hartung-Knapp standard-error decile and GENCODE v49 gene-length decile. At the study level, absolute effects were matched by sampling-variance and gene-length deciles. Exact sign-flip inference, exact sign count, t sensitivity and leave-one-study-out analyses followed the study-unit framework.')
    ])
    insert_before(d,'Software and reproducibility',[
      ('Heading 2','Statistics and reproducibility'),
      ('Normal','All hypothesis tests, sidedness, sample units and multiplicity corrections are specified above and in the Supplementary Information. The nominal alpha level was 0.05. Empirical tests used 20,000 random sets and reported (1 + exceedances)/(20,001). Exact sign-flip and sign-count tests enumerated their complete null spaces. Gene tests and each defined family of cluster or pathway tests used Benjamini-Hochberg FDR; the two post hoc gene-length union tests used Bonferroni correction. P values are reported exactly where space permits and otherwise to at least two significant digits. Biological replicate counts were inherited from source datasets and no observations were excluded after outcome inspection.')
    ])
    replace_start(d,'Analyses used Python','Analyses used Python 3.14.2 with pandas, NumPy, SciPy, statsmodels, openpyxl and pyBigWig, and R 4.5.2 with ggplot2, patchwork, scales and jsonlite. Random-set procedures used dated fixed seeds recorded in each script. The principal pipeline was rerun after dependency repair; protocols, amendments, search exports, screening decisions, environment records, derived results and executable scripts accompany the submission. Public source files are referenced by accession rather than redistributed.')
    replace_start(d,'Primary transcriptomic data are available','Primary transcriptomic data are available from NCBI GEO under GSE193021, GSE181512, GSE226411, GSE288678, GSE55867, GSE255574 and GSE310514. Independent osteoblast validation used GSE147390. External datasets are GSE166824, GSE22011, GSE152708, GSE35958, GSE230665, GSE276529, GSE317531 and GSE287556. Human femoral single-cell data are GSE317069, spatial bone data GSE299207 and the eBMD GWAS GWAS Catalog GCST90726625 [9]. Search exports, screening log, derived data and figure source tables accompany this submission.')
    replace_start(d,'All custom scripts required','All custom scripts needed to regenerate the reported analyses and figures, dated protocols, derived numerical results and environment records are archived in Zenodo release v1.0.0 at https://doi.org/10.5281/zenodo.22726636 and maintained at https://github.com/ouyangzhengxiao-bit/matrix-stiffness-osteoporosis. The version DOI identifies the exact code and data snapshot used for this manuscript; the concept DOI https://doi.org/10.5281/zenodo.22726635 resolves to the latest archived release.')
    insert_before(d,'Acknowledgements',[
      ('Normal','17. Gong, Y. et al. Aging of human bone marrow stromal cells during osteoblast differentiation induces an osteoporotic phenotype. Aging 13, 17646-17666 (2021).'),
      ('Normal','18. Page, M. J. et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ 372, n71 (2021).'),
      ('Normal','19. Knapp, G. & Hartung, J. Improved tests for a random effects meta-regression with a single covariate. Stat. Med. 22, 2693-2710 (2003).')
    ])
    # Put disclosure in Methods as required by target journal and retain author verification.
    insert_before(d,'Data availability',[
      ('Heading 2','Use of generative AI'),
      ('Normal','OpenAI Codex assisted with code review, language editing, reproducibility auditing and document preparation. The authors reviewed and verified all analyses, numerical results, citations, figures and manuscript content and take full responsibility for the work. No generative AI imagery is included.')
    ])
    # Remove duplicate terminal disclosure section.
    pars=d.paragraphs
    for i,p0 in enumerate(pars):
        if p0.text=='Generative AI disclosure':
            p0._element.getparent().remove(p0._element)
            if i+1<len(pars) and pars[i+1].text.startswith('[VERIFY AND EDIT:'):
                pars[i+1]._element.getparent().remove(pars[i+1]._element)
            break
    d.save(p)
    print('manuscript abstract words',len(re.findall(r"\b[\w-]+\b",abstract)))

def supplement():
    p=OUT/'Supplementary_Information_Communications_Biology.docx';d=Document(p);style_doc(d)
    d.add_page_break();d.add_heading('Supplementary systematic-search record',1)
    q=d.add_paragraph();q.alignment=WD_ALIGN_PARAGRAPH.CENTER;q.add_run().add_picture(str(OUT/'Supplementary_Figure_1_PRISMA.png'),width=Inches(5.7))
    cp=d.add_paragraph('Supplementary Figure 1  Public-data identification and eligibility flow. Searches were performed on 9 September 2026. No automation or machine-learning prioritization was used for screening.');cp.runs[0].bold=True
    add_table(d,['stage','n'],[['GEO records',42],['PubMed records',3],['citation/accession records',1],['screened',46],['full assessment',12],['full-record exclusions',5],['quantitative studies',7]],7.5)
    d.add_paragraph('The review-stage search was not prospectively registered. Exact query strings, machine-readable NCBI exports and a decision with reason for every record are included in search_audit/. GSE288423 and GSE288465 were excluded because stiff and soft conditions were deposited as separate, independently normalized series, making condition inseparable from series processing.')
    d.add_heading('Supplementary validation and sensitivity analyses',1)
    ind=pd.read_csv(V4/'independent_osteoblast_validation/independent_osteoblast_study_effects.csv')
    rows=[[r.study,int(r.n_program_genes),f'{r.raw_matched_effect:.4f}',f'{r.standardized_effect:.3f}'] for r in ind.itertuples()]
    add_table(d,['study','mapped genes','raw matched effect','standardized effect'],rows,7)
    d.add_paragraph('Supplementary Table 11  Independent GSE147390 osteoblast-marker validation by stiffness study. Six of seven effects were positive; exact one-sided sign-flip P=0.0391, sign-count P=0.0625. The source marker atlas contained one donor.')
    cov=pd.read_csv(V4/'human_bone_localization/covariate_matched_study_summary.csv')
    rows=[[r.program,int(r.n_studies),int(r.positive_studies),f'{r.mean_raw_effect:.4f}',f'{r.exact_signflip_p:.5f}',f'{r.exact_sign_p:.5f}',f'{r.sensitivity_t_p:.5f}'] for r in cov.itertuples()]
    add_table(d,['program','k','positive','mean','sign-flip P','sign-count P','t sensitivity P'],rows,6.6)
    d.add_paragraph('Supplementary Table 12  Post hoc study-unit robustness after simultaneous variance and gene-length matching. All leave-one-study-out mean effects remained positive.')
    rep=json.loads((V4/'representative_contrast_sensitivity/representative_contrast_summary.json').read_text())
    add_table(d,['analysis','result'],[
      ['Studies',rep['independent_studies']],['Genes analyzed',rep['genes_meta_analyzed']],['FDR-significant genes',rep['genes_fdr']],
      ['Minimum q',f"{rep['minimum_q']:.6f}"],['Spearman rho vs all-contrast analysis',f"{rep['pooled_effect_spearman_vs_all_contrasts']:.4f}"]],7)
    d.add_paragraph('Supplementary Table 13  Representative-contrast sensitivity. The first ordered eligible comparison was retained for multi-contrast series; choices are recorded in the JSON summary.')
    report=pd.read_csv(V4/'search_audit/source_study_reporting_assessment.csv')
    d.add_heading('Source-study reporting assessment',2)
    add_table(d,report.columns,report.fillna('').astype(str).values.tolist(),5.8)
    d.add_paragraph('Supplementary Table 14  Reporting assessment based on public records. Unclear denotes information not recoverable from the analyzed public files and was not imputed.')
    d.add_heading('Interpretive safeguards',2)
    d.add_paragraph('The spatial and single-cell programs from Chai et al. were treated as separate tissue representations within one multiscale publication, not independent replication. GSE147390 supplied the independent biological source but had one donor. Exact sign-flip inference was primary after review because seven studies provide limited support for normal approximations. The stricter gene-length analysis was post hoc and weakened the single-cell study-unit result to P=0.0547; this is reported as borderline rather than positive.')
    d.save(p)

def cover():
    p=OUT/'Cover_Letter_Communications_Biology.docx';d=Document(p);style_doc(d)
    replace_start(d,'Matrix stiffness is widely','Matrix stiffness is widely treated as a source of transferable mesenchymal stromal-cell programs. Across seven public studies, no gene or pathway survived independent-study inference. Response magnitude nevertheless localized to early osteoblast-lineage programs. A separate spatial dataset from the same multiscale study supported this localization, and a truly independent human osteoblast program from GSE147390 provided partial external support under exact sign-flip inference. We explicitly report its borderline sign-count result and one-donor limitation.')
    replace_start(d,'We challenged that interpretation','We challenged that interpretation using stricter gene-length matching, one-representative-contrast analysis, leave-one-study-out checks, donor-resolved ITGA11 osteogenic cultures, bidirectional mechanical-intervention ATAC-seq and eBMD gene statistics from 448,010 participants. These tests preserve the lineage-level signal while rejecting a universal 200-gene signature or osteoporosis biomarker claim.')
    replace_start(d,'The manuscript is suited','The manuscript is suited to Communications Biology because it identifies a reproducible biological organization while defining its limits. A revision-stage systematic search, full screening log, PRISMA 2020 checklist, exact study-unit inference, dated protocols, executable code, derived data and figure source tables accompany the submission.')
    d.save(p)

def prisma():
    d=Document();style_doc(d)
    p=d.add_paragraph(style='Title');p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('PRISMA 2020 Checklist')
    d.add_paragraph('Manuscript: Matrix stiffness responses localize to human osteoblast lineage programs without a universal gene signature')
    d.add_paragraph('Adapted for a revision-stage systematic search and secondary quantitative synthesis of public functional-genomics datasets. Completed 9 September 2026. PRISMA 2020 source: Page et al., BMJ 372, n71 (2021); official checklist at https://www.prisma-statement.org/prisma-2020-checklist (CC BY 4.0).')
    items=[
    (1,'Title','Identify the report as a systematic review.','Title and Abstract','Partial: title states the biological synthesis; systematic identification is explicit in Abstract/Methods.'),
    (2,'Abstract','Provide a structured summary.','Abstract','Objectives, eligibility scale, methods, main estimates and limitations reported in compact journal format.'),
    (3,'Rationale','Describe rationale in context.','Introduction','Reported.'),
    (4,'Objectives','State explicit objectives.','Introduction, final paragraph','Reported.'),
    (5,'Eligibility criteria','Specify inclusion/exclusion criteria and grouping.','Methods: Protocol and study eligibility','Reported.'),
    (6,'Information sources','Specify databases and date last searched.','Methods: Systematic search and screening','GEO and PubMed; searched 9 September 2026.'),
    (7,'Search strategy','Present full strategies.','Methods; archive search_audit/search_methods_and_flow.md','Exact strings reported.'),
    (8,'Selection process','Describe screening process.','Methods; Supplementary Fig. 1; screening log','Single-reviewer computational project; every decision and reason retained.'),
    (9,'Data collection process','Describe data extraction.','Methods: Expression processing','Scripts and accession-level provenance supplied.'),
    (10,'Data items','List outcomes and other variables.','Methods; Supplementary Tables 1 and 14','Reported.'),
    (11,'Risk of bias','Describe assessment of study limitations.','Supplementary Table 14; Discussion','Reporting assessment used; no validated tool exists for this mechanobiology omics design.'),
    (12,'Effect measures','Specify effect measures.','Methods: Expression processing','Hedges-corrected standardized mean differences; paired standardized differences.'),
    (13,'Synthesis methods','Describe eligibility, preparation, displays, synthesis and heterogeneity.','Methods: Gene and pathway meta-analysis','REML, modified Hartung-Knapp, study clustering and sensitivity analyses reported.'),
    (14,'Reporting bias assessment','Describe methods.','Discussion','Formal funnel tests were not used with seven heterogeneous studies; stated as a limitation.'),
    (15,'Certainty assessment','Describe certainty methods.','Discussion; Interpretive safeguards','Triangulation, exact inference and sensitivity concordance used; no GRADE rating applied.'),
    (16,'Study selection','Report search and selection results.','Supplementary Fig. 1; screening log','46 screened, 12 full assessed, 7 included.'),
    (17,'Study characteristics','Cite and describe included studies.','Supplementary Table 1; Data availability','Reported by accession and design.'),
    (18,'Risk of bias in studies','Present assessments.','Supplementary Table 14','Reported as public-record reporting assessment.'),
    (19,'Results of individual studies','Present summary statistics/effect estimates.','Figure 2c; Supplementary Tables 6, 11 and machine tables','Reported.'),
    (20,'Results of syntheses','Summarize contributors, heterogeneity and sensitivity.','Results; Figures 1-3; Supplementary Tables','Reported with exact P values and sensitivity analyses.'),
    (21,'Reporting biases','Present assessment.','Discussion','Not estimable reliably at k=7; limitation stated.'),
    (22,'Certainty of evidence','Present certainty.','Discussion','Claims calibrated as supported, partial, borderline or unsupported.'),
    (23,'Discussion','Interpret results, limitations and implications.','Discussion','Reported.'),
    (24,'Registration and protocol','Provide registration/protocol information.','Methods; protocol archive','Not registered; revision-stage search stated; dated analysis protocols supplied.'),
    (25,'Support','Describe financial/nonfinancial support.','Funding','Funding sources and grant numbers supplied.'),
    (26,'Competing interests','Declare competing interests.','Competing interests','No competing interests declared.'),
    (27,'Availability','Report public availability of materials.','Data and Code availability; accompanying archive','Accessions and archived release supplied; version DOI 10.5281/zenodo.22726636.')
    ]
    add_table(d,['item','topic','PRISMA requirement','location','status/comment'],items,6.1)
    d.add_paragraph('Author names, affiliations, contributions, funding and competing-interest information are supplied in the manuscript and cover letter.')
    d.save(OUT/'PRISMA_2020_Checklist.docx')

if __name__=='__main__':
    flow_figure();manuscript();supplement();cover();prisma()
