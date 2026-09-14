from pathlib import Path
import json, re
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'outputs/研究升级'
META=BASE/'meta_analysis_v2'
V4=BASE/'high_impact_v4'
OUT=BASE/'Communications_Biology_v4'
OUT.mkdir(parents=True,exist_ok=True)
BLACK=RGBColor(0,0,0); RED=RGBColor(170,0,0); BLUE='244A66'


def shade(cell,fill):
    p=cell._tc.get_or_add_tcPr(); x=OxmlElement('w:shd'); x.set(qn('w:fill'),fill); p.append(x)
def borders(cell,color='D9D9D9'):
    p=cell._tc.get_or_add_tcPr(); b=p.first_child_found_in('w:tcBorders')
    if b is None: b=OxmlElement('w:tcBorders'); p.append(b)
    for edge in ('top','left','bottom','right','insideH','insideV'):
        e=OxmlElement('w:'+edge); e.set(qn('w:val'),'single'); e.set(qn('w:sz'),'4'); e.set(qn('w:color'),color); b.append(e)
def setup(d):
    sec=d.sections[0]; sec.top_margin=Inches(.72); sec.bottom_margin=Inches(.72); sec.left_margin=Inches(.82); sec.right_margin=Inches(.82)
    st=d.styles['Normal']; st.font.name='Arial'; st.font.size=Pt(9.5); st.font.color.rgb=BLACK
    st.paragraph_format.space_after=Pt(4); st.paragraph_format.line_spacing=1.08
    for name,size in [('Title',17),('Heading 1',13),('Heading 2',10.5),('Heading 3',9.7)]:
        s=d.styles[name]; s.font.name='Arial'; s.font.size=Pt(size); s.font.bold=True; s.font.color.rgb=BLACK
        s.paragraph_format.space_before=Pt(8); s.paragraph_format.space_after=Pt(4); s.paragraph_format.keep_with_next=True
        if name=='Title':
            ppr=s.element.get_or_add_pPr(); pb=ppr.find(qn('w:pBdr'))
            if pb is not None:ppr.remove(pb)
    footer=sec.footer.paragraphs[0]; footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
    fld=OxmlElement('w:fldSimple'); fld.set(qn('w:instr'),'PAGE'); footer._p.append(fld)
def addp(d,text,boldlead=None):
    p=d.add_paragraph()
    if boldlead and text.startswith(boldlead): p.add_run(boldlead).bold=True; p.add_run(text[len(boldlead):])
    else:p.add_run(text)
    return p
def placeholder(d,text):
    p=d.add_paragraph(); r=p.add_run(text); r.bold=True; r.font.color.rgb=RED; return p
def caption(d,text):
    p=d.add_paragraph(); p.paragraph_format.keep_with_next=True; r=p.add_run(text); r.bold=True; r.font.size=Pt(8.7)
def fig(d,name,legend,width=6.65):
    p=d.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.keep_with_next=True
    p.add_run().add_picture(str(OUT/name),width=Inches(width)); caption(d,legend)
def table(d,headers,rows,font=6.8,widths=None):
    t=d.add_table(rows=1,cols=len(headers)); t.alignment=WD_TABLE_ALIGNMENT.CENTER; t.autofit=True
    tr=t.rows[0]._tr.get_or_add_trPr(); rep=OxmlElement('w:tblHeader'); rep.set(qn('w:val'),'true'); tr.append(rep)
    for j,h in enumerate(headers):
        c=t.rows[0].cells[j]; c.text=str(h); shade(c,BLUE); borders(c); c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for r in c.paragraphs[0].runs:r.font.name='Arial';r.font.size=Pt(font);r.font.bold=True;r.font.color.rgb=RGBColor(255,255,255)
    for i,row in enumerate(rows):
        cells=t.add_row().cells
        for j,v in enumerate(row):
            c=cells[j]; c.text='' if pd.isna(v) else str(v); borders(c); c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if i%2:shade(c,'F2F6F8')
            for p in c.paragraphs:
                p.paragraph_format.space_after=Pt(0)
                for r in p.runs:r.font.name='Arial';r.font.size=Pt(font)
    return t


TITLE='Matrix stiffness responses localize to human osteoblast lineage programs without a universal gene signature'
AUTHORS=('Tingting Tan, MD¹,²,†; Xia Chen, MD, PhD¹,†; Zhengxiao Ouyang, MD, PhD¹,³,*\n'
'¹ Department of Orthopedics, The Second Xiangya Hospital, Central South University, Changsha, Hunan 410011, P.R. China\n'
'² Department of Immunology, School of Basic Medical Science, Central South University, Changsha, Hunan 410013, P.R. China\n'
'³ Osteopathy Laboratory of Surgical Department, The Second Xiangya Hospital, Central South University, Changsha, Hunan 410011, P.R. China\n'
'† These authors contributed equally.\n'
'Author email addresses: Tingting Tan: tinatan@126.com; Xia Chen (Associate Chief Physician): chenxiacxv99@csu.edu.cn; Zhengxiao Ouyang: ouyangzhengxiao@csu.edu.cn.')
CORRESPONDING='*Correspondence: Zhengxiao Ouyang, MD, PhD, Department of Orthopedics, The Second Xiangya Hospital, Central South University, 139 Renmin Middle Road, Changsha, Hunan 410011, P.R. China; Email: ouyangzhengxiao@csu.edu.cn; Telephone: +86-13548560675; ORCID: https://orcid.org/0000-0002-8997-0446'
FUNDING='This work was supported by the National Natural Science Foundation of China (grant nos. 82371600 and 82571826) and the Hunan Provincial Natural Science Foundation of China (grant nos. 2026JJ20072 and 2025JJ50614). The funders had no role in the design, conduct, analysis, interpretation, or reporting of this study.'
CONTRIBUTIONS='Tingting Tan: Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization, Writing – original draft. Xia Chen: Conceptualization, Methodology, Validation, Investigation, Resources, Writing – original draft, Writing – review and editing. Zhengxiao Ouyang: Conceptualization, Methodology, Supervision, Project administration, Funding acquisition, Writing – review and editing. Tingting Tan and Xia Chen contributed equally to this work. All authors reviewed and approved the final manuscript.'
AI_DISCLOSURE='OpenAI Codex assisted with code review, language editing, reproducibility auditing and document preparation. The authors reviewed and verified all analyses, numerical results, citations, figures and manuscript content and take full responsibility for the work. No generative AI imagery is included.'


def manuscript():
    d=Document();setup(d)
    p=d.add_paragraph(style='Title');p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run(TITLE)
    addp(d,AUTHORS)
    addp(d,CORRESPONDING)
    d.add_heading('Abstract',1)
    abstract=('Matrix stiffness regulates mesenchymal stromal cell behavior, yet the degree to which its transcriptomic effects generalize across experimental contexts and connect to osteoporosis is unclear. We reanalyzed seven independent public stiffness studies comprising 13 contrasts and 12,241 genes, while preserving each GEO series as the unit of inference. No gene passed false discovery rate correction, study-level gene-effect correlations ranged from -0.253 to 0.136, and 109 pooled-gene pathway signals fell to zero when inference was repeated across independent studies. Despite this directional heterogeneity, stiffness-response magnitude was enriched in osteoblast-lineage programs from 125,063 human femoral cells (matched permutation P=0.00040) and in an independent human bone spatial program (P=0.00260). The spatial osteoblast program showed positive enrichment in all seven stiffness studies (study-unit t-test P=0.0083; sign-test P=0.0078). A fixed 200-gene signature did not enrich estimated bone mineral density associations in 448,010 participants, identify an ITGA11-high osteogenic state, or show specific promoter-accessibility reversal after mechanical stretch or soft-matrix exposure. Public data therefore support distributed, lineage-localized mechanosensitivity rather than a universal stiffness gene program or osteoporosis biomarker.')
    addp(d,abstract)
    addp(d,'Keywords: matrix stiffness; mechanotransduction; mesenchymal stromal cells; osteoblast lineage; osteoporosis; transcriptomic meta-analysis; single-cell transcriptomics')
    d.add_heading('Introduction',1)
    addp(d,'Extracellular-matrix mechanics can alter adhesion, cytoskeletal tension, nuclear signaling, lineage commitment and inflammatory function in mesenchymal stromal cells (MSCs) [1-3]. These effects have encouraged attempts to connect engineered substrate stiffness with skeletal disease. Osteoporosis, however, is an organ-level loss of bone strength arising from cell composition, remodeling, architecture and material quality. A transcriptomic response to an in-vitro substrate cannot be assumed to measure bone stiffness or disease status.')
    addp(d,'The central reproducibility problem is experimental context. Published MSC studies vary in species, donor source, ligand chemistry, dimensionality, stiffness range, viscoelasticity, duration, passage and unit of biological replication [4-8]. Pooling dependent contrasts or treating cells and technical libraries as independent observations can produce precise gene lists while concealing disagreement among experiments. Conversely, an absence of conserved individual genes does not exclude a distributed response concentrated in a biologically relevant cell lineage.')
    addp(d,'We tested these alternatives with a study-clustered analysis of seven public MSC stiffness datasets. We first asked whether any gene or pathway transferred across independent studies. We then localized response magnitude in newly available human femoral single-cell and spatial gene programs, audited localization at the independent-study level, and challenged a fixed stiffness signature using donor-resolved osteogenic states, mechanical-intervention ATAC-seq and large-scale human bone-density genetics. The analysis distinguishes lineage localization from gene-level transfer and retains negative validation results as constraints on interpretation.')
    d.add_heading('Results',1)
    d.add_heading('Independent stiffness studies show little directional concordance',2)
    addp(d,'Seven eligible GEO series contributed 13 controlled stiffness contrasts. They covered human and rat bone-marrow MSCs, two- and three-dimensional systems, approximately 150 Pa to 900 kPa, exposure from 24 h to serial passage, and both paired and unpaired designs (Supplementary Table 1). Dependent times, coatings, passages and library preparations were collapsed within a GEO series before synthesis (Fig. 1a).')
    addp(d,'Among 12,241 genes measured in at least four independent studies, none passed Benjamini-Hochberg correction. The smallest nominal Hartung-Knapp P value was 0.0283 and its q value was 1.000. Pairwise Spearman correlations between complete study-level gene effects ranged from -0.253 to 0.136 (Fig. 1b). Thus, increasing the nominal sample count by combining contrasts would have hidden substantial cross-study disagreement.')
    d.add_heading('Pooled-gene pathway signals do not survive study-unit inference',2)
    addp(d,'Rank enrichment of pooled gene effects identified 109 false-discovery-rate-significant Hallmark or Reactome pathways, with strong signals for RNA processing, MYC targets and cell-cycle functions. Genes within a pathway are correlated, so we repeated each test by deriving one enrichment score per independent GEO series. None of 1,008 pathways passed false-discovery-rate correction; the smallest study-unit q value was 0.749 (Fig. 1c). No pathway met the joint requirement of pooled-gene significance, study-unit significance and leave-one-study direction stability.')
    fig(d,'Figure_1.png','Figure 1  Study-clustered analysis reveals cross-context heterogeneity. a, Analysis architecture. b, Pairwise Spearman correlations of study-level gene effects. c, False-discovery-rate evidence for the six strongest pooled-gene pathways compared with inference using seven independent studies. The dashed line marks q=0.05.')
    d.add_heading('Stiffness-response magnitude localizes to human osteoblast-lineage programs',2)
    addp(d,'We next tested whether heterogeneous gene directions could nevertheless localize to a shared cell lineage. Author-defined gene programs from 125,063 human femoral cells were obtained independently of the stiffness meta-analysis [9]. The prespecified union of two MSC states, pre-osteoblasts and mature osteoblasts contained 1,007 genes with meta-analysis estimates. Their mean absolute meta z score exceeded variance- and study-coverage-matched background sets (standardized matched effect 3.28, permutation P=0.00040). At the individual-cluster level, both MSC states remained significant after correction across 20 cell types (q=0.0010 and 0.0015), and pre-osteoblasts met q=0.0453 (Fig. 2a). Mature osteoblasts did not show significant response-magnitude enrichment (q=1.000), placing the signal earlier in the lineage.')
    addp(d,'A prespecified replication used human bone spatial transcriptomic programs from a separate cohort [9]. The union of MSC, Osteo-MSC and osteoblast programs was enriched (184 mapped genes; standardized matched effect 2.91; permutation P=0.00260; Fig. 2b). The three individual spatial states did not survive correction when tested separately. The single-cell and spatial unions shared 125 genes. After those genes were removed, the remaining 64 mapped spatial genes retained a positive effect but were not significant (P=0.107), indicating that the replication partly depends on a shared osteoblast-lineage core.')
    d.add_heading('Localization persists when each stiffness study is the inferential unit',2)
    addp(d,'To prevent correlated genes from serving as replicates, we computed one matched enrichment effect per GEO series. The single-cell osteoblast program was enriched in five of seven studies; its mean raw matched effect was 0.132 and the one-sided study-unit t-test gave P=0.0408, whereas the exact sign test was not significant (P=0.227). The spatial osteoblast program was enriched in all seven studies, with a mean matched effect of 0.256 (t-test P=0.0083; sign-test P=0.0078; Fig. 2c). This result supports reproducible localization of response magnitude even though the directions of individual genes are not conserved.')
    fig(d,'Figure_2.png','Figure 2  Stiffness-response magnitude is concentrated in human osteoblast-lineage programs. a, Matched differences in mean absolute meta z score for 20 human femoral single-cell programs. b, Corresponding results for 11 human bone spatial programs. Orange marks prespecified osteoblast-lineage states. c, One matched enrichment effect per independent stiffness study for the single-cell and spatial osteoblast-lineage unions. Values above zero indicate greater response magnitude than variance-matched background genes.',width=5.45)
    d.add_heading('A fixed gene signature does not transfer across orthogonal biological layers',2)
    addp(d,'We fixed a 200-gene signature from the 100 largest positive and 100 largest negative pooled effects before orthogonal validation. In a donor-resolved study of 136,014 cultured human stromal cells, the signature did not distinguish ITGA11-high from ITGA11-low osteogenic cultures (directional score -0.093, permutation P=0.117) and the genome-wide correlation with the stiffness meta-effect was -0.0158. Competitive analysis of the study-defined ITGA11-high osteogenic cluster showed only nominal signed enrichment (P=0.0541, q=0.0901) and no enrichment in absolute response magnitude.')
    addp(d,'Mechanical-intervention ATAC-seq provided a direct perturbational test [10]. Replicate promoter signals were highly concordant within senescent control, mechanically stretched, young control and soft-matrix groups (rho=0.980-0.990). The prespecified rank-based signature contrast was not restored by stretch (directional score -0.00215, matched P=0.792). Soft-matrix exposure in young cells changed the score in the expected negative direction but did not reach significance (score -0.00871, P=0.235). After confirming that the deposited BigWigs were TPM-normalized, a labelled post hoc analysis reproduced the authors’ global promoter shift after stretch (+0.103 mean log1p TPM) and soft matrix (-0.126), but still found no signature-specific separation (P=0.633 and 0.238).')
    addp(d,'The genetic test used MAGMA gene statistics from an eBMD GWAS of 448,010 UK Biobank participants [9]. The fixed signature contained 191 genes with both estimates and was not enriched for stronger eBMD association (standardized matched effect -0.379, P=0.643). A regression adjusting for gene length, SNP number, model parameters and chromosome agreed (beta=-0.141, HC3 one-sided P=0.716). Thirty-four signature genes passed the published eBMD gene threshold, compared with 30.4 expected; the overlap was not significant (hypergeometric P=0.266).')
    fig(d,'Figure_3.png','Figure 3  Orthogonal datasets constrain gene-level transfer. a, Standardized matched effects for prespecified localization and transfer tests. For the soft-matrix analysis, the sign is oriented so that a positive value supports the expected reverse direction. b, Author-normalized ATAC signals reproduce global promoter shifts without separating the fixed positive and negative stiffness gene sets. c, Observed and expected numbers of eBMD-significant genes in the fixed signature.')
    d.add_heading('Discussion',1)
    addp(d,'This analysis supports a distributed model of MSC mechanosensitivity. Individual gene directions and pathway claims did not transfer across experimental systems, yet genes defining early human osteoblast-lineage states repeatedly carried larger stiffness responses. The distinction matters: cell-lineage localization can be reproducible even when no universal gene panel exists. Treating correlated genes as independent observations obscured that structure by producing 109 pathway discoveries that disappeared at the study level.')
    addp(d,'The localization was strongest in two MSC states and pre-osteoblasts, whereas mature osteoblasts were not enriched in response magnitude. This pattern is compatible with greater mechanical plasticity before terminal differentiation. The spatial program also showed positive enrichment in every stiffness study, but overlap removal weakened its gene-level significance. The strongest interpretation is therefore a shared lineage core supported in separate tissue representations, rather than fully independent replication by disjoint genes.')
    addp(d,'The orthogonal tests define the limits of translation. Moderate stretch and substrate softening produced global, directionally appropriate changes in promoter accessibility, but not selective reversal of the fixed stiffness signature. The signature also failed to identify an ITGA11-high osteogenic state and was not enriched for eBMD association. These results argue against presenting a small in-vitro signature as an osteoporosis biomarker. Human genetic association and tissue localization answer different questions: osteoblast-lineage genes can be mechanically responsive as a distributed class without the most extreme transcript changes being the genes most strongly associated with eBMD.')
    addp(d,'This study is limited by small primary studies, heterogeneous processed data and incomplete reporting of matrix chemistry and biological replication. Only seven series were available for study-unit inference, making moderator analysis underpowered. The ATAC perturbation contained two biological samples per group and supported gene-set rather than gene-level inference. The single-cell and spatial gene programs were published in the same multiscale study, although they derived from distinct human datasets. The localization protocols were frozen before gene-set testing, while the TPM and overlap-removal analyses were explicitly post hoc. None of the integrated analyses establish causality.')
    addp(d,'Public mechanobiology data can therefore support a high-resolution boundary-condition study without new wet-laboratory experiments. The reproducible feature is the concentration of response magnitude in early osteoblast-lineage programs. Future studies should retain donor identifiers, deposit raw counts and peak matrices, report modulus, ligand chemistry, dimensionality and viscoelasticity, and pair controlled mechanics with donor-resolved bone phenotypes. These design changes are needed before a transferable mechanical biomarker can be evaluated.')
    d.add_heading('Methods',1)
    d.add_heading('Protocol and study eligibility',2)
    addp(d,'A protocol was frozen before the cross-study gene meta-analysis. Primary studies imposed at least two extracellular-matrix stiffness levels on MSCs or mesenchymal progenitors, included at least two biological samples per condition and provided processed transcriptomic data. One GEO series was one independent study. Times, coatings, passages and library preparations within a series were treated as dependent contrasts. GSE166824, GSE22011, GSE152708, GSE35958, GSE230665 and GSE276529 were reserved for external analyses. Later human-bone, ATAC and genetics protocols were dated and frozen before their corresponding gene-set tests; amendments are retained with the code release.')
    d.add_heading('Expression processing and within-study effects',2)
    addp(d,'RNA-sequencing counts were converted to log2 counts per million after low-expression filtering. FPKM or RPKM values were log2-transformed; deposited log-scale arrays were not transformed again. Multiple features for one official symbol were collapsed by the median. Rat genes were restricted to one-to-one human orthologues. Unpaired two-level contrasts used Hedges-corrected standardized mean differences for stiff minus soft conditions and their sampling variances. Paired donors used within-donor differences. Multiple contrasts within a series were averaged to one study-level gene effect; variance was the mean within-contrast variance plus between-contrast variance.')
    d.add_heading('Gene and pathway meta-analysis',2)
    addp(d,'Genes measured in at least four studies entered random-effects meta-analysis. Between-study variance used restricted maximum likelihood. Confidence intervals and P values used modified Hartung-Knapp uncertainty with scale not below one and a t distribution with k-1 degrees of freedom. Benjamini-Hochberg correction covered 12,241 genes. Hallmark and Reactome sets from MSigDB 2025.1 were restricted to 15-500 genes. Pooled-gene rank tests were followed by one enrichment score per GEO series and a two-sided one-sample t-test across the seven scores; false-discovery rates were controlled across 1,008 pathways.')
    d.add_heading('Human bone cell-program localization',2)
    addp(d,'Author-defined human femoral single-cell and bone spatial gene programs were read from Supplementary Table 2 of Chai et al. [9]. Before merging, the single-cell osteoblast-lineage union was fixed as both MSC clusters, pre-osteoblasts and mature osteoblasts; the spatial union was MSCs, Osteo-MSCs and osteoblasts. The outcome was absolute gene meta z score. For 20,000 permutations, program genes were matched to non-program genes by number of contributing stiffness studies and Hartung-Knapp standard-error deciles. The observed-minus-null mean and a one-sided empirical P value were reported. Cluster tests controlled FDR separately for signed and absolute endpoints.')
    addp(d,'For study-unit audit, each GEO series contributed one difference between the mean absolute study effect in a fixed program and its variance-decile-matched null mean. The seven raw matched effects and their permutation-standardized counterparts were tested against zero with one-sided one-sample t-tests. An exact one-sided sign test was also reported. A post hoc sensitivity removed all genes shared by the single-cell and spatial unions before retesting the spatial program.')
    d.add_heading('Mechanical-intervention ATAC-seq',2)
    addp(d,'GSE287556 BigWig tracks represented senescent BMSCs, senescent BMSCs exposed to 5% strain at 0.02 Hz for 2 h, young BMSCs, and young BMSCs cultured on 1.5 kPa substrate for 24 h, with two samples per group [10]. GENCODE v49 GRCh38 protein-coding gene features defined one promoter per gene as TSS +/-2 kb. Mean BigWig signal was log1p-transformed and ranked within sample; group effects were mean treated rank minus mean control rank. The fixed positive-minus-negative signature contrast used 20,000 permutations matched by baseline accessibility and gene length. The bidirectional test expected the young soft-matrix effect to oppose stretch. BigWig files and downloads were checked for exact byte length and gzip integrity. After the primary result, the published TPM normalization was incorporated in a labelled post hoc log1p TPM sensitivity analysis.')
    d.add_heading('Osteogenic-state and eBMD analyses',2)
    addp(d,'GSE317531 provided author-defined marker genes for five donor-resolved human stromal-cell states and bulk contrasts of ITGA11-high versus ITGA11-low cultures. Competitive cluster tests used the complete meta z-score background. A fixed 200-gene signature comprised the 100 largest positive and 100 largest negative pooled stiffness effects; its directional score in each bulk contrast was tested using 20,000 random gene sets of identical sizes.')
    addp(d,'MAGMA eBMD gene statistics for 448,010 UK Biobank participants were read from Supplementary Table 7e of Chai et al. [9]. The primary test compared mean ZSTAT for the fixed signature with 20,000 non-signature sets matched by gene-length and NSNPS deciles. A prespecified sensitivity regressed ZSTAT on signature membership, log gene length, log NSNPS, log NPARAM and chromosome, using HC3 standard errors. Gene significance used the source threshold P<2.5x10^-6. Because ZSTAT measures association strength without an effect direction comparable with expression, no high-versus-low BMD direction was assigned.')
    d.add_heading('Software and reproducibility',2)
    addp(d,'Analyses used Python 3 with pandas, NumPy, SciPy and pyBigWig, and R with ggplot2. All random procedures used seed 20260908. Frozen protocols, amendments, input provenance, gene-level and study-level results, figure source tables and executable scripts accompany the submission. Large public source files are referenced by accession rather than redistributed.')
    d.add_heading('Data availability',1)
    addp(d,'Primary transcriptomic data are available from NCBI GEO under GSE193021, GSE181512, GSE226411, GSE288678, GSE55867, GSE255574 and GSE310514. External datasets are GSE166824, GSE22011, GSE152708, GSE35958, GSE230665, GSE276529, GSE317531 and GSE287556. Human femoral single-cell data are available under GSE317069, spatial bone data under GSE299207 and the eBMD GWAS under GWAS Catalog accession GCST90726625 [9]. Derived data needed to reproduce figures and tests are supplied with this submission.')
    d.add_heading('Code availability',1)
    addp(d,'All custom scripts required for the reported analyses and figures are included in the accompanying code and derived-data archive. Public source downloads are identified by stable accession or official URL.')
    d.add_heading('References',1)
    refs=[
    'Engler, A. J., Sen, S., Sweeney, H. L. & Discher, D. E. Matrix elasticity directs stem cell lineage specification. Cell 126, 677-689 (2006).',
    'Vining, K. H. & Mooney, D. J. Mechanical forces direct stem cell behaviour in development and regeneration. Nat. Rev. Mol. Cell Biol. 18, 728-742 (2017).',
    'Kechagia, J. Z., Ivaska, J. & Roca-Cusachs, P. Integrins as biomechanical sensors of the microenvironment. Nat. Rev. Mol. Cell Biol. 20, 457-473 (2019).',
    'Brielle, S. et al. Delineating the heterogeneity of matrix-directed differentiation toward soft and stiff tissue lineages via single-cell profiling. Proc. Natl Acad. Sci. USA 118, e2016322118 (2021).',
    'Lim, J. J., Vining, K. H., Mooney, D. J. & Blencowe, B. J. Matrix stiffness-dependent regulation of immunomodulatory genes in human MSCs is associated with the lncRNA CYTOR. Proc. Natl Acad. Sci. USA 121, e2404146121 (2024).',
    'Kersey, A. L. et al. Stiffness assisted cell-matrix remodeling trigger 3D mechanotransduction regulatory programs. Biomaterials 306, 122473 (2024).',
    'Schellenberg, A. et al. Matrix elasticity, replicative senescence and DNA methylation patterns of mesenchymal stem cells. Biomaterials 35, 6351-6358 (2014).',
    'Younesi, F. S. et al. A circuit of mechanically regulated transcription factors balances regenerative and fibrotic memory of mesenchymal stromal cells. Adv. Sci. 13, e22056 (2026).',
    'Chai, R. C. et al. Multiscale analysis and functional validation of the cellular and genetic determinants of skeletal disease. Nat. Genet. 58, 1530-1545 (2026).',
    'Liu, X., Ye, Y., Li, Z., Liao, L. & Wei, Q. Mechanical rejuvenation of senescent stem cells and aged bone via chromatin remodeling. Nat. Commun. 17, 1684 (2026).',
    'Al-Barghouthi, B. M. et al. Systems genetics in diversity outbred mice inform BMD GWAS and identify determinants of bone strength. Nat. Commun. 12, 3408 (2021).',
    'Morris, J. A. et al. An atlas of genetic influences on osteoporosis in humans and mice. Nat. Genet. 51, 258-266 (2019).',
    'Benisch, P. et al. The transcriptional profile of mesenchymal stem cell populations in primary osteoporosis is distinct and shows overexpression of osteogenic inhibitors. PLoS ONE 7, e45142 (2012).',
    'Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl Acad. Sci. USA 102, 15545-15550 (2005).',
    'Frankish, A. et al. GENCODE 2021. Nucleic Acids Res. 49, D916-D923 (2021).',
    'National Center for Biotechnology Information. Gene Expression Omnibus. https://www.ncbi.nlm.nih.gov/geo/ (accessed 8 September 2026).']
    for i,r in enumerate(refs,1):
        p=d.add_paragraph(f'{i}. {r}');p.paragraph_format.left_indent=Inches(.2);p.paragraph_format.first_line_indent=Inches(-.2);p.paragraph_format.space_after=Pt(2)
        for x in p.runs:x.font.size=Pt(8.2)
    d.add_heading('Acknowledgements',1);addp(d,'We thank the investigators and participants of the public datasets analyzed in this study for making these resources available.')
    d.add_heading('Funding',1);addp(d,FUNDING)
    d.add_heading('Author contributions',1);addp(d,CONTRIBUTIONS)
    d.add_heading('Competing interests',1);addp(d,'The authors declare no competing interests.')
    d.add_heading('Ethics statement',1);addp(d,'This study reanalyzed de-identified data from public repositories and performed no new human or animal experiments. Ethical approvals for source studies are reported in the original publications.')
    d.add_heading('Generative AI disclosure',1);addp(d,AI_DISCLOSURE)
    d.save(OUT/'Manuscript_Communications_Biology.docx')
    print('abstract_words',len(re.findall(r"\b[\w-]+\b",abstract)))


def supplement():
    d=Document();setup(d)
    p=d.add_paragraph(style='Title');p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.add_run('Supplementary Information')
    addp(d,TITLE);addp(d,'Tingting Tan†, Xia Chen† and Zhengxiao Ouyang*\n†These authors contributed equally. *Correspondence: ouyangzhengxiao@csu.edu.cn')
    addp(d,'This document contains the design audit and principal numerical results. Complete machine-readable tables, dated protocols and executable code are provided in the accompanying archive.')
    c=pd.read_csv(META/'cohort_characteristics_v3.csv')
    caption(d,'Supplementary Table 1  Primary and external dataset roles')
    table(d,c.columns,c.values,font=6.2)
    g=pd.read_csv(META/'gene_meta_analysis.csv').sort_values('p_hk').head(30)
    g=g[['gene','k','pooled_effect','ci_low','ci_high','p_hk','BH_q','I2_percent','max_weight']].copy()
    for col in g.columns[2:]:g[col]=g[col].map(lambda x:f'{x:.4g}')
    caption(d,'Supplementary Table 2  Top 30 gene-level meta-analysis results')
    table(d,g.columns,g.values,font=6.3)
    q=pd.read_csv(META/'pathway_study_unit_meta.csv').sort_values('p_study_unit').head(30)
    q=q[['collection','pathway','k_studies','mean_study_z','p_study_unit','BH_q_study_unit','BH_q_pooled_gene_rank','robust_pathway']]
    for col in ['mean_study_z','p_study_unit','BH_q_study_unit','BH_q_pooled_gene_rank']:q[col]=q[col].map(lambda x:f'{x:.4g}')
    caption(d,'Supplementary Table 3  Top 30 pathways using independent studies as the inference unit')
    table(d,q.columns,q.values,font=5.8)
    sc=pd.read_csv(V4/'human_bone_localization/human_scRNA_cluster_tests.csv')
    sc=sc[['cluster','endpoint','n_genes','effect','standardized_effect','permutation_p','BH_q']].copy()
    for col in sc.columns[3:]:sc[col]=sc[col].map(lambda x:f'{x:.4g}')
    caption(d,'Supplementary Table 4  Human femoral single-cell program localization')
    table(d,sc.columns,sc.values,font=6.1)
    sp=pd.read_csv(V4/'human_bone_localization/human_spatial_cluster_tests.csv')
    sp=sp[['cluster','endpoint','n_genes','effect','standardized_effect','permutation_p','BH_q']].copy()
    for col in sp.columns[3:]:sp[col]=sp[col].map(lambda x:f'{x:.4g}')
    caption(d,'Supplementary Table 5  Human bone spatial program localization')
    table(d,sp.columns,sp.values,font=6.2)
    su=pd.read_csv(V4/'human_bone_localization/study_unit_localization_effects.csv')
    for col in su.columns[3:]:su[col]=su[col].map(lambda x:f'{x:.4g}')
    caption(d,'Supplementary Table 6  Study-unit osteoblast-lineage localization audit')
    table(d,su.columns,su.values,font=6.2)
    at1=json.loads((V4/'atac_validation/atac_validation_summary.json').read_text())
    at2=json.loads((V4/'atac_validation/bidirectional_validation_summary.json').read_text())
    ats=[['Force+ vs senescent','rank directional',at1['primary']['directional_score'],at1['primary']['one_sided_permutation_p'],at1['within_group_qc']['force_plus_rho']],['Force- vs young','rank directional',at2['primary_like_directional_score'],at2['one_sided_matched_permutation_p'],at2['within_group_qc']['force_minus_1__force_minus_2']]]
    caption(d,'Supplementary Table 7  Prespecified ATAC signature tests')
    table(d,['comparison','endpoint','score','P','within-group rho'],[[a,b,f'{c:.5g}',f'{e:.5g}',f'{f:.4f}'] for a,b,c,e,f in ats],font=7)
    ge=json.loads((V4/'genetics_validation/genetics_validation_summary.json').read_text())
    genrows=[['Matched permutation',ge['primary']['n_genes'],f"{ge['primary']['standardized_effect']:.4g}",f"{ge['primary']['one_sided_permutation_p']:.4g}"],['Adjusted regression',ge['adjusted_regression_sensitivity']['n_genes'],f"{ge['adjusted_regression_sensitivity']['signature_beta']:.4g}",f"{ge['adjusted_regression_sensitivity']['one_sided_p']:.4g}"],['Significant-gene overlap',ge['significant_eBMD_overlap']['signature'],ge['significant_eBMD_overlap']['overlap'],f"{ge['significant_eBMD_overlap']['hypergeometric_p']:.4g}"]]
    caption(d,'Supplementary Table 8  eBMD genetics validation')
    table(d,['analysis','n','effect or overlap','P'],genrows,font=7)
    it=pd.read_csv(V4/'stiffness_signature_in_itga11_bulk.csv')
    caption(d,'Supplementary Table 9  ITGA11 bulk-expression validation')
    table(d,it.columns,it.round(5).values,font=6.5)
    d.add_heading('Reporting and protocol record',1)
    addp(d,'The study follows PRISMA principles where applicable to public-data identification and quantitative synthesis. Eligibility criteria, series-level inclusion, dependent-contrast handling, exclusions, full accession list, synthesis model, heterogeneity, multiple-testing correction and all prespecified sensitivity analyses are reported. Because GEO functional-genomics searches do not provide a bibliographic record count comparable with a conventional intervention review, the accession-level audit is the reproducible search record.')
    caption(d,'Supplementary Table 10  Dated protocol and amendment registry')
    records=[
        ['Cross-study meta-analysis v2','Before gene-level synthesis','Primary protocol','Eligibility, study clustering, random-effects inference and external validation'],
        ['ATAC validation v1','Before promoter extraction','Primary protocol','TSS window, rank endpoint, matched permutation and bidirectional test'],
        ['eBMD genetics validation v1','Before GWAS merge','Primary protocol','Fixed signature, gene-size and SNP-count matching, regression sensitivity'],
        ['Human bone localization v1','Before cell-program merge','Primary protocol','Single-cell osteoblast union, spatial replication and cluster multiplicity'],
        ['Study-unit localization audit','Before study-level enrichment','Primary robustness protocol','Seven GEO series as inference units with t and sign tests'],
        ['ATAC TPM normalization record','After primary ATAC test','Post hoc amendment','Added log1p TPM sensitivity after source-method normalization was confirmed'],
        ['Spatial overlap-removal record','After spatial replication','Post hoc sensitivity','Removed all genes shared with the single-cell osteoblast program']]
    table(d,['record','timing','status','scope'],records,font=6.6)
    d.add_heading('Protocol deviations and post hoc analyses',2)
    addp(d,'The prespecified ATAC rank result remained primary. A log1p TPM sensitivity was added only after source normalization was confirmed. A second post hoc analysis removed 125 genes shared by the single-cell and spatial osteoblast programs; the remaining 64 mapped genes retained a positive effect but were not significant (P=0.107). No primary endpoint or threshold was changed. Complete dated protocols are included in the archive.')
    d.save(OUT/'Supplementary_Information_Communications_Biology.docx')


def cover():
    d=Document();setup(d);p=d.add_paragraph(style='Title');p.add_run('Cover letter for Communications Biology')
    addp(d,'15 September 2026');addp(d,'Editors\nCommunications Biology');addp(d,'Dear Editors,')
    addp(d,f'Please consider our Article, “{TITLE},” for publication in Communications Biology.')
    addp(d,'Matrix stiffness is widely treated as a source of transferable mesenchymal stromal-cell programs. We find a more specific biological organization. Across seven public studies, no gene or pathway survived independent-study inference, yet the magnitude of stiffness responses was concentrated in early osteoblast-lineage programs from 125,063 human femoral cells and reproduced in human bone spatial programs. Crucially, the spatial program was enriched in all seven independent stiffness studies. The result identifies lineage-localized, distributed mechanosensitivity while explaining why small universal gene panels fail across experimental contexts.')
    addp(d,'We challenged that interpretation rather than relying on one positive enrichment. A frozen 200-gene signature failed in donor-resolved ITGA11 osteogenic cultures, bidirectional mechanical-intervention ATAC-seq and eBMD gene statistics from 448,010 participants. The ATAC data reproduced global force-dependent accessibility changes but not signature-specific reversal. These orthogonal negative results delimit the claim and distinguish cell-state localization from an osteoporosis biomarker.')
    addp(d,'The manuscript is suited to Communications Biology because it uses transparent secondary analysis to provide a biological advance of broad relevance to mechanobiology and skeletal research: response magnitude can be reproducibly lineage-localized even when individual transcript directions are context dependent. Dated protocols, study-unit analyses, code, derived data and figure source tables are included.')
    addp(d,'No new human or animal experiments were performed. All participant-level source data were de-identified and publicly released by the original investigators.')
    addp(d,'The author team comprises Tingting Tan, Xia Chen and Zhengxiao Ouyang; Tingting Tan and Xia Chen contributed equally. All authors have reviewed and approved the manuscript and agree with its submission to Communications Biology. The manuscript is original, has not been published previously and is not under consideration elsewhere. '+FUNDING+' The authors declare no competing interests.')
    addp(d,'Thank you for your consideration.')
    addp(d,'Sincerely,\nZhengxiao Ouyang, MD, PhD\nCorresponding author\nDepartment of Orthopedics, The Second Xiangya Hospital, Central South University\n139 Renmin Middle Road, Changsha, Hunan 410011, P.R. China\nEmail: ouyangzhengxiao@csu.edu.cn\nORCID: https://orcid.org/0000-0002-8997-0446')
    d.save(OUT/'Cover_Letter_Communications_Biology.docx')


if __name__=='__main__':manuscript();supplement();cover()
