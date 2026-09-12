# Search methods and PRISMA-style flow

Search date: 9 September 2026. No date or language restriction was applied.

## Sources and exact queries

**NCBI GEO DataSets (Entrez E-utilities, series records only)**

`gse[ETYP] AND (stiffness[All Fields] OR "matrix elasticity"[All Fields] OR "substrate stiffness"[All Fields] OR "mechanical memory"[All Fields]) AND (mesenchymal[All Fields] OR BMSC[All Fields])`

This returned 42 GEO series records. The raw ESearch response and all 42 ESummary records are archived.

**PubMed (Entrez E-utilities, focused omics search)**

`(("substrate stiffness"[Title/Abstract] OR "matrix stiffness"[Title/Abstract] OR "matrix elasticity"[Title/Abstract]) AND ("mesenchymal stromal cell"[Title/Abstract] OR "mesenchymal stem cell"[Title/Abstract] OR BMSC[Title/Abstract]) AND (transcriptom*[Title/Abstract] OR RNA-seq[Title/Abstract] OR RNAseq[Title/Abstract] OR microarray[Title/Abstract]))`

This returned three publication records. Citation and known-accession checking contributed one additional eligible GEO series (GSE193021), which the terminology-dependent GEO query missed.

## Flow counts

- Database records: GEO 42 + PubMed 3 = 45
- Records from citation/accession checking: 1
- Duplicates removed: 0
- Records screened: 46
- Excluded from title/summary or structured GEO summary: 34
- Full GEO record/reports assessed: 12
- Excluded after full eligibility assessment: 5
  - one donor with non-independent cells: 1 (retained only as descriptive external data)
  - stiffness conditions separated across independently normalized GEO series: 2
  - glass-only series without soft-versus-stiff contrast: 1
  - strain/vibration rather than matrix-stiffness perturbation: 1
- Studies included in quantitative synthesis: 7

The accession-level decision and reason for every record are in `study_screening_log_2026-09-09.csv`. The search was performed during manuscript revision and was not prospectively registered; this limitation must remain explicit.
