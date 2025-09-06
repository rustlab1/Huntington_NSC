# RNA-seq Analysis: HTT-Mutant vs. Control (GSE270472)

Reproducible workflow for **HTT Mutant (Mutant)** vs. **Control** samples from [GSE270472](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270472).  
Paper: *Huntingtin loss-of-function contributes to transcriptional deregulation in Huntington’s disease*, Cell & Bioscience (2025), doi: 10.1186/s13578-025-01443-5.

## Workflow
### 1. Data processing (bash)
1. **Download FASTQ** (ENA SRRs for HTT-Mutant & Control):contentReference[oaicite:0]{index=0}  
2. **Trim & Subsample** with `fastp` + `reformat.sh` (~300 MB/sample):contentReference[oaicite:1]{index=1}  
3. **Align & Count** with **HISAT2** → BAM → **featureCounts** → `counts/final_counts_symbols.tsv`:contentReference[oaicite:2]{index=2}

### 2. Downstream analysis (R)
1. **Load counts & annotation** → filter low-expressed genes  
2. **Sample grouping**: HTT-Mutant vs. Control  
3. **Differential expression** with **DESeq2** (Wald test, FDR correction)  
4. **Top table export**: padj, log2FC, gene symbols, description  
5. **Volcano plot**: visualize significant DEGs  
6. **GO enrichment (BP)** with **clusterProfiler** and `org.Hs.eg.db`  
   - Export enrichment results  
   - Barplot of top enriched categories

## Output
- `counts/final_counts_symbols.tsv` (counts matrix)  
- `top250.tsv` (top DEGs with annotation)  
- `normalized_counts.tsv` (optional)  
- `enrichGO_BP_results.tsv` (GO enrichment results)  
- Plots: volcano plot, GO barplot

## Data
- GEO: **GSE270472**  
- Organism: *Homo sapiens*  
- Groups: HTT-KO, Control
