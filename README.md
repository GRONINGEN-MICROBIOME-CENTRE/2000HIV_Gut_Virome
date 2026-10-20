# 2000HIV Gut Virome
Gut virome characterization, strain sharing, and clinical associations in the **2000HIV study**.

<p align="center">
  <img src="images/Overview LLNEXT.png" alt="Overview of pregnancy and early-life gut virome analysis workflow" width="900">
</p>

<table align="center">
<tr>

<td align="center" valign="top" style="padding-right:40px;">
<div style="margin-bottom:6px;"><b>Publication</b></div>
<a href="https://doi.org/10.64898/2026.03.05.709809"> 
<img src="images/bioRxiv_logo.png" height="39">
</a>
</td>

<td align="center" valign="top" style="padding-left:40px;">
<div style="margin-bottom:6px;"><b>Viral catalogs</b></div>
<a href="https://doi.org/10.5281/zenodo.18981021">
<img src="images/zenodo_logo.png" height="45">
</a>
</td>

</tr>
</table>

## Overview
Here, we present a comprehensive characterization of the **gut DNA virome in people with HIV (PWH)**, investigating its composition, strain sharing, and associations with behavioural, immunological, and clinical features in the 2000HIV study.

- 🧬 **Large-scale virome profiling:** We analyzed faecal shotgun metagenomes from 1,661 participants in the 2000HIV cohort, generating a comprehensive viral catalog that substantially expands known gut viral diversity in PWH.

- 🔗 **Sexual behaviour and viral sharing:** Men who have sex with men (MSM) showed increased viral diversity and a higher *Prevotella*/*Bacteroides* phage ratio. Viral strain sharing increased with the number of sexual partners and receptive anal sex, but did not mirror HIV phylogenetic relationships.

- 🛡️ **Immunological status and HIV history:** Specific viral populations were associated with immunological status, duration of HIV infection, and markers of viral control and persistence, including HIV controller status and HIV reservoir size.

- 💉 **COVID-19 vaccination and immune responses:** Viral community profiles differed by COVID-19 vaccination status, with these differences partly mediated by cytokine responses measured in ex vivo stimulation assays.
  
---

## Table of Contents
- [Cohort](#cohort)
- [Methods](#methods)
  - [Assembly](#assembly)
  - [Virus Identification and Abundance Profiling](#virus-identification-and-abundance-profiling)
  - [Virus Characterization](#virus-characterization)
  - [Strain-Level Profiling and Characterization of Shared Viruses](#strain-level-profiling-and-characterization-of-shared-viruses)
- [Data Availability](#data-availability)
- [Code Structure](#code-structure)
- [Citation](#citation)


---

## Cohort

This study is based on the 2000HIV cohort, a longitudinal study of PWH recruited between 2019 and 2021 from four HIV treatment centres and followed using extensive multi-omics profiling, immunological assays, and clinical assessments.

- **Original cohort description**:
Vos W, Groenendijk A, et al. The 2000HIV study: Design, multi-omics methods and participant characteristics. Frontiers in Immunology, 2022.

---

## Methods

### Virus Identification and Abundance Profiling

Viral sequences were identified from assembled contigs generated using **metaSPAdes**.

- **Viral identification:** [geNomad](https://github.com/apcamargo/genomad) was used to identify viral sequences in the assembled contigs.
- **Prophage pruning and quality control:** Predicted prophage regions were initially pruned using geNomad. [CheckV](https://bitbucket.org/berkeleylab/checkv) was then used to remove host contamination, perform final prophage pruning, estimate viral genome completeness, and filter viral sequences.
- **Dereplication:** CheckV scripts (`anicalc.py` and `aniclust.py`) were used for dereplication against public viral genome databases: 95% ANI / 85% AF against the Gut Human Viral Catalogue (UHGV) and RefSeq.
- **Abundance estimation:** Quality-controlled reads were mapped to the dereplicated viral genome catalogue using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/).

---

### Virus Characterization and Strain-Level Profiling

- **Taxonomic classification:** Viral taxonomy was assigned using a hierarchical approach combining RefSeq annotations, **geNomad**, and **VITAP**.
- **Lifestyle prediction:** Bacteriophage lifestyles were predicted using **BACPHLIP**.
- **Host prediction:** **iPHoP** was used to predict viral hosts, with host assignments summarized at the genus level.
- **Strain-level viral profiling:** **inStrain** was used for strain-level comparisons, with strain sharing defined as population ANI (popANI) > 99.999%. Strain-sharing patterns were assessed in relation to sexual behaviour, including the number of sexual partners and receptive anal sex, and HIV phylogenetic clusters.

---

### *Caudoviricetes* and HIV Phylogenetic Analysis

- ***Caudoviricetes* phylogeny:** A maximum-likelihood phylogeny was reconstructed from concatenated alignments of conserved viral marker proteins, including the large terminase subunit, major capsid protein, and portal protein. Marker detection and alignment were performed using **HMMER**, **MUSCLE**, and **ClipKIT**, followed by tree inference with **FastTree**.
- **HIV phylogeny:** HIV *pol* sequences were amplified from CD4⁺ T-cell genomic DNA and sequenced using Sanger sequencing. Consensus sequences were aligned using **MAFFT**, and maximum-likelihood phylogenies were reconstructed using **IQ-TREE**, with branch support assessed using 1,000 ultrafast bootstrap replicates and the approximate likelihood-ratio test. Putative transmission clusters were identified using **ClusterPicker**.

### Association and Mediation Analyses

- **Association analysis:** Associations between gut viral profiles and host phenotypes were assessed using **MaAsLin3**, modelling viral abundance and prevalence separately and adjusting for relevant host and technical covariates. Models were additionally adjusted for the abundance of the predicted bacterial host genus where available. 
- **Mediation analysis:** Causal mediation analysis was performed using the R package **mediation** to assess whether cytokine levels mediated associations between host phenotypes and phage presence. Models accounted for relevant covariates, including predicted bacterial host abundance. Mediation effects were estimated using 1,000 quasi-Bayesian simulations.

---

## Data Availability

The HIV gut viral catalog (HGVC), associated metadata, and viral abundance table are available via [**Zenodo**](https://doi.org/10.5281/zenodo.18981021):

| File | Description | Link |
|------|-------------|------|
| `HGVC.fa.gz` | HGVC vOTU representative viral genomes (n=23,996) | [Download](https://zenodo.org/record/18981021/files/LLNEXT_GV_catalog.fa.gz) |
| `HGVC_metadata.tsv` | Metadata for all HGVC species-level vOTUs | [Download](https://zenodo.org/record/18981021/files/LLNEXT_GV_catalog_metadata.tsv) |
| `HGVC_gut_viral_abundance_table.tsv` | Per-sample RPKM abundances of HGVC vOTUs | [Download](https://zenodo.org/record/18981021/files/LLNEXT_gut_viral_abundance_table.tsv) |

The datasets provide the inputs required to run the downstream analysis scripts in `scripts/F_Downstream_scripts/`

---
## Code Structure
All scripts used in this project can be found within scripts and are organized as follows:

```text
scripts/
├── A_Virus_Identification/
│   ├── Virus_detection/                  # geNomad
│   └── Contamination_quality_check/      # Quality control
│
├── B_Dereplication/                      # BLAST DB generation and ANI-based dereplication (aniclust.py, anicalc.py)
│
├── C_Abundance_estimation/               # Viral abundance estimation (read mapping)
│
├── D_Virus_Characterization/
│   ├── Genome_quality/                   # Viral genome completeness estimation
│   ├── Taxonomy/                         # Viral taxonomy (geNomad, ViTAP)
│   └── Host_and_lifestyle_prediction/    # Host prediction (iPHoP) and lifestyle (BACPHLIP)
│
├── E_Strain_sharing/
│   ├── Strain_level_profiling/            # Strain-level viral profiling (inStrain) 
│   └── X
│
└── F_Downstream_scripts/                   # Statistical analyses and figure generation

```
---

## Citation

If you use this resource in your research, please cite both the associated publication and the data resource.

### Publication

> **Gut virome alterations in people living with HIV associate with host phenotypes and immune function**  
> Asier Fernández-Pato*, Haoran Peng*, Ángela Del Castillo-Izquierdo, Victoria Rios-Vazquez, Kavita Mehta, Rob J. W. Arts, Jéssica C. dos Santos, Mareva Delporte, Janneke E. Stalenhoef, Casper Rokx, Marvin A. H. Berrevoets, Leo A. B. Joosten, Mihai G. Netea, Alexander Kurilshikov, Linos Vandekerckhove, Andre J. A. M. van der Ven†, Alexandra Zhernakova† & Jingyuan Fu†
> bioRxiv (2026).
> [https://doi.org/10.64898/2026.03.05.709809](https://doi.org/10.64898/2026.03.05.709809)

---

### Data resource

> Fernández Pato, A. (2026). **HIV gut viral catalog (HGVC)(v1.0)** [Data set]. Zenodo.  
> [https://doi.org/10.5281/zenodo.18981021](https://doi.org/10.5281/zenodo.18981021)

---

## Contact

For questions regarding the analysis workflow or datasets, please feel free to [open an issue](https://github.com/GRONINGEN-MICROBIOME-CENTRE/2000HIV_Gut_Virome/issues) or contact:

