
---
16S rRNA Microbiome Analysis in R

This project focuses on analyzing microbiome data based on sequencing of the V3–V4 regions of the 16S rDNA gene. The analysis was conducted using **R**, with key packages including `dada2`, `phyloseq`, and `ggplot2`. The main goal is to build a reproducible pipeline for raw data processing and taxonomy assignment.

Project Description
The input data comes from **Illumina NovaSeq 6000** sequencing and includes fecal samples collected from guinea pigs. The project consists of the following steps:

- Filtering and merging of paired-end reads
- Dereplication and detection of ASVs (Amplicon Sequence Variants)
- Taxonomic assignment using the SILVA database

---
Analiza mikrobiomu 16S rDNA w R
Projekt analizy danych mikrobiomu na podstawie sekwencjonowania regionów V3–V4 genu 16S rDNA. Analiza została przeprowadzona w języku **R** z wykorzystaniem bibliotek takich jak `dada2`, `phyloseq` i `ggplot2`. Głównym celem projektu jest opracowanie pipeline'u do przetwarzania danych surowych oraz przypisania taksonomii.


Opis projektu
Dane wejściowe pochodzą z sekwencjonowania platformą **Illumina NovaSeq 6000** i obejmują próbki kału pobrane od świnek morskich. Projekt składa się z następujących etapów:

- filtrowanie i łączenie odczytów parzowanych
- dereplikacja i wykrycie ASV (Amplicon Sequence Variants)
- przypisanie taksonomii na podstawie bazy **SILVA**
