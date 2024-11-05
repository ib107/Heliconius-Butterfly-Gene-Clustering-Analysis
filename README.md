# Heliconius Butterfly Gene Clustering Analysis

This project examines clustering patterns in *Heliconius* butterfly species by analyzing DNA sequences from two genes, COI (mitochondrial) and EF1a (nuclear), to understand potential ecological adaptations or behavioral traits, such as pollination behavior. Using clustering techniques and clustering metrics, the project investigates differences in clustering patterns across genes and their relation to ecological factors.

## Table of Contents
- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Data Sources](#data-sources)
- [Requirements](#requirements)
- [Usage](#usage)
- [Analysis Workflow](#analysis-workflow)
- [References](#references)

## Introduction
*Heliconius* butterflies exhibit various ecological adaptations, especially in relation to pollination. This project explores whether COI (cytochrome oxidase I) and EF1a (elongation factor 1-alpha) gene sequences cluster differently within *Heliconius* species. The analysis focuses on clustering strength using measures such as the Silhouette index, with ecological traits potentially incorporated to provide further context.

## Project Structure
- `Data/`: Contains raw DNA sequence data and processed alignment files. 
- `R Script/`: Script for sequence preprocessing, clustering analysis, and visualization.
- `Plots/`:  Significant figures, and visualizations.
- Assignment 2 PDF: Written evaluation and analysis.
- `README.md`: Project overview and instructions.

## Data Sources
1. **NCBI GenBank**: COI and EF1a gene sequences for *Heliconius* species.

## Requirements
- R (version 4.0 or higher)
- R packages: `rentrez`,`Biostrings`, `seqinr`, `stringr`, `dplyr`, `ggplot2`, `DECIPHER`, `ape`, `cluster`, `fpc`, `ggdendro`, `ggfortify`

## Usage
1. **Data Import**: Load COI and EF1a sequences along with any trait data.
2. **Data Exploration and Quality Control**: Filter sequences based on quality metrics such as sequence length and proportion of N bases.
3. **Clustering Analysis**: Perform k-means and hierarchical clustering and calculate Silhouette scores to compare. 
4. **Visualization**: Generate clustering visualizations and dendrograms to illustrate gene clustering patterns.

## Analysis Workflow
1. **Data Preprocessing**:
   - Import sequences and trait data.
   - Filter and align sequences, handling missing data as needed.

2. **Clustering and Metrics**:
   - Cluster sequences by gene using k-means and hierarchical clustering. 
   - Calculate internal clustering measures like Silhouette indices for pattern assessment.

3. **Visualization**:
   - Create dendrograms and other visualizations to compare clustering patterns between COI and EF1a genes.

## References
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
