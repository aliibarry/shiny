---
title: "datasetsummary"
author: "amb"
date: "22/08/2021"
output: html_document
---

Deep RNA-seq of mouse DRG subpopulations after injury was performed to interrogate subtype-specific and shared injury signatures. All work was done in accordance with the UK Home Office and the University of Oxford Policy on the Use of Animals in Scientific Research. These studies conform to ARRIVE guidelines. Specific transgenic details and methodologies will be available in the published reports, and are currently available by request to [Ali Barry](https://www.ndcn.ox.ac.uk/team/allison-barry). All work presented here is currently unpublished, and is the work of Ali Barry, Giorgos Baskozos, and Dave Bennett at the University of Oxford (NDCN).

---


**DRG subpopulation RNAseq**

Briefly, spared nerve injury (SNI) was performed on transgenic animals, each labelling one of five sensory neuron subtypes of interest (general nociceptors (*Scn10a*), peptidergic (PEP, *Calca*) and non-peptidergic (NP, *Mrgprd*) nociceptors, C-LTMRs (*Th*), and AB-RA + Ad-LTMRs (*Ntrk2*)). Paired ipsilateral and contralateral lumbar DRG samples were then collected 3 days and 4 weeks after surgery. In each condition, both male and female samples were collected to enable sex analyses. Together, 160 samples were collected, with 154 passing QC.

Spared nerve injury (SNI) is one of many neuropathic injury models used. Here, it was selected as a localized model of neuropathic pain with an internal control (contralateral, uninjured tissue). SNI is a behaviourally relevant animal model, with sensory dysfunction similar to that in neuropathic pain patients. 

Reads were mapped to the GRCm38 genome using STAR alignment and Samtools was used to sort, index, and merge BAM files. Quality control (QC) was performed with both FastQC and Samtools prior to gene counting with HTSeq. Count corrections for effective library sizes were performed in R using DESeq2, and normalized gene counts were fitted to a negative binomial distribution. Count transformations were performed using VST.

---
