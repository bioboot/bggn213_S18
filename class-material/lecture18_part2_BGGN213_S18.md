---
layout: page
title: BGGN-213, Lecture 18 (Part 2)
---


Designing a personalized cancer vaccine
================

**BGGN-213 Lecture 18:**  
Barry Grant &lt; <http://thegrantlab.org> &gt;  
Date: 2018-05-31   (15:19:17 PDT on Thu, May 31)  
{:.message}  


**Notes**: To identify somatic mutations in a tumor, DNA from the tumor is sequenced and compared to DNA from normal tissue in the same individual using *variant calling algorithms*.  

Comparison of tumor sequences to those from normal tissue (rather than ‘the human genome’) is important to ensure that the detected differences are not germline mutations.    

To identify which of the somatic mutations leads to the production of aberrant proteins, the location of the mutation in the genome is inspected to identify non-synonymous mutations (i.e. those that fall into protein coding regions and change the encoded amino acid).  


### Section 1. Protein sequences from healthy and tumor tissue
The following sequences resulted from such an NGS analysis of patient healthy and tumor tissue. You can also download these sequences (named`P53_wt` and `P53_mutant`) as the following FASTA format sequence file [lecture18_sequences.fa]({{ site.baseurl }}//class-material/lecture18_sequences.fa){:.no-push-state}.

~~~
>P53_wt Cellular tumor antigen p53 - Healthy Tissue
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD

>P53_mutant Cellular tumor antigen p53 - Tumor Tissue
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
SCMGGMNRRPILTIITLEV
~~~

Your task is to identify tumor specific mutations that could potentially be used for vaccine development. 

> **Q1**: **Identify sequence regions that contain all 9-mer peptides that are only found in the tumor**.  Hint: You will need to first identify the sites of mutation in the above sequences and then extract the surrounding subsequence region. This subsequence should encompass all possible 9-mers in the tumor derived sequence. In other words extract the subsequence from 8 residues before and 8 residues after all point mutations in the tumor sequence.

You could simply align the above sequences using your favorite alignment method and then visually identify differences and extract the required subsequences. However, this apprach is error-prone and obviously will not scale to further examples (i.e. many additional proteins and 100s of mutations). 

If you are motivated to use your new R coding skills to address this question (and you should be!) then you may be interested in the following [example analysis]({{ site.baseurl }}//class-material/lecture18_part2_example/){:.no-push-state}. However, we encourage you to **try on your own first** before consulting this worked example.  


### Section 2. Patient HLA typing results and HLA binding prediction: 
The **HLA molecules** expressed in an individual can be identified based on DNA or RNA sequencing of normal tissue or by targeted HLA sequencing assays, and by comparison of the results to specialized databases listing all known human HLA alleles.

The following HLA typing results were obtained for our patient of interest.

~~~
HLA-A*02:01   
HLA-A*68:01 
HLA-B*07:02 
HLA-B*35:01
~~~


To prioritize which of the mutations in a tumor should be included in a vaccine, they can be scanned for those resulting in mutated peptides that bind HLA molecules of the patient with high affinity. This is done using **HLA binding algorithms** generated using machine learning algorithms trained on large sets of experimentally determined peptide:HLA binding data. We will here use algorithms developed as part of the **Immune Epitope Database (IEDB)** project hosted at the La Jolla Institute for Allergy and Immunology.  

See: IEDB HLA binding prediction website [http://tools.iedb.org/mhci/](http://tools.iedb.org/mhci/){:.no-push-state}.  

> **Q2**: **Identify 9-mer peptides in the identified sequence regions unique to the tumor that can be potentially presented to T cells.** Hint: Use the IEDB HLA binding prediction server above to identify the top ranked 9-mer peptides for each patient HLA (see above for HLA typing results).  

 
> **Q3**: **Identify the top peptide for each patient HLA allele (see above)**. Hint: You can download a CSV formated result file for all predictions and use R or a spreadsheet application to identify the top ranked peptides for each allele. Which approach would you rather use if you were going to scale to 100s of analysis tasks?


### Section 3. Identifying tumor specific peptides
Any peptide resulting from a mutation should be checked if it is present in other proteins of the human genome, which would mean that it is not specific to the tumor. 

> **Q4**: **Are any of your top ranked peptides unique to the tumor?** Hint: Use [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins){:.no-push-state} to check your top identified peptides against Human sequences in the NR database.

 
> **Q5**: **Which peptide would you choose for potential vaccine design?** 


