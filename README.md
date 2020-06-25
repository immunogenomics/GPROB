
# G-PROB

Multiple diseases can present with similar initial symptoms, making it difficult to clinically differentiate between these conditions. G-PROB uses patients' genetic profile to expidite the differentiation between multiple diseases. G-Prob is described by Knevel et al in Science Translation Medicine May 2020 https://stm.sciencemag.org/content/12/545/eaay1548/tab-e-letters for the situation of inflammatory arthritis. This genetic diagnostic tool can be applied to any situation of phenotypically similar diseases with different underlying genetics.

## This function requires the following input
- pheno_prev = a data frame with on each row of column 1 the  disease name and the prevalence of the diseases in column 2 
-  df_SNPs = a data frame with column 1: PT_id, columns 2 till k: headers are SNP_names, each field contains the number of minor allele a patient has (0,1 or 2) for that particular SNP. 
-  df_ORs = data frame column 1: SNP_names, column 2 till k: each columns representing the ORs for the minor allele of SNPs for one disease 
-  N = the number of iteration to find newtons theta (see equatino 3 in the paper's Materials and Methods)

### Calling the function 
run:
G_PROB(pheno_prev,df_SNPs, df_ORs, 20)


## Citation
If you use the G-Prob tool, please cite:

Knevel R, le Cessie S, Terao CC, et al. Using genetics to prioritize diagnoses for rheumatology outpatients with inflammatory arthritis. Sci Transl Med. 2020;12(545):eaay1548. doi:10.1126/scitranslmed.aay1548
