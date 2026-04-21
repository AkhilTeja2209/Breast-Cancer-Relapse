# Breast-Cancer-Relapse
Analysing breast cancer relapse free survival of patients.

This reository allows you to analyse the relapse free survival of breast cancer patients, tested on patients tracked for upto 160 months after their first case of breast cancer. The data has been taken from [GSE2034](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034) databse.

Place the .tsv file and the R file in the same folder and run the R file using Rstudio or any other IDE that supports R code. The code should make a DEGs file for the data (GSE2034_DEGs.csv) and the final output should have 4 files made in the same folder; Volcano plot, KEGG enrichment, GO enrichment and Heatmap of top 50 samples, all in .png format. 

The sample outputs are also attached in the repository.

#Attaching the inference for more understanding: [BreastCancer_GSE2034_Analysis_final.pdf](https://github.com/user-attachments/files/26942868/BreastCancer_GSE2034_Analysis_final.pdf)

#NOTE:
For a better result, it is recommended that the .tsv file is created by the user itself. Steps to do so:

1. Go to the GEO Series page for GSE2034. (Breast cancer relapse free survival with relapse, brain relapse and no relapse). 
2. Click “Analyze with GEO2R”.
3. In GEO2R: Define groups (e.g., relapse vs no relapse).
4. Keep default log2 transform (auto-detect). Use Benjamini & Hochberg for multiple testing correction.
5. Click “Top 250” or “Analyze”.
6. Download full results (all genes) as a table (.tsv).

This .tsv file can replace the given .tsv file in the repository, which only analyses the top 250 values.
