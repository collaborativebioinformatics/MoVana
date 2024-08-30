import gseapy as gp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your list of Ensembl gene IDs that were produced by a bash script 
# Assuming your gene IDs are in a text file, one per line
gene_list_file = "/g/korbel/olisov/hackathon/del_genes.txt"
with open(gene_list_file) as f:
    gene_list = [line.strip() for line in f]

# Convert the list to a pandas DataFrame
genes_df = pd.DataFrame(gene_list, columns=['gene_id'])

# Perform Gene Set Enrichment Analysis using gseapy
# You can choose the library you want to use for enrichment, e.g., "KEGG_2019_Human"
enrichment_results = gp.enrichr(gene_list=genes_df['gene_id'].tolist(),
                                gene_sets='KEGG_2019_Human',
                                organism='Human',  # specify the organism, it can be 'Human', 'Mouse', etc.
                                outdir='/g/korbel/olisov/hackathon/enrichment_results',  # the output directory
                                cutoff=0.05)

# Extract the results DataFrame
results_df = enrichment_results.results

# Extract the number of genes from the 'Overlap' column
results_df['Num_Genes'] = results_df['Overlap'].apply(lambda x: int(x.split('/')[0]))

# Plot 1: Barplot with number of genes on x-axis and Adjusted P-value as color
plt.figure(figsize=(8, 6))
sns.barplot(
    x='Num_Genes', 
    y='Term', 
    data=results_df.head(10), 
    dodge=False, 
    palette="coolwarm", 
    hue='Adjusted P-value'
)
plt.title('Top Enriched Terms')
plt.xlabel('Number of Genes')
plt.legend(title='Adjusted P-value')
plt.savefig(f'/g/korbel/olisov/hackathon/GSEA.png', bbox_inches='tight', dpi=300)
plt.show()
