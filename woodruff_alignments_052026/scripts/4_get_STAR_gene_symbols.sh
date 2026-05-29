## this is a very slow and stupid way to do this, 
## i am only doing this because i don't want to instal R or python in this docker

bash helper_scripts/ensembl_to_symbol.sh \
      ../abundances/STAR/counts_matrix.tsv > ../abundances/STAR/counts_matrix_genesymbols.tsv