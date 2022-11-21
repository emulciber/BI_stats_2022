# Homework 5 program

This program creates a table with CI intervals intersect and z-test results.

How to launching the program:
```
python mean_tests.py <path_to_first_table> <path_to_second_table> <output_path> --adjust <method>
```

Example:
```
python mean_tests.py data/b_cells_expression_data.csv data/nk_cells_expression_data.csv tests_results.csv --adjust bonferroni
```
The program has 3 required and 1 optional arguments:
* Positional argument 1    
Path to the table with genes expressions for first cell type.
* Positional argument 2   
Path to the table with genes expressions for second cell type.
* Positional argument 3   
The results table name.
* _--adjust_   
Set method to adjust z-test p-values (methods are used from function statsmodels.stats.multitest.multipletests).

The program creates a csv table with next columns:
* __gene_name__: names of genes from input tables
* __mean_diff__: the difference in means between two cell types (from first to second)
* __ci_test_results__: contains 'True' or 'False'. 'True' if 95% confidence intervals do not intersect, so the difference in gene expression between two cell types is significant
* __z_test_p_values__: p values of z-test
* __z_test_adj_p_values__: adjusted p values of z-test with adjust method set by user
* __z_test_results__: contains 'True' or 'False'. 'True' if z-test showed the difference in gene expression between two cell types is significant (p-value <0.05)

