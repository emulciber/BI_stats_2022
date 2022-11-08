# Homework 5 program

This program creates a table with CI intervals intersect and z-test results.

Example of launching the program:
```
python mean_tests.py data/b_cells_expression_data.csv data/nk_cells_expression_data.csv tests_results.csv
```
The program requires 3 arguments:
* Argument 1: Path to the table with genes expressions for first cell type.
* Argument 2: Path to the table with genes expressions for second cell type.
* Argument 3: The results table name.

The program creates a csv table with next columns:
* __gene_name__: names of genes from input tables
* __ci_test_results__: contains 'True' or 'False'. 'True' if 95% confidence intervals do not intersect, so the difference in gene expression between two cell types is significant
* __z_test_results__: contains 'True' or 'False'. 'True' if z-test showed the difference in gene expression between two cell types is significant (p-value <0.05)
* __z_test_p_values__: p values of z-test
* __mean_diff__: the difference in means between two cell types (from first to second)