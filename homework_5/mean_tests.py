import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests
import scipy.stats as st
import argparse


def check_intervals_intersect(first_ci, second_ci):
    if min(first_ci[1], second_ci[1]) - max(first_ci[0], second_ci[0]) <= 0:
        return False
    else:
        return True


def check_dge_with_ci(first_table, second_table):
    ci_test_results = []
    gene_names = []
    for gene_name in first_table.iloc[:, :-1].columns:
        gene_names.append(gene_name)
        first_intervals = st.t.interval(confidence=0.95, 
                                        df=len(first_table[gene_name]) - 1, 
                                        loc=np.mean(first_table[gene_name]),
                                        scale=st.sem(first_table[gene_name]))
        second_intervals = st.t.interval(confidence=0.95, 
                                df=len(second_table[gene_name]) - 1, 
                                loc=np.mean(second_table[gene_name]),
                                scale=st.sem(second_table[gene_name]))
        ci_test_results.append(not check_intervals_intersect(first_intervals, second_intervals))
    
    return gene_names, ci_test_results


def check_dge_with_ztest(first_table, second_table):
    z_test_results = []
    z_test_p_values = []
    for gene_name in first_table.iloc[:, :-1].columns:
        z_result = ztest(first_table[gene_name], second_table[gene_name])
        z_test_results.append(z_result[1] < 0.05)
        z_test_p_values.append(z_result[1])

    return z_test_results, z_test_p_values


def get_means_difference(first_table, second_table):
    mean_diff = []
    for gene_name in first_table.iloc[:, :-1].columns:
        md = np.mean(second_table[gene_name]) - np.mean(first_table[gene_name])
        mean_diff.append(md)
    return mean_diff


def correct_p_values(results, method):
    test_result, adj_p_values = multipletests(results['z_test_p_values'], method = method)[0:2]
    results['z_test_adj_p_values'] = adj_p_values
    results['z_test_results'] = test_result
    results = results[['gene_name', 'mean_diff', 'ci_test_results', 'z_test_p_values', 'z_test_adj_p_values', 'z_test_results']]
    return results


def get_pandas_dataframe(first_table, second_table):
    gene_names, ci_test_results = check_dge_with_ci(first_table, second_table)
    z_test_results, z_test_p_values = check_dge_with_ztest(first_table, second_table)
    mean_diff = get_means_difference(first_table, second_table)

    results = {
        "gene_name": gene_names,
        "mean_diff": mean_diff,
        "ci_test_results": ci_test_results,
        "z_test_p_values": z_test_p_values,
        "z_test_results": z_test_results,
    }

    results = pd.DataFrame(results)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('first_cell_type_expressions_path', type=str, help='Path to table with genes expressions for first cell type')
    parser.add_argument('second_cell_type_expressions_path', type=str, help='Path to table with genes expressions for second cell type')
    parser.add_argument('save_results_table', type=str, help='Results table name')
    parser.add_argument('--adjust', type=str, help='Set method to adjust p values')
    args = parser.parse_args()
    

    first_table = pd.read_csv(args.first_cell_type_expressions_path)
    second_table = pd.read_csv(args.second_cell_type_expressions_path)

    results = get_pandas_dataframe(first_table, second_table)

    if args.adjust:
        results = correct_p_values(results, method=args.adjust)
        
    results.to_csv(args.save_results_table, index=False)