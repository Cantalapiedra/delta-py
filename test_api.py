from delta import compute_delta_pvalue

tree="../tree_trait_examples_claudia/cluster_138_tree.nw"
data="../tree_trait_examples_claudia/cluster_138_trait.csv"
columns=["cluster"]
threads=10
pval_reps=10
data_sep=","
id_index=0

delta, p_value = compute_delta_pvalue(tree, data, columns, threads, pval_reps, data_sep, id_index, verbose = True)

print(delta)
print(p_value)
