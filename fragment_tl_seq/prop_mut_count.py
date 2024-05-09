def prop_mutation_count(woc, c, woc_total, c_total, treatment_kind=""):

    # function to get the proprotion of all the mutations in the reads

    """
    The functions takes in the filtered reads from sample without treatment(woc), with treatment
    (c), total read without treatment (woc_total) and with treatment (c_total)

    """
    woc_mut = woc["ref|query"].to_list()
    woc_mut = [mut.split(",") for mut in woc_mut]
    woc_mut = sum(woc_mut, [])  # flattent the list of lists

    c_mut = c["ref|query"].to_list()
    c_mut = [t.split(",") for t in c_mut]
    c_mut = sum(c_mut, [])  # flattent the list of lists

    def count_mutation_types(x):
        df = pd.DataFrame({"mutations": x, "S/N": np.arange(0, len(x))})
        df_mut_counts = df["mutations"].value_counts().to_frame()
        return df_mut_counts

    woc_mut_counts = count_mutation_types(woc_mut)
    c_mut_counts = count_mutation_types(c_mut)

    # calculate the proportion
    woc_mut_counts = woc_mut_counts["count"] / woc_total
    c_mut_counts = c_mut_counts["count"] / c_total

    # merge the woc and the c_mut_counts
    merged_counts = pd.concat([woc_mut_counts, c_mut_counts], axis=1, keys=["woc", "c"])
    merged_counts = merged_counts.reset_index()
    merged_counts = merged_counts.dropna()

    return merged_counts
