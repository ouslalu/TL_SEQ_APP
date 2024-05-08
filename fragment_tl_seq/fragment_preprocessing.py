import os
import pandas as pd
import numpy as np
from Bio import pairwise2
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import itertools
import pysam
import argparse


def file_preparation(pear_merged_loc, metadata, folder_location, base_quality):
    """
    preapares all the needed directory

    Arguments:
        pear_merged_loc (str): location of the read input
        metadata_location (csv): metadata with the information about the inputs and treatments
        folder_location (dir): location
        base_quality (dir): location of the base quality txt files.

    return metadata pandas dataframe


    """
    if not os.path.exists(pear_merged_loc):
        raise ValueError("Input folder not found")

    # change working directory to the location where the files are
    os.chdir(folder_location)

    for directory in ("plots", "files", "summary"):
        if not os.path.exists(directory):
            os.mkdir(directory)

    # labelling_summary to get the information on the samples.
    global labelling_summary_df
    labelling_summary_df = metadata
    # print(labelling_summary_df)


def process_pos_ref_query(x):

    # Remove numbers using regex
    x_1 = re.sub(r"\d+\|", "", x)

    return x_1


def get_reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",
    }
    reverse_complement = ""

    for nucleotide in sequence[::-1]:
        if nucleotide in complement:
            reverse_complement += complement[nucleotide]
        else:
            reverse_complement += nucleotide

    return reverse_complement


def check_primers(x, f_primer, r_primer):
    """
    This function takes the nucleotide fragment x and check if the forward primer and the reverse primer is in the sequence
    """
    # print(f_primer, x[0:len(f_primer)], r_primer, x[len(x)-len(r_primer):len(x)])
    if x[0 : len(f_primer)] == f_primer:
        if x[len(x) - len(r_primer) : len(x)] == r_primer:
            return True
    else:
        return False


def get_mismatches(read, product_sequence):
    # This functions get the list and the mismatches between the product sequence and individual read.
    # pos = position, product_nucelotide, read_nucleotide
    mismatch_list = []
    for i in range(0, len(product_sequence), 1):
        # print(i)
        mismatch = {}
        if product_sequence[i] != read[i]:
            mismatch["pos"] = i
            mismatch["product_nucleotide"] = product_sequence[i]
            mismatch["read_nucleotide"] = read[i]
            mismatch_list.append(mismatch)
    return mismatch_list


def read_length_plot(read_length_prop_df, sna):
    """
    plots the distribution of read length for the sequence.
    read_length_prop_df: pandas dataframe with read lenght and proportion columns
    sna: Sample name

    """

    try:
        if "read_length" and "%_reads" in read_length_prop_df.columns.to_list():
            # check if the read_length and % reads in the dataframe
            plt.figure(figsize=(4, 4))

            read_length_prop_df = read_length_prop_df.sort_values(by=["read_length"])
            # print(read_length_prop_df)

            # Calculate the maximum read length proportion and determine the index
            max_length = max(read_length_prop_df["%_reads"])
            max_length_index = (
                read_length_prop_df["%_reads"].to_list().index(max_length)
            )

            # print(max_length_index)
            # Create a color list, default color is 'grey', highlight color is 'blue'
            colors = ["grey"] * len(read_length_prop_df)
            colors[max_length_index] = "#296EB4"

            # print(colors)
            read_length_prop_df["read_length"] = read_length_prop_df[
                "read_length"
            ].astype("object")
            ax = sns.barplot(
                x="read_length",
                y="%_reads",
                data=read_length_prop_df,
                palette=colors,
                edgecolor="black",
                lw=2,
            )
            # print(read_length_prop_df.dtypes)

            for axis in ["bottom", "left"]:
                ax.spines[axis].set_linewidth(2.5)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            plt.xticks(rotation=30)
            plt.xlabel("read length")
            plt.ylabel("Proportion (%)")
            plt.title(sna.replace("_", " "))
            plt.savefig(f"plots/{sna}_trial.png", dpi=250, bbox_inches="tight")
            plt.close()

    except KeyError:
        print(f"dataframe does not have the read_length and the %_reads columns")


def analysis_for_nacent_nucleotides(sna):
    """
    Add docstring to this function


    """
    # call the function that return

    sequenced_df = pd.read_csv(
        f"{pear_merged_loc}/{sna}_merged.txt",
        sep="\t",
        header=None,
        encoding="unicode_escape",
    )

    # get the the base quality files
    base_qual = pd.read_csv(
        f"{base_quality}/{sna}_merged_base_quality.txt",
        sep="\t",
        header=None,
    )

    base_qual.columns = ["base_quality"]
    sequenced_df.columns = ["sequence"]

    try:
        # This is to catch error resulting from the total reads not equal to total lines of base qualities

        t_read = sequenced_df.shape[0]
        qual = base_qual.shape[0]

        print(t_read, qual)
        if t_read == qual:
            sequenced_df["base_quality"] = base_qual["base_quality"]

    except ValueError:
        return f"Error: The number of reads({t_read}) and base qualities({qual}) are different.Check the fastq files again"

    # Number of reads
    total_reads = sequenced_df.shape[0]

    # because the sample name in the labbeling summary df does not have the S.* behind out fatsq, I replaced it with ""
    # sna_e = re.sub("_S.*", "", sna)
    sna_e = sna

    product_sequence = (
        labelling_summary_df[labelling_summary_df["Sample name"] == sna_e][
            "Product sequence"
        ]
        .to_list()[0]
        .upper()
    )

    # length of the reads.
    df_read = sequenced_df["sequence"].apply(lambda x: len(x))
    df_read = df_read.value_counts(normalize=True).to_frame()
    df_read["%_reads"] = df_read["proportion"] * 100
    df_read["read_length"] = df_read.index
    df_read = df_read[["read_length", "%_reads"]]
    df_read = df_read.reset_index(drop=True)
    df_read["read_length"] = df_read["read_length"].astype("category")
    # plot the top 5
    df_top_read = df_read.head()

    read_length_plot(df_top_read, sna)

    # return df_top_read

    # filter the reads with read length = 148 and has both the forward and the reverse primer
    sequenced_df["read_length"] = sequenced_df["sequence"].apply(lambda x: len(x))
    sequenced_df_filtered = sequenced_df[
        sequenced_df["read_length"] == len(product_sequence)
    ]
    # sequenced_df_filtered

    # function to check if the primers are in the sequence
    f_primer = (
        labelling_summary_df.loc[labelling_summary_df["Sample name"] == sna_e, "F"]
        .values[0]
        .upper()
    )
    r_primer = (
        labelling_summary_df.loc[labelling_summary_df["Sample name"] == sna_e, "R"]
        .values[0]
        .upper()
    )
    r_primer = get_reverse_complement(r_primer)
    # print(f_primer,r_primer)

    sequenced_df_filtered["fr_primers"] = sequenced_df_filtered["sequence"].apply(
        lambda x: check_primers(x, f_primer, r_primer)
    )
    sequenced_df_filtered = sequenced_df_filtered[
        sequenced_df_filtered["fr_primers"] == True
    ]

    sequenced_df_filtered["mismatches"] = sequenced_df_filtered["sequence"].apply(
        lambda read: get_mismatches(read, product_sequence)
    )

    mismatches = sequenced_df_filtered["mismatches"].to_list()

    processed_mismatches = []
    for mis in mismatches:
        if len(mis) == 0:
            processed_mismatches.append("")
        else:
            mismatch_string_list = []
            for mis_dict in mis:  # get the mismatches in the dictionary
                mismatch_string = f"{mis_dict['pos']}|{mis_dict['product_nucleotide']}|{mis_dict['read_nucleotide']}"
                mismatch_string_list.append(mismatch_string)
            mismatch_sting_list = ",".join(mismatch_string_list)
            processed_mismatches.append(mismatch_sting_list)

    processed_mismatches = processed_mismatches

    # return processed_mismatches

    # selecting reads all mutations
    sequenced_df_filtered["pos|ref|query"] = processed_mismatches
    c_con_mutation = sequenced_df_filtered[
        ["sequence", "base_quality", "pos|ref|query"]
    ]

    # filter to retain mutations that pas. phred score

    # def filter_by_phred_score
    def phred_to_quality(ascii_char):
        return ord(ascii_char) - 33

    ascii_to_quality = {chr(i): phred_to_quality(chr(i)) for i in range(33, 127)}

    # Convert the dictionary to a DataFrame
    df_ascii = pd.DataFrame(
        list(ascii_to_quality.items()), columns=["ASCII Character", "Quality Score"]
    )

    sequenced_df_mutation_mut = sequenced_df_filtered[
        sequenced_df_filtered["pos|ref|query"] != ""
    ]
    mutations = sequenced_df_mutation_mut["pos|ref|query"].to_list()
    # print("muatations")
    # print(mutations[1:10])
    sequenced_df_mutation_mut.head()
    # mutations = [m for m in mutations if m != ''] #keep only mutations

    # make list of base quality
    base_quality_list = sequenced_df_mutation_mut[
        "base_quality"
    ].to_list()  # base_qualtiy

    # len(mutations)

    mut_count = []
    phred_score_test = []
    mutation_q_score = []

    for list_index in range(0, len(mutations), 1):
        # getting the position of the mutation and count the number of the mutations
        string_mut = mutations[list_index]
        mut_slash = string_mut.count("|")  # one mutation should have 2 slashes
        mut_count.append(mut_slash / 2)
        string_mut_list = string_mut.split(",")
        # print(mut_count)
        # print(string_mut_list)
        bq_score_list = []
        q_score = []
        for i in string_mut_list:
            pos = int(i.split("|")[0])  # get the position of the mutation
            # print(pos)
            bq = base_quality_list[list_index][
                pos
            ]  # get the base quality at the position
            bq_score = df_ascii[df_ascii["ASCII Character"] == bq][
                "Quality Score"
            ].to_list()[
                0
            ]  # get the score of the base quality from the ascii dictionary
            q_score.append(str(bq_score))
            if bq_score > 27:
                bq_score_list.append("Y")
            else:
                bq_score_list.append("N")
        bq_score_st = "".join(bq_score_list)  # convert the list to string
        q_score_st = "|".join(q_score)
        phred_score_test.append(
            bq_score_st
        )  # append the phred test to list phred score tend list
        mutation_q_score.append(q_score_st)

    sequenced_df_mutation_mut["mutation_phred_score_test"] = phred_score_test
    sequenced_df_mutation_mut["Num_of_mutations"] = mut_count
    sequenced_df_mutation_mut["Q_score"] = mutation_q_score

    # filter TC mutations with phred score greater than 27. i.e mutation_phred_score_test > 27
    phred_filtered_mutations = sequenced_df_mutation_mut[
        sequenced_df_mutation_mut["mutation_phred_score_test"].str.contains("Y")
    ]
    phred_filtered_mutations = phred_filtered_mutations.reset_index(drop=True)

    phred_filtered_mutations["ref|query"] = phred_filtered_mutations[
        "pos|ref|query"
    ].apply(process_pos_ref_query)

    return total_reads, phred_filtered_mutations, df_top_read


# t_reads, phred_filtered_mutations, df_top_read = analysis_for_nacent_nucleotides("c_24h_exon2_fr1_S33")


def count_mutation_types(x):
    # function counts different mutations counts and take in a list
    # check if argument is list and it is not empty
    assert isinstance(x, list), f"{x} is not a list"
    assert len(x) > 1, f"{x} is an empty list"

    df = pd.DataFrame({"mutations": x, "S/N": np.arange(0, len(x))})
    df_mut_counts = df["mutations"].value_counts().to_frame()
    return df_mut_counts


def prop_mutation_count(sna):

    # function to get the proprotion of all the mutations in the reads

    """
    The functions takes in sample name and calculate the proportion of all mutations in its read
    """

    t_reads, phred_filtered_mutations, df_top_read = analysis_for_nacent_nucleotides(
        sna
    )
    sample_mut = phred_filtered_mutations["ref|query"].to_list()
    sample_mut = [mut.split(",") for mut in sample_mut]
    sample_mut = sum(sample_mut, [])  # flatten the list of lists

    # count the number of different mutations
    sample_mut_counts = count_mutation_types(sample_mut)

    # calculate the proportion
    sample_mut_prop = sample_mut_counts["count"] / t_reads

    return sample_mut_prop


# sample_mut_prop = prop_mutation_count("c_24h_exon2_fr1_S33")
# print(sample_mut_prop)


def query(x):
    # mutation wrangling
    x = x.split(",")
    tc_count = 0

    for i in x:
        if "T|C" in i:
            tc_count += 1

    return tc_count


def count_number_of_tc_mutations(sna):
    """
    This function count the number of tc muations
    """
    t_reads, phred_filtered_mutations, df_top_read = analysis_for_nacent_nucleotides(
        sna
    )

    phred_filtered_mutations["TC_counts"] = phred_filtered_mutations["ref|query"].apply(
        query
    )
    df_tc_mutations = phred_filtered_mutations[
        phred_filtered_mutations["TC_counts"] > 0
    ]
    df_tc_mutations = (
        df_tc_mutations["TC_counts"].value_counts().to_frame(name="num_tc_mutations")
    )

    # divide the number of reads with TC_mutations by total reads after droping mutations with less than 10 reads on read count

    df_prop_tc_mutations = df_tc_mutations[df_tc_mutations["num_tc_mutations"] > 10]
    df_prop_tc_mutations = df_tc_mutations["num_tc_mutations"] / t_reads
    df_prop_tc_mutations = df_prop_tc_mutations.to_frame()

    df_prop_tc_mutations = df_prop_tc_mutations.reset_index()
    df_prop_tc_mutations.columns = ["Number of T|C counts", "Proportion of reads"]

    return df_prop_tc_mutations


def count_number_of_different_pos_tc_mutations(region, time):
    """
    Count all the positions T|C mutations in sample

    Argument:
        region: sample name
        time
    return:
        None
        make the plots for the T|C mutations

    """

    c_sample = labelling_summary_df[
        (labelling_summary_df["Region"] == region)
        & (labelling_summary_df["Time"] == time)
        & (labelling_summary_df["Treatment"] == "c")
    ]["Sample name"].iloc[0]

    woc_sample = labelling_summary_df[
        (labelling_summary_df["Region"] == region)
        & (labelling_summary_df["Time"] == time)
        & (labelling_summary_df["Treatment"] == "woc")
    ]["Sample name"].iloc[0]

    t_reads_c, phred_filtered_mutations_c, df_top_read_c = (
        analysis_for_nacent_nucleotides(c_sample)
    )

    t_reads_woc, phred_filtered_mutations_woc, df_top_read_woc = (
        analysis_for_nacent_nucleotides(woc_sample)
    )

    phred_filtered_mutations_c["pos|ref|query_count"] = phred_filtered_mutations_c[
        "pos|ref|query"
    ].apply(query)

    phred_filtered_mutations_c = phred_filtered_mutations_c[
        phred_filtered_mutations_c["pos|ref|query_count"] > 0
    ]
    process_pos_ref_query_c = phred_filtered_mutations_c["pos|ref|query"].to_list()

    phred_filtered_mutations_woc["pos|ref|query_count"] = phred_filtered_mutations_woc[
        "pos|ref|query"
    ].apply(query)
    phred_filtered_mutations_woc = phred_filtered_mutations_woc[
        phred_filtered_mutations_woc["pos|ref|query_count"] > 0
    ]
    process_pos_ref_query_woc = phred_filtered_mutations_woc["pos|ref|query"].to_list()

    process_pos_ref_query_c = [t.split(",") for t in process_pos_ref_query_c]
    process_pos_ref_query_c = sum(process_pos_ref_query_c, [])

    process_pos_ref_query_woc = [t.split(",") for t in process_pos_ref_query_woc]
    process_pos_ref_query_woc = sum(process_pos_ref_query_woc, [])

    process_pos_ref_query_c_df = pd.DataFrame(
        {"pos|ref|query": process_pos_ref_query_c}
    )
    process_pos_ref_query_woc_df = pd.DataFrame(
        {"pos|ref|query": process_pos_ref_query_woc}
    )

    process_pos_ref_query_c_df = (
        process_pos_ref_query_c_df["pos|ref|query"]
        .value_counts()
        .to_frame()
        .reset_index()
        .head(10)
    )
    process_pos_ref_query_c_df["prop"] = process_pos_ref_query_c_df["count"] / t_reads_c

    process_pos_ref_query_woc_df = (
        process_pos_ref_query_woc_df["pos|ref|query"]
        .value_counts()
        .to_frame()
        .reset_index()
        .head(10)
    )
    process_pos_ref_query_woc_df["prop"] = (
        process_pos_ref_query_woc_df["count"] / t_reads_woc
    )
    process_pos_ref_query_c_df["treatment"] = ["c"] * process_pos_ref_query_c_df.shape[
        0
    ]
    process_pos_ref_query_woc_df["treatment"] = [
        "woc"
    ] * process_pos_ref_query_woc_df.shape[0]

    c_woc_df = pd.concat(
        [process_pos_ref_query_c_df, process_pos_ref_query_woc_df], axis=0
    )

    plt.figure(figsize=(10, 4))
    ax = sns.barplot(
        data=c_woc_df,
        x="pos|ref|query",
        y="prop",
        hue="treatment",
        edgecolor="black",
        lw=2,
    )

    ax.spines["left"].set_linewidth(2.5)
    ax.spines["bottom"].set_linewidth(2.5)

    ax.spines["top"].set_edgecolor("white")
    ax.spines["right"].set_edgecolor("white")
    for i in ["top", "right"]:
        ax.spines[i].set_visible("False")

    # plt.yscale("log")
    plt.title(f"{region}_{time}")
    plt.savefig(f"plots/pos_tc_{region}_{time}.png", dpi=250, bbox_inches="tight")
    # print(c_woc_df)
    plt.close()


def tc_mutations_plot(df, **kwargs):
    """
    plots the proportion of TC mutations in the control(woc) and treated samples (c)

    Arguments:
        df: dataframe with 'number of T|C counts', 'Proportion of reads' and treatment.


    returns None
    """
    segment = kwargs.get("segment", "")
    time = kwargs.get("time", "")

    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame.")

    required_columns = ["Number of T|C counts", "Proportion of reads", "treatment"]
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"DataFrame missing required columns: {', '.join(missing_columns)}"
        )

    plt.figure(figsize=(10, 4))
    ax = sns.barplot(
        data=df,
        x="Number of T|C counts",
        y="Proportion of reads",
        hue="treatment",
        edgecolor="black",
        lw=2,
    )

    ax.spines["left"].set_linewidth(2.5)
    ax.spines["bottom"].set_linewidth(2.5)

    ax.spines["top"].set_edgecolor("white")
    ax.spines["right"].set_edgecolor("white")
    for i in ["top", "right"]:
        ax.spines[i].set_visible("False")

    # plt.yscale("log")
    plt.title(f"{segment}_{time}")

    # create the plot directory if it does no exist
    if not os.path.exists("plots"):
        os.mkdir("plots")

    plt.savefig(f"plots/{segment}_{time}.png", dpi=250, bbox_inches="tight")
    plt.close()


def number_of_tc_mutations_plot(region, time):
    """
    function plots the proportions of different number of T|C mutations between converted (c) and uconverted (woc)

    Arguments:
        region (str): the name of the region in comparison
        time (str): hour of treatment

    returns None
    """

    c_sample = labelling_summary_df[
        (labelling_summary_df["Region"] == region)
        & (labelling_summary_df["Time"] == time)
        & (labelling_summary_df["Treatment"] == "c")
    ]["Sample name"].iloc[0]

    woc_sample = labelling_summary_df[
        (labelling_summary_df["Region"] == region)
        & (labelling_summary_df["Time"] == time)
        & (labelling_summary_df["Treatment"] == "woc")
    ]["Sample name"].iloc[0]

    df_prop_tc_mutation_c = count_number_of_tc_mutations(c_sample)
    df_prop_tc_mutation_woc = count_number_of_tc_mutations(woc_sample)

    df_prop_tc_mutation_c = df_prop_tc_mutation_c.reset_index()
    df_prop_tc_mutation_woc = df_prop_tc_mutation_woc.reset_index()

    df_prop_tc_mutation_c["treatment"] = ["c"] * df_prop_tc_mutation_c.shape[0]
    # print(["c"] * df_prop_tc_mutation_c.shape[0])
    # print(df_prop_tc_mutation_c)
    df_prop_tc_mutation_woc["treatment"] = ["woc"] * df_prop_tc_mutation_woc.shape[0]

    df_prop_tc_mutation_c_woc = pd.concat(
        [df_prop_tc_mutation_c, df_prop_tc_mutation_woc], axis=0
    )
    # print(df_prop_tc_mutation_c_woc)
    tc_mutations_plot(df_prop_tc_mutation_c_woc, segment=region, time=time)


# number_of_tc_mutations_plot("exon2", "24h")


def get_region_from_snps_vcf(vcf_file, chromosome, start, end):
    """
    This function identifies SNPs in the reference sequence.

    Argunments:
        vcf_loc (str): location of the vcf file with the SNPs
        chromosome (str): chromosome which the gene is located
        start (int): starting position of the gene on the chromosome
        end (int): end position of the gene on the chromosome

    returns a dataframe with SNPs in that region of genes.

    """

    chromosomes = [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        "X",
        "x",
        "Y",
        "y",
        "MT",
        "mt",
    ]

    # check if the inputs are valid
    if chromosome not in chromosomes:
        raise ValueError(f"{chromosome} not a valid chromosome: [1-22, X, Y]")

    for inp in [start, end]:
        if not isinstance(inp, int):
            raise ValueError("start and end arguments can only be an integer")

    f = open("SNP/snp_in_region.txt", "w")

    # use pysam to work with the vcf files
    vcf_file = pysam.VariantFile(vcf_file)

    # header
    for line in str(vcf_file.header).split("\n"):
        if line.startswith("#CHROM"):
            f.write(line)
            f.write("\n")
            break

    # get the region of interest
    chromosome = str(chromosome)
    for reg in vcf_file.fetch(chromosome, start - 1, end):
        f.write(str(reg))

    f.close()

    # put the
    snp_df = pd.read_csv("SNP/snp_in_region.txt", sep="\t")

    return snp_df


def reference_reformating(reference_loc, gene_direction="+"):
    """
    This function reformat the reference to have the location same location designation as
    the reference
    NB:  The chromosome numbering follows the 5' to 3'. Therefore, if the gene_direction is 3' - 5',
    the function does reverse complentary before numbering

    Arguments:
        reference_loc (str):  the txt_location of the reference sequences
        gene_direction (str): direction of the gene expression. The function expects "+" or "-"

    returns:
        formatted reference

    """
    if gene_direction not in ["+", "-"]:
        raise ValueError("gene_direction can only be '+' or '-'")

    reference = open(reference_loc, "r")

    for line in reference:
        if gene_direction == "-":
            # print(line)
            reference_ = get_reverse_complement(line)

        if gene_direction == "+":
            reference_ = line

    return reference_.upper()


def reference_coordinate(
    vcf="", chromosome="", start=0, end=0, reference_loc="", gene_direction="+"
):
    """
    This function labels give position to every nucleotide in the
    reference sequence that it can be compared to the nucleotides
    in the snp from the get_region_from_snps_vcfb functions

    Arguments:
        vcf (str): the directory containing the vcf file for the SNPs.
        reference_loc (str) =  the text file with the reference sequence
        start (int) = starting position of the gene on the chromosome
        end (int) = end position of the gene on the chromosome
        gene_direction (str): ["+", "-"]. Gives the direction in which the gene is expressed

    returns:
        dataframe of SNPs found in the reference gene

    """

    # get the formatted reference and snp_df from the vcf
    reference_ = reference_reformating(reference_loc, gene_direction)
    snp_df = get_region_from_snps_vcf(vcf, chromosome, start, end)

    # get snp positions
    positions = snp_df["POS"].to_list()

    pos = [i for i in range(start, end + 1, 1)]
    nuc = list(reference_)

    # print(len(pos))
    # print(len(nuc))

    ref_df = pd.DataFrame({"pos": pos, "nucleotide": nuc})

    # nucleotides positions in reference that snps have been identified for
    ref_df_2 = ref_df.set_index("pos")

    ref_snp = ref_df_2.loc[positions, :]

    ref_snp.reset_index(inplace=True)

    ref_snp_df = pd.concat([ref_snp, snp_df], axis=1)

    # print(ref_snp_df)

    return ref_snp_df


def check_the_snp_position_in_samples(
    vcf="",
    reference_loc="",
    region="",
    chromosome=None,
    start=0,
    end=0,
    gene_direction="+",
    labelling_summary_df=None,
):
    """
    Identifies SNP in the sample

    arguments:
        vcf (str): the directory containing the vcf file for the SNPs.
        reference_loc =  the text file with the reference sequence
        region (str): region name
        chromosome (int or string): chromosome number [1-22, X, Y]
        start (int): the starting position of the reference on the human genome
        end (int): the ending position of the reference on the human genome
        gene_direction (str): ["+", "-"]. Gives the direction in which the gene is expressed
        labbeling_summary_df (panda DataFrame): dataframe of the metadata.

    returns:
        product sequence with snps in lower case
        list of the SNPs in the sample sequence
    """

    # get the treated sample name to use to check the SNP for the region
    sample = labelling_summary_df.loc[
        (labelling_summary_df["Region"] == region)
        & (labelling_summary_df["Treatment"] == "c"),
        "Sample name",
    ].values[0]

    # get the forward prime
    f_primer = (
        labelling_summary_df.loc[labelling_summary_df["Sample name"] == sample, "F"]
        .values[0]
        .upper()
    )

    # get the reverse primer
    r_primer = (
        labelling_summary_df.loc[labelling_summary_df["Sample name"] == sample, "R"]
        .values[0]
        .upper()
    )

    # the product sequence (i.e amplified region)
    product_sequence = (
        labelling_summary_df.loc[
            labelling_summary_df["Sample name"] == sample, "Product sequence"
        ]
        .values[0]
        .upper()
    )

    # get reference reformatted
    reference = reference_reformating(reference_loc, gene_direction)

    # change the direction of the sequence if the gene is expressed from 3' to 5' end
    if gene_direction == "-":
        f_primer = get_reverse_complement(f_primer)
        r_primer = get_reverse_complement(r_primer)
        product_sequence = get_reverse_complement(product_sequence)

    # get the position of the sample sequence in the reference
    f_primer_index = reference.find(f_primer)

    # This list contains the beginning and the end position of the sample sequence on the chromosome.
    first_nuc_pos = start + f_primer_index
    last_nuc_pos = start + f_primer_index + len(product_sequence)
    sample_region_on_the_chromosome = [i for i in range(first_nuc_pos, last_nuc_pos, 1)]

    # print(sample_region_on_the_chromosome)

    # get SNPs
    ref_snp_df = reference_coordinate(
        vcf, chromosome, start, end, reference_loc, gene_direction
    )

    # save ref_snp_df
    if not os.path.exists("SNP"):
        os.mkdir("SNP")
    ref_snp_df.to_csv(f"SNP/snp_in_gene_region.csv")

    # make dataframe with the SNP positions and the reference nucleotide (from the VCF file)
    snp_positions = ref_snp_df["pos"].to_list()
    snp_reference = ref_snp_df["REF"].to_list()

    snp_pos_ref_dict = {"snp_positions": snp_positions, "snp_reference": snp_reference}
    snp_pos_ref_dict_df = pd.DataFrame(snp_pos_ref_dict)

    # list of SNPs found in sample and their coressponding position on human reference genome
    pos_on_sequence_list = []

    # check if the first position of the sample is within the SNP positions identified in the human genome reference
    if sample_region_on_the_chromosome[0] >= snp_positions[0]:
        for reg in sample_region_on_the_chromosome:
            if reg in snp_positions:
                # reference from the snp
                ref_from_snp = snp_pos_ref_dict_df[
                    snp_pos_ref_dict_df["snp_positions"] == reg
                ]["snp_reference"].to_list()[0]
                # print(ref_from_snp)

                # get the position of the SNP on the sample and not on the human genome reference
                pos_on_sequence = reg - sample_region_on_the_chromosome[0]

                # format SNP identified in the sample as :
                # 'SNP_position_on_human_reference_genome|nucelotide_on_human_genome_reference|SNP_position_in_sample_sequence|nucelotide_in_sample
                ref_snp_sample = f"{reg}|{ref_from_snp}|{pos_on_sequence}|{product_sequence[pos_on_sequence]}"

                # add all identified SNP to a list
                pos_on_sequence_list.append(ref_snp_sample)

                # turn SNP to lowercase to distinguish from other positions.
                product_sequence_list = list(product_sequence)

                product_sequence_list[pos_on_sequence] = product_sequence_list[
                    pos_on_sequence
                ].lower()

                product_sequence = "".join(product_sequence_list)

    # print(product_sequence)
    # print(pos_on_sequence_list)

    return product_sequence, pos_on_sequence_list


def summary_report(
    vcf,
    reference_loc,
    region,
    chromosome,
    start,
    end,
    gene_direction,
    labelling_summary_df,
):
    """
    This reports the result of the analysis

    """

    # get the product sequence showing SNPs in lowercase and the positions of the SNPs
    product_sequence, pos_on_sequence_list = check_the_snp_position_in_samples(
        vcf,
        reference_loc,
        region,
        chromosome,
        start,
        end,
        gene_direction,
        labelling_summary_df,
    )

    if not os.path.exists("summary"):
        os.mkdir("summary")

    with open(f"summary/{region}_snp_report.txt", "w") as summary_file:
        # header
        summary_file.write(
            f"###.................{region}.......................####\n\n\n"
        )

        # sequence with SNP
        summary_file.write(
            "###.................Sequence (with SNP in lowercase).......................####\n"
        )
        summary_file.write(product_sequence + "\n\n\n")

        # SNPs positions
        summary_file.write(
            "SNP_position_on_human_reference_genome|nucelotide_on_human_genome_reference|SNP_position_in_sample_sequence|nucelotide_in_sample\n"
        )
        for SNPs in pos_on_sequence_list:
            summary_file.write(SNPs + "\n")


# change working directory to the location where the files are
# os.chdir(folder_location)


# for directory in ("plots", "files", "summary"):
# if not os.path.exists(directory):
#  os.mkdir(directory)


# labelling_summary to get the information on the samples.
# labelling_summary_df = pd.read_csv(f"{pear_merged_loc}/sample_metadata.csv")

# get the sample names that I'm doing analysis for
# samples = labelling_summary_df["Sample name"].to_list()

"""
# folder containing all files
folder_location = "/Users/olalekan/Desktop/TL_seq_App"
pear_merged_loc = "inputs/pear_merged"
base_quality = "inputs/base_quality"
metadata = "inputs/sample_metadata.csv"

vcf_file = "SNP/EndoC.all.chromosomes.eva.vcf.gz"
# get_region_from_snps_vcf(vcf_file, 11, 2181009, 2182439)
# number_of_tc_mutations_plot("exon2", "24h")


ref_loc = "inputs/reference_gene/gene_sequence_from_refseq.txt"
# reference_coordinate(vcf_file, 11, 2181009, 2182439, ref_loc, "-")
# check_the_snp_position_in_samples(vcf_file, 11, "c_24h_exon2_fr1_S33", ref_loc, 2181009, 2182439, gene_direction="-")

chromosome = 11

gene_pos_start = 2181009
gene_pos_end = 2182439

strandedness = "-"

# analyse all samples
df = pd.read_csv(metadata)
regions = list(set(df["Region"].to_list()))
times = list(set(df["Time"].to_list()))


for region in regions[0:1]:
    file_preparation(
        pear_merged_loc,
        metadata=metadata,
        folder_location=folder_location,
        base_quality=base_quality,
    )

    summary_report(
        vcf_file,
        ref_loc,
        region,
        chromosome,
        gene_pos_start,
        gene_pos_end,
        strandedness,
        labelling_summary_df,
    )
    for time in times:
        number_of_tc_mutations_plot(region, time)
        count_number_of_different_pos_tc_mutations(region, time)



TO-DO:
    make pear merged reads an input: pear_merged_dir = "/Users/olalekan/Desktop/TL_seq_App/inputs/pear_merged"
    make base quality directory an input: base_qual_dir= "/Users/olalekan/Desktop/TL_seq_App/inputs/base_quality"

 """


def main():
    parser = argparse.ArgumentParser(description="Identify newly transcribed reads")
    parser.add_argument(
        "--folder_location",
        type=str,
        default="",
        help="Directory to save the output for the analysis",
    )
    parser.add_argument(
        "--pear_merged_loc",
        type=str,
        default="",
        help="Directory for merged output files.",
    )
    parser.add_argument(
        "--base_quality",
        type=str,
        default="",
        help="directory of the base_quality",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default="",
        help="Path to the sample metadata csv file.",
    )
    parser.add_argument(
        "--vcf_file",
        type=str,
        default="",
        help="Path to the VCF file for SNPs",
    )
    parser.add_argument(
        "--ref_loc",
        type=str,
        default="",
        help="Path to the reference sequence file for the sample",
    )
    parser.add_argument(
        "--chromosome",
        type=int,
        default=None,
        choices=[
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            "X",
            "Y",
            "x",
            "y",
        ],
        help="chromosome in which the gene of interest is located.",
    )
    parser.add_argument(
        "--gene_pos_start",
        type=int,
        default=0,
        help="Start position of the gene of interest.",
    )
    parser.add_argument(
        "--gene_pos_end",
        type=int,
        default=0,
        help="End position of the gene of interest.",
    )
    parser.add_argument(
        "--strandedness",
        type=str,
        default="-",
        choices=["+", "-"],
        help="Strandedness of the gene, either '+' or '-'.",
    )

    args = parser.parse_args()

    # Read the metadata to determine regions and times
    df = pd.read_csv(args.metadata)
    regions = list(set(df["Region"].to_list()))
    times = list(set(df["Time"].to_list()))

    # print(df)

    global pear_merged_loc
    global base_quality

    pear_merged_loc = args.pear_merged_loc
    base_quality = args.base_quality

    for region in regions:
        file_preparation(
            pear_merged_loc=args.pear_merged_loc,
            metadata=df,
            folder_location=args.folder_location,
            base_quality=args.base_quality,
        )

        summary_report(
            args.vcf_file,
            args.ref_loc,
            region,
            args.chromosome,
            args.gene_pos_start,
            args.gene_pos_end,
            args.strandedness,
            labelling_summary_df=df,  # Define or modify this according to actual use
        )

        for time in times:
            number_of_tc_mutations_plot(region, time)
            count_number_of_different_pos_tc_mutations(region, time)


if __name__ == "__main__":
    main()


# The to fix in the above is that it runs multipe samples at once - such that it does a sample at a time.
