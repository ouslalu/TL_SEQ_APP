# cli file to set up the terminal app
import argparse
import pandas as pd
from .fragment_preprocessing import file_preparation
from .fragment_preprocessing import process_pos_ref_query
from .fragment_preprocessing import get_reverse_complement
from .fragment_preprocessing import check_primers
from .fragment_preprocessing import get_mismatches
from .fragment_preprocessing import read_length_plot
from .fragment_preprocessing import analysis_for_nacent_nucleotides
from .fragment_preprocessing import count_mutation_types
from .fragment_preprocessing import prop_mutation_count
from .fragment_preprocessing import query
from .fragment_preprocessing import count_number_of_tc_mutations
from .fragment_preprocessing import count_number_of_different_pos_tc_mutations
from .fragment_preprocessing import tc_mutations_plot
from .fragment_preprocessing import number_of_tc_mutations_plot
from .fragment_preprocessing import get_region_from_snps_vcf
from .fragment_preprocessing import reference_reformating
from .fragment_preprocessing import reference_coordinate
from .fragment_preprocessing import check_the_snp_position_in_samples
from .fragment_preprocessing import summary_report


def main():
    global pear_merged_loc
    global base_quality

    parser = argparse.ArgumentParser(description="Identify newly transcribed reads")
    parser.add_argument(
        "--folder_location",
        "-f",
        type=str,
        required=True,
        help="Directory to save the output for the analysis",
    )
    parser.add_argument(
        "--pear_merged_loc",
        "-p",
        type=str,
        required=True,
        help="Directory for merged output files.",
    )
    parser.add_argument(
        "--base_quality",
        "-b",
        type=str,
        required=True,
        help="directory of the base_quality",
    )
    parser.add_argument(
        "--metadata",
        "-m",
        type=str,
        required=True,
        help="Path to the sample metadata csv file.",
    )
    parser.add_argument(
        "--vcf_file",
        "-v",
        type=str,
        required=True,
        help="Path to the VCF file for SNPs",
    )
    parser.add_argument(
        "--ref_loc",
        "-r",
        type=str,
        required=True,
        help="Path to the reference sequence file for the sample",
    )
    parser.add_argument(
        "--chromosome",
        "-chr",
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
        required=True,
        help="chromosome in which the gene of interest is located.",
    )
    parser.add_argument(
        "--gene_pos_start",
        "-gs",
        type=int,
        required=True,
        help="Start position of the gene of interest.",
    )
    parser.add_argument(
        "--gene_pos_end",
        "-ge",
        type=int,
        required=True,
        help="End position of the gene of interest.",
    )
    parser.add_argument(
        "--strandedness",
        "-s",
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
