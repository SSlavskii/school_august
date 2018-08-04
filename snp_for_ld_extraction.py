import subprocess
import pandas as pd
import argparse
import os

plink_path = '/home/ubuntu/gwas/old_gwas/tools/plink/plink'
ref_path = '/home/ubuntu/gwas/old_gwas/sega/ld_real/by_chr/'
table_path = '/home/ubuntu/gwas/old_gwas/sega/school_august/core_snps_filtered.csv'

RS_ID = 'ref_rs_id'
CHR = 'ref_chr'
R_threshold = 0.1
MAF_threshold = 0.05

def get_snp_list(chr_num, table_path, bim_path):
    df = pd.read_table(table_path)
    df = df[df[CHR] == chr_num]
    bim = pd.read_table(bim_path, header=None)  # rs_id column has index 1
    snp_list = list(df[df[RS_ID].isin(bim[1])][RS_ID])
    return snp_list


def get_table_for_chr(chr_num, keep_path=None, output_name=""):
    bfile_path = f"{ref_path}ALL.filt.chr{chr_num}"
    snp_list = get_snp_list(chr_num, table_path, bfile_path+'.bim')
    snps = ' '.join(snp_list)
    chr_num = str(chr_num)
    output_name = output_name + chr_num

    if keep_path is None:
        query = f"{plink_path}  --bfile {bfile_path} --r with-freqs -out {output_name} --ld-snps {snps}" \
                f" --ld-window-kb 250 --ld-window 10000000"
    else:
        query = f"{plink_path}  --bfile {bfile_path} --r with-freqs -out {output_name} --ld-snps {snps} " \
                f"--ld-window-kb 250  --ld-window 10000000 --keep {keep_path}"

    subprocess.call(query, shell=True)
    df = pd.read_csv(f"{output_name}.ld", sep='\s+')
    df.dropna(inplace=True)
    df.drop(df[abs(df["R"]) < R_threshold].index, inplace=True)
    df.drop(df[(df["MAF_A"] < MAF_threshold) | (df["MAF_B"] < MAF_threshold)].index, inplace=True)


    df.to_csv(f"{output_name}.ld", index=False, sep='\t')
    os.remove(f"{output_name}.log")
    os.remove(f"{output_name}.nosex")


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-chr', type=int, default=0, dest="chr_num")
    parser.add_argument('-keep', type=str, default=None, dest="keep_path")
    parser.add_argument('-out', type=str, default="", dest="output_name")
    return parser


if __name__ == "__main__":
    args, _ = create_parser().parse_known_args()
    if args.chr_num > 0:
        print(f"Calculation for chromosome {args.chr_num}")
        get_table_for_chr(args.chr_num, keep_path=args.keep_path, output_name=args.output_name)
    else:
        print('Calculation for all chromosomes')
        for i in range(1,23):
            get_table_for_chr(i, keep_path=args.keep_path, output_name=args.output_name)
        frames = [pd.read_csv(f"{args.output_name + str(i)}.ld", sep='\t') for i in range(1, 23)]
        result = pd.concat(frames)
        result.to_csv(f"{args.output_name}_all_chr.ld", index=False, sep='\t')
        for i in range(1,23):
            os.remove(f"{args.output_name + str(i)}.ld")