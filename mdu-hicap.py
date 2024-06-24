import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
import os
import glob
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG)

def run_hicap(id, path_contigs, output_dir):
    # Ensure that each ID has its own subdirectory
    id_output_dir = os.path.join(output_dir, id)
    os.makedirs(id_output_dir, exist_ok=True)
    
    temp_fa = f"{id}.fa"
    os.system(f"cp {path_contigs} {temp_fa}")
    os.system(f"hicap --query_fp {temp_fa} --output_dir {id_output_dir}")
    
    tsv_path = os.path.join(id_output_dir, f"{id}.tsv")
    if not os.path.isfile(tsv_path):
        with open(tsv_path, 'w') as f:
            f.write("#isolate\tpredicted_serotype\tattributes\tgenes_identified\tlocus_location\tregion_I_genes\tregion_II_genes\tregion_III_genes\tIS1016_hits\n")
            f.write(f"{id}\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\n")
    
    os.remove(temp_fa)
    logging.debug(f"Finished processing {id}, output in {id_output_dir}")

def process_entry(entry, output_dir):
    try:
        run_hicap(entry[0], entry[1], output_dir)
    except Exception as e:
        logging.error(f"Error processing {entry[0]}: {e}")

def process_file(file_path, output_dir, parallel):
    with open(file_path) as f:
        lines = [line.strip().split('\t') for line in f]
    
    if parallel:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_entry, line, output_dir) for line in lines]
            for future in as_completed(futures):
                try:
                    future.result(timeout=300)  # 5-minute timeout
                except TimeoutError:
                    logging.error("Processing timed out")
                except Exception as e:
                    logging.error(f"Error occurred during parallel execution: {e}")
    else:
        for id, path_contigs in lines:
            try:
                run_hicap(id, path_contigs, output_dir)
            except Exception as e:
                logging.error(f"Error processing {id}: {e}")

def concat_hicap(output_dir, to_csv=False):
    if not os.path.isdir(output_dir):
        print(f"Directory {output_dir} does not exist.")
        return
    
    tsv_files = glob.glob(f"{output_dir}/*/*.tsv")
    if not tsv_files:
        print(f"No .tsv files found in {output_dir}.")
        return

    all_dfs = []
    for i, tsv_file in enumerate(tsv_files):
        df = pd.read_csv(tsv_file, delimiter="\t")
        if i == 0:
            all_dfs.append(df)
        else:
            all_dfs.append(df)

    if all_dfs:
        # Concatenate all dataframes while keeping the header from the first file only
        summary_df = pd.concat(all_dfs, ignore_index=True)
        summary_path = os.path.join(output_dir, "hicap_summary.tab")
        summary_df.to_csv(summary_path, sep='\t', index=False, header=True)
        print(f"Created concatenated file: {summary_path}")

        if to_csv:
            csv_path = summary_path.replace(".tab", ".csv")
            summary_df.to_csv(csv_path, index=False)
            print(f"Converted {summary_path} to {csv_path}")
    else:
        print("No data to concatenate.")

def print_hicap_version():
    try:
        version_output = subprocess.check_output(['hicap', '--version']).decode('utf-8').strip()
        print(f"Hicap version: {version_output}")
    except Exception as e:
        print(f"Failed to get hicap version: {e}")
        
def main():
    parser = argparse.ArgumentParser(
        description="Runs the mdu-hicap function with the sample ID and path to contigs.",
        usage=(
            "## Recommended:\n"
            "python mdu-hicap.py --sample_id <ID> --path_contigs <path/to/contigs> [--parallel] [--concat] [--output_dir <dir>] [--csv]\n"
            "python mdu-hicap.py --contigs_file <path/to/contigs_file> [--parallel] [--concat] [--output_dir <dir>] [--csv]\n\n"
            "## Options:\n"
            "mdu-hicap.py [--sample_id <ID> --path_contigs <path/to/contigs>] [--contigs_file <file>]\n"
            "                      [--parallel] [--concat] [--output_dir <dir>] [--csv]\n"
        )
    )

    parser.add_argument('-i', '--sample_id', help="Sample ID for the hicap analysis")
    parser.add_argument('-c', '--path_contigs', help="Path to contigs file")
    parser.add_argument('-f', '--contigs_file', help="Path to tab-delimited file containing sample IDs and contigs paths", required=False)

    parser.add_argument('-p', '--parallel', action='store_true', help="Run parallel jobs")
    parser.add_argument('-t', '--concat', action='store_true', help="Concatenate .tsv files and convert to .csv")
    parser.add_argument('-o', '--output_dir', default='hicap_output', help="Directory for hicap outputs and concatenation")
    parser.add_argument('--csv', action='store_true', help="Convert concatenated summary to CSV format. NOTE: Use it with --concat option.")

    args = parser.parse_args()

    if args.sample_id and args.path_contigs:
        run_hicap(args.sample_id, args.path_contigs, args.output_dir)
    elif args.contigs_file:
        process_file(args.contigs_file, args.output_dir, args.parallel)
    else:
        print("Error: Provide either sample_id and path_contigs, or contigs_file.")
        parser.print_help()
        exit(1)
    
    if args.concat:
        concat_hicap(args.output_dir, to_csv=args.csv)
    
    print_hicap_version()

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        parser = argparse.ArgumentParser(
            description="Runs the mdu-hicap function with the sample ID and path to contigs.",
            usage=(
                "## Recommended:\n"
                "python mdu-hicap.py --contigs_file <path/to/contigs_file> [--parallel] [--concat] [--output_dir <dir>] [--csv]\n"
                "python mdu-hicap.py --sample_id <ID> --path_contigs <path/to/contigs>\n\n"
                "## Options:\n"
                "mdu-hicap.py [--sample_id <ID> --path_contigs <path/to/contigs>] [--contigs_file <file>]\n"
                "                      [--parallel] [--concat] [--output_dir <dir>] [--csv]\n"
            )
        )
        parser.print_help(sys.stderr)
        sys.exit(1)
    main()
