import glob
import os
import argparse
import subprocess
import pandas as pd
from io import StringIO
import re


def parser_args():
    parser = argparse.ArgumentParser(
        description='Compares the results between two sets of paired VCFs')
    parser.add_argument('--truth_folder', '-T',
                        help='Folder containing the truth VCF files used as a benchmark',
                        required=True)
    parser.add_argument('--truth_suffix', '-t',
                        nargs='+',
                        help='Suffix to remove from truth VCF files in the benchmark folder',
                        required=True)
    parser.add_argument('--query_folder', '-Q',
                        help='Folder containing the query VCF files requiring verification',
                        required=True)
    parser.add_argument('--query_suffix', '-q',
                        nargs='+',
                        help='Suffix to remove from query VCF files',
                        required=True)
    parser.add_argument('--output_name', '-o',
                        type=str,
                        help='File path to output directory',
                        required=True)
    parser.add_argument('--reference_genome', '-r',
                        type=str,
                        help='FASTA genome reference (Must have index in same directory)',
                        required=True)
    return parser.parse_args()


# Import list of VCF files
class compareVCFs(object):
    # TODO Split functions into classes
    def __init__(self, truth_folder,
                 query_folder,
                 truth_pattern,
                 query_pattern,
                 output_name,
                 reference_genome):
        self.truth_folder = truth_folder
        self.query_folder = query_folder
        self.truth_pattern = truth_pattern[0]
        self.query_pattern = query_pattern[0]
        self.output_name = output_name
        self.reference_genome = reference_genome
        self.truth_VCFs = None
        self.truth_samples = None
        self.query_VCFs = None
        self.query_samples = None
        self.sample_overlap = None
        self.truth_mismatch = None
        self.query_mismatch = None
        self.all_samples = None
        self.truth_dict = None
        self.query_dict = None

    def get_filenames(self):
        # Get list of all VCFs in truth set
        search_pattern = self.truth_folder + "*" + self.truth_pattern
        self.truth_VCFs = glob.glob(search_pattern)
        regex_pattern = '\.' + self.truth_pattern.strip('.') + '$'
        self.truth_samples = [re.sub(regex_pattern, '', os.path.basename(i)) for i in self.truth_VCFs]

        # Get list of all VCFs in query set
        search_pattern = self.query_folder + "*" + self.query_pattern
        self.query_VCFs = glob.glob(search_pattern)
        regex_pattern = '\.' + self.query_pattern.strip('.') + '$'
        self.query_samples = [re.sub(regex_pattern, '', os.path.basename(i)) for i in self.query_VCFs]
        # Create a dictionary to reiterate through
        self.truth_dict = dict(zip(self.truth_samples, self.truth_VCFs))
        self.query_dict = dict(zip(self.query_samples, self.query_VCFs))

    def check_file_matches(self):
        # TODO Split out print function into new function
        # Provide metrics on matching files
        self.sample_overlap = set(self.truth_samples).intersection(self.query_samples)
        self.truth_mismatch = set(self.truth_samples).difference(self.query_samples)
        self.query_mismatch = set(self.query_samples).difference(self.truth_samples)
        self.all_samples = set(self.truth_samples).union(self.query_samples)

        # Print import stats
        print("Paired Samples: " + str(len(self.sample_overlap)))
        # print(self.sample_overlap)
        print("Total Samples: " + str(len(self.all_samples)))
        # print(self.all_samples)
        print("Unpaired Samples: In Query Set, absent Truth Set: " + str(len(self.truth_mismatch)))
        if len(self.truth_mismatch) == 0:
            print("\tNo Mismatched samples from Truth Set")
        else:
            print("\tSamples present in the Truth Set but not the Query Set:")
            print(self.truth_mismatch)
        print("Unpaired Samples: In Truth Set, absent Query Set: " + str(len(self.query_mismatch)))
        if len(self.query_mismatch) == 0:
            print("\tNo Mismatched samples from Query Set")
        else:
            print("\tSamples present in the Query Set but not the Truth Set:")
            print(self.query_mismatch)

        # Run vcftools, bcftools, or som.py for each pair

    def run_comparison(self):
        # Uses som.py - https://github.com/Illumina/hap.py/blob/master/doc/sompy.md
        # See link for describtion of output.
        # Pull down docker image fro bcftools
        # cmd = "docker pull quay.io/biocontainers/bcftools@sha256:e262ac00f3e4e67bd1627ababc296e990393f45857f54eaa5f2e223dfbb5cc02"
        # returned_value = os.system(cmd)  # returns the exit code in unix

        # Initialise dataframe

        count = 0

        for sample_name in self.sample_overlap:
            truth_file = self.truth_dict.get(sample_name)
            query_file = self.query_dict.get(sample_name)


            # Create command to pass to som.py
            cmd = "docker run -v " \
                  + self.output_name + ":/home/results" \
                  + " -v " + os.path.dirname(self.reference_genome) + ":/home/reference/" \
                  + " -v " + query_file + ":/home/query/" + os.path.basename(query_file) \
                  + " -v " + truth_file + ":/home/truth/" + os.path.basename(truth_file) \
                  + " --rm 7593fa57027b som.py " \
                  + "/home/truth/" + os.path.basename(truth_file) + " " \
                  + "/home/query/" + os.path.basename(query_file) \
                  + " -r /home/reference/" + os.path.basename(self.reference_genome) + " -o temp.txt"

            output = subprocess.Popen(cmd,
                                      shell=True,
                                      stdout=subprocess.PIPE)

            output = StringIO(output.communicate()[0].decode('utf-8'))

            # Remove multiple spaces in string
            b = re.sub(" +", " ", output.getvalue())
            # Remove trailing newline
            b = re.sub("^\n", "", b)
            # Read sanitised string into pandas dataframe using StringIO object
            b = StringIO(b)
            df = pd.read_csv(b, sep="\s")

            # Add sample name to each row of the dataframe
            df.insert(0, "sample_name", sample_name, True)

            # One first pass
            # TODO check that df is correctly formatted, if not print to log
            if count == 0:
                results_df = df
            else:
                results_df = results_df.append(df, ignore_index=True)
            count = 1
        # TODO Split out into separate function
        results_df.to_csv(self.output_name + ".csv",
                          index=False)


        return results_df


    def sort_vcfs(self):
        # TODO Sort all VCFs
        '''
        # Sort and compress files
        cmd = "vcf-sort " + truth_file + " | bgzip -c > " + truth_file + ".gz"
        print(cmd)
        returned_value = os.system(cmd)
        cmd = "vcf-sort " + query_file + " | bgzip -c > " + query_file + ".gz"
        returned_value = os.system(cmd)
        '''

        return 1

    def compress_vcfs(self):
        # TODO use bgzip to compress files (required for some tools)
        '''
        # Sort and compress files
        cmd = "vcf-sort " + truth_file + " | bgzip -c > " + truth_file + ".gz"
        print(cmd)
        returned_value = os.system(cmd)
        cmd = "vcf-sort " + query_file + " | bgzip -c > " + query_file + ".gz"
        returned_value = os.system(cmd)
        '''
        return 1

    def index_vcfs(self):
        # TODO produce index files for all input (required for some tools)
        '''
        # Index files
        # cmd = "tabix -p vcf ONCW_02_EK6037_SWIFT57_Pan2684_S2_R1_001.refined.primerclipped.vardict_out.vcf.gz"
        cmd = "tabix -p vcf " + truth_file + ".gz"
        returned_value = os.system(cmd)
        cmd = "tabix -p vcf " + query_file + ".gz"
        returned_value = os.system(cmd)
        '''
        return 1

    def run_bcftools_isec(self):
        # TODO produce union and Venn diagrams between 2 vcfs
        '''
        # Use bcftools to compare pairs of samples
        # cmd = "bcftools isec -n =2 -p isec " + truth_file + ".gz " + query_file + ".gz"
        # cmd = "bedtools intersect -u -a " + truth_file + ".gz -b " + query_file + ".gz"
        # returned_value = os.system(cmd)
        # cmd = "bedtools jaccard -a " + query_file + ".gz -b " + truth_file + ".gz >> my_results.txt"
        returned_value = os.system(cmd)
        '''
        return 1

    def run_snpEff(self):
        # TODO use snpEff to check for concordance between
        '''
        cmd = "java -jar /home/graeme/snpEff/SnpSift.jar concordance -v " + truth_file + " " + query_file + " >> snp_concordance.txt"
        # print(cmd)
        returned_value = os.system(cmd)
        # TODO PRODUCE report
        # docker pull quay.io/biocontainers/hap.py@sha256:c4d4631dd1050209702edbd78deb2c133a1a97d6c783139966e41182f960d3f2
        '''
        return 1

# Summarise data
# Visualise data
# Heatmap
# venn diagram

def main():
    # Parse commandline args

    args = parser_args()

    # Initialise compareVCFs object
    imported_VCFs = compareVCFs(args.truth_folder,
                                args.query_folder,
                                args.truth_suffix,
                                args.query_suffix,
                                args.output_name,
                                args.reference_genome, )
    imported_VCFs.get_filenames()
    imported_VCFs.check_file_matches()
    imported_VCFs.run_comparison()


if __name__ == '__main__':  # if this .py script is executing as the main function, the run main()
    main()
