# file stream script to integrate read count ID with full lineage
import csv
from pathlib import Path
import argparse
from collections import Counter


def f(seq): # Order preserving
  seen = set()
  return [x for x in seq if x not in seen and not seen.add(x)]


def read_metascope_out(args):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    output_name = args.outDir / args.file.name.replace('.sorted', '')

    write_header(output_name)

    # open csv files
    with open(args.file, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')

        # get all read_counts to setup proportion
        # read_count_list = []
        id_list = []
        for line in reader:
            id_list.append(line.get('TaxonomyID'))
            # read_count_list.append(int(line.get('read_count')))

        # total_reads = sum(read_count_list)

        unique_id_list = f(id_list)
        counts = Counter(id_list)

    sliding_window = []

    i = 0

    # open csv files
    with open(args.file, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')

        # use ID occurrence to prime the sliding window
        for id in unique_id_list:
            print(counts[id])

            for line in reader:
                # building a sliding window
                sliding_window.append(line)

                if len(sliding_window) >= counts[id]:
                    write_file(output_name, sliding_window, args)
                    break


def write_header(output_name):
    # open output file in append mode
    with open(output_name, 'w') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['TaxonomyID', 'Genome', 'read_count', 'Proportion', 'EMreads', 'EMProportion']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter=',')

        writer.writeheader()


def write_file(output_name, sliding_window, line):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['TaxonomyID', 'Genome', 'read_count', 'Proportion', 'EMreads', 'EMProportion']

        line = sliding_window[0]

        if line['TaxonomyID'] == 'NA':
            print(line['TaxonomyID'])
            sliding_window.clear()
            return

        if len(sliding_window) > 1:
            new_read_count_list = []

            for line in sliding_window:
                new_read_count_list.append(int(line.get('read_count')))

            new_read_count = sum(new_read_count_list)

            print(new_read_count_list, new_read_count)
            new_read_count_list.clear()

            line['read_count'] = new_read_count

        # remove the duplicates
        sliding_window.clear()

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter=',')
        writer.writerow(line)
        
        return


def main():
    # setup command line arguments to be passed to the script
    parser = argparse.ArgumentParser(
        prog='merge NCBI ids',
        description='A program to convert the output of merged MetaScope IDs into a non-redundant file'
    )
    parser.add_argument('-f', '--file', nargs='?', required=True, help='MetaScope ID\'s output csv file')
    parser.add_argument('-o', '--outDir', nargs='?', required=True, help='directory for the output file(s)')
    args = parser.parse_args()

    # directory to store the results
    args.outDir = Path(args.outDir)

    # directory containing the raw CHH files
    args.file = Path(args.file)

    read_metascope_out(args)


if __name__ == "__main__":
    main()