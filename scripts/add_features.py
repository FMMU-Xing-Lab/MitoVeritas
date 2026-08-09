#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Add annotation features to the feature table (parameterized version of 添加10个特征.py).

Added columns (9):
    repeat-region, Population-freq, NAV, type, if_Trans, if_version,
    region, VAF_mitomap, dbSNP

Usage:
    python3 add_features.py --input feature.txt \
                            --output feature_add10.txt \
                            --data-dir ../data
"""

import argparse
import os
import sys
import time
from collections import Counter


def read_lines(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.readlines()


def parse_args():
    parser = argparse.ArgumentParser(description="Add annotation features to the feature table.")
    parser.add_argument("--input", required=True, help="feature table (after delete_features.py)")
    parser.add_argument("--output", required=True, help="output feature table (e.g. feature_add10.txt)")
    parser.add_argument(
        "--data-dir",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"),
        help="directory containing dbSNP.txt, mitomap.txt, mitomap-snp.txt, region.txt, mtDNA_region_ge5.txt",
    )
    return parser.parse_args()


def main():
    start_time = time.time()
    args = parse_args()
    data_dir = args.data_dir

    def data_file(name):
        path = os.path.join(data_dir, name)
        if not os.path.exists(path):
            sys.exit("ERROR: required data file not found: %s" % path)
        return path

    # ==================== read input feature table ====================
    site_file = read_lines(args.input)
    print("Read %d lines from %s" % (len(site_file), args.input))

    # remove repeated header lines (rows whose first column is 'id')
    site_file = [
        line for idx, line in enumerate(site_file)
        if idx == 0 or line.split("\t")[0] != "id"
    ]
    if len(site_file) < 2:
        sys.exit("ERROR: no variant rows found in %s" % args.input)
    print("Feature table: %d variants + header" % (len(site_file) - 1))

    # ==================== add repeat-region (position in a >=5 bp repeat region) ====
    region_file = read_lines(data_file("mtDNA_region_ge5.txt"))
    region_list = []
    for i in region_file:
        cols = i.strip("\n").split("\t")
        if len(cols) < 4:
            continue
        region_list.append([cols[1], cols[2], cols[-1]])

    for i in range(1, len(site_file)):
        cols = site_file[i].strip("\n").split("\t")
        try:
            site = int(cols[0].split("~")[2])
        except (IndexError, ValueError):
            print("WARNING: cannot parse position from row %d, setting repeat-region=0" % i,
                  file=sys.stderr)
            site = -1
        is_in_region = 0
        for j in region_list:
            if int(j[0]) <= site <= int(j[1]):
                is_in_region = 1
                break
        site_file[i] = site_file[i].strip("\n") + "\t" + str(is_in_region) + "\n"

    print("Added repeat-region feature (%.1f s)" % (time.time() - start_time))

    # ==================== add Population-freq (position frequency across samples) ====
    freq_list = []
    sample_list = []
    for i in range(1, len(site_file)):
        parts = site_file[i].split("\t")[0].split("~")
        if len(parts) >= 3:
            freq_list.append(parts[2])
            sample_list.append(parts[0])
    sample_set = set(sample_list)
    freq_counter = Counter(freq_list)
    n_samples = len(sample_set)

    for i in range(1, len(site_file)):
        parts = site_file[i].split("\t")[0].split("~")
        site = parts[2] if len(parts) >= 3 else "0"
        count_freq = freq_counter[site]
        site_file[i] = site_file[i].strip("\n") + "\t" + str(count_freq / n_samples) + "\n"

    print("Added Population-freq feature (%.1f s)" % (time.time() - start_time))

    # ==================== add NAV (known polymorphism in mitomap-snp) ====
    f_ref = read_lines(data_file("mitomap-snp.txt"))
    nav_index = {}
    for j in range(1, len(f_ref)):
        cols = f_ref[j].split("\t")
        if len(cols) < 4:
            continue
        nav_index.setdefault(cols[0], []).append((cols[2], cols[3]))

    for i in range(1, len(site_file)):
        name = site_file[i].strip("\n").split("\t")[0]
        parts = name.split("~")
        if len(parts) < 5:
            print("WARNING: cannot parse variant id '%s', NAV=0" % name, file=sys.stderr)
            site_file[i] = site_file[i].strip("\n") + "\t0\n"
            continue
        site, bg1, bg2 = parts[2], parts[3], parts[4]
        nav_value = 0
        for refer, alter in nav_index.get(site, []):
            if (bg1 == refer and bg2 == alter) or (bg1 == alter and bg2 == refer):
                nav_value = 1
                break
        site_file[i] = site_file[i].strip("\n") + "\t" + str(nav_value) + "\n"

    print("Added NAV feature (%.1f s)" % (time.time() - start_time))

    # ==================== add type / if_Trans / if_version / region / VAF_mitomap / dbSNP ====
    file_region = read_lines(data_file("region.txt"))
    file_mitomap = read_lines(data_file("mitomap.txt"))
    ref_dbSNP = set(line.strip() for line in read_lines(data_file("dbSNP.txt")))

    ref_region = {}
    for i in file_region:
        cols = i.strip("\n").split("\t")
        if len(cols) >= 2:
            ref_region[cols[0]] = cols[1]

    ref_mitomap = {}
    for i in file_mitomap:
        cols = i.strip("\n").split("\t")
        if len(cols) >= 2:
            ref_mitomap[cols[0]] = cols[1]

    trans = ["A>G", "G>A", "C>T", "T>C"]

    for i in range(1, len(site_file)):
        parts = site_file[i].strip("\n").split("\t")[0].split("~")
        if len(parts) < 5:
            print("WARNING: cannot parse variant id on row %d, skipping annotations" % i,
                  file=sys.stderr)
            site_file[i] = site_file[i].strip("\n") + "\tNA\t0\t1\tNA\t0\t0\n"
            continue
        site = parts[2]
        base1 = parts[3]
        base2 = parts[4]
        base_type = "%s>%s" % (base1, base2)
        sitetype_index = "%s%s>%s" % (site, base1, base2)

        if base_type in trans:
            if_trans = "1"
            if_version = "0"
        else:
            if_trans = "0"
            if_version = "1"

        region = ref_region.get(site, "NA")
        mitomap = ref_mitomap.get(sitetype_index, "0")
        dbSNP = "1" if sitetype_index in ref_dbSNP else "0"

        site_file[i] = site_file[i].strip("\n") + "\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            base_type, if_trans, if_version, region, mitomap, dbSNP,
        )

    # ==================== write output ====================
    title = (site_file[0].strip("\n") +
             "\trepeat-region\tPopulation-freq\tNAV\ttype\tif_Trans\tif_version\tregion\tVAF_mitomap\tdbSNP\n")

    with open(args.output, "w", encoding="utf-8", newline="\n") as new_file:
        new_file.write(title)
        for i in range(1, len(site_file)):
            new_file.write(site_file[i])

    print("Added all annotation features and wrote %d variants to %s (%.1f s total)" %
          (len(site_file) - 1, args.output, time.time() - start_time))


if __name__ == "__main__":
    main()
