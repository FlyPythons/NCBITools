#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

from ncbitools import __author__, __email__, __version__
from ncbitools.FastaReader import open_fasta


LOG = logging.getLogger(__name__)


def read_tbl(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if line.startswith("#"):
            continue

        if line:
            yield line.split(sep)


def get_seq_by_acc(seqobj, accessions):
    """
    get the seq of accession in accessions list from seqobj object

    :param seqobj:
    :param accessions:
    :return:
    """

    for record in seqobj:
        id = record.id
        acc = ".".join(id.split(".")[:-1])

        if acc in accessions:
            print(str(record))

    return 0


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
author:  %s
email: %s
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("accession", help="")
    args.add_argument("fasta")

    return args.parse_args()


def main():
    args = set_args()

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    seqobj = open_fasta(args.fasta)

    LOG.info("Parsing accession from %r" % args.accession)
    accessions = [".".join(a[0].split(".")[:-1]) for a in read_tbl(args.accession)]

    LOG.info("Getting sequences")
    get_seq_by_acc(seqobj, accessions)

    LOG.info("All done!")


if __name__ == "__main__":
    main()
