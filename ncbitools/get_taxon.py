#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

from ncbitools import __author__, __email__, __version__

__all__ = []

LOG = logging.getLogger(__name__)


class Tax(object):

    __slots__ = ["id", "name", "parent", "rank"]

    def __init__(self):
        self.id = None
        self.name = None
        self.parent = None
        self.rank = {}


def read_tbl(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if line.startswith("#"):
            continue

        if line:
            yield line.split(sep)


def read_names_dmp(file, name_class_filter="scientific name"):
    """
    read names.dmp and return a dict of tax_id and it's scientific name

    names.dmp contains:
    tax_id					-- the id of node associated with this name
    name_txt				-- name itself
    unique name				-- the unique variant of this name if name not unique
    name class				-- (synonym, common name, ...)

    :param file:
    :param name_class_filter:
    :return:
    """

    r = {}
    LOG.info("Parsing tax names from %r" % file)

    for line in read_tbl(file):
        tax_id, _, name_txt, _, unique_name, _, name_class, _ = line

        if name_class != name_class_filter:
            continue

        r[tax_id] = name_txt

    return r


def read_nodes_dmp(file):
    """
    tax_id					-- node id in GenBank taxonomy database
    parent tax_id				-- parent node id in GenBank taxonomy database
    rank					-- rank of this node (superkingdom, kingdom, ...)
    embl code				-- locus-name prefix; not unique
    division id				-- see division.dmp file
    inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
    genetic code id				-- see gencode.dmp file
    inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
    mitochondrial genetic code id		-- see gencode.dmp file
    inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
    GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
    hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
    comments				-- free-text comments and citations

    :param file:
    :return:
    """

    r = {}
    LOG.info("Parsing tax nodes from %r" % file)

    for line in read_tbl(file):
        tax_id = line[0]
        parent = line[2]
        rank = line[4]

        r[tax_id] = [parent, rank]

    return r


def merge_dmp(names_dmp, nodes_dmp):
    """

    :param names_dmp:
    :param nodes_dmp:
    :return: a dict contains {tax_id: [name_txt, parent, rank]}
    """

    r = {}

    names = read_names_dmp(names_dmp)
    nodes = read_nodes_dmp(nodes_dmp)

    for n in nodes:
        r[n] = [names[n]] + nodes[n]

    return r


def get_rank(database, tax_id, tax):

    if tax_id not in database:
        return 0

    name, parent, rank = database[tax_id]

    if parent == tax_id:
        return tax

    tax.rank[rank] = name

    get_rank(database, parent, tax)


def get_taxons(names_dmp, nodes_dmp, tax_ids, ranks):

    database = merge_dmp(names_dmp, nodes_dmp)
    print("\t".join(["tax_id"] + ranks))

    for i in read_tbl(tax_ids):
        LOG.info("Parsing tax of %s" % i[0])
        tax = Tax()
        tax.id = i[0]
        get_rank(database, i[0], tax)

        r = [i[0]]

        for j in ranks:

            if j in tax.rank:
                r.append(tax.rank[j])
            else:
                r.append("")

        print("\t".join(r))


def set_args():

    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
author:  %s
email: %s
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("file", help="file contains tax id per line")
    args.add_argument("--names", default="names.dmp", help="path of names.dmp download from NCBI")
    args.add_argument("--nodes", default="nodes.dmp", help="path of nodes.dmp download from NCBI")
    args.add_argument("--ranks", default="superkingdom,kingdom,phylum,class,order,family,genus,species",
                      help="ranks to show in the result, "
                           "default: 'superkingdom,kingdom,phylum,class,order,family,genus,species' ")

    return args.parse_args()


def main():
    args = set_args()

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    get_taxons(args.names, args.nodes, args.file, args.ranks.split(","))


if __name__ == "__main__":
    main()

