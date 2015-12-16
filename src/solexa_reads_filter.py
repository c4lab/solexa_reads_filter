#!/usr/bin/env python
#
# Solexa Reads Filter
#
# Filters:
# * s35: From the 5' end of the read, the first 25 of 35 bases must have
#   quality scores of at least 30, otherwise the read is discarded.
#   (If the read length is lower than 35, the read will be discarded.)
# * Ns: If there exists an ambiguous base (N call) then the read is
#   discarded.
# * polyN: If 85% of the read contains one type of base call (A, T, G, C),
#   then the read deemed to have low complexity and is discarded.

__version__ = '1.2'

import sys
import argparse
import logging
import gzip
import bz2
from logging import StreamHandler, Formatter
from os import makedirs
from os.path import join, exists, splitext
from itertools import izip


def parse_fastq(path):
    ext = splitext(path)[1]

    if ext == '.gz':
        fi = gzip.open(path, 'r')
    elif ext == '.bz2':
        fi = bz2.open(path, 'r')
    else:
        fi = open(path, 'r')

    n = 0
    for i in fi:
        t = n % 4
        i = i.strip()

        if t == 0:
            name = i
        elif t == 1:
            seq = i
        elif t == 2:
            spacer = i
        else:
            qua = i
            yield name, seq, spacer, qua

        n += 1

    fi.close()


def solexa_reads_filter(r1_path, r2_path, quafmt, min_len=1,
                        s35=True, ns=True, polyn=True):
    logger = logging.getLogger('filter')
    n_total_seq = 0
    n_total_base = 0
    n_mlen_drop_seq = 0
    n_mlen_drop_base = 0
    n_s35_drop_seq = 0
    n_s35_drop_base = 0
    n_ns_drop_seq = 0
    n_ns_drop_base = 0
    n_polyn_drop_seq = 0
    n_polyn_drop_base = 0
    n_report = 2000

    for r1, r2 in izip(parse_fastq(r1_path), parse_fastq(r2_path)):
        name_1, seq_1, spacer_1, qua_1 = r1
        name_2, seq_2, spacer_2, qua_2 = r2
        n_total_seq += 2
        n_total_base += len(seq_1) + len(seq_2)

        if n_total_seq % n_report == 0:
            logger.info('Processed: {0}'.format(n_total_seq))
            n_report *= 2

        if len(seq_1) < min_len or len(seq_2) < min_len:
            n_mlen_drop_seq += 2
            n_mlen_drop_base += len(seq_1) + len(seq_2)
            continue

        if s35 and (is_s35_bad(qua_1, quafmt) or is_s35_bad(qua_2, quafmt)):
            n_s35_drop_seq += 2
            n_s35_drop_base += len(seq_1) + len(seq_2)
            continue

        if ns and (is_ns(seq_1) or is_ns(seq_2)):
            n_ns_drop_seq += 2
            n_ns_drop_base += len(seq_1) + len(seq_2)
            continue

        if polyn and (is_polyn(seq_1) or is_polyn(seq_2)):
            n_polyn_drop_seq += 2
            n_polyn_drop_base += len(seq_1) + len(seq_2)
            continue

        yield r1, r2

    n_retained_seq = n_total_seq - n_mlen_drop_seq - n_s35_drop_seq
    n_retained_seq -= n_ns_drop_seq + n_polyn_drop_seq
    n_retained_seq_pct = round(float(n_retained_seq) / n_total_seq * 100, 2)
    n_retained_base = n_total_base - n_mlen_drop_base - n_s35_drop_base
    n_retained_base -= n_ns_drop_base + n_polyn_drop_base
    n_retained_base_pct = round(float(n_retained_base) / n_total_base * 100, 2)

    logger.info('Min len dropped reads: {0}'.format(n_mlen_drop_seq))
    logger.info('Min len dropped bases: {0}'.format(n_mlen_drop_base))
    logger.info('s35 dropped reads: {0}'.format(n_s35_drop_seq))
    logger.info('s35 dropeed bases: {0}'.format(n_s35_drop_base))
    logger.info('Ns dropped reads: {0}'.format(n_ns_drop_seq))
    logger.info('Ns dropeed bases: {0}'.format(n_ns_drop_base))
    logger.info('polyN dropped reads: {0}'.format(n_polyn_drop_seq))
    logger.info('polyN dropeed bases: {0}'.format(n_polyn_drop_base))
    logger.info('Total reads: {0}'.format(n_total_seq))
    logger.info('Total bases: {0}'.format(n_total_base))
    logger.info('Retained reads: {0} ({1} %)'.format(n_retained_seq, n_retained_seq_pct))
    logger.info('Retained bases: {0} ({1} %)'.format(n_retained_base, n_retained_base_pct))


def is_s35_bad(qua, quafmt):
    """s35 filter

    From the 5' end of the read, the first 25 of 35 bases must
    have quality scores of at least 30, otherwise the read is discarded.
    If the read length is lower than 35, the read is discarded.
    """
    assert quafmt in [33, 64]

    if len(qua) < 35:
        return True

    q30 = 0
    for i in range(35):
        if ord(qua[i]) - quafmt >= 30:
            q30 += 1

    if q30 < 25:
        return True
    else:
        return False


def is_ns(seq):
    """Ns filter

    If there exists an ambiguous base (N call) then the read is discarded.
    """

    if 'N' in seq:
        return True
    else:
        return False


def is_polyn(seq, threshold=0.85):
    """polyN filter

    If n% (n = threshold * 100) of the read contains one type of base call
    (A, T, G, C), then the read deemed to have low complexity and is discarded.
    """

    if float(seq.count('A')) / len(seq) >= threshold:
        return True

    if float(seq.count('T')) / len(seq) >= threshold:
        return True

    if float(seq.count('C')) / len(seq) >= threshold:
        return True

    if float(seq.count('G')) / len(seq) >= threshold:
        return True

    return False


def main():
    parser = argparse.ArgumentParser(prog='solexa_reads_filter')
    parser.add_argument('-r1', dest='r1', metavar='<file>',
                        required=True,
                        help='read 1 in FASTQ format (also support gzip and bzip2 format)')
    parser.add_argument('-r2', dest='r2', metavar='<file>',
                        required=True,
                        help='read 2 in FASTQ format (also support gzip and bzip2 format)')
    parser.add_argument('-Q', dest='quafmt', type=int,
                        metavar='(33|64)', choices=[33, 64],
                        required=True,
                        help='format of quality scores (33 or 64)')
    parser.add_argument('-r', dest='min_len', metavar='<int>',
                        type=int, default=1,
                        help='minimum read length to be retained after trimming (default: 1)')
    parser.add_argument('-o', dest='outdir', metavar='<dir>',
                        required=True,
                        help='output directory')
    parser.add_argument('-m', dest='merge', action='store_true',
                        help='merge read 1 and read 2 after filtering')
    parser.add_argument('-z', dest='fr_s35', action='store_false',
                        help='turn off s35 filtering')
    parser.add_argument('-x', dest='fr_ns', action='store_false',
                        help='turn off Ns filtering')
    parser.add_argument('-v', dest='fr_polyn', action='store_false',
                        help='turen off polyN filtering')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()

    if not exists(args.outdir):
        makedirs(args.outdir)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=join(args.outdir, 'log'),
        filemode='w',
    )
    stream_handler = StreamHandler(sys.stdout)
    formatter = Formatter('%(asctime)s [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(stream_handler)

    if args.merge:
        with open(join(args.outdir, 'pe.merged.fq'), 'w') as fo:
            for r1, r2 in solexa_reads_filter(args.r1, args.r2, args.quafmt,
                                              min_len=args.min_len,
                                              s35=args.fr_s35,
                                              ns=args.fr_ns,
                                              polyn=args.fr_polyn):
                fo.write('\n'.join(r1))
                fo.write('\n')
                fo.write('\n'.join(r2))
                fo.write('\n')
                fo.flush()
    else:
        with open(join(args.outdir, 'pe.r1.fq'), 'w') as fo_r1, \
                open(join(args.outdir, 'pe.r2.fq'), 'w') as fo_r2:
            for r1, r2 in solexa_reads_filter(args.r1, args.r2, args.quafmt,
                                              min_len=args.min_len,
                                              s35=args.fr_s35,
                                              ns=args.fr_ns,
                                              polyn=args.fr_polyn):
                fo_r1.write('\n'.join(r1))
                fo_r1.write('\n')
                fo_r1.flush()

                fo_r2.write('\n'.join(r2))
                fo_r2.write('\n')
                fo_r2.flush()

if __name__ == '__main__':
    main()
