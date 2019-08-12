#!/usr/bin/env python

import sys
import argparse
import logging
from Bio import SeqIO
from collections import defaultdict
#by jmxia

def n50(lengths):

    sum_length = sum(lengths)
    accu_length = 0

    for i in sorted(lengths, reverse=True):
        accu_length += i

        if accu_length >= sum_length*0.5:
            return i

def stat_genome(genome,contig_out):

    lengths = []
    base_dict = {}

    for record in SeqIO.parse(genome, "fasta"):
        seq = record.seq.upper()
        a_number = seq.count('A')
        t_number = seq.count('T')
        g_number = seq.count('G')
        c_number = seq.count('C')
        n_number = seq.count('N')
    
        lengths.append(len(seq))
        base_dict[record.id] = [a_number,t_number,g_number,c_number,n_number]

    with open(contig_out, "w") as fh:
        fh.write("""#stat result\n
Total_bases\t{0:,}
Contig_number\t{1:,}
Longest_Contig\t{2:,}
Shortest_Contig\t{3:,}
Contig_N50\t{4:,}
""".format(sum(lengths), len(lengths), max(lengths), min(lengths), n50(lengths)))
    
    return base_dict

def stat_gc(gc_dict, gc_out):

    sum_a = sum_t = sum_g = sum_c = sum_n = 0

    for a, t, g, c, n in gc_dict.values():
        sum_a += a
        sum_t += t
        sum_g += g
        sum_c += c
        sum_n += n

    sum_all = sum_a + sum_t + sum_g + sum_c + sum_n

    with open(gc_out, "w") as fh:
        fh.write("""\
#Base\tNumber\t%of_total
A\t{0:,}\t{1:.2f}
T\t{2:,}\t{3:.2f}
G\t{4:,}\t{5:.2f}
C\t{6:,}\t{7:.2f}
N\t{8:,}\t{9:.2f}
GC\t{10:,}\t{11:.2f}
Total\t{12:,}\t{13:.2f}
""".format(
    sum_a, sum_a*100.0/sum_all,
    sum_t, sum_t*100.0/sum_all,
    sum_g, sum_g*100.0/sum_all,
    sum_c, sum_c*100.0/sum_all,
    sum_n, sum_n*100.0/sum_all,
    sum_g+sum_c, (sum_g+sum_c)*100.0/sum_all,
    sum_all, 100.0
    ))

def set_args(parser):
  parser.add_argument('-i','--input',metavar='FILE',type=str,required=True,help='Input genome file')
  parser.add_argument('-c','--contig',metavar='STR',type=str,default='contig_stat.out',help='Output contig stat')
  parser.add_argument('-g','--gc',metavar='STR',type=str,default='gc_stat.out',help='Output gc stat')
  return parser

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
Description:
	gc_contig_stat.py  stat gc and length of contig in assembled file
''')

    args = set_args(parser).parse_args()
    dict_base = stat_genome(args.input,args.contig)
    stat_gc(dict_base,args.gc)

if __name__ == "__main__":
  main()

