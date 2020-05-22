import gg_primer
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import sys
import csv
# import argparse
#
# parser = argparse.ArgumentParser(
#         description = 'Program takes a csv file with a name, part type for gg cloning, and the sequence you want to design primers for in order to domesticate. The output is a tab delimited file'
# )
# parser.add_argument('-i', help='a csv file')
# parser.add_argument('-o', help='a tab delimited file')
# parser.add_argument('-t', type=int, help='primer melting temp')
# args = parser.parse_args()


inFile = sys.argv[1]
outFile = sys.argv[2]
set_tmp = sys.argv[3]

output = {}
with open(inFile,'r') as i:
        for line in i:
                item = line.split(',')
                seq_name = item[0].strip()
                seq_type = item[1].strip()
                seq = item[2].strip()

                fwd_primer_seq = gg_primer.fwd_tm(s=seq, temp=set_tmp)
                rev_primer_seq = gg_primer.rev_tm(s=seq, temp=set_tmp)

                primers = [fwd_primer_seq, rev_primer_seq]

                final_primers = gg_primer.append_overhangs(primers, type=seq_type)
                output[seq_name] = final_primers


with open (outFile, 'w') as o:
        for key, value in output.items():

                for line in value:
                        s = "\t".join(line)
                        o.writelines(key + '\t' + "%s\n" % s)
i.close()
o.close()
