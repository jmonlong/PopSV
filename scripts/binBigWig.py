import pyBigWig
import argparse

parser = argparse.ArgumentParser(description='Compute bin counts from a BigWig file')
parser.add_argument('-bw', dest='bw', help='the BigWig file', required=True)
parser.add_argument('-bin', dest='bin', help='the BED file with the bin definition', required=True)
parser.add_argument('-o', dest='out', help='the output file', required=True)
parser.add_argument('-readL', dest='rl', help='the read length', default=100)
parser.add_argument('-approx', dest='approx', default=False,
                    help='should the pre-computed summaries be used (faster but not exact)')

args = parser.parse_args()

bw = pyBigWig.open(args.bw)
out = open(args.out, 'w')
out.write('chr\tstart\tend\tbc\n')

bins = open(args.bin, 'r')
for bin in bins:
    bin = bin.rstrip()
    bint = bin.split('\t')
    bc = bw.stats(bint[0], int(bint[1]), int(bint[2]), exact=not args.approx)
    if(bc[0] is None):
        bc = 0
    else:
        bc = float(bc[0]) * (float(bint[2]) - float(bint[1])) / args.rl
    out.write(bin + '\t' + str(bc) + '\n')

out.close()
