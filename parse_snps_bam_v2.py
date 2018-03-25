import argparse
import pysam
from numpy import mean


def parameters(snp_data):
    to_print = []
    allele2_freq = 40
    mean_coverage = snp_data[1]
    parse_data = snp_data[0]

    for snp, data in parse_data.items():
        if int(data[7]) > int(allele2_freq):
            if int(data[0]) > int(mean_coverage):
                to_print.append([snp, data])
        else:
            if int(data[0]) > int(mean_coverage):
                to_print.append([snp, data])

    return to_print


def snp_or_not(snp_data, map_file):
    """Expect Data"""
    # Position
    # Number of reads aligned
    # "A" count
    # "C" count
    # "G" count
    # "T" count
    # Reference base
    # Allele freq
    # Allele freq2
    # Piled reads
    with open(map_file, 'r') as g:
        g.readline()
        dbSNP = {}
        for line in g:
            line = line.split('\t')
            if not line[3] in [i for i, _ in dbSNP.items()]:
                dbSNP[line[3]] = [line[0], line[5].strip()]

    for position, info in snp_data:
        if int(position) in [int(i) for i, yo in dbSNP.items()]:
            print('{}\t{}\t{}\t{}'.format(str(dbSNP[str(position)][0]),
                                          str(dbSNP[str(position)][1]),
                                          str(position),
                                          '\t'.join([str(i) for i in info])
                                          ))

        else:
            print('not_dbSNP\t-\t{}\t{}'.format(
                str(position),
                '\t'.join([str(i) for i in info])))


def parse_bam(args):
    """Parse all positions from Bam pileup"""
    header = ['RS',
              'REF',
              'POS',
              'ReadsNum',
              'A',
              'C',
              'G',
              'T',
              'RefAllele',
              'AlleleFreq',
              'AlleleFreq2',
              'ReadNum_samtools']
    if args.dbsnp:
        print('\t'.join(header))
    else:
        print('\t'.join(header[2:]))

    with open(args.ref_seq, 'r') as f:
        f.readline()
        ref_seq = f.readline()

    samfile = pysam.AlignmentFile(args.bam, 'rb')
    data = {}
    coverage = []

    for column in samfile.pileup():
        position = {'a': 0,
                    'c': 0,
                    'g': 0,
                    't': 0
                    }

        for read in column.pileups:
            if not read.is_del and not read.is_refskip:
                base = read.alignment.query_sequence[read.query_position]
                position[base.lower()] += 1

        reads_num = sum([i for _, i in position.items()])
        coverage.append(reads_num)
        sort = sorted(position.items(), key=lambda x: x[1])

        if reads_num:
            allele_freq = (sort[-1][1] * 100.0) / reads_num
            allele2_freq = (sort[-2][1] * 100.0) / reads_num
        else:
            allele_freq = 0
            allele2_freq = 0

        if sort[-1][0] != ref_seq[column.pos].lower() or int(allele2_freq) > 20:
            data[column.reference_pos + 1] = [reads_num,
                                              position['a'],
                                              position['c'],
                                              position['g'],
                                              position['t'],
                                              ref_seq[column.pos],
                                              round(allele_freq, 2),
                                              round(allele2_freq, 2),
                                              column.n]

    samfile.close()

    return (data, mean(coverage))


def args_parse():
    """Parse the args duhhh"""
    parser = argparse.ArgumentParser(description="Get SNPs for amplicons")
    parser.add_argument('bam', help="File name")
    parser.add_argument('ref_seq', help="One fasta record/file")
    parser.add_argument('--dbsnp', help="DBsnp info")

    return parser.parse_args()


if __name__ == '__main__':
    args = args_parse()
    parse = parse_bam(args)
    snps = parameters(parse)
    if args.dbsnp:
        snp_or_not(snps, args.dbsnp)
    else:
        for snp in snps:
            print('\t'.join(snp))
