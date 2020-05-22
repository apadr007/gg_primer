import sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# using benchling parameters
    # dna_conc = 1000 # nM
    # NA_conc = 50 # mM
    # Mg_conc = 0 # mM
    # dNTPs_conc = 0.2 # mM

# myseq = 'ATGCATGC'
# myseq = 'CAGTCAGTACGTACGTGTACTGCCGT'
# print(melting.temp(myseq, DNA_c=dna_conc, Na_c=NA_conc, Mg_c=Mg_conc, dNTPs_c=dNTPs_conc))

def parts_type():
    ''' {NoneType} -> dict

    this function returns a dict that contains the sequence overhangs
    for gg cloning based on part names.

    Following guidelines from Dueber YTK paper

    dict values are a list. the first element in the list is the 5' end
    the second item in the list is the sequence that belongs in the 3' end

    both sequences are 5' -> 3' ## VERIFY THIS


    >>> d = parts_type()
    >>> isinstance(d, dict)
    True


    '''

    parts = {'2': ['GCATCGTCTCATCGGTCTCAAACG', 'TATGTGAGACCTGAGACGGCAT'],
    '3a': ['GCATCGTCTCATCGGTCTCATATG', 'TTCTTGAGACCTGAGACGGCAT'],
    '3b': ['GCATCGTCTCATCGGTCTCATTCT', 'ATCCTGAGACCTGAGACGGCAT'],
    '4a': ['GCATCGTCTCATCGGTCTCAATCC', 'TGGCTGAGACCTGAGACGGCAT']}

    return parts


def fwd_tm(s, temp):
    ''' {str, float} -> str

    Design a primer for a sequence  "s" at a desired tm "temp".
    I am using the mt.Tm_NN() function for primer melting calculation.

    Specifically, I found that the nn_table=mt.DNA_NN4 parameter to be
    what Benchling uses for their own primer design.

    # test for str output
    >>> myseq = 'GCATCGTCTCATCGGTCTCAAACGTTCCCGTAAATGCATCAG'
    >>> t = 55
    >>> output = fwd_tm(s=myseq, temp=t)
    >>> type(output) == str
    True

    # test for expected output (cross referenced with Benchling tm = ~56C)
    >>> output
    'GCATCGTCTCATCGGTCTCAAAC'
    '''
    s = str(s)
    tm = 0
    fwd_pri = ''
    for letter in s:
        fwd_pri = ''.join([fwd_pri, letter])
        if tm >= int(temp):
            break
        tm = mt.Tm_NN(fwd_pri, nn_table=mt.DNA_NN4)
        tm = round(tm, 2)

    return fwd_pri

def rev_tm(s, temp):
    ''' {str, float} -> str

    Design a primer for a sequence  "s" at a desired tm "temp".
    I am using the mt.Tm_NN() function for primer melting calculation.

    Specifically, I found that the nn_table=mt.DNA_NN4 parameter to be
    what Benchling uses for their own primer design.

    # test for str output
    >>> myseq = 'GCATCGTCTCATCGGTCTCAAACGttcccgtaaatgcatcagttcccgtaaatgcatcagttcccgtaaatgcatcagttcccgtaaatgGCGATTGGAGTGGATAAATTCACTAGTCTAGAGGGT'
    >>> t = 55
    >>> output = rev_tm(s=myseq, temp=t)
    >>> type(output) == str
    True

    # test for expected output (cross referenced with Benchling tm = ~56C)
    >>> output
    'GTGGATAAATTCACTAGTCTAGAGGGT'
    '''

    tm = 0
    rev_pri = ''

    # iterate through string in reverse
    for letter in reversed(s):
        # maintain 5' -> 3' direction
        rev_pri = ''.join([letter, rev_pri])
        if tm >= int(temp):
            break
        tm = mt.Tm_NN(rev_pri, nn_table=mt.DNA_NN4)
        tm = round(tm, 2)

    return rev_pri

def append_overhangs(myseqs, type):
    ''' {list of Str, Str} -> list of Str

    Attach the appropriate overhang for GG cloning.

    type is a PART TYPE sequence in a 5' -> 3' manner,
    which contains a BsmBI and a BsaI site for cloning
    and domestication.

    NOTE: for the reverse primer, returns the reverse complement
    in a manner ready for ordering


    >>> myseq = ['TTATATATATATAT', 'CGCGCCGCGCGCG']
    >>> myparts = parts_type()
    >>> output = append_overhangs(myseqs=myseq, type='4a')

    #test FWD primer with overhangs
    >>> output[0][2]
    'GCATCGTCTCATCGGTCTCAATCCTTATATATATATAT'
    >>> output[0][1]
    '4a'

    #reverse complement of REV primer & overhangs
    >>> output[1][2]
    'ATGCCGTCTCAGGTCTCAGCCACGCGCGCGGCGCG'

    '''

    myparts = parts_type()

    fwd_primer = myparts[type][0] + myseqs[0]

    rev_primer = myseqs[1] + myparts[type][1]

    # generate reverse complement primer
    rev_primer = Seq(rev_primer)
    rev_primer_rc = rev_primer.reverse_complement()

    seq = [ ['fwd', type, fwd_primer], ['rev', type, str(rev_primer_rc)] ]

    return seq



if __name__ == '__main__':
    import doctest
    doctest.testmod()
