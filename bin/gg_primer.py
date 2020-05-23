import sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Restriction import BsmBI
from pydna.dseq import Dseq


# using benchling parameters
    # dna_conc = 1000 # nM
    # NA_conc = 50 # mM
    # Mg_conc = 0 # mM
    # dNTPs_conc = 0.2 # mM

# myseq = 'ATGCATGC'
# myseq = 'CAGTCAGTACGTACGTGTACTGCCGT'
# print(melting.temp(myseq, DNA_c=dna_conc, Na_c=NA_conc, Mg_c=Mg_conc, dNTPs_c=dNTPs_conc))
def pTMP139():
    ''' {NoneType} -> str

    outputs the full pTMP139 plasmid sequence


    >>> my_plas = pTMP139()
    >>> my_plas == 'ccatatgaagcagcatgacttcttcaagtccgccatgccggaaggctatgtgcaggaacgcacgatttcctttaaggatgacggcacgtacaaaacgcgtgcggaagtgaaatttgaaggcgataccctggtaaaccgcattgagctgaaaggcattgactttaaagaagacggcaatatcctgggccataagctggaatacaattttaacagccacaatgtttacatcaccgccgataaacaaaaaaatggcattaaagcgaattttaaaattcgccacaacgtggaggatggcagcgtgcagctggctgatcactaccagcaaaacactccaatcggtgatggtcctgttctgctgccagacaatcactatctgagcacgcaaagcgttctgtctaaagatccgaacgagaaacgcgatcatatggttctgctggagttcgtaaccgcagcgggcatcacgcatggtatggatgaactgtacaaatgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatacgtctctgaccagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtagtcggcgagacggaaagtgaaacgtgatttcatgcgtcattttgaacattttgtaaatcttatttaataatgtgtgcggcaattcacatttaatttatgaatgttttcttaacatcgcggcaactcaagaaacggcaggttcggatcttagctactagagaaagaggagaaatactagatgcgtaaaggcgaagagctgttcactggtgtcgtccctattctggtggaactggatggtgatgtcaacggtcataagttttccgtgcgtggcgagggtgaaggtgacgcaactaatggtaaactgacgctgaagttcatctgtactactggtaaactgccggttccttggccgactctggtaacgacgctgacttatggtgttcagtgctttgctcgttatccgga'
    True

    '''

    output = 'ccatatgaagcagcatgacttcttcaagtccgccatgccggaaggctatgtgcaggaacgcacgatttcctttaaggatgacggcacgtacaaaacgcgtgcggaagtgaaatttgaaggcgataccctggtaaaccgcattgagctgaaaggcattgactttaaagaagacggcaatatcctgggccataagctggaatacaattttaacagccacaatgtttacatcaccgccgataaacaaaaaaatggcattaaagcgaattttaaaattcgccacaacgtggaggatggcagcgtgcagctggctgatcactaccagcaaaacactccaatcggtgatggtcctgttctgctgccagacaatcactatctgagcacgcaaagcgttctgtctaaagatccgaacgagaaacgcgatcatatggttctgctggagttcgtaaccgcagcgggcatcacgcatggtatggatgaactgtacaaatgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatacgtctctgaccagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtagtcggcgagacggaaagtgaaacgtgatttcatgcgtcattttgaacattttgtaaatcttatttaataatgtgtgcggcaattcacatttaatttatgaatgttttcttaacatcgcggcaactcaagaaacggcaggttcggatcttagctactagagaaagaggagaaatactagatgcgtaaaggcgaagagctgttcactggtgtcgtccctattctggtggaactggatggtgatgtcaacggtcataagttttccgtgcgtggcgagggtgaaggtgacgcaactaatggtaaactgacgctgaagttcatctgtactactggtaaactgccggttccttggccgactctggtaacgacgctgacttatggtgttcagtgctttgctcgttatccgga'

    return output

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

def gen_expected_plasmid(s, type):
    ''' {str, str} -> str

    s    = the sequence to domesticate
    type = part type

    The function does the following:
        1) attaches overhangs in the domestication sequence
        2) Digests the sequence to domesticate with BsmBI
        3) Digests pTMP139 with BsmBI
        4) Inserts the sequence to domesticate into the pTMP139 backbone
        5) Outputs the full plasmid

    NOTE: rev sequence is NOT reverse complemented!

    # test 1)
    >>> myp = parts_type()
    >>> seq = 'gggggggggg'
    >>> digest_s = myp['3a']
    >>> output = digest_s[0] + seq + digest_s[1]
    >>> output
    'GCATCGTCTCATCGGTCTCATATGggggggggggTTCTTGAGACCTGAGACGGCAT'

    # test 2)
    >>> test_seq = 'GCATCGTCTCATCGGTCTCATATGggggggggggTTCTTGAGACCTGAGACGGCAT'
    >>> test_seq = Dseq(test_seq)
    >>> cut_seq = test_seq.cut(BsmBI)[1]

    # There is an extra "GACC" at the 3' end of this sequence. this is the complement of the overhang enzyme makes.
    >>> str(cut_seq)
    'TCGGTCTCATATGggggggggggTTCTTGAGACC'

    # test 3)
    >>> bb = pTMP139()
    >>> bb = Dseq(bb)
    >>> cut_bb = bb.cut(BsmBI)[1]
    >>> str(cut_bb)
    'gaccagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtagtcgg'

    # test 4)
    # generated expected sequence here based on Benchling assembly of "assembly_test_bb"
    >>> theoretic_plas = (cut_bb + cut_seq).looped()
    >>> str(theoretic_plas)
    'gaccagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtagTCGGTCTCATATGggggggggggTTCTTGA'

    '''

    # 1)
    myparts = parts_type()
    digest_seq = myparts[type]
    amplicon = digest_seq[0] + s + digest_seq[1]

    # 2)
    amplicon = Dseq(amplicon)
    # digest with BsmBI
    cut_amplicon = amplicon.cut(BsmBI)[1]

    # 3)
    backbone = pTMP139()
    backbone = Dseq(backbone)
    cut_backbone = backbone.cut(BsmBI)[1]

    # 4)
    theoretic_plasmid = (cut_backbone + cut_amplicon).looped()

    # 5)
    return theoretic_plasmid







if __name__ == '__main__':
    import doctest
    doctest.testmod()
