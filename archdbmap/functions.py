# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>
"""
# Standard Libraries
import os
import sys
import itertools

# External Libraries
import wget
import pandas as pd

# This Library


def fasta_input( options ):
    """Convert PDB inputs to FASTA if needed.
    """
    from SBI.sequence import Sequence
    from SBI.sequence.Fasta import Fasta
    from SBI.structure import PDB

    # We already have the fasta sequences
    if options.fasta is not None:
        options.fasta = Fasta(options.fasta)
        return options

    # Process PDB list
    pdblist = [x.replace('.pdb', '').strip().split('_') for x in open(options.pdblist).readlines()]

    def get_sequence(pdb):
        seqs = []
        sys.stdout.write(pdb[0] + '\n')
        if not os.path.isdir(pdb[0]):
            os.mkdir(pdb[0])
        pdbf = os.path.join(pdb[0], '{}.pdb'.format(pdb[0]))
        if not os.path.isfile(pdbf):
            wget.download('http://files.rcsb.org/view/{}.pdb'.format(pdb[0]), out=pdbf)
        pdbstr = PDB(pdbf)
        qchains = pdb[1] if len(pdb) > 1 else pdbstr.chain_identifiers
        for chain in qchains:
            print('\t--' + chain + '--')
            seqs.append(Sequence('{}_{}'.format(pdb[0], chain), pdbstr.get_chain_by_id(chain).protein_sequence))
        return seqs

    sys.stdout.write('Processing PDB inputs\n')
    options.fasta = Fasta.build_multifasta('{}.fa'.format(options.pdblist), itertools.chain(*map(get_sequence, pdblist)), True)
    return options


def run_blast(query_id, query_sequence, db, idxf, blastexe, blastpath, blastdb, dofilter):
    from SBI.external.blast import BlastExe
    from SBI.sequence import Fasta

    # Run Blast
    fasta = Fasta.build('{}.fa'.format(query_id), query_id, query_sequence, force=True)
    print(os.path.abspath(str(fasta.file)), os.path.abspath(db))
    BlastExe.dynamic(executable=blastexe, path=blastpath, dbformater=blastdb)
    blst = BlastExe(db, clean=False)
    print(blst._EXE.full_executable)

    blastResult = blst.execute_query_seq(
        sequenceID=query_id,
        sequence=query_sequence,
        blast_input_file=os.path.join(os.getcwd(), 'blast.fa.tmp'),
        blast_output_file=os.path.join(os.getcwd(), 'blast.out.tmp'))
    blastResult.correct_hit_count(count_hit_file=idxf)
    if dofilter:
        blastResult.set_hit_filter(tz_type='ID', tz_parameter=0)

    idx_data = {}
    fd = open(idxf)
    for line in fd:
        if len(line.strip()) > 0:
            k = line.split('\t')
            idx_data[k[0].lstrip('>')] = k[1].strip().split(';')
    fd.close()

    return blastResult, idx_data


def process_hits(dbcon, blastResult, query_id, query_sequence, idx_data):

    alignment = {'id': [], 'seq': [], 'src': [], 'eval': [], 'loopid': []}
    alignment['id'].append(query_id)
    alignment['seq'].append(query_sequence)
    alignment['src'].append("TRG")
    alignment['loopid'].append("-")
    alignment['eval'].append(0)
    blastdata = []
    print('Checking a total of {} raw hits'.format(len(blastResult)))
    print('Checking a total of {} hits over the threshold'.format(len(blastResult.hits)))
    for blastHit in blastResult.hits:
        print('Analyzing hit {0} .......................'.format(blastHit.sequenceID))
        print(blastHit)
        print('.........................................')
        blastdata.append(str(blastHit).split())
        blastdata[-1].append('global')
        pLoops = dbcon.get_loops_for_protein(blastHit.sequenceID, 'MCL')
        for loop in pLoops:
            loop = {'name':   loop[0], 'ss1l': loop[1], 'ss2l': loop[2],
                    'length': loop[3], 'ini':  loop[4], 'end':  loop[5],
                    'subclass':     loop[6].split(','),
                    'subclass_nid': loop[7].split(','), 'sse': loop[8], 'seq': loop[9]}
            try:
                loopali = blastHit.get_section_from_sequence_position('hit', loop['ini'], loop['end'])
                if loopali.gaps != 0:
                    if "-" in loopali.hit_seq[:int(loop['ss1l'])]:
                        raise AttributeError
                    if "-" in loopali.hit_seq[-int(loop['ss2l']):]:
                        raise AttributeError
                    if "-" in loopali.query_seq[:int(loop['ss1l'])]:
                        raise AttributeError
                    if "-" in loopali.query_seq[-int(loop['ss2l']):]:
                        raise AttributeError
                qp = loopali.query_pos
                if len(qp) == 0:
                    continue
                blastdata.append(str(loopali).split())
                blastdata[-1].append(loop['subclass'][0].replace("_", ""))
                alignment['id'].append(loopali.sequenceID.replace('_', ''))
                alignment['eval'].append(blastHit.evalue)
                alignment['src'].append(str(loop['subclass'][0]))
                alignment['loopid'].append(str(loop['name']))
                hitseq = [".", ] * len(query_sequence)
                seqtype = 'seq'
                if loopali.gaps == 0:
                    hitseq[int(qp[0]) - 1:int(qp[1])] = list(loop[seqtype])
                else:
                    in1 = [int(qp[0]) - 1, int(qp[0]) + int(loop['ss1l']) - 1]
                    in2 = [int(qp[-1]) - int(loop['ss2l']), int(qp[-1])]
                    hitseq[int(qp[0]) - 1:int(qp[-1])] = ["-", ] * len(hitseq[int(qp[0]) - 1:int(qp[-1])])
                    hitseq[in1[0]:in1[1]] = list(loop[seqtype][:int(loop['ss1l'])])
                    hitseq[in2[0]:in2[1]] = list(loop[seqtype][-int(loop['ss2l']):])
                alignment['seq'].append(str("".join(hitseq)))
            except IndexError:
                print('{} Hit segment does not contain classified Smotif assignation'.format(blastHit.sequenceID))
            except AttributeError:
                print('{} Hit contains gaps in the secondary structures'.format(blastHit.sequenceID))

    wi = []
    for i, k in enumerate(alignment['seq']):
        if len(k) != len(query_sequence):
            wi.append(i)
    for i in reversed(wi):
        alignment['id'].pop(i)
        alignment['eval'].pop(i)
        alignment['src'].pop(i)
        alignment['seq'].pop(i)
        alignment['loopid'].pop(i)

    return alignment, blastdata


def make_consensus(alidf):
    consensus = {'id': [alidf.iloc[0].id, ], 'seq': [alidf.iloc[0].seq, ], 'eval': [alidf.iloc[0].eval, ]}

    def unify(col):
        ops = list(set(col.unique()).difference(set(['.'])))
        if len(ops) == 0:
            return '.'
        else:
            return ops[0]

    for seqid, grp in alidf.groupby('id'):
        if seqid == consensus['id'][0]:
            continue
        dfg = grp.seq.apply(lambda x: pd.Series(list(x)))
        consensus['id'].append(seqid)
        consensus['seq'].append(''.join(list(dfg.apply(unify))))
        consensus['eval'].append(grp.iloc[0].eval)

    return pd.DataFrame(consensus)[['id', 'seq', 'eval']]
