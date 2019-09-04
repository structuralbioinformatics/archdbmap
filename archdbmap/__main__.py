# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>
"""
# Standard Libraries
import argparse
import os
import sys

# External Libraries
from SBI import SBIglobals
from SBI.external.blast import BlastExe
from archdbsql import ArchDBsql
import pandas as pd

# This Library
import functions


SBIglobals.verbose = True


def get_options():
    # Parse Arguments
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Inputs
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-pdblist', dest='pdblist', action='store',
                       help='List of PDB identifiers. Accepts formats: xxxx, xxxx.pdb, xxxx_x, xxxx_x.pdb to pick or not a chain of interest.')
    group.add_argument('-fasta', dest='fasta', action='store', help='FASTA file for the query, accepts multiFASTA.')

    # Blast Configuration
    group = parser.add_argument_group('blast')
    parser.add_argument('-blastexe', dest='blastexe', action='store', default='psiblast', help='Blast executable.')
    parser.add_argument('-blastpath', dest='blastpath', action='store', default='/usr/local/bin/', help='Blast path.')
    parser.add_argument('-makeblastdb', dest='makeblastdb', action='store', default='makeblastdb', help='Blast fasta to db formater.')
    parser.add_argument('-blastnofilter', dest='blastnofilter', action='store_true', default=False,
                        help='Do not filter matches by Rost\'s twilight curve.')

    # ArchDB Configuration
    parser.add_argument('-archdbFA', dest='archdbFA', action='store', help='Fasta of the ArchDB database. Includes the IDX file.')
    parser.add_argument('-archdbHST', dest='host', action='store', default='lpdisrv1.epfl.ch', help='Host for ArchDB SQL server.')
    parser.add_argument('-archdbDB', dest='dbname', action='store', default='archdb14mini', help='Database name for ArchDB SQL server.')
    parser.add_argument('-archdbUSR', dest='user', action='store', default='archdbuser', help='User for ArchDB SQL server.')
    parser.add_argument('-archdbPWD', dest='pswd', action='store', default=None, help='User password for ArchDB SQL server.')

    options = parser.parse_args()

    # Check ArchDB14 files and connections are present.
    archdbdownURL = 'http://sbi.upf.edu/archdb/download'
    # - FASTA files
    sys.stdout.write('Checking the existence of the ArchDB14 FASTA files.\n')
    options.archdbFA = os.path.abspath(options.archdbFA)
    options.archdbFAidx = options.archdbFA + '.idx'
    if not os.path.isfile(options.archdbFA) or not os.path.isfile(options.archdbFAidx):
        msg = ''
        if not os.path.isfile(options.archdbFA):
            msg += 'Missing file {}.\n'.format(options.archdbFA)
        if not os.path.isfile(options.archdbFAidx):
            msg += 'Missing file {}.\n'.format(options.archdbFAidx)
        msg += 'Files can be downloaded from {}\n'.format(archdbdownURL)
        raise IOError(msg)

    # - SQL connection
    sys.stdout.write('Checking the connection to the ArchDB14 SQL.\n')
    options.conn = ArchDBsql(dbhost=options.host, dbname=options.dbname, dbuser=options.user, dbpass=options.pswd)
    if not options.conn._db.ping():
        msg = 'Database is not accessible.\n'
        msg += 'SQL can be downloaded from {}\n'.format(archdbdownURL)
        raise IOError(msg)

    # Check BLAST availability
    sys.stdout.write('Linking Blast\n')
    options.blast = BlastExe.dynamic(executable=options.blastexe, path=options.blastpath, dbformater=options.makeblastdb)

    # Fix input
    options = functions.fasta_input( options )

    return options


def main():
    """Map Smotifs from the ArchDB14 database over queries by sequence.
    """
    options = get_options()

    cwd = os.getcwd()

    for query in options.fasta.live_show():
        print('PDB: {}; CHAIN: {}'.format(*query.id.split('_')))
        wdir = os.path.join(*query.id.split('_'))
        if not os.path.isdir(wdir):
            os.makedirs(wdir)
        os.chdir(wdir)
        blastR, IDXdata = functions.run_blast(query.id, query.sequence, options.archdbFA, options.archdbFAidx,
                                              options.blastexe, options.blastpath, options.makeblastdb,
                                              not options.blastnofilter)

        alignment, blastdata = functions.process_hits(options.conn, blastR, query.id, query.sequence, IDXdata)

        alidf = pd.DataFrame(alignment)[['id', 'seq', 'eval', 'src', 'loopid']]

        consensus = functions.make_consensus(alidf)

        with open('{}.archsearch'.format(query.id), 'w') as fd:
            fd.write('#ARCHMATCH\n')
            fd.write(alidf.to_csv(sep='\t', index=False, header=False) + '\n')
            fd.write('#MATCHCONSENSUS\n')
            fd.write(consensus.to_csv(sep='\t', index=False, header=False) + '\n')

        print('\n' + '---' * 50 + '\n')
        os.chdir(cwd)


if __name__ == '__main__':
    main()
