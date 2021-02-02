package DeconSeqConfig;

use strict;

use constant DEBUG => 0;
use constant PRINT_STUFF => 1;
use constant VERSION => '0.4.3';
use constant VERSION_INFO => 'DeconSeq version '.VERSION;

use constant ALPHABET => 'ACGTN';

use constant DB_DIR => '/N/dc2/projects/solanumgenome/library/deconseq_DB/';
use constant TMP_DIR => '/N/dc2/projects/solanumgenome/library/deconseq_DB/temp/';
use constant OUTPUT_DIR => ;

use constant PROG_NAME => 'bwa64';  # should be either bwa64 or bwaMAC (based on your system architecture)
use constant PROG_DIR => ;      # should be the location of the PROG_NAME file (use './' if in the same location at the perl script)

use constant DBS => {human => {name => 'Human Reference GRCh38',  #database name used for display and used as input for -dbs and -dbs_retain
                               db => 'human'},            #database name as defined with -p for "bwa index -p ..." (specify multiple database chunks separated with commas without space; e.g. hs_ref_s1,hs_ref_s2,hs_ref_s3)
                     bacteria => {name => 'Bacterial genomes',
                              db => 'bacteria_c1,bacteria_c2,bacteria_c3,bacteria_c4,bacteria_c5,bacteria_c6,bacteria_c7,bacteria_c8,bacteria_c9,bacteria_c10'},
                     viral => {name => 'Viral genomes',
                             db => 'viral'},
                     solanum => {name => 'Solanum genomes',
                             db => 'lyco,tube,penn'}};
use constant DB_DEFAULT => 'human';

#######################################################################

use base qw(Exporter);

use vars qw(@EXPORT);

@EXPORT = qw(
             DEBUG
             PRINT_STUFF
             VERSION
             VERSION_INFO
             ALPHABET
             PROG_NAME
             PROG_DIR
             DB_DIR
             TMP_DIR
             OUTPUT_DIR
             DBS
             DB_DEFAULT
             );

1;
