#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

bam_prep_for_lgt_seq.pl

=head1 SYNOPSIS

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 DESCRIPTION

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 AUTHOR - Karsten Sieber & David R. Riley

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use warnings;
use strict;
use Carp;
$Carp::MaxArgLen = 0;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions(
    \%options,         'input_list|I=s',  'input|i=s',       'name_sort_input=s', 'sort_mem=s',          'split_bam=i',   'encrypt=i',     'key=s',
    'prelim_filter=i', 'seqs_per_file=i', 'keep_softclip=i', 'Qsub|Q=i',          'excl=i',              'sub_mem=s',     'sub_name=s',    'threads|t=i',
    'projects=s',      'output_dir|o=s',  'subdirs=i',       'overwrite=s',       'samtools_bin=s',      'ergatis_dir=s', 'output_list=s', 'bin_dir=s',
    'fs=s',            'clovr=s',         'diag',            'verbose=i',         'print_hostname|ph=i', 'config_file=s', 'help|h',        'help_full',
    'tcga_dirs=i',     'aln_human=i',     'hg19_ref=s',      'no_gal=i',          'hostname=s',          'sub_mail=s',
) or die "Error: Unrecognized command line option. Please try again.\n";
use print_call;

# print_hostname(\%options);                                    ## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc.
use Scalar::Util qw(reftype);
use POSIX;
use run_cmd;
use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use LGTSeek;
use File::Basename;
use setup_input;
use setup_output;

if ( $options{help} ) {
    die "Help Basic Info: This script will remove M_M reads and split output bams into smaller bams. 
        --input=            <BAM>
        --name_sort_input=  <0|1> [0] 1= Resort the input bam by read names.  
        --output_dir=       Directory for output. 
        --help_full         Full Help Info\n";
}

if ( $options{help_full} ) {
    die "Help Full Info: This script will remove M_M reads, keeping M_U, U_U, M_M with Softclip. The output bams are split into smaller chunks. 
        ----------------------------------------------------------------------------------------
        --input|i=          Full path to input BAM. Fastq_1 is acceptable to map at human then prelim_filter. 
        --input_list|I=     <LIST of BAMS> 1 bam per line.
        ----------------------------------------------------------------------------------------
        --aln_human =       Align the input bam or fastq's to hg19 before filtering. **MUST** be name sorted input. 
        ----------------------------------------------------------------------------------------
        --paired_end=       <0|1> [1] 1= Paired End Sequencing reads.
        --prelim_filter=    <0|1> [1] 1= Filter out M_M reads, keeping MU,UU,and SC. 
          --keep_softclip=  <0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT) 
        --name_sort_input=  <0|1> [0] 1= Resort the input bam by read names.  
          --sort_mem=       [1G] Mem per thread to sort with. 
          --threads=        [1] # of threads 
        --split_bam=        <0|1> [1] 1= Split bam(s)
          --seqs_per_file=  <lines per bam> [50,000,000]
        --encrypt=          <0|1> [0] 1= encrypt ** untested **
          --key=            GPG key to use for encryption. ** untested **
        ----------------------------------------------------------------------------------------
        --Qsub|Q=           <0|1> [0] 1= Qsub this script for each input. 
          --project=        <project> [jdhotopp-lab] SGE group project to submit command under.
          --sub_mem=        [5G] --sub_mem MUST >= (--threads * --sort_mem)
          --sub_name=       < > Name of the SGE Job.
          --sub_mail=       [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
        ----------------------------------------------------------------------------------------
        --output_dir|o=     Directory for all output. Example: /path/to/{output_dir}/{tcga_dirs}/{subdirs}/ || /path/to/{output_dir}/{subdirs}/
         --tcga_dirs=       <0|1> [0] 1= Make the sub-dir prefix the input's last folder in path (Maintain TCGA analysis_id directory structure)
          --subdirs=        <0|1> [0] 1= Make the sub-dir prefix in output_dir based on input name.
        --output_list=      <0|1> [1] 1= Make a list of the output created.
        ----------------------------------------------------------------------------------------
        --overwrite=        <0|1> [0] 1= Overwrite previous files.
        ----------------------------------------------------------------------------------------
        --bin_dir=          [/local/projects-t3/HLGT/scripts/lgtseek/bin/] Directory where LGTSeek.pm is stored.
        ----------------------------------------------------------------------------------------
        --verbose           <0|1> [0] 1= Verbose reporting of progress. 0 =Turns off reports. 
        --help              Basic Help Info
        --help_full         Full Help Info
        --conf_file=        [~/.lgtseek.conf]
        ----------------------------------------------------------------------------------------\n";
}

if ( !$options{input} && !$options{input_list} ) { confess "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n"; }

## Set default values
$options{prelim_filter} = defined $options{prelim_filter} ? "$options{prelim_filter}" : "1";
$options{sub_mem}       = defined $options{sub_mem}       ? "$options{sub_mem}"       : "5G";
$options{output_list}   = defined $options{output_list}   ? "$options{output_list}"   : "1";

my $lgtseek = LGTSeek->new2( \%options );

my $input = setup_input( \%options );    ## $input is a array ref.
my $x     = 0;

my $original_output_dir = $lgtseek->{output_dir};

foreach my $input (@$input) {
    $x++;
    my ( $fn, $path, $suf ) = fileparse( $input, ( '_resorted.bam', '.bam' ) );
    my $subdir     = $fn;
    my @split_path = split( /\//, $path );
    my $tcga_dir   = $split_path[-1];
    $lgtseek->{output_dir} = $original_output_dir;
    if ( $lgtseek->{tcga_dirs} == 1 ) { $lgtseek->{output_dir} = $lgtseek->{output_dir} . "$tcga_dir\/"; }
    if ( $lgtseek->{subdirs} == 1 )   { $lgtseek->{output_dir} = $lgtseek->{output_dir} . "$subdir\/"; }
    $options{output_dir} = $lgtseek->{output_dir};
    $lgtseek->_run_cmd("mkdir -p -m u=rwx,g=rwx,o= $lgtseek->{output_dir}");

    ## Qsub this script foreach input and any of the options passed
    if ( $lgtseek->{Qsub} == 1 ) {
        ## If we are in the orignal call, change input from list to a single file
        if ( $options{input_list} ) { $options{input} = $input; }
        ## Check $sub_mem is enough for sorting
        my $original_sub_mem;
        my $original_sort_mem;
        if ( $options{sub_mem} ) {
            if ( $options{sub_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sub_mem = $1; }
        }
        if ( $options{sort_mem} ) {
            if ( $options{sort_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sort_mem = $1; }
        }
        if ( $original_sub_mem < ( $original_sort_mem * $options{threads} ) ) {
            $options{sub_mem} = ( ceil( ( $original_sort_mem * $options{threads} ) * 1.1 ) ) + 1 . "G";
        }

        ## Build qsub command
        my $cmd = "$^X $0";
        foreach my $key ( keys %options ) {
            next if ( $key =~ /Qsub/ );
            next if ( $key =~ /subdirs/ );
            next if ( $key =~ /tcga_dirs/ );
            next if ( $options{input_list} && $key =~ /input_list/ );    ## If we are in the orignal call, we don't want to qsub more lists
            if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
        }
        my $sub_name = $options{sub_name} ? $options{sub_name} : "prelim$x";
        ## submit command to grid
        Qsub(
            {   cmd      => $cmd,
                wd       => $lgtseek->{output_dir},
                sub_name => $sub_name,
                sub_mem  => $options{sub_mem},
                threads  => $lgtseek->{threads},
                project  => $lgtseek->{project},
                hostname => $options{hostname},
                excl     => $options{excl},
            }
        );

        ## Skip to next input for qsub
        next;
    }

    $options{output_dir} = $lgtseek->{output_dir};
    print_notebook( \%options );

    # Align to Human.
    my $human_aln = defined $options{aln_human} ? $options{aln_human} : "0";
    if ( $human_aln == 1 ) {
        print STDERR "======== RUNBWA-Human ========\n";
        my $human_bam1;
        ## Map input bam @ hg19
        if ( $input =~ /\.bam$/ ) {
            $human_bam1 = $lgtseek->runBWA(
                {    ## &runBWA returns an array
                    input_bam   => $input,
                    output_bam  => 1,
                    threads     => $lgtseek->{threads},
                    output_dir  => "$lgtseek->{output_dir}/aln_human/",
                    reference   => $lgtseek->{hg19_ref},
                    overwrite   => $lgtseek->{overwrite},
                    cleanup_sai => 1,
                }
            );
        }
        ## Map fastqs @ hg19
        elsif ( $input =~ /.fq$/ || $input =~ /.fastq$/ || $input =~ /.fastq.gz$/ ) {
            my ( $in1, $in2 ) = split( /,/, $input );    ## split input if needed
            my ( $input_base, $input_dir, $input_suffix ) = fileparse( $in1, @{ $lgtseek->{fastq_suffix_list} } );
            ## &runBWA returns an array
            $human_bam1 = $lgtseek->runBWA(
                {   input_dir   => $input_dir,
                    input_base  => $input_base,
                    output_bam  => 1,
                    threads     => $lgtseek->{threads},
                    output_dir  => "$lgtseek->{output_dir}/aln_human/",
                    reference   => $lgtseek->{hg19_ref},
                    overwrite   => $lgtseek->{overwrite},
                    cleanup_sai => 1,
                }
            );
        }
        $input = @$human_bam1[0];
    }

    ## PrelimFilter to remove M_M reads.
    my $bams;
    if ( $lgtseek->{prelim_filter} == 1 || $lgtseek->{split_bam} == 1 || $lgtseek->{name_sort_input} == 1 ) {
        $bams = $lgtseek->prelim_filter(
            {   input_bam       => $input,
                output_dir      => $lgtseek->{output_dir},
                name_sort_input => $lgtseek->{name_sort_input},    ## Default = 0
                sort_mem        => $lgtseek->{sort_mem},           ## Default = 1G lgtseek default.
                split_bam       => $lgtseek->{split_bam},          ## Default = 1
                seqs_per_file   => $lgtseek->{seqs_per_file},      ## Default = 50M
                keep_softclip   => $lgtseek->{keep_softclip},      ## Default = 1
                overwrite       => $lgtseek->{overwrite},          ## Default = 0
            }
        );
    }

    ## Encrypt the file(s)
    my @encrypted;
    if ( $lgtseek->{encrypt} == 1 ) {
        foreach my $files (@$bams) {
            my $cmd = "gpg -o $files\.gpg --default-key $options{key} -r $options{key} -e $files";
            if   ( $lgtseek->{Qsub} == 1 ) { Qsub($cmd); }
            else                           { $lgtseek->_run_cmd($cmd); }
            push( @encrypted, "$files\.gpg" );
        }
    }

    ## Print out the output list
    if ( $lgtseek->{output_list} == 1 ) {

        # Open a list file to write the output bams to
        open( my $olistfh, ">$lgtseek->{output_dir}/output.list" ) or confess "Unable to open: $lgtseek->{output_dir}/output.list because: $!\n";
        if ( $lgtseek->{encrypt} == 1 ) {
            foreach my $out (@encrypted) { print $olistfh "$out\n"; }
        }
        elsif ( $lgtseek->{split_bam} == 1 || $lgtseek->{prelim_filter} == 1 ) {
            foreach my $out2 (@$bams) { print $olistfh "$out2\n"; }
        }
    }

    print STDERR "+++ Complete: $0 on: $input\t+++\n";
}

__END__

# my $hostname = ( defined $options{hostname} ) ? "$options{hostname}" : "*";
# if   ( defined $options{hostname} ) { $no_gal = "0"; }
# else                                { $no_gal = defined $options{no_gal} ? "$options{no_gal}" : "0"; }
