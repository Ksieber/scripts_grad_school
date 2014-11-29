#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

lgtseq_analysis.pl

=head1 SYNOPSIS

Search an bam for bacterial-human LGT.

=head1 DESCRIPTION


=head1 AUTHOR - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut

my $LGTSEQ_ANALYSIS = '1.00';

use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,        'input|i=s',           'input_list|I=s',   'Qsub|q=i',        'sub_mem|mf=s',      'sub_name=s',        'excl=i',           'no_gal=i',
    'hostname=s',     'decrypt=i',           'url=s',            'prelim_filter=i', 'name_sort_input=i', 'keep_softclip=i',   'split_bam=i',      'seqs_per_file=i',
    'aln1_human=i',   'aln2_human=i',        'split_bac_list=s', 'hg19_ref=s',      'refseq_list=s',     'output_dir|o=s',    'subdirs=i',        'tcga_dirs=i',
    'lgt_coverage=i', 'max_overlap=i',       'min_length=i',     'bin_dir=s',       'samtools_bin=s',    'ergatis_bin=s',     'prinseq_bin=s',    'donor_lineage=s',
    'host_lineage=s', 'threads|t=i',         'taxon_host=s',     'taxon_dir=s',     'taxon_idx_dir=s',   'path_to_blastdb=s', 'best_hits_only=i', 'evalue_cutoff=s',
    'verbose|V=i',    'print_hostname|ph=i', 'conf_file=s',      'help|h',          'help_full|?',       'workflow_help',     'conf_help',        'sub_mail=s',
    'Qsub_iterate=i', 'overwrite=i',         'project=s',
) or die "Error: Unrecognized command line option. Please try again.\n";

## Check if the user needs help information
if ( $options{help} )          { &help; }             ## @ end of script
if ( $options{help_full} )     { &help_full; }        ## @ end of script
if ( $options{workflow_help} ) { &workflow_help; }    ## @ end of script
if ( $options{conf_help} )     { &conf_help; }        ## @ end of script

use print_call;
$options{print_hostname} = $options{print_hostname} ? $options{print_hostname} : "0";
print_hostname( \%options );                          ## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc.
### May need to change this depending on where the script is being run
use LGTSeek;
use Time::SoFar;
use run_cmd;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;

## Check we have the necessary inputs
if ( !$options{input} and !$options{input_list} ) {
    confess "Error: Please give an input.bam with --input=<FILE> or --input_list=<LIST>. Try again or use --help_full.\n";
}
if ( !$options{output_dir} ) {
    print "It is HIGHLY recommended you STOP, restart, and use a --output_dir=</some/path/>.\n";
    sleep 60;
}

## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2( \%options );

## Qsub the job instead of running it
if ( $options{Qsub} or $options{Qsub_iterate} ) {
    if ( !$options{sub_name} ) { $options{sub_name} = "lgtseq"; }
    if ( !$options{project} )  { $options{project}  = $lgtseek->{project}; }
    Qsub_script( \%options );
}

## Print the script call
print_call( \%options, "LGTSEQ_ANALYSIS_VERSION=$LGTSEQ_ANALYSIS" );

## Setup array ref of inputs
my $inputs = setup_input( \%options );

foreach my $input (@$inputs) {
    ## Setup output directory
    my ( $name, $path, $suf ) = fileparse( $input, $lgtseek->{suffix_regex} );
    my $output_dir = $lgtseek->{output_dir} ? $lgtseek->{output_dir} : "$path/lgtseq/";
    chomp $name;
    my $subdir     = $name;
    my @split_path = split( /\//, $path );
    my $tcga_dir   = $split_path[-1];
    if ( $lgtseek->{tcga_dirs} == 1 ) { $output_dir = $output_dir . "/$tcga_dir\/" unless ( $output_dir =~ /$tcga_dir\/*$/ ); }
    if ( $lgtseek->{subdirs} == 1 )   { $output_dir = $output_dir . "/$subdir\/"   unless ( $lgtseek->{output_dir} =~ /$subdir\/*$/ ); }
    $output_dir =~ s/\/{2,}/\//g;
    run_cmd("mkdir -p -m u=rwx,g=rwx,o= $output_dir");

    print_notebook( \%options );

    # Primary aln to Human.
    if ( $lgtseek->{aln1_human} ) {
        print STDERR "======== RUNBWA-Human1 ========\n";
        my $human_bam1;
        ## Map input bam @ hg19
        if ( $input =~ /\.bam$/ ) {
            $human_bam1 = $lgtseek->runBWA(
                {    ## &runBWA returns an array
                    input_bam   => $input,
                    output_bam  => 1,
                    threads     => $lgtseek->{threads},
                    output_dir  => "$output_dir/human_aln1/",
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
                    output_dir  => "$output_dir/human_aln1/",
                    reference   => $lgtseek->{hg19_ref},
                    overwrite   => $lgtseek->{overwrite},
                    cleanup_sai => 1,
                }
            );
        }
        $input = $human_bam1->[0];
    }

    # Prelim_filter
    if ( $lgtseek->{prelim_filter} or $lgtseek->{name_sort_input} ) {
        print STDERR "======== Prelim_filter ========\n";
        my $prelim_filtered_bam = $lgtseek->prelim_filter(
            {    ## &prelim_filter returns an array
                input_bam       => $input,
                output_dir      => "$output_dir/prelim_filter/",    ## vv lgtseek &prelim_filter defaults vv
                name_sort_input => $lgtseek->{name_sort_input},     ## Default = 0
                sort_mem        => $lgtseek->{sort_mem},            ## Default = 1G lgtseek default. lgt_prep overides to 40G.
                threads         => $lgtseek->{threads},             ## Default = 1
                split_bam       => "0",                             ## Default = 1
                keep_softclip   => $lgtseek->{keep_softclip},       ## Default = 1, better to split with lgt_seq_prelim
                overwrite       => $lgtseek->{overwrite},           ## Default = 0
            }
        );                                                          ## ^^                              ^^
        $input = $prelim_filtered_bam->[0];                         ## This works because we are not splitting bams
    }

    # Secondary aln human
    if ( $lgtseek->{aln2_human} ) {
        print STDERR "======== RUNBWA-Human2 ========\n";
        my $human_bam2 = $lgtseek->runBWA(
            {                                                       ## &runBWA returns an array
                input_bam   => $input,
                output_bam  => 1,
                threads     => $lgtseek->{threads},
                output_dir  => "$output_dir/human_aln2/",
                reference   => $lgtseek->{hg19_ref},
                overwrite   => $lgtseek->{overwrite},
                cleanup_sai => 1,
            }
        );
        $input = $human_bam2->[0];
    }

    # Align to the Bacteria.
    print STDERR "======== RUNBWA-Bacteria ========\n";
    my $bacterial_bams = $lgtseek->runBWA(
        {    ## &runBWA returns an array
            input_bam      => $input,
            output_bam     => 1,
            threads        => $lgtseek->{threads},
            output_dir     => "$output_dir/bacterial_alignments/",
            reference_list => $lgtseek->{split_bac_list},
            overwrite      => $lgtseek->{overwrite},
            cleanup_sai    => 1,
        }
    );

    # Postprocess the BWA mappings to human and bacteria
    print STDERR "======= POSTPROCESS =======\n";
    my $host_bams = ["$input"];                 ## Put human input bam into array ref for post-proc.
    my $pp_data   = $lgtseek->bwaPostProcess(
        {   donor_bams    => $bacterial_bams,
            host_bams     => $host_bams,
            output_dir    => $output_dir,
            output_prefix => $name,
            overwrite     => $lgtseek->{overwrite},
        }
    );

    # Create file with number of counts
    my @header = ('run_id');
    my @vals   = ($name);
    open OUT, ">$output_dir/$name\_post_processing.tab" or die;
    map {
        push( @header, $_ );
        my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
        push( @vals, $foo );
    } ( 'total', 'host', 'no_map', 'all_map', 'single_map', 'integration_site_human', 'integration_site_bac', 'microbiome', 'lgt' );
    &print_tab( "$output_dir/$name\_post_processing.tab", \@header, \@vals );

    # Clean up output we don't need anymore; This is IMPORTANT on nodes. Not so much on filesystem.
    $lgtseek->_run_cmd("rm -rf $output_dir/human_aln1/");
    $lgtseek->_run_cmd("rm -rf $output_dir/human_aln2/");

    # $lgtseek->_run_cmd("rm -rf $output_dir/prelim_filter/"); #FOOBAR
    $lgtseek->_run_cmd("rm -rf $output_dir/bacterial_alignments/");

    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## LGT Analysis.
    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## Check to make sure we found LGT.
    print STDERR "======== LGT ========\n";
    if ( $lgtseek->empty_chk( { input => $pp_data->{files}->{lgt_donor} } ) ) {
        print STDERR "***Warning*** No LGT in: $pp_data->{files}->{lgt_donor}\. Skipping further LGT analysis.\n";
    }
    else {
        # Prinseq filter the putative lgts
        # This is for manual curation later. Use all potential LGT for further analysis.
        print STDERR "======== LGT-PRINSEQ =======\n";
        my $filtered_bam = $lgtseek->prinseqFilterBam(
            {   input_bam  => $pp_data->{files}->{lgt_host},
                output_dir => $output_dir,
            }
        );

        # Add filtered count to counts.
        push( @header, 'lgt_pass_prinseq' );
        push( @vals,   $filtered_bam->{count} );
        &print_tab( "$output_dir/$name\_post_processing.tab", \@header, \@vals );

        # Calculate BWA LCA's for LGTs
        print STDERR "========== LGT-BWA-LCA  ==========\n";
        $lgtseek->runBWA(
            {   input_bam      => $filtered_bam->{bam},
                output_dir     => "$output_dir\/lgt_bwa-lca/",
                out_file       => "$output_dir\/lgt_bwa-lca\/$name\_lgt_lca-bwa.txt",
                reference_list => $lgtseek->{refseq_list},
                threads        => $lgtseek->{threads},
                run_lca        => 1,
                overwrite      => $lgtseek->{overwrite},
                cleanup_sai    => 1,
            }
        );

        # Blast & get best hits
        print STDERR "========== LGT-BESTBLAST2 ==========\n";
        my $best_blasts = $lgtseek->bestBlast2(
            {   bam        => $filtered_bam->{bam},
                db         => $lgtseek->{path_to_blastdb},
                threads    => $lgtseek->{threads},
                lineage1   => $lgtseek->{donor_lineage},
                lineage2   => $lgtseek->{host_lineage},
                output_dir => "$output_dir/blast_validation/"
            }
        );

        print STDERR "======== LGT-Blast-LCA =========\n";
        my $lgt_blast_lca = $lgtseek->blast2lca(
            {   blast          => $best_blasts->{overall_blast},
                output_dir     => "$output_dir\/lgt_blast-lca/",
                evalue_cutoff  => $options{evalue_cutoff},         # Default = 1
                best_hits_only => 1,                               # Default = 0
            }
        );

        # LGTFinder
        print STDERR "======== LGT-FINDER =========\n";
        my $valid_lgts = $lgtseek->runLgtFinder(
            {   input_file_list => $best_blasts->{list_file},
                lineage1        => $lgtseek->{donor_lineage},
                lineage2        => $lgtseek->{host_lineage},
                output_prefix   => $name,
                max_overlap     => $lgtseek->{max_overlap},
                min_length      => $lgtseek->{min_length},
                output_dir      => "$output_dir/lgt_finder/",
            }
        );

        print STDERR "======== VALID-BLAST-PP =========\n";
        ## Create a new bam from blast validation & lgtfinder results and add #'s to post_processing.tab
        my $blast_validated_lgt = $lgtseek->validated_bam(
            {   input    => $filtered_bam->{bam},
                by_clone => $valid_lgts->{by_clone},
                output   => "$output_dir/$name\_lgt_host_filtered_validated.bam"
            }
        );

        ## Add numbers for validated-LGT to post_processing.tab
        push( @header, 'lgt_valid_blast' );
        push( @vals,   "$blast_validated_lgt->{count}" );
        &print_tab( "$output_dir/$name\_post_processing.tab", \@header, \@vals );

        print STDERR "======== LGT-ALN-HG19 =======\n";
        my $hg19_lgt_bam = $lgtseek->runBWA(
            {   input_bam   => $blast_validated_lgt->{bam},
                output_bam  => 1,
                threads     => $lgtseek->{threads},
                output_dir  => $output_dir,
                reference   => $lgtseek->{hg19_ref},
                overwrite   => $lgtseek->{overwrite},
                cleanup_sai => 1,
            }
        );

        # Fix the header for the final LGT bam.
        my $Picard        = "$lgtseek->{java_bin} \-$lgtseek->{java_opts} -jar $lgtseek->{Picard_jar}";
        my $cmd           = "$Picard AddCommentsToBam I=$hg19_lgt_bam->[0] O=$output_dir/$name\_lgt_final.bam";
        my $header_string = $lgtseek->_run_cmd("samtools view -H $blast_validated_lgt->{bam}");
        print STDERR "FOO1: $header_string\n";
        my @pg_list = grep( /^\@PG|^\@CO/, split( /\n/, $header_string ) );
        print STDERR "Foo2: @pg_list\n";
        map {
            if ( $_ !~ /ID:bwa/ ) {
                my @split_pg_line = split( /\t/, $_ );
                shift(@split_pg_line);    # remove BWA \@PG line
                my $comment_line = join( "\t", @split_pg_line );
                $cmd = $cmd . " C=\"$comment_line\"";
            }
        } @pg_list;
        $lgtseek->_run_cmd("$cmd");
        $lgtseek->_run_cmd("rm $hg19_lgt_bam->[0]");

        # Run blast and keep raw output
        #my $blast_ret = $lgtseek->_run_cmd("blastall -p blastn -a $lgtseek->{threads} -e 10e-5 -T F -d $lgtseek->{path_to_blastdb} -i $lgt_fasta > $output_dir/blast_validation/$name\_blast.raw");
    }

    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## Microbiome Analysis
    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## Calculate BWA LCA's for Microbiome Reads
    print STDERR "====== Microbiome ======\n";
    ## Check to make sure we found Microbiome Reads. If no microbiome reads skip this step.
    if ( $lgtseek->empty_chk( { input => "$output_dir\/$name\_microbiome.bam" } ) ) {
        print STDERR "***Warning*** No Microbiome reads in: $output_dir\/$name\_microbiome.bam. Skipping microbiome LCA calculation.\n";
        print_complete( \%options, "LGTSEQ_ANALYSIS_VERSION=$LGTSEQ_ANALYSIS" );
    }
    else {

        # Prinseq filter the putative microbiome reads
        # Use only prinseq quality reads for microbiome analysis
        print STDERR "===== Microbiome-PRINSEQ =====\n";
        my $filtered_bam = $lgtseek->prinseqFilterBam(
            {   input_bam  => $pp_data->{files}->{microbiome_donor},
                output_dir => $output_dir,
            }
        );

        # Add filtered count to counts.
        push( @header, 'microbiome_pass_prinseq' );
        push( @vals,   $filtered_bam->{count} );
        &print_tab( "$output_dir/$name\_post_processing.tab", \@header, \@vals );

        print STDERR "==== Microbiome-BWA-LCA ====\n";
        $lgtseek->runBWA(
            {   input_bam      => $filtered_bam->{bam},
                output_dir     => "$output_dir\/microbiome_bwa-lca/",
                out_file       => "$output_dir\/microbiome_bwa-lca/$name\_microbiome_lca-bwa.txt",    ## KBS 01.07.14
                reference_list => $lgtseek->{refseq_list},
                threads        => $lgtseek->{threads},
                run_lca        => 1,
                overwrite      => $lgtseek->{overwrite},
                cleanup_sai    => 1,
            }
        );

        ## Cleanup Intermediate Files
        print STDERR "RM: $output_dir/microbiome_prinseq_filtering/\n";
        $lgtseek->_run_cmd("rm -rf $output_dir/microbiome_prinseq_filtering/");

    }

    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## Calculate coverage of LGT on human side. (Not apropriate most the time)
    ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( $lgtseek->{lgt_coverage} == 1 ) {
        print STDERR "====== Calculating Coverage of Hg19 LGT ======\n";
        $lgtseek->mpileup(
            {   input      => "$output_dir\/$name\_lgt_host_filtered.bam",
                output_dir => $output_dir,
                ref        => $lgtseek->{hg19_ref},
                cleanup    => 1,
                overwrite  => $lgtseek->{overwrite},
            }
        );

    }

    $lgtseek->time_check;
    print_complete( \%options, "LGTSEQ_ANALYSIS_VERSION=$LGTSEQ_ANALYSIS" );
}

## Subroutines
sub print_tab {
    my ( $file, $header, $vals ) = @_;
    open OUT, ">$file" or $lgtseek->fail("Couldn't open $file\n");
    print OUT join( "\t", @$header );
    print OUT "\n";
    print OUT join( "\t", @$vals );
    print OUT "\n";
}

sub help {
    die "Help: This script will identify bacterial human LGT.
    --input=                <Input> Accepts bams, fastq, and fastq.gz. With fatsq's only use 1 of the pair for input. (ie: foo_1.fastq.gz)
    --output_dir=           Directory for all output. Will be created if it doesn't exist. 
    --help_full|?           Full help info on options.
    --workflow_help         Help setting up an efficient lgtseq workflow with the optional steps.\n";
}

sub help_full {
    die "   LGTSEQ_ANALYSIS_VERSION    $LGTSEQ_ANALYSIS
    Help: This script takes a bam and identifies bacterial human LGT.
         _____
    ____/Input\\_________________________________________________________________________________
    --input|i=              <Input BAM or fastq> **MUST** be name sorted. Also use --prelim_filter & --name_sort_input for position sorted inputs.
    --input_list|I=         <List of inputs> 1 per line.
         ______
    ____/Output\\________________________________________________________________________________
    --output_dir|o=         Directory for all output. Example: /path/to/{output_dir}/{tcga_dirs}/{subdirs}/ || /path/to/{output_dir}/{subdirs}/
     --tcga_dirs=           <0|1> [0] 1= Make the sub-dir prefix the input's last folder in path (Maintain TCGA analysis_id directory structure)
      --subdirs=            <0|1> [0] 1= Make the sub-dir prefix in output_dir based on input name.
    --overwrite=            <0|1> [0] 1= Overwrite previously data. 
         _______________________
    ____/Primary Human Alignment\\_______________________________________________________________
    --aln1_human=           <0|1> [0] 1= Primary aln to hg19. MUST be used if input=fastq's
         _____________________
    ____/Preliminary Filtering\\__________________________________________________________________
    --prelim_filter=        <0|1> [0] 1= Filter a human mapped bam, keeping potential LGT & Microbiome reads.
    --name_sort_input=      <0|1> [0] 1= Resort the input bam by read names.
    --keep_softclip=        <0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT)
    --sort_mem=             [1G] Mem per thread to sort with. Careful this corresponds with --threads. 
         _________________________
    ____/Secondary Human Alignment\\______________________________________________________________
    --aln2_human=           <0|1> [0] 1= Secondary aln to hg19. Useful for mapping to standarized hg19 after prelim_filter.
         __________________
    ____/Submit to SGE-grid\\_____________________________________________________________________
    --Qsub|q=               <0|1> [0] 1= qsub the job to the grid.
    --Qsub_iterate=         <0|1> [0] 1= Iterate over --input_list instead of qsub'ing each file from a list. 
      --threads=|t          [1] # of CPU's to use for multithreading BWA sampe
      --sub_mem|mf=         [4G] Min mem to qsub for on grid
      --sub_name=           [lgtseq] 
      --sub_mail=           [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
      --project=            Grid project to use. 
      --print_hostname|ph=  <0|1> [0] Print hostname. Defaults to 1 if verbose=1.
         ________________
    ____/Help Information\\_______________________________________________________________________
    --verbose|V             <0|1> [0] 1= Verbose reporting of progress. 0 =Turns off reports. 
    --help|h                Help Basic Info
    --help_full|?           Help Full Info
    --workflow_help         Examples how to use optional portions of lgt_seq to increase efficiency.
    --conf_file=            [~/.lgtseek.conf]
    --conf_help             Help Information on lgtseek.conf requirements.
    _____________________________________________________________________________________________\n";
}

sub workflow_help {
    die "----------------------------------------------------------------------------------------
    Workflow of modules:
    input    -> Primary aln @ hg19 -> Prelim_filter     -> Secondary aln @ hg19 -> map_bacteria -> LGTFINDER
    Options: -> aln1_human=0|1     -> prelim_filter=0|1 -> aln2_human=0|1       -> Always 1     -> Always 1
    Ex1 - Input bacteria mapped:                --aln1_human=1 --prelim_filter=1 --aln2_human=0
    Ex2 - Input human mapped 2 non-hg19:        --aln1_human=0 --prelim_filter=1 --aln2_human=1
    Ex3 - Input hg19 position sorted            --aln1_human=0 --prelim_filter=1 --aln2_human=0
    Ex4 - Post prelim-filter.pl                 --aln1_human=0 --prelim_filter=0 --aln2_human=0
    ----------------------------------------------------------------------------------------\n";
}

## This section isn't complete yet and not check to make sure these param's ARE in the .conf file
sub conf_help {
    die "== The .lgtseek.conf file options. Options with \"\*\" are mandatory. This list is NOT an exhaustive list of mandatory options (yet).
== The format of the ~/.lgtseek.conf is \"tab\" (white space) delimited \$option_name \\t \$option_value. ex: threads    4
             __________
        ____/References\\_____________________________________________________________________________
        hg19_ref*            Path to hg19 reference
        split_bac_list*      Path to the list of split bacterial references (A2D, E2P, R2Z)
        refseq_list*         Path to all bacterial references in refseq. 
             ____________
        ____/SGE-Options\\____________________________________________________________________________ 
        project*             SGE Project name
        sub_mem            
        sub_name
        threads
             ______________
        ____/Bin Directories\\________________________________________________________________________
        bin_dir*              [/local/projects-t3/HLGT/scripts/lgtseek/bin/]
        taxon_host*           [mongotest1-lx.igs.umaryland.edu:10001]
        taxon_dir*            [/local/db/repository/ncbi/blast/20120414_001321/taxonomy/]
        taxon_idx_dir*        [/local/projects-t3/HLGT/idx_dir/20120414]
        path_to_blastdb*      [/local/db/repository/ncbi/blast/20120414_001321/nt/nt]\n";
}

