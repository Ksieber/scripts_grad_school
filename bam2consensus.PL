#!/usr/bin/perl
use lib ( 
  '/home/ksieber/perl5/lib/perl5/', 
  '/home/ksieber/scripts/', 
  '/local/projects-t3/HLGT/scripts/lgtseek/lib/',
  '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/'
  );
use warnings;
use strict;
use File::Basename;
use print_call;
use read_in_list;
use run_cmd;
use read_bam;

if ( !@ARGV and -t STDIN ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,    'input|i=s',      'sort=i',            'reference|ref|r=s', 'region=s', 'region_list=s', 'smart_regions=i', 'read_list=s',
    'output|O=s', 'output_dir|o=s', 'output_prefix|p=s', 'cleanup=i',         'help|?',   'more_info',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }
if ( $options{more_info} ) { &more_info; }

if ( !$options{input} )        { die "\n\tError: Must pass --input=<position sorted bam>. Please try again\n\n"; }
if ( !-e $options{input} )     { die "\n\tError: The input file doesn't exist: $options{input}\n\n"; }
if ( !$options{reference} )    { die "\n\tError: Must pass --reference=<ref.fa> Please try again\n\n"; }
if ( !-e $options{reference} ) { die "\n\tError: The reference file doesn't exist: $options{reference}\n\n"; }
if ( ( !$options{region} and !$options{region_list} ) and !$options{read_list} ) { die "\n\tError: Must pass a region or read list. See --help for more info.\n\n"; }
if ( defined $options{read_list}   and !-e $options{read_list} )   { die "\n\tError: The --read_list file doesn't exist: $options{read_list}\n\n"; }
if ( defined $options{region_list} and !-e $options{region_list} ) { die "\n\tError: The --region_list file doesn't exist: $options{region_list}\n\n"; }

my $bam = $options{input};
my ( $in_fn, $in_dir, $in_suf ) = fileparse( $bam, qr/\.[^\.]+/ );
my $sort_check = run_cmd("/home/ksieber/bin/samtools view -H $bam | head -n 1");
if ( $sort_check !~ /^\@HD.+SO\:coordinate/ ) { print STDERR "Warning: This bam does not look position sorted. Sorting it now ...\n"; $options{sort} = 1; }
if ( defined $options{sort} and $options{sort} == 1 ) {
    run_cmd("/home/ksieber/bin/samtools sort $bam $in_dir$in_fn\.psort");
    run_cmd("/home/ksieber/bin/samtools index $in_dir$in_fn\.psort.bam");
    $bam = "$in_dir$in_fn\.psort.bam";
}
my %ids;
if ( defined $options{read_list} ) {
    open( IN, "<", "$options{read_list}" ) or die "Error: Unable to open read_list: $options{read_list}\n";
    while (<IN>) {
        chomp( my $line = $_ );
        $ids{$line}++;
    }
    close IN;
}
my $out_dir       = defined $options{output_dir}    ? $options{output_dir}    : $in_dir;
my $out_fn        = defined $options{output_prefix} ? $options{output_prefix} : $in_fn;
my $smart_regions = defined $options{smart_regions} ? $options{smart_regions} : 1;
my $cleanup       = defined $options{cleanup}       ? $options{cleanup}       : 1;
my $reference     = $options{reference};
my $out_pref      = "$out_dir$out_fn";

if ( defined $options{region} or defined $options{region_list} ) {
    my @region_array;
    if ( defined $options{region} ) {
        push( @region_array, $options{region} );
    }
    if ( defined $options{region_list} ) {
        open( IN, "<", "$options{region_list}" ) or die "Error: Unable to open the region_list: $options{region_list}\n";
        while (<IN>) {
            chomp( my $line = $_ );
            next if ( $line =~ /^[\#]+/ );
            my @split = split( /\s+/, $line );
            push( @region_array, $split[0] );
        }
        close IN;
    }

    foreach my $region (@region_array) {
        my $clean_region = $region;
        $clean_region =~ s/([\|\:])+/__/g;
        if ( $smart_regions == 1 ) {
            ## Open input FH
            my ( $header, $BAM_FH ) = open_bam( $bam, $region );
            ## Open output FH
            open(
                my $OUT_FH,
                "|-",
                "/home/ksieber/bin/samtools view - -u | /home/ksieber/bin/samtools sort -O bam -T tmp_bam2_sort - | /home/ksieber/bin/samtools mpileup -uAf $reference - 2>/dev/null | /home/ksieber/bin/bcftools call -m -O z > $out_pref\_$clean_region.vcf.gz"
            );
            print $OUT_FH "$header\n";

            my $bam_chr;
            my $bam_min;
            my $bam_max;

            while ( my $read = read_bam($BAM_FH) ) {
                next if ( $read->{flag}->{qunmapped} );
                next if ( defined $options{read_list} and !defined $ids{ $read->{id} } );

                my $read_5_pos = $read->{position};
                my $read_3_pos = $read->{position} + length( $read->{sequence} ) - 1;

                if ( !defined $bam_chr ) { $bam_chr = $read->{chr}; }
                if ( !defined $bam_min or ( $read_5_pos <= $bam_min ) ) { $bam_min = $read_5_pos; }
                if ( !defined $bam_max or ( $read_3_pos >= $bam_max ) ) { $bam_max = $read_3_pos; }

                print $OUT_FH "$read->{line}\n";
            }

            close $OUT_FH;
            close $BAM_FH;

            my $actual_region       = "$bam_chr\:$bam_min\-$bam_max";
            my $actual_clean_region = $actual_region;
            $actual_clean_region =~ s/([\|\:])+/__/g;

            run_cmd("/home/ksieber/bin/tabix $out_pref\_$clean_region.vcf.gz");
            if ( defined $options{output_dir} or defined $options{output_prefix} ) {
                run_cmd("/home/ksieber/bin/samtools faidx $reference \'$actual_region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz > $out_pref\_$actual_clean_region\.fa");
            }
            elsif ( defined $options{output} ) {
                run_cmd("/home/ksieber/bin/samtools faidx $reference \'$actual_region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz >> $options{output}");
            }
            else {
                my $consensus_sequence = run_cmd("/home/ksieber/bin/samtools faidx $reference \'$actual_region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz");
                print STDOUT "$consensus_sequence\n";
            }
        }
        else {
            run_cmd(
                "/home/ksieber/bin/samtools view -u $bam \'$region\' | /home/ksieber/bin/samtools mpileup -uAf $reference - 2>/dev/null | /home/ksieber/bin/bcftools call -m -O z > $out_pref\_$clean_region.vcf.gz"
            );
            run_cmd("/home/ksieber/bin/tabix $out_pref\_$clean_region.vcf.gz");
            if ( defined $options{output_dir} or defined $options{output_prefix} ) {
                run_cmd("/home/ksieber/bin/samtools faidx $reference \'$region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz > $out_pref\_$clean_region\.fa");
            }
            elsif ( defined $options{output} ) {
                run_cmd("/home/ksieber/bin/samtools faidx $reference \'$region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz >> $options{output}");
            }
            else {
                my $consensus_sequence = run_cmd("/home/ksieber/bin/samtools faidx $reference \'$region\' | /home/ksieber/bin/bcftools consensus $out_pref\_$clean_region.vcf.gz");
                print STDOUT "$consensus_sequence\n";
            }
        }
        if ( $cleanup == 1 ) { run_cmd("rm $out_pref\_$clean_region.vcf.gz $out_pref\_$clean_region.vcf.gz.tbi"); }
    }
}
elsif ( defined $options{read_list} ) {
    my ( $header, $BAM_FH ) = open_bam($bam);
    open(
        my $OUT_FH,
        "|-",
        "/home/ksieber/bin/samtools view - -u | /home/ksieber/bin/samtools sort -O bam -T tmp_bam2_sort - | /home/ksieber/bin/samtools mpileup -uAf $reference - 2>/dev/null | /home/ksieber/bin/bcftools call -m -O z > $out_pref\_tmp.vcf.gz"
    );
    print $OUT_FH "$header\n";

    my $bam_chr;
    my $bam_min;
    my $bam_max;

    while ( my $read = read_bam($BAM_FH) ) {
        next if ( $read->{flag}->{qunmapped} );
        next if ( !defined $ids{ $read->{id} } );
        next if ( defined $bam_chr and ( $read->{chr} ne $bam_chr ) );

        my $read_5_pos = $read->{position};
        my $read_3_pos = $read->{position} + length( $read->{sequence} ) - 1;

        if ( !defined $bam_chr ) { $bam_chr = $read->{chr}; }
        if ( !defined $bam_min or ( $read_5_pos < $bam_min ) ) { $bam_min = $read_5_pos; }
        if ( !defined $bam_max or ( $read_3_pos > $bam_max ) ) { $bam_max = $read_3_pos; }

        print $OUT_FH "$read->{line}\n";
    }

    close $OUT_FH;
    close $BAM_FH;

    run_cmd("/home/ksieber/bin/tabix $out_pref\_tmp.vcf.gz");
    if ( defined $options{output_dir} or defined $options{output_prefix} ) {
        run_cmd("/home/ksieber/bin/samtools faidx $reference \'$bam_chr\:$bam_min\-$bam_max\' | /home/ksieber/bin/bcftools consensus $out_pref\_tmp.vcf.gz > $out_pref\.fa");
    }
    elsif ( defined $options{output} ) {
        run_cmd("/home/ksieber/bin/samtools faidx $reference \'$bam_chr\:$bam_min\-$bam_max\' | /home/ksieber/bin/bcftools consensus $out_pref\_tmp.vcf.gz >> $options{output}");
    }
    else {
        my $consensus_sequence = run_cmd("/home/ksieber/bin/samtools faidx $reference \'$bam_chr\:$bam_min\-$bam_max\' | /home/ksieber/bin/bcftools consensus $out_pref\_tmp.vcf.gz");
        print STDOUT "$consensus_sequence\n";
    }
    if ( $cleanup == 1 ) { run_cmd("rm $out_pref\_tmp.vcf.gz $out_pref\_tmp.vcf.gz.tbi"); }
}
else {
    die "Error: Unsure how to process the data. Please use --help for proper use.\n";
}
if ( $cleanup == 1 and ( defined $options{sort} and $options{sort} == 1 ) ) { run_cmd("rm $in_dir$in_fn\.psort.bam $in_dir$in_fn\.psort.bam.bai"); }
if ( defined $options{output} or defined $options{output_dir} or defined $options{output_prefix} ) { print_complete( \%options ); }

sub help {
    die "This script will create a consensus sequence for a subset of reads from a bam.
    --input=                < Position sorted bam >  Create consensus from reads aligned to reference.
      --sort=               <0|1> [0] 1= Position sort the input bam.
    --reference|ref|r=      < /path/to/ref.fa > That the bam is aligned against. This will be used generate the consensus sequence.
    --region=               < /home/ksieber/bin/samtools style region > (chr:1-10) region to pull reads from to create the consensus sequence.
    --region_list=          List of regions, 1 per line. Parses for first column, other data in alt coloumns ok. 
      --smart_regions=      <0|1> [1] 1= *May be different than --region(s). **See --more_info**
    --reads_list=           < /path/to/list.txt > List of read-ids to create a consensus sequence from. 
                            All --read_list reads must be from the same chromosome, otherwise only the first chr will be analyzed.
    --output=               < /full/path/for/output.fa >    [Standard Output]
    --output_dir=           < /directory/for/output/ >      [/Input-bam/dir/]
    --output_prefix=        < {output_prefix}.fa >          [Input-bam-name_region] 
    --cleanup=              <0|1> [1] 1= Delete intermediate .psort.bam and .vcf files. 0=Keeps files.\n";
}

sub more_info {
    die
        "\n  --smart_regions=1 Will pull all the reads that map to the --region and create a new \"smart region\" based on the coordinates of the reads aligning to the --region.\n  --smart_regions=0 Will use ONLY the --region coordinates.\n\n";
}
