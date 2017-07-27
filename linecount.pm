package linecount;
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use run_cmd;
use strict;
use warnings;
use File::Basename;
use Cwd;
use POSIX;
use File::Util;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( wc );

sub wc {
    my $file = shift;
    my $count;
    if ( $file =~ /.*\.bam$/ ) {
        my $indexed = 0;

        my $FileUtil           = File::Util->new();
        my $file_modified_last = scalar localtime $FileUtil->last_modified("$file");
        $file_modified_last =~ s/\s+/_/g;

        my ( $bam_fn, $path, $suffix ) = fileparse( $file, ".bam" );
        $path =~ s/(\/)$//;
        if ( $path =~ /^\.{1}$/ ) {
            my $cwd = getcwd;
            if ( -e "$cwd\/$file" ) {
                $path = $cwd;
                $file = "$cwd\/$file";
            }
        }
        foreach my $header_line ( split( /\n/, `samtools view -H $file` ) ) {
            next if ( $header_line !~ /^\@CO\s+LineCount/ );
            if ( $header_line =~ /^\@CO\s+LineCount:(\d+)\s+BamName:.*$bam_fn.*\.bam\s+TimeStamp:$file_modified_last/ ) {
                $count   = $1;
                $indexed = 1;
            }
        }
        if ( $indexed != 1 ) {
            $count = `samtools view $file | wc -l`;
            $count =~ /^(\d+)\s+/;
            my $wc             = $1;
            my $random_int     = ceil( rand(1) * 1000000 );
            my $raw_time_stamp = localtime();
            my $time_stamp     = $raw_time_stamp;
            $time_stamp =~ s/\s+/_/g;

            run_cmd("cp $file $path\/tmp_bam_for_counting_$random_int\.bam");

            # Remove potentially outdated linecount header_line and print the new one
            open( HEADER, ">", "$path\/tmp_$random_int\.header" ) or die "Error: Unable to open the output header: $path\/tmp_$random_int\.header because: $!\n";
            my $header = run_cmd("samtools view -H $file");
            my @split_header = split( /\n/, $header );
            foreach my $header_line (@split_header) { next if ( $header_line =~ /^\@CO\s+Reads:/ ); print HEADER "$header_line\n" if ( $header_line !~ /^\@CO\s+LineCount:/ ); }
            print HEADER "\@CO\tLineCount:$wc\tBamName:$file\tTimeStamp:$time_stamp";
            close HEADER;

            # Adjust the new header
            run_cmd("samtools reheader -P $path\/tmp_$random_int\.header $path\/tmp_bam_for_counting_$random_int\.bam > $file");

            # Cleanup intermediate files
            run_cmd("rm $path\/tmp_bam_for_counting_$random_int.bam $path\/tmp_$random_int.header");

            # Adjust the bam's modify time to match the $time_stamp
            run_cmd("touch -d \"$raw_time_stamp\" $file");
        }
    }
    elsif ( $file =~ /.*\.sam$/ ) {
        $count = `samtools view -S $file | wc -l`;
    }
    elsif ( $file =~ /.*\.gz$/ ) {
        $count = `zcat $file | wc -l`;
    }
    else {
        $count = `wc -l $file`;
    }
    $count =~ /^(\d+)\s*/;
    my $num = $1;
    if ( $num == 0 ) {
        print STDERR "* Warning * | This file is empty: $file\n";
    }
    return "$num";
}
1;

__END__
