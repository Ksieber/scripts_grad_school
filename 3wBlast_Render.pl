#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Bio::Index::Blast;
use Bio::SearchIO::Writer::TextResultWriter;


my $file = shift or die "Usage: render4.pl <blast.filei> <ReadID>\n";
my $id = shift or die; 

my $indexed_file = "$file.idx";
my $indexed = -e $indexed_file;

my $index = Bio::Index::Blast->new(-filename => $indexed_file,
                                   -write_flag => 1,
                                   -verbose => 1);
unless ($indexed) {
   print STDERR "Creating .idx\n";
                              #$index->id_parser(\&get_id);
                              # sub get_id {
                              # my $line = shift;
                              #      $line =~ /Query= (\S+)/;
                              ##      print $1;
                              #      $1;
                              #  }
   $index->make_index($file);
}

my $fh = $index->get_stream($id);
my $searchio = Bio::SearchIO->new(-noclose => 1,
                                  -format  => 'blast',
                                  -fh      => $fh); 

my $result = $searchio->next_result() or die "No Result.\n";

my $panel = Bio::Graphics::Panel->new(
        -length  => $result->query_length,
         -width   => 800,
         -pad_left   => 10,
         -pad_right  => 10,
      );

my $full_length = Bio::SeqFeature::Generic->new(
         -start   => 1,
         -end           => $result->query_length,
         -display_name  => $result->query_name,
      );

$panel->add_track($full_length,
         -glyph   => 'arrow',
         -tick    => 2,
         -fgcolor => 'orange',
         -double  => 1,
         -label   => 1,
      );

my $track = $panel->add_track(
         -glyph   => 'graded_segments',
         -label   => 1,
         -connector  => 'dashed',
         -bgcolor => 'blue',
         -font2color => 'orange',
         -sort_order => 'high_score',
         -description => sub {
            my $feature = shift;
            return unless $feature->has_tag('description');
            my ($description) = $feature->each_tag_value('description');
            my $score = $feature->score;
            "$description, score=$score";
            },
      );

while (my $hit = $result->next_hit) {
 ##next unless $hit->significance < 1E-20;
   my $feature = Bio::SeqFeature::Generic->new(
         -display_name  => $hit->name,
         -score         => $hit->raw_score,
         -tag           => {
                  description => $hit->description
                  },
         );
   print STDERR "$hit->name\n";
   while(my $hsp = $hit->next_hsp ) {
      $feature->add_sub_SeqFeature($hsp,'EXPAND');
   }
   $track->add_feature($feature);
}

my ($url,$map,$mapname)= 
   $panel->image_and_map(-root => '/',
                         -url => '/home/ksieber/tmp/images/',
                         -link => '#$name');

my $stream = $index->get_stream($id);
open (WH, ">", "$file.html");
print WH "<html>\n<body>\n<pre>";
my $last=0;
while (<$stream>) {
   if (/^.?BLAST/) {
      $last++;
      if ($last >1) {last};
   }
   my $line = $_;
   if ($line =~ m/^\>(\S+)/) {
     print WH "<a name=$1> $line</a>\n";
   }
   else {print WH $line;}
   if ($line =~ m/Searching.*done/) {
      print WH "</pre>\n"; 
      print WH "$map\n";
      print WH  qq(<img src="$url" usemap="#$mapname">),"\n";
      print WH "<pre>\n";
   }
} 
print WH "</pre>\n</body>\n</html>\n";

close $stream || die "error\n";
close WH || die "error\n";


