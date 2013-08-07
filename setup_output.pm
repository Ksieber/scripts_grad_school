package setup_output;
use strict;
use File::Basename;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( setup_output );
## Returns hash(ref) of output dirs. Key = input, val = output_dir
## Suggested use: my $out_dir=setup_output($input);
## Looks for main script's $options{output_dir}
## Also looks for main $options{subdirs} 

sub setup_output {
    my %outdir;
    our %options;
    our $results;
    my $in = shift;
    if(!$in){die "Must pass an ref array of input into setup_output.pm\n";}
    foreach my $file (@$in){
        my ($fn,$path,$suf)=fileparse($file);
        if($main::options{subdirs}==1){
            $outdir{$file}=$main::options{output_dir} ? "$main::options{output_dir}/$fn/" : "$path/$fn/";
        } else {
            $outdir{$file}=$main::options{output_dir} ? $main::options{output_dir} : $path; 
        }   
        `mkdir -p $outdir{$file}`;
    }
    return \%outdir;
}
1;
