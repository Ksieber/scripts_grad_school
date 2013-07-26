package setup_output;
use strict;
use warnings;
use File::Basename;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( setup_output );
## Returns hash(ref) of output dirs. Key = input, val = output_dir
## Suggested use: my $out_dir=setup_output($input);
## Looks for main script's $options{output_dir}
## Also looks for main $options{test_print} &| $options{subdirs} 

sub setup_output {
    my %outdir;
    our %options;
    our $results;
    if($main::options{output_dir} && !$main::options{test_print}){
        `mkdir -p $main::options{output_dir}`;
    }
    
    my $in = shift;
    if(!$in){die "Must pass an ref array of input into setup_output.pm\n";}
    foreach my $file (@$in){
        my ($fn,$path,$suf)=fileparse($file);
        if($main::options{subdirs}){
            $outdir{$file}=$main::options{output_dir} ? "$main::options{output_dir}/$fn/" : "$path/$fn/";
            if(!$main::options{test_print}){
               `mkdir -p $outdir{$file}`; 
            }
        } else {
            $outdir{$file}=$main::options{output_dir} ? $main::options{output_dir} : $path; 
        }   
    }
    return \%outdir;
}
1;
