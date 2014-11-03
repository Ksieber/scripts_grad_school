#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
use run_cmd;
use setup_input;
use setup_output;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions (\%options,
    'input=s',
    'input_list=s',
    'output_dir=s',
    'subdirs=s',
    'gpg_key=s',
    'recipient=s',
    'test_print=s',
    'help',
);

my $key = $options{gpg_key} ? $options{gpg_key} : "Portable";
my $recipient = $options{recipient} ? $options{recipient} : "Portable";
if(!$options{subdirs}){$options{subdirs}=0};
if(!$options{test_print}){$options{test_print}=0};

if($options{help}){die
"Help: This script will encrypt a file(s) with gpg. 
    --input=
    --input_list=
    --gpg_key=          GPG Key to be used for encryption
    --recipient=        GPG recipient user
    --output_dir=		[cwd]
    --subdirs=          <0|1> [0] If 1, a new sub-directory will be made in the output_dir for each output file.
    --test_print=       <0|1> [0] If 1, the encryption command will be printed to STDERR instead of being executed. 
    --help\n";
}

if(!$options{input} && !$options{input_list}){die
	"Error: Must give an input. use --input= or --input_list=.\n";
}

my $input=setup_input();           ## Returns array(ref) of input(s)
my $out_dir=setup_output($input);   ## Returns hash(ref) of output dirs: input=>output_dir;
foreach my $file (@$input){
    my($fn,$path)=fileparse($file);
    my $cmd = "gpg -o $out_dir->{$file}/$fn\.gpg --default-key $key -r $recipient -e $file";
    if($options{test_print}==1){
        print STDERR "$cmd\n";
        next;
    } else { 
        run_cmd($cmd);
    }
}

