package run_cmd;
use strict;
use File::Basename;
use warnings;
use Carp;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( run_cmd setup_logs );
## run_cmd records a unix cmd, executes it, checks it did not fail
## run_cmd can also take a fh to print the commands to the fh
## setup_logs takes a hash ref. of input=>output_dirs and returns a hash reference (input==>fh);
## if setup_logs isn't run, run_cmd will print the cmd to STDERR
## Suggested use: my $cmd_log_fhs=setup_logs($out_dirs);
## Suggested use: run_cmd("bwa aln foo.bam ...");
## Suggested use: run_cmd("bwa aln foo.bam",$cmd_log_fhs->{$file});

sub run_cmd {
	my ($cmd,$log)=@_;
	my $fh;
	if($log){
		open($fh,">>","$log") or die "Can't open: $fh because: $!\n";
	} else { 
		$fh = *STDERR;
	}
	
	print $fh "$cmd\n";
	my $res = `$cmd`;
    if($?){
        print $fh "$?\n";
        confess "Error!\n";
        die "$cmd died with message:\n$res\n\n";
    }
    return $res;
}

sub setup_logs {
	my $out_dirs = shift;
	my %logs;
	foreach my $file (keys %$out_dirs){
		my ($fn,$path,$suf)=fileparse($file);
		$logs{$file}="$out_dirs->{$file}$fn\.log";
	}
	return \%logs;
}
1;