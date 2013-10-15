package run_cmd;
#use warnings;
use strict;
use File::Basename;
use Carp;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( run_cmd setup_logs Qsub Qsub2);
## &run_cmd records a unix cmd, executes it, checks it did not fail
## &run_cmd can also take a fh to log the commands to the fh
## &setup_logs takes a hash ref. of input=>output_dirs and returns a hash reference (input=>fh);
## if &setup_logs isn't run, &run_cmd will log the cmd to STDERR
## &Qsub will take a shell script || command and qsub it to the grid
## &Qsub prints qsub info to STDERR 
## Suggested use: my $cmd = "perl foo.pl";
## Suggested use: my $cmd_log=setup_logs($out_dirs);
## Suggested use: my $log = $cmd_log->{$file};
## Suggested use: Qsub("echo foo",$log);
## Suggested use: run_cmd("bwa aln foo.bam ...");
## Suggested use: run_cmd($cmd,$cmd_log->{$file});

=head2 run_cmd
	Title		: run_cmd
	Usage		: run_cmd($cmd,$log);
	Function	: qsub command with good error logging if it fails
	Returns 	: return value for executing unix command
=cut 

sub run_cmd {
	my ($cmd,$log)=@_;
	my $fh;
	if($log){
		open($fh,">>","$log") or die "Can't open: $fh because: $!\n";
	} else { 
		$fh = *STDERR;
	}
	
	print $fh "$cmd\n";
	chomp(my $res = `$cmd`);
    if($?){
        print $fh "$?\n";
        confess "Error!\n";
        die "$cmd died with message:\n$res\n\n";
    }
    return $res;
}

=head2 setup_logs
	Title		: setup_logs
	Usage		: setup_logs($output_dirs)	
				  ## $output_dirs is a hash_ref: $hash{$file}=$output_directory
	Function	: setup log name based on file and output dirs.
	Returns 	: Hash_ref for log file based on file name
					$return_value : $hash{$file}=$log_file

=cut 

sub setup_logs {
	my $out_dirs = shift;
	my %logs;
	foreach my $file (keys %$out_dirs){
		my ($fn,$path,$suf)=fileparse($file);
		$logs{$file}="$out_dirs->{$file}$fn\.log";
	}
	return \%logs;
}

=head2 Qsub
	Title		: Qsub
	Usage		: Qsub($cmd,$log);
	Function	: qsub commands
=cut 

sub Qsub {
	my ($cmd,$log) = @_;
	my $fh;
   	if($log){
   		open($fh,">>","$log") or die "Can't open: $fh because: $!\n";
	} else { 
		$fh = *STDERR;
	}
	print $fh "Going to qsub: $cmd\n";
	my $qsub;
	if($cmd=~/\.sh$/){
		$qsub="qsub -V -P jdhotopp-lab $cmd";
	} else {
		$qsub="echo \"$cmd\" | qsub -V -P jdhotopp-lab";
   	}
   	chomp(my $report = `$qsub`);
   	if($?){
   		print $fh "$?\n";
   		confess "Error!\n";
   		die "$qsub died with message: $report\n";
   	}
   	print $fh "$report\n";
   	return $report;
 }

=head2 Qsub2 
	Title		: Qsub2
	Usage		: Qsub2({ cmd => $cmd, threads => $threads });
	Function	: More verbose version of qsub.
	Args	:
		cmd 	 	=> command/shell script to be qsub'ed
		log 	 	=> file to log the qsub command in.
		name	 	=> Job Name.
		threads	 	=> # cpu threads to use.
		mem		 	=> Memory free @ start.
	 	project	 	=> grid project to use.
		cwd		 	=> <0|1> [0] 1= Use current working directory to work from.
		wd 		 	=> Directory for grid to work from.
		hostname 	=> Hostname to run on. ex: => "grid*"
		no_gal	 	=> <0|1> [0] 1= Don't run on galactus. ie hostname="magneto|juggernaut|grid*"
		excl		=> <0|1> [0] 1= Run exclusively on a node.
=cut 

sub Qsub2 {
 	my $config=$_[0];
 	if(!$config->{cmd}){die "Must use cmd => <Command to qsub>. Please Fix.\n";}
 	my $cmd = $config->{cmd};
 	my $log = $config->{log};
 	my $wd = defined $config->{wd} ? " -wd $config->{wd}" : undef;
 	my $threads=undef;
 	my $t = defined $config->{threads} ? $config->{threads} : 1;
 	if($t>=2){$threads = " -pe thread $t -q threaded.q";}
 	my $mem = defined $config->{mem} ? " -l mf=$config->{mem}" : undef;
 	my $project = defined $config->{project} ? " -P $config->{project}" : " -P jdhotopp-lab";
 	my $cwd = defined $config->{cwd} ? " -cwd" : undef;
    my $name = defined $config->{name} ? " -N $config->{name}" : undef;
    my $hostname = defined $config->{hostname} ? " -l hostname=\"$config->{hostname}\"" : undef;
    if($config->{no_gal}==1){$hostname= " -l hostname=\"magneto|juggernaut|grid*\"";}
    my $exclusive;
    if($config->{excl}==1){
    	$exclusive = " -l excl=true";
    } else {
    	$exclusive = undef;
    }
    print STDERR "exclusive\t$config->{excl},\t$exclusive.\n";
 	my $fh;
 	if($log){
   		open($fh,">>","$log") or die "Can't open: $fh because: $!\n";
	} else { 
		$fh = *STDERR;
	}
	my $qsub = "qsub -V$name$project$wd$cwd$mem$threads$hostname$exclusive";
	print $fh "QSUB: echo \"$cmd\" | $qsub\n";
	my $que;
	if($cmd=~/\.sh$/){
		$que = "$qsub $cmd";
	} else {
		$que = "echo \"$cmd\" | $qsub";
	}
	chomp(my $report = `$que`);
   	if($?){
   		print $fh "$?\n";
   		confess "Error!\n";
   		die "$que died with message: $report\n";
   	}
   	print $fh "QSUB: $report\n";
   	sleep 2;
   	return $report;
 }
  

1;
