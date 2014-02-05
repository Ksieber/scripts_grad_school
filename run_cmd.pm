package run_cmd;
#use warnings;
use strict;
use File::Basename;
use Carp;
$Carp::MaxArgLen = 0; ## Report full length error 
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( run_cmd setup_logs Qsub Qsub2 Qsub3);
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
		open($fh,">>","$log") or confess "Can't open: $fh because: $!\n";
	} else { 
		$fh = *STDERR;
	}
	
	print $fh "CMD: $cmd\n";
	chomp(my $res = `$cmd`);
    if($?){
        print $fh "$?\n";
        confess "FAIL_CMD: $cmd died with message:\n$res\n\n";
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
   		open($fh,">>","$log") or confess "Can't open: $fh because: $!\n";
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
   		confess "Error: $qsub died with message: $report\n";
   	}
   	print $fh "$report\n";
   	return $report;
 }

=head2 Qsub2 
	Title		: Qsub2
	Usage		: Qsub2({ cmd => $cmd, threads => $threads });
	Function	: qsub but with greater control than &Qsub
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
 	if(!$config->{cmd}){confess "Must use cmd => <Command to qsub>. Please Fix.\n";}
 	my $cmd = $config->{cmd};
 	my $log = $config->{log};
 	my $wd = defined $config->{wd} ? " -wd $config->{wd}" : undef;
 	my $threads=undef;
 	my $t = defined $config->{threads} ? $config->{threads} : 1;
 	if($t>=2){$threads = " -pe thread $t -q threaded.q";}
 	my $mem = defined $config->{mem} ? " -l mf=$config->{mem}" : undef;
 	my $project = defined $config->{project} ? " -P $config->{project}" : " -P jdhotopp-lab";
 	my $cwd = defined $config->{cwd} ? " -cwd" : undef;
    my $name = defined $config->{name} ? " -N $config->{name}" : " -N Foobar";
    my $hostname = defined $config->{hostname} ? " -l hostname=\"$config->{hostname}\"" : undef;
    if($config->{no_gal}==1){
    	if($config->{excl}==1){$hostname = " -l hostname=\"juggernaut|grid*\""}
    	else{$hostname= " -l hostname=\"magneto|juggernaut|grid*\"";}
    }
    my $exclusive;
    if($config->{excl}==1){$exclusive = " -l excl=true";} else {$exclusive = undef;}
 	my $fh;
 	if($log){
   		open($fh,">>","$log") or confess "Can't open: $fh because: $!\n";
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
   		confess "Error: $que died with message: $report\n";
   	}
   	print $fh "QSUB: $report\n";
   	sleep 2;
   	return $report;
 }
 
=head2 Qsub3
	Title		: Qsub3
	Usage		: Qsub3(\%options);
	Function	: This Qsub is used ONLY to resubmit the same script call back to the grid.
	Args		: \%hash with:
			--input=s
			--input_list=s
			--output_dir=s
			--subdirs=i
			--Qsub=i
			--threads=i
			--sub_mem=s
			--sub_name=s
			--project=s
			--wd=s
			--cwd=i
			--hostname=s
			--excl=i
			--no_gal=i
=cut 
sub Qsub3 {
  	my $options=$_[0] || confess "Error Qsub3 didn't get receive a hash ref properly: $!\n";

  	## Setup Qsub command:
  	my $opts;
  	# $opts->{wd}= defined $options->{wd} ? " -wd $options->{wd}" : undef;
 	$opts->{threads}=undef;
 	my $t = defined $options->{threads} ? $options->{threads} : 1;
 	if($t>=2){$opts->{threads} = " -pe thread $t -q threaded.q";}
 	$opts->{mem} = defined $options->{sub_mem} ? " -l mf=$options->{sub_mem}" : undef;
 	$opts->{project} = defined $options->{project} ? " -P $options->{project}" : " -P jdhotopp-lab";
 	$opts->{cwd} = defined $options->{cwd} ? " -cwd" : undef;
    $opts->{name} = defined $options->{sub_name} ? " -N $options->{sub_name}" : undef;			## KBS 01.07.14 " -N ksieber";
    $opts->{hostname} = defined $options->{hostname} ? " -l hostname=\"$options->{hostname}\"" : undef;
    if($options->{no_gal}==1){$opts->{hostname}= " -l hostname=\"magneto|grid*\"";}
    $opts->{exclusive} = undef;
    if($options->{excl}){$opts->{exclusive} = " -l excl=true";} else {$opts->{exclusive} = undef;}
    my $qsub = "qsub -V";
    foreach my $op (keys %$opts){
    	if(defined $opts->{$op}){
    		$qsub = "$qsub"."$opts->{$op}";
    	}
    }

    ## Setup input and output dirs
  	if(!$options->{input} && !$options->{input_list}){confess "Error: No input given. Use --input or --input_list.\n";}
  	if(!$options->{output_dir}){confess "Error: No output directory given. Use --output_dir.\n";}
  	my @inputList;
  	my %out_dirs;
  	if($options->{input_list}){
     	open(LIST,"<","$options->{input_list}") or confess "Can't open --input_list: $options->{input_list} because: $!\n";
   	    while(<LIST>){
   		    chomp;
  		    push(@inputList,$_);
  		    my($fn,$dir,$suff)=fileparse($_,qr/\.[^\.]+/);
  		    $fn=~/^([A-Za-z0-9]+).*/;
  		    $out_dirs{$_}="$options->{output_dir}/$1";
  	    }
  	    close LIST;
    }
    if($options->{input}){
       push(@inputList,$options->{input});
	}
	my $original_output_dir = $options->{output_dir};
	run_cmd("mkdir -p $options->{output_dir}"); 

	## Build and submit commands
	foreach my $input (@inputList){
		my($name,$path,$suf)=fileparse($input,qr/\.[^\.]+/);
		if($options->{input_list}){
			$options->{output_dir}=$out_dirs{$input};
		} else {
			$options->{output_dir} = $original_output_dir;
		}
		if($options->{subdirs}==1){
			$options->{output_dir} = "$options->{output_dir}\/"."$name\/"; 
			run_cmd("mkdir -p $options->{output_dir}");
		}
		my $qsub_dir = defined $options->{cwd} ? undef : "-wd $options->{output_dir}";
		if($options->{input_list}){ $options->{input} = $input; }										## If we are in the orignal call, we need to make sure to qsub a single input
		my $cmd = "$^X $0";
		$cmd = $cmd." --print_hostname=1";
		foreach my $key (keys %$options){
			next if($key eq 'input_list');
			next if($key eq 'Qsub');
			next if($key eq 'subdirs');
			next if($key eq 'excl');
			next if($key eq 'no_gal');
			next if($key eq 'sub_mem');
			next if($key eq 'sub_name');
			if($options->{$key}){$cmd = $cmd." --$key=$options->{$key}"};
		}
		my $que;
		if($cmd=~/\.sh$/){
			$que = "$qsub $qsub_dir $cmd";
		} else {
			$que = "echo \"$cmd\" | $qsub $qsub_dir";
		}
		print STDERR "QSUB: $que\n";
		chomp(my $report = `$que`);
   		if($?){
   			print STDERR "$?\n";
   			confess "Error: $que died with message: $report\n";
   		}
   		print STDERR "QSUB: $report\n";
   		sleep 2;
	}
	die "+++ Finished submiting all jobs for: $0 +++\n";
}

1;
