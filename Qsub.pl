#!/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use lib('/home/kbs14104/scripts/grad_school/');
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options, 'cmd|c=s',       'log=s',  'sub_name|n=s', 'threads|t=s', 'sub_mem|mem=s', 'project|p=s', 'cwd=s',
    'wd=s',    'sub_mail|ml=s', 'help|?', 'hostname=s',   'sub_queue=s'
) or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "This script will Qsub a command. The following options can be passed.
    --cmd=          command/shell script to be qsub'ed [STDIN].
    --threads=      # cpu threads to use [1].
    --sub_name=     Job name [STDIN].
    --sub_mem=      min amount of RAM to use.
    --sub_mail=     [0|1]  1=email ~/.run_cmd.conf email when the job is complete. 
    --sub_queue=    String for the SGE queues (can be comma sep, no spaces)
    --cwd=          <0|1> [0] 1= Use Current Working Directory. 
    --wd=           Directory for grid to work from [~/].
    --hostname=     Specifiy a hosename\n";
}

if ( defined $options{cmd} ) {
    my $report = Qsub(
        {   cmd       => $options{cmd},
            log       => $options{log},
            sub_name  => $options{sub_name},
            sub_mem   => $options{sub_mem},
            sub_mail  => $options{sub_mail},
            sub_queue => $options{sub_queue},
            threads   => $options{threads},
            project   => $options{project},
            cwd       => $options{cwd},
            wd        => $options{wd},
            hostname  => $options{hostname}
        }
    );
}
else {
    while (<>) {
        chomp( my $cmd = $_ );
        my $report = Qsub(
            {   cmd       => $cmd,
                log       => $options{log},
                sub_name  => $options{sub_name},
                sub_mem   => $options{sub_mem},
                sub_mail  => $options{sub_mail},
                sub_queue => $options{sub_queue},
                threads   => $options{threads},
                project   => $options{project},
                cwd       => $options{cwd},
                wd        => $options{wd},
                hostname  => $options{hostname}
            }
        );
    }
}

