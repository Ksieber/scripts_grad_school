#!/usr/bin/perl -I /home/ksieber/scripts/
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'opt1|foo=s', 'opt2=s', );
use Data::Dumper;
use Statistics::R;
my $R = Statistics::R->new( r_bin => "/usr/local/bin/R" );
$R->run(
    'calc_JSD <- function(inMatrix, pseudocount=0.0000001, ...) {
                    KLD <- function(x,y) { sum(x *log(x/y)) }
                    JSD<- function(x,y) { sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) }
                    matrixColSize <- length( colnames( inMatrix ) )
                    matrixRowSize <- length( rownames( inMatrix ) )
                    colnames <- colnames( inMatrix )
                    resultsMatrix <- matrix( 0, matrixColSize, matrixColSize )

                    inMatrix = apply( inMatrix, 1:2, function(x) ifelse ( x==0, pseudocount, x ) )
                    for ( i in 1:matrixColSize ) {
                        for ( j in 1:matrixColSize ) {
                            resultsMatrix[ i, j ] = JSD( as.vector( inMatrix[ , i ] ), as.vector( inMatrix[ , j ] ) )
                        }
                    }
                    colnames -> colnames( resultsMatrix ) -> rownames( resultsMatrix )
                    as.dist( resultsMatrix ) -> resultsMatrix
                    attr( resultsMatrix, "method" ) <- "dist"
                    return( resultsMatrix )
                }'
);
$R->run(
    'calc_JSD_boot_fxn = function (x_df, index) {
                tmp_df <- data.frame( x_df[,1], x_df[index,2] )
                return ( calc_JSD(tmp_df) ) 
                }'
);

$R->run("pop1=c(1,2,3,4)");
$R->run("pop2=c(1,2,3,4)");
$R->run("pop3=c(9,10,11,12)");
$R->run("pop4=c(9,10,11,12)");
$R->run("pop5=c(9,10,11,12)");
$R->run("counts = data.frame(pop1,pop5)");

$R->run('ct=prop.table(as.matrix(counts), margin=2)');

my $JSD_lines = $R->run('calc_JSD(ct)');
my @JSD_split = split(/\n/,$JSD_lines);
my $calc_JSD = (split/\s+/,$JSD_split[1])[1];
printf STDERR ("%-15s %-15.3f %-15.3f %-15.3f","$calc_JSD","$calc_JSD","$calc_JSD","$calc_JSD");

$R->run("library(boot)");
$R->run("JSD_boot <- boot(ct, calc_JSD_boot_fxn, R=1000, parallel=\"multicore\", ncpus=4)");
my $JSdist_ci = $R->run('boot.ci(JSD_boot, type="norm")');
my $ci_data_line = ( split /\n/, $JSdist_ci )[8];
my $jsd_ci_lower;
my $jsd_ci_upper;
if ( defined $ci_data_line ) {
    $ci_data_line =~ /\s+\((.+)\,\s+(.+)\)/;
    $jsd_ci_lower = $1;
    $jsd_ci_upper = $2;
}
else {
    $jsd_ci_lower = "NULL";
    $jsd_ci_upper = "NULL";
}
printf STDERR ( "\t%20s", "$jsd_ci_lower -> $jsd_ci_upper" );
print STDERR "\n";


__END__

