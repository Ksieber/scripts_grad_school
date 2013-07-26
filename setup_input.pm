package setup_input;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( setup_input );
### Returns array(ref) of input files
## Suggested use: my $input=setup_input();
## Relies on main script using $options{input} &| $options{input_list}

sub setup_input {
   my @inputList;
   our %options;
   our $results;
   if(!$main::options{input} && !$main::options{input_list}){ die
       "Error: setup_input.pm requires either $options{input} or $options{input_list}.\n";
   }
   if($main::options{input}){
       push(@inputList,$main::options{input});
   }	
   if($main::options{input_list}){
     	open(LIST,"<","$main::options{input_list}") or die "Can't open --input_list.\n";
   	    while(<LIST>){
   		    chomp;
  		    push(@inputList,$_);
  	    }
  	    close LIST;
    } 
    return \@inputList;      
}
1;
