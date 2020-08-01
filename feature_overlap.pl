#!/usr/bin/perl                                                                                                  
#load some preliminary things - not important
BEGIN{
    $VERSION = "1.0";
    $home = '/mnt/iusers01/bk01/mobcdsw2/scripts';
    unshift (@INC,"$home/perl_modules/");

}

use strict;
use Getopt::Long;
use Data::Dumper;
#####


#get the file information from 'config' - basically the command line options. These are set further down
my($config);
$config = &configure(scalar @ARGV);


main($config);

sub main{
    my($config)=@_;


    #open the feature region list file
    open FEAT, $config->{feature} or die "Cannot open $config->{feature}\n";
    my %feat;
    while (<FEAT>){
	chomp;
	#split into chr, start,end
	my ($chr,$st,$end) = (/(.+):(.+)-(.+)/);
	$chr =~ s/chr//;
	#save these values in a hash
	$feat{$chr}{$st}{$end}++;
    }
    close FEAT;

    
    my %pir_count;
    my $hit_count = '0';
    my @overlap;

    #open the query interactions file
    open INT, $config->{query} or die "Cannot open $config->{query}\n";
    while (<INT>){
	chomp;
	#these are separated by comma - split up based on this
	my @ln = split(',');
	$ln[0] =~ s/chr//;
	#save values in a hash
	$pir_count{$ln[0]}{$ln[1]}{$ln[2]}++;
    }
    close INT;

    my $pir_c = '0';
    my $v_hit;

    


    #print output file header
    open OUT, ">$config->{out}" or die "Cannot open $config->{out}\n";
    print OUT "PIR_CHR,PIR_START,PIR_STOP,FEAT_CHR,FEAT_START,FEAT_STOP\n";

    #the next bit basically goes through all PIRs and for each, loops through all features (in this case CTCF sites) and assesses whether they overlap. Overlap can be in 3 ways: 1) CTCF fully within PIR. 2) PIR fully within CTCF (unlikely but accounted for) 3) CTCF straddling either side of the PIR boundary. If any of these criteria matched, we print out the PIR and CTCF. 

    #loop through PIRs...
    foreach my $pir_chr (keys %pir_count){
	foreach my $pir_st (keys %{$pir_count{$pir_chr}}){
	    foreach my $pir_end (keys %{$pir_count{$pir_chr}{$pir_st}}){

		$v_hit = 'N';
		$pir_c++;
		
		#loop through features...
		foreach my $feat_st (keys %{$feat{$pir_chr}}){
		    foreach my $feat_end (keys %{$feat{$pir_chr}{$feat_st}}){
			
			#full overlap - feature within ROI boundary
			if ($feat_st >= $pir_st && $feat_end <= $pir_end){
			    $v_hit = 'Y';
			    my $match = $pir_chr.",".$pir_st.",".$pir_end.",".$pir_chr.",".$feat_st.",".$feat_end;
			    print OUT "$match\n";
			    push(@overlap, $match);
			}
			
			#full overlap - feature exceeds ROI boundary
			if ($feat_st < $pir_st && $feat_end > $pir_end){
			    $v_hit = 'Y';
			    my $match = $pir_chr.",".$pir_st.",".$pir_end.",".$pir_chr.",".$feat_st.",".$feat_end;
			    print OUT "$match\n";                                                              
			    push(@overlap, $match);
			}
			
			#partial overlap
			if (($feat_st >= $pir_st && $feat_st <= $pir_end && $feat_end > $pir_end) ||
			    ($feat_st < $pir_st && $feat_end >= $pir_st && $feat_end <= $pir_end)){
			    
			    $v_hit = 'Y';
			    my $match = $pir_chr.",".$pir_st.",".$pir_end.",".$pir_chr.",".$feat_st.",".$feat_end;
			    print OUT "$match\n";
			    push(@overlap, $match);
			}
		    }
		}
		
		#if this PIR overlaps a feature...
		if ($v_hit eq 'Y'){
		    $hit_count++;
		}
		
	    }
	}
    }
    
    my $total = $hit_count/$pir_c;
    my $percent = $total*100;
    #open OUT, ">$config->{out}" or die "Cannot open $config->{out}\n";
    #print OUT "$hit_count,$pir_c,$total,$percent\n";
    close OUT;
    
}



sub configure{
    my $args = shift;

    my $config = {};
    my @samples;

#    $config = {'samples' => \@samples};                                                                         

    GetOptions($config,
               "feature=s",
               "query=s",
               "random",)
        || warn "error : $!\n";


    if ((!defined $config->{feature}) || (!defined $config->{query})){
        warn "Please enter file details\n";
        die;
    }

    my $n1 = $config->{feature};
    $n1 =~ s/^.+\///;
    $n1 =~ s/_stripped//;
    my $n2 = $config->{query};
    $n2=~ s/^.+\///;
    
    $config->{out} = $n1."_".$n2."_overlap.txt";


    return ($config);

}


