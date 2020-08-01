#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;

my($config);
$config = &configure(scalar @ARGV);


main($config);

sub main{
    my($config)=@_;
    
    my %pir;

    #read in Promoter-interacting DNA regions (PIRs)
    open PIRS, $config->{DNA_PIR} or die "Cannot open $config->{DNA_PIR}\n";
    while (<PIRS>){
	chomp;
	unless (/\#/){
	    my ($chr, $start, $stop) = (/.+\t(.+?),(.+?),(.+)/);
	    $pir{$chr}{$start}{$stop}++;
	}
    }
    
    my %rir;

    #read in RNA-interacting DNA regions (RIRs)
    open RIRS, $config->{DNA_RNA} or die "Cannot open $config->{DNA_RNA}\n";
    while (<RIRS>){
	chomp;
	unless (/^#/){
	    my ($chr, $start, $stop) = (/(.+?)\t(.+?)\t(.+?)\t.+/);
	    $rir{$chr}{$start}{$stop}++;
	    #print "#$chr# #$start# #$stop#\n";

	}

    }
    
    my %pir_with_rna;

    #overlap between the two 
    foreach my $chr (keys %pir){
	foreach my $p_st (keys %{$pir{$chr}}){
	    foreach my $p_end (keys %{$pir{$chr}{$p_st}}){
		
		#compare PIR regions to RIR regions
		foreach my $r_st (keys %{$rir{$chr}}){
		    foreach my $r_end (keys %{$rir{$chr}{$r_st}}){

			my $complex = $chr."_".$p_st."_".$p_end."_".$chr."_".$r_st."_".$r_end;

			#RIR within the PIR
			if ($r_st >= $p_st && $r_end <= $p_end){
			    $pir_with_rna{$complex}++;
			}

			#RIR overlapping either end
			if ($r_st >= $p_st && $r_st <= $p_end && $r_end > $p_end){
                            $pir_with_rna{$complex}++;
			}
		    
			if ($r_st < $p_st && $r_end >= $p_st && $r_end <= $p_end){
                            $pir_with_rna{$complex}++;
			}

			#PIR is within RIR
			if ($r_st < $p_st && $r_end > $p_end){
                            $pir_with_rna{$complex}++;
                        }
		    }
		}
	    }
	}
    }
    
    
    open OUT, ">$config->{out}";

    foreach my $int (keys %pir_with_rna){
	print OUT "$int\n";
    }
    
    close OUT;

}



sub configure{
    my $args = shift;

    my $config = {};
    my @samples;

#    $config = {'samples' => \@samples};                                                            

    GetOptions($config,
               "DNA_PIR=s",
               "DNA_RNA=s",
	       "random")
        || warn "error : $!\n";


    if ((!defined $config->{DNA_RNA}) || (!defined $config->{DNA_PIR})){
        warn "Please enter file details\n";
        die;
    }

    #if a random input file, take the number and use in output
    if (defined $config->{random}){    
	$config->{out} = $config->{DNA_PIR};
	$config->{out} =~ s/\.txt/_PIR_RIR_overlap\.txt/;
    }
    else {
	$config->{out} = 'PIR_RIR_overlap.txt';
    }

    return ($config);

}



=head




#open the ROI file
my $ = $ARGV[1];

open EQTL, $eQTL or die "Cannot open $eQTL\n";
my %c;
while (<EQTL>){
    chomp;
    my @ln = split('\t');
    my ($chr,$pos) = ($ln[0] =~ /(\d+?)_(\d+?)_.+/);
    my $slope = $ln[4];

    $c{$chr}{$pos} = $slope;
}
close EQTL;


my %match1;
my %match2;

my %uniq_counter;
my $total = '0';
my @hit_slope;

my $int = $ARGV[0];
open INT, $int or die "Cannot open $int\n";
while (<INT>){
    chomp;
    my @ln = split('\t');
    
    my ($pir_chr, $pir_st, $pir_end) = ($ln[1] =~ /(.+?),(.+?),(.+)/);
    
    $pir_chr =~ s/chr//;
    my $v_hit = 'N';

    #match eQTL chr and loop through matches
    foreach my $v_pos (keys %{$c{$pir_chr}}){
	#if pos falls within the PIR
	if ($v_pos >= $pir_st && $v_pos <= $pir_end){
	    $v_hit = 'Y';
	    push(@hit_slope, $c{$pir_chr}{$v_pos});
	}
    }
}

=head
    $total++;

    my $p = $pir_chr.",".$pir_st.",".$pir_end;
    $uniq_counter{$p}++;

    if ($v_hit eq 'Y'){

	#add whole promoter-PIR
	$match1{$_}++;
	#add unique PIR
	$match2{$p}++;
     
    }



close INT;

my $total2 = keys %uniq_counter;


my $c1 = keys %match1;
my $c2 = keys %match2;

print "OVERLAPPING unique PIR = $c1 out of $total\n";
print "OVERLAPPING unique PIR = $c2 out of $total2\n";
=cut
