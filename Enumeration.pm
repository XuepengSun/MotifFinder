#Version 2.0 (May-16-2016)

package Enumeration;

use Exporter qw(import);

our @ISA = qw(Exporter);
our @EXPORT_OK= qw(seed);

use Algorithm::Combinatorics qw(variations_with_repetition);
use Data::Dumper;
use threads;
use threads::shared;
use File::Find;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Descriptive;

my @bases = ('A','T','G','C');
my $count :shared;	# motif counts
my %seed_allmotif :shared;				# Alignment for conserved loci for all motifs, loci are separated by '##', '+' for each species
my %seed_MCS :shared;					# MCS for all motifs
my $nfile :shared;					
my %data :shared;
my %seed_detail :shared;
my $extension_process :shared;
my %seed_extended :shared;
our (%conv,%seedout,%seed_sorted);		# seed_sorted contains all seed sequences filtered by MCS

sub seed{
	my ($dir,$gap,$window,$cpu,$ref,$out,$extension,$mcs,$backgroud)=@_;	
	my @files = glob("$dir/*.aln");
	
# calculate how many alignments are valid
# Criteria for valid alignments : 3' and 5' sequences > 50bp; coding region > 100bp; informative site >= 0.5 for each species
# coding, 5' and 3' sequences are stored in alignment separately

	unless(-d 'MotifTmpDir2'){mkdir 'MotifTmpDir2'}
	for(my $j=0;$j<=scalar (@files) -1;$j++){
                foreach my $t(threads->list(threads::joinable)){$t->join()}
                if(scalar threads->list(threads::running) <= $cpu){
                        my $thr=threads->new(\&AlignmentProcessing,$j,$files[$j],$ref);
                }
                else{$j--}
	}
    while(scalar threads->list(threads::running) > 0){ sleep(1) }
    foreach my $t(threads->list(threads::joinable)){$t->join()}

# combine gene alignments into one, which is separated by 40 Ns

	my $tmpseq = 'MotifTmpDir2/'.$out.'.tmpseq';
	open OUT,">$tmpseq"; 
	my (%nal,%all,@motifs);
	foreach $ky(keys %data){
		my @tp=split(/\+/,$ky);
		$nal{$tp[0]}++;
		$all{$tp[0]}{$tp[1]}{$tp[2]}=$data{$ky};
	}
	foreach $kk(keys %all){	
		foreach $pp(keys %{$all{$kk}}){
			push @{$conv{'promoter'}{$pp}},$all{$kk}{$pp}{'p'};       # for promoter regions
			push @{$conv{'coding'}{$pp}},$all{$kk}{$pp}{'c'};		  # for coding regions
			push @{$conv{'terminator'}{$pp}},$all{$kk}{$pp}{'t'};	  # for terminator regions
		}
	}

# Here, I put complementary motifs into one, and merge the promoter, coding, and terminator sequencing with their reverse complementary sequences
# So, for each motif, both strands of genome sequneces are scanned

	foreach $qw(keys %conv){
		foreach $sn(keys %{$conv{$qw}}){
			$conv{$qw}{$sn}=join('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',@{$conv{$qw}{$sn}});
		#Here, we conjunct the sequence with its reverse complement, so the seed motif will be half-numbered	
			my $revseq = join("",reverse split("",$conv{$qw}{$sn}));
			$revseq =~tr/ATGC/TACG/;
			$conv{$qw}{$sn}=$conv{$qw}{$sn}.'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'.$revseq;   # data structure: $conv{'promoter'}{'scer'} = seq + seq_revcomp
			print OUT '>',$qw,'_',$sn,"\n",$conv{$qw}{$sn},"\n";
		}
	}
	close OUT;
	print STDERR "\n\nTotal alignments used: ",scalar keys %nal,"\n\n";

# Split motifs to the number equal to threads 
	
	my %check_redudancy=(); 	# remove reverse complementary motifs
	$comb = variations_with_repetition(\@bases,6);
	while($m = $comb->next){
		my $motig_full = join("",@$m);
		my $m_revcomp = join("",reverse @$m);
		$m_revcomp =~tr/ATGC/TACG/;
		next if exists $check_redudancy{$motig_full} || exists $check_redudancy{$m_revcomp};
		$check_redudancy{$motig_full}++;
		my $fthree = join("",@$m[0..2]);
		my $lthree = join("",@$m[3..5]);
		push @motifs,join("",@$m);
		for my $gp(1..$gap){
			my $gapped = $fthree.'.{'.$gp.'}'.$lthree;
			push @motifs,$gapped;
		}
	}
	print STDERR "Total Motifs: ",scalar @motifs,"\n\n";
	my $nsplit = int(scalar @motifs / $cpu) + 1;
	while(my @submotif = splice @motifs,0,$nsplit){
		my $thr=threads->new(\&motifsearch,\@submotif,$window,$backgroud,$extension,$ref);	 # initial search by seed
	}	
# wait and join all remaining threads

	while(scalar threads->list(threads::running) > 0){
		sleep(1);
	}
	foreach my $t(threads->list(threads::joinable)){$t->join()}
	
	#print "\n",scalar keys %seed_MCS,"\n";
	#print "\n",scalar keys %seed_allmotif,"\n";
	print STDERR "\n\nMotif Searching by Seed:  Finished!\n\n";
	print STDERR "Motif sorting by MCS:  $mcs\n\n";
	
	MotifSortByMCS($mcs);			# Sort Motif by Motif conservation score
	
	print STDERR "Number of motif >= $mcs: ",scalar keys %seed_sorted,"\n\n";
	
	my @sortted_motif = keys %seed_sorted;
	my $split_motif = int((scalar @sortted_motif) / $cpu) + 1;        # Parallel processing for motif extension
		
	while(my @submotif = splice @sortted_motif,0,$split_motif){
		my $thr=threads->new(\&MotifExtension,\@submotif,$extension);	 # initial search by seed
	}
	
	while(scalar threads->list(threads::running) > 0){
		sleep(1);
	}
	foreach my $t(threads->list(threads::joinable)){$t->join()}
	
	open OUT,'>Motif_tmp';
	foreach(sort keys %seed_extended){print OUT $_,"\t",$seed_extended{$_},"\n"}
	close OUT;
}


sub MotifExtension{
	
#	Principle for Motif extension is to enumerate all possible extensions iteratively one base at a time.
#	The control set consists of a set of close motits, e.g. casual motif: ABC-m-XYZ, the control sets are (notA)BC-m-(notX)YZ, (notA)BC-m-X(notY)Z and so on
#	For each Iteration, the best chi-squared base were added, extension is ceased if chisq does not excess the 3sd of z score (chi), this maybe changable
	my ($ConMotif,$bp_extension) = @_;
	my %IUB=(
				'A' => 'A',
				'T' => 'T',
				'G' => 'G',
				'C' => 'C',
				'R' => '[AG]',
				'Y' => '[CT]',
				'K' => '[GT]',
				'M' => '[AC]',
				'S' => '[GC]',
				'W' => '[AT]',
				'B' => '[CGT]',
				'D' => '[AGT]',
				'H' => '[ACT]',
				'V' => '[ACG]',
				'.' => '.'
	);
	
	foreach my $seedsorted(@$ConMotif){
		$extension_process++;
		print STDERR "Extension process: motif $extension_process \n" if $extension_process % 10 == 0;
		my ($gap_size)=$seedsorted=~/(\d+)/;
		$gap_size //=0;
		my ($three_left,$three_right)=$seedsorted=~/(^.{3}).*(.{3}$)/;
		my $motif_new = '.' x  $bp_extension . $three_left . '.' x $gap_size . $three_right . '.' x $bp_extension;   # ....ABC......XYZ....
		my @seed_control=();	# store the seed variaties: [^A]BC[^X]YZ ...
		my %control_sets=();
		
		for my $i(0..2){
			my @split_left = split("",$three_left);
			$split_left[$i]='[^'.$split_left[$i].']';
			my $new_left = join("",@split_left);
			for my $j(0..2){
				my @split_right = split("",$three_right);
				$split_right[$j]='[^'.$split_right[$j].']';
				my $new_tmp = $new_left.join("",@split_right);
				push @seed_control,$new_tmp;
			}
		}

		foreach my $allseed(keys %seed_MCS){
			my ($seed_gap)=$allseed=~/(\d+)/;	
			$seed_gap //=0;
			next if $seed_gap != $gap_size;
			$gap_remove = $allseed;
			$gap_remove =~s/[^A-Z]//g;
			foreach (@seed_control){
				if($gap_remove=~/$_/){
					my ($first,$last)=$gap_remove=~/(^.{3}).*(.{3}$)/;
					$control_sets{$allseed}='.' x  $bp_extension . $first . '.' x $seed_gap . $last . '.' x $bp_extension;
					my @tpp = split("",$control_sets{$allseed});
					$control_sets{$allseed}=\@tpp;
					last;
				}
			}
		}
		
		my @causal_original = split("",$motif_new);
		
		for(1..$gap_size+$bp_extension*2){
			my @extension_site=();
			while($motif_new=~/\./g){push @extension_site,length($`)}
			my %ext_causal_result=();
			my %ext_control_result=();
			my @causal_ConAlignments = split(/\#\#/,$seed_allmotif{$seedsorted});
			foreach(@causal_ConAlignments){
				my @subsplit = split(/\+/,$_);
				foreach my $site(@extension_site){
					foreach my $base(keys %IUB){
						my @tmp = @causal_original;
						$tmp[$site]=$base;
						my $tmp_motif;
						foreach(@tmp){ $tmp_motif.=$IUB{$_}	}
						my $align_count=0;
						foreach(@subsplit){$align_count++ if $_=~/$tmp_motif/i}
						$ext_causal_result{$site}{$base}++ if $align_count == (scalar @subsplit);
					}
				}
			}
			
			my $num_controlAlignments = 0;
			foreach my $control_motif(keys %control_sets){
				my @control_ConAlignments = split(/\#\#/,$seed_allmotif{$control_motif});
				foreach(@control_ConAlignments){
					$num_controlAlignments ++ ;
					my @subsplit = split(/\+/,$_);
					foreach my $site(@extension_site){
						foreach my $base(keys %IUB){
							my @tmp = @{$control_sets{$control_motif}};
							$tmp[$site]=$base;
							my $tmp_motif;
							foreach(@tmp){ $tmp_motif.=$IUB{$_}	}
							my $align_count=0;
							foreach(@subsplit){$align_count++ if $_=~/$tmp_motif/i}
							$ext_control_result{$site}{$base}++ if $align_count == (scalar @subsplit);
						}
					}
				}
			}

# χ2 = [n(ad – bc)2] / [(a + b) (c + d) (a + c) (b + d)]
			my $ext_site = '';
			my $ext_base = '';
			my $chisq = -999;
					
			foreach my $site(keys %ext_causal_result){
				foreach my $base(keys %{$ext_causal_result{$site}}){
					$ext_causal_result{$site}{$base} //=0;
					$ext_control_result{$site}{$base} //=0;
					my $causal_conserved = $ext_causal_result{$site}{$base};
					my $causal_nonconserved = (scalar @causal_ConAlignments) - $causal_conserved;
					my $control_conserved = $ext_control_result{$site}{$base};
					my $control_nonconserved = $num_controlAlignments - $control_conserved;
					my $n = $num_controlAlignments + (scalar @causal_ConAlignments);
					
					# Here, chi square is two tailed, but we only need over-represented motif, so we should remove the lower-represented ones
					next if scalar @causal_ConAlignments == 0 || $num_controlAlignments ==0;
					next if $causal_conserved /(scalar @causal_ConAlignments) <= $control_conserved / $num_controlAlignments;
					
					my $chi_square = ($n*($causal_conserved*$control_nonconserved - $control_conserved*$causal_nonconserved)**2)/(($causal_conserved+$control_conserved)*($causal_nonconserved+$control_nonconserved)*($causal_conserved+$causal_nonconserved)*($control_conserved+$control_nonconserved));
					
#					print STDERR $causal_conserved,"\t",$causal_nonconserved ,"\t",$control_conserved,"\t",$control_nonconserved,"\t";
#					print STDERR $chi_square,"\n";
					
					if($chi_square >= $chisq){
						$ext_site = $site;
						$ext_base = $base;
						$chisq = $chi_square;
					}
				}
			}
			# Set criteria to continue extension : z**2 = χ2; χ2 >= 3*z;
			if($chisq < 0 || $chisq <= 3 * sqrt($chisq) ){			# For unknown situations
				my $motif_final = join("",@causal_original);
				$motif_final=~s/^\.+|\.+$//g;
				$seed_extended{$seedsorted} = $motif_final;
				last;
			}
			else{
				$causal_original[$ext_site] = $ext_base;
				$motif_new = join("",@causal_original);
				foreach (keys %control_sets){
					$control_sets{$_}->[$ext_site] = $ext_base;
				}
			}
		}
		
		my $motif_final2 = join("",@causal_original);
		$motif_final2 =~s/^\.+|\.+$//g;
		$seed_extended{$seedsorted} = $motif_final2;
	}
}


sub MotifSortByMCS{
	my $mcs_cutoff = shift;
	my %tempdata=();
	foreach $seed(keys %seed_MCS){
		my ($gapsize)=$seed=~/(\d+)/;
		$gapsize //=0;
		$tempdata{$gapsize}{$seed}=$seed_MCS{$seed};
	}

# To caculate MCS foreach motif, the geometric mean for each group of motifs were caculated, and the stadard deviation from this mean was generated
# The motifs outside the 3sd were removed before caculating the final mean for each group. MCS is the Z score under binomial model
# Alternative method for removing outliers is using quartile, but the results may be similar.	
	
	foreach $gs(sort{$a<=>$b} keys %tempdata){
		my @tempvalues = values %{$tempdata{$gs}};
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@tempvalues);
		my $geometric_mean = $stat->geometric_mean();   		# also named as log-average
		my @temp = map{($_ - $geometric_mean)**2} @tempvalues;
		my $sd = sqrt((sum @temp)/(scalar @temp -1));
		my @sort = ();
		foreach (@tempvalues){push @sort,$_ if $_ >= $geometric_mean-3*$sd && $_<=$geometric_mean+3*$sd}
		my $mean = (sum @sort)/(scalar @sort);
		foreach(keys %{$tempdata{$gs}}){
			my ($motif_con,$motif_total)= $seed_detail{$_}=~/(\d+)_(\d+)/;
			next if $motif_total == 0;
			$seedMCS = ($motif_con - $motif_total*$mean)/sqrt($motif_total*$mean*(1-$mean));     # Z score under binomial model
			if($seedMCS >= $mcs_cutoff){
				$seed_sorted{$_}=$seedMCS;
			}
		}	
	}
}


sub motifsearch{

# For each motif, record its position in Reference sequence, then extract the sequences from all speceies by windows centered at matched positions.
# Record the occurance and conservation, for the conserved motif, extract the 'extension' bp of sequences from both end, this will be used for motif extension

	my ($sub_motifs,$window_size,$region_type,$motif_extension,$ref_species) = @_;
	$|=1;
	foreach my $motif(@$sub_motifs){
		$count++;
		print STDERR "\rMotif Processed: $count" if $count % 20 ==0;
		my ($dot)=$motif=~/(\d+)/;
		$dot //=0;
		my %data=();
		my %pos_conserved=();
		while($conv{$region_type}{$ref_species}=~/$motif/ig){
			my $match_pos = length($`)+1;
			$data{$match_pos}++;
		}
		foreach my $pos(keys %data){
			my @seq_with_extension=();
			my $start_pos = max(1,$pos-$window_size);
			my $extract_len = $pos+6+$dot+$window_size-$start_pos;
			foreach my $species(keys %{$conv{$region_type}}){
				my $temp_seq = substr($conv{$region_type}{$species},$start_pos-1,$extract_len);
				if($temp_seq=~/$motif/ig){
					my $sub_start = length($`)+$start_pos-$motif_extension;
					if($sub_start < 1){@seq_with_extension=(); last}
					else{
						my $sub_len = $motif_extension*2 + 6 + $dot;
						push @seq_with_extension,substr($conv{$region_type}{$species},$sub_start-1,$sub_len);
					}
				}
				else{
					@seq_with_extension=();
					last;
				}
			}
			if(@seq_with_extension){
				$pos_conserved{$pos} = join('+',@seq_with_extension);
			}
		}
		my $mcsore = (scalar keys %pos_conserved) / (scalar keys %data);
		$mcsore //=0;
		$seed_MCS{$motif} = $mcsore;
		my @values = values %pos_conserved;
		$seed_allmotif{$motif}=join('##',@values);
		$seed_detail{$motif}=(scalar keys %pos_conserved).'_'.(scalar keys %data);
	}
}


sub AlignmentProcessing{

# Alignments should be in the standard form: lower case for promoter and terminator, upper case for coding sequences

	my ($id,$fname,$refsp)= @_;
	$nfile++;
	$| =1;
	print STDERR "\rProcessing Alignments: $nfile";
    my (%sep,$strand,%tmp,%align);
	open IN,$fname || die "cannot open $fname\n";
    while(<IN>){
      	next if /^\s+|^CLUSTAL/;
		chomp;
		my @s=split /\s+/;
		$tmp{$s[0]}.=$s[1];
	}
    close IN;
    $tmp{$refsp}=~/[A-Z].*[A-Z]/;
    my $match = $&;
	my $mlen = length $match;
    my $llen = length $`;
    my $rlen = length $';
    $match =~s/[^A-Z]//g;
    my $lw = substr($match,0,3);
    my $rw = substr($match,length($match)-3,3);
    if($lw eq 'ATG'){$strand = 'pos'}
    elsif($rw eq 'CAT'){$strand = 'neg'}
    my $check = 0;
	foreach $sp(keys %tmp){         # set criterion for valid alignment
        my $lseq = substr($tmp{$sp},0,$llen);
        my $mseq = substr($tmp{$sp},$llen,$mlen);
        my $rseq = substr($tmp{$sp},$mlen+$llen,$rlen);
		my $pro = $id.'+'.$sp.'+'.'p';
		my $cds = $id.'+'.$sp.'+'.'c';
		my $trm = $id.'+'.$sp.'+'.'t'; 
        if($strand eq 'pos'){
            $align{$pro}=$lseq;
            $align{$cds}=$mseq;
            $align{$trm}=$rseq;
        }
        elsif($strand eq 'neg'){
            $align{$pro}=$rseq;
            $align{$cds}=$mseq;
            $align{$trm}=$lseq;
        }
        else{next}
        $lseq=~s/[^A-Za-z]//g;
        $mseq=~s/[^A-Za-z]//g;
        $rseq=~s/[^A-Za-z]//g;
        if (length($lseq) < 50 || length($rseq) < 50 || length($mseq) < 100){
			$check++;
        }
	}
	if($check == 0){
		foreach (keys %align){$data{$_}=$align{$_}}
	}
}

return 1;





