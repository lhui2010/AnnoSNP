#!/usr/bin/env perl

use 5.12.0;
use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::FeatureIO;#THE MODULE GFF.pm was modified to bypass irregular Terms
use Bio::Align::RadicalChanges;#THE MODULE RadicalChanges was adopted from DNAStatistics, adding radical mutation detection to get_syn_changes()


use Bio::SeqFeature::Annotated;
use Bio::SeqIO;

use warnings;
#use Data::Dumper;

use Bio::RangeI;
use Bio::Coordinate::GeneMapper;

#use Time::HiRes;

#for bugs, please write to lhui2010@gmail.com


#compare the derived allele frequency based on concensus genotype of the four populations
#Tue Jan  8 12:48:22 CST 2013



our %CodonStat = Bio::Align::RadicalChanges ->new()->get_syn_changes;#$CodonStat{'ATG'}{'ATG'} returns 0, after Description becomes Synony..
our @DescriptionForCodonStat = (qw/Synonymous NonSynonymous/, ('Radical')x15);


#use Bio::Tools::CodonTable;
#use Bio::Align::DNAStatistics;


#my $myCodonTable   = Bio::Tools::CodonTable->new();


if(@ARGV <3)
{
	die "usage:  $0.pl  gff fasta snp \n";
}


my $gff=Bio::FeatureIO->new(-file => "$ARGV[0]", -format => 'GFF', -version => 3);
my $fa = Bio::SeqIO->new(-file => "$ARGV[1]", -format => 'Fasta');
open my $snp_file, $ARGV[2] or die;
open my $out_file, ">$ARGV[2].snp_stat" or die;
open my $spec_file, ">$ARGV[2].snp_spec" or die;


warn "Reading fasta...";
#open $snp_file, "tmp.snp" or die;
warn "Done\n";



#######################read fasta##########################
my %seq;
while(my $tmp = $fa->next_seq())
{
	        $seq{$tmp->display_id()} = $tmp;
}


warn "Reading gff...\n";



######################read gff into memory##################


##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
our %array_of_gene_start;#used for determin the loci's nearest gene, sorted by ascending;
our %start2end;#return end of the gene given its start and chr
our %start2gene;#return the gene of this start

my %gene; #containing total genes from this gff file
while ( my $feature = $gff->next_feature() ) {


#       no strict 'refs';
#       print "Instance METHOD IS  " . Dumper( \%{ref ($seq)."::" }) ;exit;
        if($feature->type ->name eq "mRNA")
        {
                $start2gene{$feature->seq_id}{$feature -> start} = $feature -> get_Annotations ('ID')->value;
                $start2end{$feature->seq_id}{$feature -> start}  = $feature -> end;
        }

        next unless($feature-> type ->name eq "CDS");

        push  @{$gene{($feature -> get_Annotations ('Parent'))->value}}, $feature;

}
for my $sid( %start2gene)
{
        $array_of_gene_start{$sid} = [sort {$a<=>$b} (keys %{$start2gene{$sid}})];
}

warn "Done\nBuilding gene mapper...";


###################Building Gene Mapper#####################

use Storable 'dclone';

my (%chr_mapper, %cds_mapper);

for my $this_gene (sort keys  %gene)
{
	$chr_mapper{$this_gene} = Bio::Coordinate::GeneMapper->new(
		-in  => 'chr',
		-out => 'cds',
		-exons => [
			map {
				Bio::Location::Simple->new(
					-start  => $_ ->start,#location->{_start},
					-end =>    $_ ->end, #location->{_end},
					-strand => $_ ->strand, #location->{_strand},
					-seq_id => $_ ->seq_id,
					)
			} @{$gene{$this_gene}}
			],
		);
	$cds_mapper{$this_gene} =  dclone($chr_mapper{$this_gene});
	$cds_mapper{$this_gene} -> swap;
}



warn "Done\nReading SNP file...\n\n\n";
##############read SNP file#############


my %gene_mutation_count;#count of three type per gene, nonsynonymous (including radical), synonymous, radical
my %gene_mutation_specific;


while(<$snp_file>)
{
#preparint SNP
        chomp;
        next if($_ eq "");
	my ($chr, $location, $ref_base, $snp_base, undef) = split /\t/, $_;
	next if ($ref_base eq '.' or $snp_base eq '.');


	my @tmp_array = &SNP_stat($chr, $location,$ref_base, $snp_base, \%chr_mapper, \%cds_mapper, \%seq);#exit;

	for my $tmp(@tmp_array)
	{
		if($$tmp{strand} ne "NA")
		{
			$gene_mutation_count{$$tmp{gene_name}}{$$tmp{mutation_type}} ++;#= $this_frequency;

			if ($$tmp{mutation_type} eq "Radical")
			{
				$gene_mutation_count{$$tmp{gene_name}}{Nonsynonymous} ++;
			}

			push @{$gene_mutation_specific{$$tmp{gene_name}}}, {%$tmp};
		}
	}

}#end reading snp




warn "Done\nProgram finished..Outputing result into $ARGV[2].stat\n";

print $out_file $_,"\t" for (qw/GeneName Nonsynonymous Synonymous Radical/);
print $out_file "\n";

for my $name( sort keys %gene_mutation_count)
{
	print $out_file "$name\t";

	for my $type(qw/Nonsynonymous Synonymous Radical/)
	{
		if(exists $gene_mutation_count{$name}{$type})
		{
			print $out_file $gene_mutation_count{$name}{$type};
		}
		else
		{
			print $out_file 0;
		}
		print $out_file "\t";
	}
	print $out_file "\n";
}

for my $name (sort keys %gene_mutation_specific)
{
	for my $ar( @{$gene_mutation_specific{$name}})
	{
		if ( $$ar{strand}<0)
		{
			$$ar{strand} = "-";
		}
		else
		{
			$$ar{strand} = "+";
		}
		for my $key (qw/gene_name chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type/)
		{
			print $spec_file $$ar{$key}, "\t";
		}
		print $spec_file "\n";
	}
}


warn "Clearing memory...\n";

#############sub routine###############


#accept: chromosome_id  coordinate_on_chr Ref-base SNP-base \%chr_mapper, \%cds_mapper
#return an a hash containing the following elements
#chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type

sub SNP_stat
{
#my $start = Time::HiRes::time();

#print @DescriptionForCodonStat;
#print %CodonStat;
#exit;
        #all results will be stored in here
        #chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type
        my @results;
        my %hash_return;
        $hash_return{$_} = "NA" for(qw/chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type ref_pep snp_pep gene_name/);


        my ($chr_mapper, $cds_mapper, $seq);
        ($hash_return{"chromosome"}, $hash_return{"location_chr"}, $hash_return{"ref_base"}, $hash_return{"snp_base"}, $chr_mapper, $cds_mapper, $seq) = @_;

        #like: chromosome01 100 A G \%chr_mapper \%_cds_mapper

#warn (Time::HiRes::time()-$start);
#$start = Time::HiRes::time();
        return @results if($hash_return{"ref_base"} eq $hash_return{"snp_base"} or !exists $hash_return{"ref_base"} or !exists $hash_return{"snp_base"});


##Sear for nearest gene##
        if($hash_return{"location_chr"} <$array_of_gene_start{$hash_return{"chromosome"}}[0])
        {
                return @results;# this SNP gets up too early, no gene is evolved yet!
        }
#from the first to the last, no one escapes!
        my (@nearest_gene);#, $nearest_gene_start);
        for my $search_for_nearest_gene( @{$array_of_gene_start{$hash_return{"chromosome"}}})
        {
                if($hash_return{"location_chr"} >= $search_for_nearest_gene and
                        $hash_return{"location_chr"}<= $start2end{$hash_return{"chromosome"}}{$search_for_nearest_gene})#found it's stage, next is future , now is most important
                {
#                       last;
                        push @nearest_gene, $start2gene{$hash_return{"chromosome"}}{$search_for_nearest_gene};
                }
        #       else
        #       {
        #               $nearest_gene_start = $search_for_nearest_gene;
        #       }
        }
#TODO 改为正链找一次，然后负链找一次


        my $loci_chr = Bio::Location::Simple->new(
                -start => $hash_return{"location_chr"},
                -end => $hash_return{"location_chr"},
                -seq_id => $hash_return{"chromosome"}
                );
        for my $this_gene(@nearest_gene)
        {
#clear cache
                $hash_return{$_} = "NA" for(qw/strand location_cds phase ref_codon snp_codon mutation_type ref_pep snp_pep gene_name/);
                $hash_return{gene_name} = $this_gene;
                next unless (exists $$chr_mapper{$this_gene});

                my $loci_cds = $$chr_mapper{$this_gene}->map($loci_chr);

#found the corresponding gene
                if(defined $loci_cds ->start)
                {
                        $hash_return{"location_cds"} = $loci_cds ->start;

                        $hash_return{"strand"} = $$chr_mapper{$this_gene} ->map($loci_chr) ->strand;

                        my @codon_phase_delay;
                        my (@ref_codon, @snp_codon);
                        my (@codon_loci, @codon_base);

                        $hash_return{"phase"}=$hash_return{"location_cds"}%3 -1;#get phase and transform from 1 2 0 to 0 1 2 ;
                        $hash_return{"phase"} +=3 if($hash_return{"phase"} <0);#0 1 2

#                       if($loci_start->strand >0) seems free of strand problem
#                               {
                                 if($hash_return{"phase"} == 0)
                                {
                                        @codon_phase_delay=(0, 1, 2);
                                }
                                elsif($hash_return{"phase"} == 1)
                                {
                                        @codon_phase_delay = (-1, 0, 1);
                                }
                                elsif($hash_return{"phase"} == 2)
                                {
                                        @codon_phase_delay = (-2, -1, 0);
                                }

 #       eval
#        {
                                @codon_loci = map {
                                        $$cds_mapper{$this_gene}->map(
                                        Bio::Location::Simple->new(
                                        -start=>$hash_return{"location_cds"}+$_,
                                        -end=>$hash_return{"location_cds"}+$_,
                                        )
                                        ) ->start
                                        }@codon_phase_delay;

                                        my $tmp_ref = $hash_return{"ref_base"};
                                        my $tmp_snp = $hash_return{"snp_base"};

                                if($hash_return{"strand"} > 0)
                                {
                                        @codon_base = map{uc($$seq{$hash_return{"chromosome"}}->subseq($_,$_))}@codon_loci;
                                }
                                else
                                {
                                        @codon_base = map{ uc($$seq{$hash_return{"chromosome"}}->subseq($_,$_))}@codon_loci;
                                        $_ =~ tr/ATCG/TAGC/ for(@codon_base);

                                        $tmp_ref =~ tr/ATCG/TAGC/;
                                        $tmp_snp =~ tr/ATCG/TAGC/;
                                }

                                @ref_codon = @codon_base;#TODO use Seq object instead
                                $ref_codon[$hash_return{"phase"}] = $tmp_ref;
                                @snp_codon = @codon_base;
                                $snp_codon[$hash_return{"phase"}] = $tmp_snp;
  #      };
#                       }

                        $hash_return{"ref_codon"} = join("", @ref_codon);
                        $hash_return{"snp_codon"} = join("", @snp_codon);

                        $hash_return{"mutation_type"} = $DescriptionForCodonStat[$CodonStat{$hash_return{"ref_codon"}}{$hash_return{"snp_codon"}}];

                        use Bio::Tools::CodonTable;
                        my $myCodonTable   = Bio::Tools::CodonTable->new();

                        $hash_return{"ref_pep"} = $myCodonTable->translate($hash_return{"ref_codon"});
                        $hash_return{"snp_pep"} = $myCodonTable->translate($hash_return{"snp_codon"});

#               return %hash_return; #ignoring overlapping genes
                        push @results, {%hash_return};
                }#endif
        }#end for
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();exit;
#       return %hash_return;#no corresponding geen found
        return @results;
}



