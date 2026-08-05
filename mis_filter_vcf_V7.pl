#!/usr/bin/perl -w
#use strict;
use PerlIO::gzip;
use Getopt::Long;


my ($vcf,$outpath,$gqa,$outvcf,$outgene,$name,$dpmin,$dpmax,$DPmin1,$DPmax1,$gqs,$miss,$maf,$pop,$dp,$gq);
GetOptions(
    "vcf=s"=>\$vcf,
    "out=s"=>\$outpath,
	"gq_all=s"=>\$gqa,
	"dp_min1=s"=>\$dpmin,
	"dp_max1=s"=>\$dpmax,
	"DP_min=s"=>\$DPmin1,
	"DP_max=s"=>\$DPmax1,
	"gq_single=s"=>\$gqs,
	"miss=s"=>\$miss,
	"maf=s"=>\$maf,
	"num=s"=>\$pop,
	"dp_pos=s"=>\$dp,
	"gq_pos=s"=>\$gq,
	"name=s"=>\$name,
);
if ($vcf=~/gz$/){
	#open IN ,"<:gzip",$vcf or die "can not open $vcf,$!.\n";
	open(IN,"gunzip -c $vcf|") or die "can't open!";
}else{
	open IN,'<',$vcf or die "can not open $vcf,$!.\n";
}
open VCF,'>:gzip',"$outpath/$name\.vcf\.gz" || die "can not open $!.\n";
open VCF1,'>:gzip',"$outpath/$name\.vcf\.all\.gz" || die "can not open $!.\n";
open GENO,'>:gzip',"$outpath/$name\.geno\.gz" || die "can not open $!.\n";
open GENO1,'>:gzip',"$outpath/$name\.geno\.gapit.gz" || die "can not open $!.\n";
open TXT,">$outpath/$name.xls" || die "can not open $!.\n";
open GENO2,'>:gzip',"$outpath/$name\.geno\.maf\.gz" || die "can not open $!.\n";
##save sample name##
while (my $line = <IN>){
	chomp $line;$line=~s/\t0\|1\:/\t0\/1\:/g;$line=~s/\t1\|0\:/\t0\/1\:/g;$line=~s/\t0\|0\:/\t0\/0\:/g;$line=~s/\t1\|1\:/\t1\/1\:/g;
	if ($line=~/^##/){
	    print VCF "$line\n";
		print VCF1 "$line\n";
	}
	elsif($line=~/^#C/){
		my @seq=split/\s+/,$line;
		my @sep=(split(/\t/,$line))[9..$pop+8];
		my $sum1=@sep;
		my $se=join "\t",@sep;
		print GENO "$seq[0]\t$seq[1]\t$seq[3]\t$se\n";
		print GENO1 "$seq[0]\t$seq[1]\t$seq[3]\t$seq[4]\t$se\n";
		print GENO2 "$seq[0]\t$seq[1]\t$seq[3]\t$seq[4]\tmaf\t$se\n";
		print VCF "$seq[0]\t$seq[1]\t$seq[2]\t$seq[3]\t$seq[4]\t$seq[5]\t$seq[6]\t$seq[7]\t$seq[8]\t$se\n";
		print VCF1 "$seq[0]\t$seq[1]\t$seq[2]\t$seq[3]\t$seq[4]\t$seq[5]\t$seq[6]\t$seq[7]\t$seq[8]\t$se\n";
		last;
	}
}
my @seq1;my @seq2;my @sep3;my $i=0;my $j=0;my $mk=0;my $mk0=0;my $mk1=0;my $sum1=0;my $total=0;my $sum=0;my $i1;my $i2;my $i3;my $sum2;my $miss01=0;my $miss02=0;my $miss03=0;my $miss04=0;my $miss05=0;my $double=0;my $sumgqa=0;
while (my $line = <IN>){
	chomp $line;$line=~s/\t0\|1\:/\t0\/1\:/g;$line=~s/\t1\|0\:/\t0\/1\:/g;$line=~s/\t0\|0\:/\t0\/0\:/g;$line=~s/\t1\|1\:/\t1\/1\:/g;
	$mk0=0;	$mk1=0;	$j=0;$total++;
	my @sep2=(split(/\t/,$line))[0..8];
	my @sep3=(split (/DP\=/,$sep2[7]));
	my @sep4=(split (/\;/,$sep3[1]));
	my @alt_list=(split (/,/,$sep2[4]));
##filter multiple alleles && low quality && total depth##	
	#if ((length ($sep2[4])> 1) || ($sep2[5]<$gqa) || $sep4[0]>$DPmax1 || $sep4[0]<$DPmin1){
	#	if (length ($sep2[4])> 1){$double++;}
	if ((@alt_list> 1) || ($sep2[5]<$gqa) || $sep4[0]>$DPmax1 || $sep4[0]<$DPmin1){
		if (@alt_list> 1){$double++;}
		elsif($sep2[5]<$gqa){$sumgqa++;}
		next;
	}
##filter individual low quality &&individual depth##	
	else {
		$sum1++;
		my @sep1=(split(/\t/,$line))[9..$pop+8];
		for ($i=0;$i<$pop;$i++){
			if($sep1[$i] =~ /\.\/\./) {    ###解决miss样品只有./.信息时的报错问题. by youqing
				$j++;
				$seq1[$i] = "./.";
				$seq2[$i] = "--";
				$seq3[$i] = "NN";
				next;
			}
			my @sep5=(split /\:/,$sep1[$i]);$sep5[$dp] = 0 if $sep5[$dp] eq "\.";$sep5[$gq] = 0 if $sep5[$gq] eq "\.";
			if (($sep5[$dp] >= $dpmin) && ($sep5[$dp] <= $dpmax) && ($sep5[$gq] >= $gqs)){
				$seq1[$i]=$sep1[$i];
				if (((split(/\:/,$sep1[$i]))[0])eq "0/0"){$seq3[$i]=$seq2[$i]="$sep2[3]$sep2[3]";$mk0=$mk0+2;}
				elsif (((split(/\:/,$sep1[$i]))[0])eq "0/1"){$seq3[$i]=$seq2[$i]="$sep2[3]$sep2[4]";$mk0++;$mk1++;}
				elsif (((split(/\:/,$sep1[$i]))[0])eq "1/1"){$seq3[$i]=$seq2[$i]="$sep2[4]$sep2[4]";$mk1=$mk1+2;}
				}
			else {
				$j++;
				$seq1[$i]="./.";
				$seq2[$i]="--";
				$seq3[$i]="NN";
			}
		}
##filter  miss && maf ##		
		$i1=$j/$pop;
		if($i1<=0.1){$miss01++;}
		if($i1<=0.2){$miss02++;}
		if($i1<=0.3){$miss03++;}
		if($i1<=0.4){$miss04++;}
		if($i1<=0.5){$miss05++;}
		if ($i1 <= $miss){
			$sum2++;
			my $mk=$mk0+$mk1;
			if ($mk0==0){$mk=1;}
			$i2=$mk0/($mk);
			$i3=$mk1/($mk);
			if ($i2>$i3){$i2=$i3;}
			if(($i2 >=$maf) && ($i3 >= $maf)){
				$sum++;
				my $line1=join("\t",@sep2)."\t".join("\t",@seq1);
				my $line2="$sep2[0]\t$sep2[1]\t$sep2[3]\t$sep2[4]\t$i2\t".join("\t",@seq2);
				my $line3="$sep2[0]\t$sep2[1]\t$sep2[3]\t".join("\t",@seq2);
				my $line4="$sep2[0]\t$sep2[1]\t$sep2[3]\t$sep2[4]\t".join("\t",@seq3);
				print VCF "$line1\n";
				print VCF1 "$line\n";
				print GENO "$line3\n";
				print GENO2 "$line2\n";
				print GENO1 "$line4\n";
			}
		}
	}
}
##print table##
my $do2=$total - $double;
my $gqsv=$total - $sumgqa - $double;
my $popfl=$gqsv - $sum1;
my $misfl=$sum1 - $sum2;
my $maffl=$sum2 - $sum;
print TXT "filter options:\n";
print TXT "GQ-single\t$gqs\n";
print TXT "DP-single\t$dpmin...$dpmax\n";
print TXT "GQ-total\t$gqa\n";
print TXT "DP-total\t$DPmin1...$DPmax1\n";
print TXT "miss\t$miss\n";
print TXT "maf\t$maf\n\n";

print TXT "SNP:\n";
print TXT "option\traw_snp_num\tdoble_alle\tGQ\tDP\tmiss\tmaf\n";
print TXT "DP-single:$dpmin...$dpmax\t \t \t$gqa\t$DPmin1...$DPmax1\t$miss\t$maf\n";
print TXT "saved\t$total\t$do2\t$gqsv\t$sum1\t$sum2\t$sum\n";
print TXT "filtered\t0\t$double\t$sumgqa\t$popfl\t$misfl\t$maffl\n\n";
print TXT "othes_miss_rate\tmis0.1\tmis0.2\tmis0.3\tmis0.4\tmis0.5\n";
print TXT "pop_miss_saved\t$miss01\t$miss02\t$miss03\t$miss04\t$miss05\n";
close IN;
close VCF;
close GENO;
close TXT;
