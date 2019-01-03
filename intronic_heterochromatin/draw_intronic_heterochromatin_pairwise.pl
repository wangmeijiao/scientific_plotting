use strict;
use warnings;
use Getopt::Long;
use SVG;
use Data::Dumper;
use POSIX qw[ _exit ];


#use methods in /home/mjwang/pwdexx/oryza_epiCompara_sixSpecies/05.promotorHetBorder/localView_big/localView_big.pl
#               /exporthome/mjwang/pwdex/oryza_epiCompara_knob1/drawKnob/others/drawDensity_orthGenesPair_polygon/drawDensityCurve_orthGenesPair.pl
##gene level, simple and quick draw
##data tracks use curve instead of barplot
##use densi_customized_norm with input substraction
##combine replications into one track with opacity and transparency

#perl this_file.pl -info genome.info -anchor "japo|Os08g16559|chr08:10130791-10131842|1051,brac|Ob08g17720|chr08:5897489-5899053|1564" -outdir test -quiet | display

my ($info_f,$anchor,$anno,$outdir);
GetOptions("info=s",\$info_f,"anchor=s",\$anchor,"anno=s",\$anno,"outdir=s",\$outdir);

my $f_norm_brac_seedling = 1.051244; # all data of mRNA_bracRS_seedling use this
my $f_norm_brac_tillering = 0.837; # all data of mRNA_bracRS_tillering use this

my $f_norm_brac_k43_seedling = 0.231273 ; # all data of k43_bracRS_seedling use this
my $f_norm_brac_k43_tillering = 0.57123 ; # all data of k43_bracRS_tillering use this


#1, initialize
if(! -d $outdir){`mkdir -p $outdir`} #perl mkdir can't make dir recursively
die "no bin dir exists\n" if(! -d "bin");
die "no all2embl.pl exists\n" if(! -f "bin/all2embl.pl");
my $info=&readInfo($info_f); #print Dumper $info;exit;
my ($order,$anchors,$regions, $maxLen)=&anchorParse($anchor);
my ($anno1,$anno2) = split/\n/,$anno;
my ($height,$width)=(600,600);
my $ppx=$width/$maxLen;
my $space=100; 
my $boxH=22; #track total height
my $n=scalar @{$order};
my $svg=SVG->new('height'=>$height+2*$space,'width'=>$width+2*$space);
#print Dumper $info,$order,$anchors,$regions, $maxLen,$ppx; exit;


my %color=(
#"k92_bw_seedling_rep1","rgb(151,3,1)","k92_bw_seedling_rep2","rgb(98,13,95)","k92_bw_tilling_60d","rgb(206,68,42)" 
           "japo" => {"gene","rgb(0,51,153)","DNA-TE","grey","DNA","grey","LTR","black","LINE","grey","SINE","grey","RC","grey","link","black","link_rev","black","k92_bw_seedling_rep1","rgb(151,3,1)","k92_bw_seedling_rep2","rgb(151,3,1)","k92_bw_tilling_60d","rgb(98,13,95)","H3K27me","red","k43_bdg_seedling","rgb(51,5,156)","k43_bdg_tilling_60d","rgb(51,5,156)","H3K56Ac","green","mRNA_bdg_seedling","rgb(51,51,255)","mRNA_bdg_tilling_60d","rgb(51,51,255)","siRNA","purple","GC_content","black","AT_content","black","DHsite","grey","CG","#121a2a","CHG","#ef4136","CHH","#FF1FDF","gap","grey","mappability","grey","gene_island","black"},

           "brac" => {"gene","rgb(0,51,153)","DNA-TE","grey","DNA","grey","LTR","black","LINE","grey","SINE","grey","RC","grey","link","black","link_rev","black","k92_bw_seedling","rgb(151,3,1)","k92_bw_tilling_60d","rgb(98,13,95)","k92_bw_tilling_vari","rgb(98,13,95)","H3K27me","red","k43_bdg_seedling","rgb(51,5,156)","k43_bdg_tilling_60d","rgb(51,5,156)","H3K56Ac","green","mRNA_bdg_seedling","rgb(51,51,255)","mRNA_bdg_tilling_60d","rgb(51,51,255)","siRNA","purple","GC_content","black","AT_content","black","DHsite","grey","CG","#121a2a","CHG","#ef4136","CHH","#FF1FDF","gap","grey","mappability","grey","gene_island","black"}

);
#"k92_bw_tilling_vari","rgb(206,68,42)"
my $color_link = "black";

#my $max_fpkm = 550;
my $max_fpkm = 300;

my %range=(

           "japo" => {"k43_bdg_seedling","0-0-15", "k43_bdg_tilling_60d","0-0-55","k56ac","1-1-55","k92_bw_seedling_rep1","0-0-2","k92_bw_seedling_rep2","0-0-2","k92_bw_tilling_60d","0-0-2","k271","2-2-55","mRNA_bdg_seedling","0-0-55","mRNA_bdg_tilling_60d","0-0-$max_fpkm","siRNA","1-1-25","DHsite","2-2-25" },

           "brac" => {"k43_bw_seedling","0-0-3","k43_bdg_seedling","0-0-15","k43_bdg_tilling_60d","0-0-55","k56ac","1-1-55","k92_bw_seedling","0-0-3","k92_bw_tilling_60d","0-0-3","k92_bw_tilling_vari","0-0-6","k271","1-1-95","mRNA_bdg_seedling","0-0-55","mRNA_bdg_tilling_60d","0-0-$max_fpkm","siRNA","1-1-45"}


);



$svg->text("x",$space+10,"y",$space-20,"-cdata",$anno1,"text-anchor","start","stroke","black");
$svg->text("x",$space+10,"y",$space-10,"-cdata",$anno2,"text-anchor","start","stroke","black");

#2, get data and draw for each species
my $y_next=$space;
my $fa_extract_file_last;
foreach my $i(0..$#$order){
   my $spec=$order->[$i];
   print STDERR "process species $spec ";
   my ($chr,$range,$len)=split/:|\|/,$regions->{$spec};
   my ($s,$e) = split/-/,$range;
   my $len_chr=&readLen($info->{$spec}->{"chrSize"})->{$chr};
   die "len_chr illegal\n" if(!defined $len_chr || $len_chr <= 0);
   if($s <0 ){die "s <0 in $spec region $regions->{$spec}"}
   if($e > $len_chr){die "e > len_chr($len_chr) in $spec region $regions->{$spec}"}
   my $region="$chr:$s-$e";
   my $anchorId=$anchors->{$spec};
   my $prefix = $spec."_".$anchorId;
   print STDERR "with region $region, anchor $anchorId\n";

   #1,functional tracks

   #k43_seedling 
   my ($k43_extract,$k43_extract_file,$peaks_k43,$peaks_k43_extract_file);
   if(exists $info->{$spec}->{"k43_bdg_seedling"} ){
      print STDERR "  k43_bdg_seedling..\n";
      ($k43_extract,$k43_extract_file)=&bdgExtract_wig($info->{$spec}->{"k43_bdg_seedling"},$region,$outdir,$prefix,"k43");
      if($spec eq "brac"){foreach my $i(0..$#$k43_extract){my ($c,$s,$e,$v) = split/[ \t]+/,$k43_extract->[$i]; $v=$v*$f_norm_brac_k43_seedling;$k43_extract->[$i]="$c\t$s\t$e\t$v"   } }
      if(exists $info->{$spec}->{"k43_peak_seedling"}){
        ($peaks_k43,$peaks_k43_extract_file)=&bedExtract($info->{$spec}->{"k43_peak_seedling"},$region,$outdir,$prefix,"k43_peaks");
         $y_next=&drawBdg_slim($k43_extract,$peaks_k43,$y_next,$region,"k43_bdg_seedling",$color{$spec}->{'k43_bdg_seedling'},$ppx,$range{$spec}->{'k43_bdg_seedling'}, "thin","barplot");
      }else{ 
             $y_next=&drawBdg_slim($k43_extract,"none",$y_next,$region,"k43_bdg_seedling",$color{$spec}->{'k43_bdg_seedling'},$ppx,$range{$spec}->{'k43_bdg_seedling'}, "thin","barplot");
           }
   }

   #k43_tilling_60d
   #$y_next += 15;
   if(exists $info->{$spec}->{"k43_bdg_tilling_60d"} ){
      print STDERR "  k43_bdg_tilling_60d..\n";
      ($k43_extract,$k43_extract_file)=&bdgExtract_wig($info->{$spec}->{"k43_bdg_tilling_60d"},$region,$outdir,$prefix,"k43");
      if($spec eq "brac"){foreach my $i(0..$#$k43_extract){my ($c,$s,$e,$v) = split/[ \t]+/,$k43_extract->[$i]; $v=$v*$f_norm_brac_k43_tillering;$k43_extract->[$i]="$c\t$s\t$e\t$v"   } }

      if(exists $info->{$spec}->{"k43_peak_tilling_60d"}){
        ($peaks_k43,$peaks_k43_extract_file)=&bedExtract($info->{$spec}->{"k43_peak_tilling_60d"},$region,$outdir,$prefix,"k43_peaks");
         $y_next=&drawBdg_slim($k43_extract,$peaks_k43,$y_next,$region,"k43_bdg_tilling_60d",$color{$spec}->{'k43_bdg_tilling_60d'},$ppx,$range{$spec}->{'k43_bdg_tilling_60d'},"thin","barplot");
      }else{ 
             $y_next=&drawBdg_slim($k43_extract,"none",$y_next,$region,"k43_bdg_tilling_60d",$color{$spec}->{'k43_bdg_tilling_60d'},$ppx,$range{$spec}->{'k43_bdg_tilling_60d'},"thin","barplot");
           }
   }


=pod
   #k56ac densi
   my ($k56ac_extract,$k56ac_extract_file,$peaks_k56ac,$peaks_k56ac_extract_file);
   if(exists $info->{$spec}->{"k56ac_bdg"} ){
      print STDERR "  k56ac..\n";
      ($k56ac_extract,$k56ac_extract_file)=&bdgExtract_wig($info->{$spec}->{"k56ac_bdg"},$region,$outdir,$prefix,"k56ac");
      if(exists $info->{$spec}->{"k56ac_peaks"}){
        ($peaks_k56ac,$peaks_k56ac_extract_file)=&bedExtract($info->{$spec}->{"k56ac_peaks"},$region,$outdir,$prefix,"k56ac_peaks");
         $y_next=&drawBdg_slim($k56ac_extract,$peaks_k56ac,$y_next+20,$region,"H3K56Ac",$ppx,$range{$spec}->{'k56ac'});
      }else{ 
             $y_next=&drawBdg_slim($k56ac_extract,"none",$y_next+20,$region,"H3K56Ac",$ppx,$range{$spec}->{'k56ac'});
           }
   }
=cut

   #k92_seedling_rep1 customizedNorm
   my ($k92_extract,$k92_extract_file,$peaks_k92,$peaks_k92_extract_file);
   my %mark1 = ("japo","k92_bw_seedling_rep1","brac","k92_bw_seedling");
   my %peak1 = ("japo","k92_peak_seedling_rep1","brac","k92_peak_seedling");
   if(exists $info->{$spec}->{$mark1{$spec}} ){
      print STDERR "  $mark1{$spec}..\n";
      ($k92_extract,$k92_extract_file)=&bdgExtract_wig($info->{$spec}->{$mark1{$spec}},$region,$outdir,$prefix,"k92");
      if(exists $info->{$spec}->{$peak1{$spec}}){
        ($peaks_k92,$peaks_k92_extract_file)=&bedExtract($info->{$spec}->{$peak1{$spec}},$region,$outdir,$prefix,"k92_peaks");
         #$y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,"H3K9me2_seedling",$ppx,$range{$spec}->{'k92_bw_seedling'},"bold","barplot");
         $y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,$mark1{$spec},$color{$spec}->{$mark1{$spec}},$ppx,$range{$spec}->{$mark1{$spec}},"NA","curve");
      }else{ 
             $y_next=&drawBdg_slim($k92_extract,"none",$y_next,$region,"k92_bw_seedling_rep1",$color{$spec}->{'k92_bw_seedling_rep1'},$ppx,$range{$spec}->{'k92_bw_seedling_rep1'},"bold","barplot");
           }
   }

   #k92_bw_seedling_rep2 customizedNorm
   my %mark2 = ("japo","k92_bw_seedling_rep2","brac","k92_bw_tilling_60d");
   my %peak2 = ("japo","k92_peak_seedling_rep2","brac","k92_peak_tilling_60d");
   if(exists $info->{$spec}->{$mark2{$spec}} ){
      print STDERR "  $mark2{$spec}..\n";
      ($k92_extract,$k92_extract_file)=&bdgExtract_wig($info->{$spec}->{$mark2{$spec}},$region,$outdir,$prefix,"k92");
      if(exists $info->{$spec}->{$peak2{$spec}}){
        ($peaks_k92,$peaks_k92_extract_file)=&bedExtract($info->{$spec}->{$peak2{$spec}},$region,$outdir,$prefix,"k92_peaks");
         #$y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,"k92_bw_tilling_60d_vari",$ppx,$range{$spec}->{'k92_bw_tilling_60d_vari'},"bold","barplot");
         $y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,$mark2{$spec},$color{$spec}->{$mark2{$spec}},$ppx,$range{$spec}->{$mark2{$spec}},"NA","curve");
      }else{ 
             $y_next=&drawBdg_slim($k92_extract,"none",$y_next,$region,"k92_bw_seedling_rep2",$color{$spec}->{'k92_bw_seedling_rep2'},$ppx,$range{$spec}->{'k92_bw_seedling_rep2'},"bold","barplot");
           }
   }


   #k92_tilling_60d customizedNorm
   my %mark3 = ("japo","k92_bw_tilling_60d","brac","k92_bw_tilling_vari");
   my %peak3 = ("japo","k92_peak_tilling_60d","brac","k92_peak_tilling_vari");
   if(exists $info->{$spec}->{$mark3{$spec}} ){
      print STDERR "  $mark3{$spec}..\n";
      ($k92_extract,$k92_extract_file)=&bdgExtract_wig($info->{$spec}->{$mark3{$spec}},$region,$outdir,$prefix,"k92");
      if(exists $info->{$spec}->{$peak3{$spec}}){
        ($peaks_k92,$peaks_k92_extract_file)=&bedExtract($info->{$spec}->{$peak3{$spec}},$region,$outdir,$prefix,"k92_peaks");
         #$y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,"k92_bw_tilling_60d",$ppx,$range{$spec}->{'k92'},"bold","barplot");
         $y_next=&drawBdg_slim($k92_extract,$peaks_k92,$y_next,$region,$mark3{$spec},$color{$spec}->{$mark3{$spec}},$ppx,$range{$spec}->{$mark3{$spec}},"NA","curve");
      }else{ 
             $y_next=&drawBdg_slim($k92_extract,"none",$y_next,$region,"k92_bw_tilling_60d",$color{$spec}->{'k92_bw_tilling_60d'},$ppx,$range{$spec}->{'k92_bw_tilling_60d'},"bold","barplot");
           }
   }



=pod
   #k271 densi
   my ($k271_extract,$k271_extract_file,$peaks_k271,$peaks_k271_extract_file);
   if(exists $info->{$spec}->{"k271_bdg"} ){
      print STDERR "  k271..\n";
      ($k271_extract,$k271_extract_file)=&densiExtract_wig($info->{$spec}->{"k271_bdg"},$region,$outdir,$prefix,"k271");
      my $k271_extract_pack = &bdgPack($k271_extract);
      if(exists $info->{$spec}->{"k271_peaks"}){
        ($peaks_k271,$peaks_k271_extract_file)=&bedExtract($info->{$spec}->{"k271_peaks"},$region,$outdir,$prefix,"k271_peaks");
         $y_next=&drawBdg_slim($k271_extract_pack,$peaks_k271,$y_next+20,$region,"H3K27me",$ppx,$range{$spec}->{'k271'},"thin");
      }else{ 
             $y_next=&drawBdg_slim($k271_extract_pack,"none",$y_next+20,$region,"H3K27me",$ppx,$range{$spec}->{'k271'});
           }
   }
=cut


   #2,annotation tracks
   print STDERR "  annotation tracks..\n";

   #mrna_seedling
   if(exists $info->{$spec}->{"mRNA_bdg_seedling"} ){
      print STDERR "  mRNA_bdg_seedling..\n";
      my ($mrna_extract,$mrna_extract_file)=&bdgExtract_wig($info->{$spec}->{"mRNA_bdg_seedling"},$region,$outdir,$prefix,"mRNA");
     if($spec eq "brac"){foreach my $i(0..$#$mrna_extract){my ($c,$s,$e,$v) = split/[ \t]+/,$mrna_extract->[$i]; $v=$v*$f_norm_brac_seedling;$mrna_extract->[$i]="$c\t$s\t$e\t$v"   } }
     $y_next=&drawBdg_slim($mrna_extract,"none",$y_next,$region,"mRNA_bdg_seedling",$color{$spec}->{'mRNA_bdg_seedling'},$ppx,$range{$spec}->{'mRNA_bdg_seedling'},"thin","barplot");   
   }

   #mrna_tilling_60d
   my ($mrna_extract,$mrna_extract_file);
   if(exists $info->{$spec}->{"mRNA_bdg_tilling_60d"} ){
      print STDERR "  mRNA_bdg_tilling_60d..\n";
      my ($mrna_extract,$mrna_extract_file)=&bdgExtract_wig($info->{$spec}->{"mRNA_bdg_tilling_60d"},$region,$outdir,$prefix,"mRNA");
     if($spec eq "brac"){foreach my $i(0..$#$mrna_extract){my ($c,$s,$e,$v) = split/[ \t]+/,$mrna_extract->[$i]; $v=$v*$f_norm_brac_tillering;$mrna_extract->[$i]="$c\t$s\t$e\t$v"   } }
     $y_next=&drawBdg_slim($mrna_extract,"none",$y_next,$region,"mRNA_bdg_tilling_60d",$color{$spec}->{'mRNA_bdg_tilling_60d'},$ppx,$range{$spec}->{'mRNA_bdg_tilling_60d'},"thin","barplot");   
   }

   ##gene
   print STDERR "    gene\n";
   my ($gff_extract,$gff_extract_file_gene)=&gffExtract($info->{$spec}->{"geneGff"},$region,$outdir,$prefix);
   $y_next=&drawGff($gff_extract,$region,$y_next,$anchorId,$color{$spec}->{'gene'});          
   ##LTR
   $y_next+=15;
   print STDERR "    LTR\n";
   my ($bed_extract_LTR,$gff_extract_file_LTR)=&bedExtract($info->{$spec}->{"teBed_LTR"},$region,$outdir,$prefix,"LTR");
   $y_next=&drawBed($bed_extract_LTR,$region,$y_next,"LTR",$color{$spec}->{'LTR'}); 
   ##DNA
   print STDERR "    DNATE\n";
   my ($bed_extract_DNA,$gff_extract_file_DNA)=&bedExtract($info->{$spec}->{"teBed_DNA"},$region,$outdir,$prefix,"repeat_region");
   $y_next=&drawBed($bed_extract_DNA,$region,$y_next,"DNA-TE",$color{$spec}->{'DNA-TE'}); 
   ##LINE
   print STDERR "    LINE\n";
   my ($bed_extract_LINE,$gff_extract_file_LINE)=&bedExtract($info->{$spec}->{"teBed_LINE"},$region,$outdir,$prefix,"repeat_region");
   $y_next=&drawBed($bed_extract_LINE,$region,$y_next,"LINE",$color{$spec}->{'LINE'}); 
   ##SINE
   print STDERR "    SINE\n";
   my ($bed_extract_SINE,$gff_extract_file_SINE)=&bedExtract($info->{$spec}->{"teBed_SINE"},$region,$outdir,$prefix,"repeat_region");
   $y_next=&drawBed($bed_extract_SINE,$region,$y_next,"SINE",$color{$spec}->{'SINE'}); 
   ##RC
   print STDERR "    RC\n";
   my ($bed_extract_RC,$gff_extract_file_RC)=&bedExtract($info->{$spec}->{"teBed_RC"},$region,$outdir,$prefix,"repeat_region");
   $y_next=&drawBed($bed_extract_RC,$region,$y_next,"RC",$color{$spec}->{'RC'}); 


   #3,links: blastm8/orthTab  (space between species, output for ACT)
   print STDERR "  draw links..\n";
   $y_next+=10;
   my $fa_extract_file_out;
   if($i != $#$order){
      my $spec_next=$order->[$i+1];
      my ($chr,$s,$e,$strand,$len)=split/:|-|\|/,$regions->{$spec_next};
      #my $len=$e-$s+1;
      my $len_chr=&readLen($info->{$spec_next}->{"chrSize"})->{$chr};
      die "len_chr illegal\n" if(!defined $len_chr || $len_chr <= 0);
      if($s <0 ){die "s <0 in $spec_next region $regions->{$spec_next}"}
      if($e > $len_chr){die "e > len_chr($len_chr) in $spec_next region $regions->{$spec_next}"}
      my $region_next="$chr:$s-$e";
      my $anchorId_next=$anchors->{$spec_next};
      my $prefix_next = $spec_next."_".$anchorId_next;
      my ($blastm8,$blastm8_file,$fa_extract_file,$fa_extract_file_next)=&getBlastm8($info->{$spec}->{"fa"},$info->{$spec_next}->{"fa"},$region,$region_next,$outdir,$prefix,$prefix_next,"blast",1e-5);
      $fa_extract_file_out=$fa_extract_file;
      if($i == $#$order-1){ $fa_extract_file_last = $fa_extract_file_next}
      $y_next=&drawBlastm8($blastm8,$y_next,$color_link);
   }


}#species foreach

#$svg->text("x",$space+10,"y",$y_next+15,"-cdata",$anno2,"text-anchor","start","stroke","black");

print STDERR "done\n\n";

my $prefix_outfile = `basename $outdir`;
chomp $prefix_outfile;
my $outfile = $outdir."/$prefix_outfile.svg";
&svgWrite($svg,$outfile);
POSIX::_exit(0); #to speed up exit by ommitting perl_destruct for big hash



###############################subs####################################
sub readInfo(){ #genome info file
  my %info;
  #my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n+/,$_;
   my $species=shift @box;
   #push @order,$species;
   foreach(@box){
     $_=~s/[\t ]+$//;
     my ($type,$file)=split/[\t ]+/,$_;
     next if($type =~/^"/);
     if(!exists $info{$species}->{$type}){ $info{$species}->{$type}=$file }else{die "dup $type in species $species\n"} 
   }
  }
  close INF;
  $/="\n";
  return \%info;
}


sub anchorParse(){
    my $string=shift @_;#japo|Os08g16559|chr08:10130791-10131842|1051,brac_1.4sub|Ob08g17720|chr08:5897489-5899053|1564
    my %anchor;
    #my %chr;
    my %coord;
    my @order;
    my $maxLen;
    my @temp=split/,/,$string;
    foreach(@temp){
      my($spec,$id,$region,$len)=split/\|/,$_;
      my ($chr, $range) = split/:/,$region;
      my ($start,$end) = split/-/,$range;
      push @order,$spec;
      if(!exists $anchor{$spec}){$anchor{$spec}=$id}else{die "species dup: $string\n"}
      #if(!exists $chr{$spec}){$chr{$spec}=$chr}else{die "species dup: $string\n"}
      if(!exists $coord{$spec}){$coord{$spec}="$chr:$start-$end|$len"}else{die "species dup: $string\n"}
      if(!defined $maxLen){$maxLen = $len}else{if($maxLen < $len){ $maxLen = $len}  }
    }
    return \@order,\%anchor,\%coord,$maxLen;
    #return \@order,\%chr,\%anchor,\%coord;
   
}


sub readBed(){
   my $file=shift;
   die "file $file not exists\n" if(! -f $file);
   my %bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my @box=split/\t/,$_;
     #die "trunct bed file $file at $_\n" if(scalar @box < 6);
     if(!exists $bed{$box[3]}){$bed{$box[3]}=join"\t",@box   }
   }
   close BED;
   return \%bed;
}

sub readLen(){
  my %len;
  my $file=shift @_;
  die "file $file not exists\n" if(! -f $file);
  open LEN, $file or die "$!";
  while(<LEN>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,$len)=split/[\t| ]+/,$_;
    if(!exists $len{$chr}){$len{$chr}=$len}else{die"dup $chr\n"}
  }
  close LEN;
  return \%len;
}

sub getLenAll(){
  my ($info,$archi)=@_;
  my @len;
  foreach my $item(@{$archi}){
   my($species,$chr)=split/:/,$item; #don't use $_
   my $len=&readLen($info->{$species}->{"faSize"})->{$chr} ;  # here $_ will be unpredictable when call another sub in a sub
   push @len, $len;
  }
  return \@len;
}

sub gffExtract(){
   my ($file,$region,$outdir,$prefix)=@_;
   die "file $file not exists\n" if(! -f $file);
   my ($chr,$s,$e)=split/:|-/,$region;
   my %gff_region;
   my $outfile="$outdir/${prefix}_${chr}_${s}_${e}.gff";
   open GFF, $file or die "$!";   
   open GFFO, ">$outfile" or die "$!";
   while(<GFF>){  #only capture mrna and CDS lines
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my ($chr_in,$foo1,$ftype,$start,$end,$foo2,$strand,$fshift,$attr)=split/\t/,$_;
     die"empty value at $_\n" if(!defined $chr_in || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
     ($start,$end)=($start<$end)?($start,$end):($end,$start);
     next if($chr ne $chr_in); #not the same chr
     next if( $start >= $e || $end <= $s); #not overlap at all
     if($start < $s ){$start = $s}; # cut pending ends
     if($end > $e){$end = $e};
     my ($id,$pid);
     if($ftype eq "mRNA"){
        $attr=~/ID=([^;]+)/;
        $id=$1;
        if(!exists $gff_region{$id}){
          $gff_region{$id}->{"range"}="$start-$end"; 
          $gff_region{$id}->{"strand"}=$strand; 
          my ($start_real,$end_real);
          $start_real=$start-$s+1; #change to real coord for gff+fa -> embl 
          $end_real=$end-$s;
          print GFFO "$chr_in\t$foo1\t$ftype\t$start_real\t$end_real\t$foo2\t$strand\t$fshift\tID=$id\n";
        }else{die "CDS comes before mrna, quit\n"}
     }elsif($ftype eq "exon"){
        $attr=~/Parent=([^;]+)/;    
        $pid=$1;
        if(exists $gff_region{$pid}){
          push @{$gff_region{$pid}->{"exon"}},"$start-$end";
          my ($start_real,$end_real);
          $start_real=$start-$s+1; #change to real coord for gff+fa -> embl 
          $end_real=$end-$s;
          print GFFO "$chr_in\t$foo1\t$ftype\t$start_real\t$end_real\t$foo2\t$strand\t$fshift\tParent=$pid\n";
        }else{die "exon comes before mrna, quit\n"}
      }elsif($ftype eq "CDS"){
        $attr=~/Parent=([^;]+)/;
        $pid=$1;
        if(exists $gff_region{$pid}){
          push @{$gff_region{$pid}->{"CDS"}},"$start-$end-$fshift";
          my ($start_real,$end_real);
          $start_real=$start-$s+1; #change to real coord for gff+fa -> embl 
          $end_real=$end-$s;
          print GFFO "$chr_in\t$foo1\t$ftype\t$start_real\t$end_real\t$foo2\t$strand\t$fshift\tParent=$pid\n";
        }else{die "CDS comes before mrna, quit\n"}
      }      
   }#while end
   close GFF;
   close GFFO;
   #print Dumper \%gff_region;     
   return (\%gff_region,$outfile);

}

sub drawGff(){
  my ($table,$region,$y_gene,$anchorId,$col)=@_;
  my ($chr,$start_seg,$end_seg)=split/:|-/,$region;
  my $x_left=$space;
  my $x_right=$space+($end_seg-$start_seg+1)*$ppx;
  #1,draw detailed mrna-cds structure
  my $flag_drawname = 0;
  if(keys %{$table} <= 5){$flag_drawname = 1}
  foreach my $id( keys %{$table}){
    #if($id eq $anchorId){$col = "red"}else{$col = "black"}
    my $strand=$table->{$id}->{"strand"};
    my ($start,$end)=split/-/,$table->{$id}->{"range"};
    my (@exon, @cds);
    if(defined $table->{$id}->{"exon"}){@exon=@{$table->{$id}->{"exon"}}}else{print STDERR "no exon for gene $id\n"}
    if(defined $table->{$id}->{"CDS"}){@cds=@{$table->{$id}->{"CDS"}}}else{print STDERR "no cds for gene $id\n"}
    my $edge=$space;

    #draw gene backbone line
    $svg->line("x1",$edge+($start-$start_seg)*$ppx,"y1",$y_gene,"x2",$edge+($end-$start_seg)*$ppx,"y2",$y_gene,"stroke",$col,"stroke-width",1);
    if($flag_drawname == 1){$svg->text("x",$edge+($start+($end-$start)*0.45-$start_seg)*$ppx,"y",$y_gene+15,"font-family","Arial-Narrow","font-size",8,"stroke",$col,"text-anchor","middle","-cdata",$id)} #only draw gene names if gene number <= 5

 #draw exon and CDS details without consider the strand
  if (!defined $table->{$id}->{"exon"}){}else{
    my $barH=5;
    foreach my $region(@exon){
      my ($start,$end)=split/-/,$region;
      $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_gene-0.5*$barH,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$col,"stroke",$col   );

    }
  }

  if (!defined $table->{$id}->{"CDS"}){}else{
    my $barH=10;
    foreach my $region(@cds){
      my ($start,$end)=split/-/,$region;
      $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_gene-0.5*$barH,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$col,"stroke",$col   );

    }
  }

 }#foreach end
 return $y_gene;
}#sub end


sub bdgExtract_wig(){
   my ($file,$region,$outdir,$prefix,$mark)=@_;
   die "file $file not exists\n" if(! -f $file);
   my ($chr,$range)=split/:/,$region;
   my ($s,$e) = split/-/,$range;
   my $outfile="$outdir/${prefix}_${mark}_${chr}_${s}_${e}.wig";
   my @data_region;
   if($file=~/bedgraph\.bgz$/){ @data_region=`tabix $file $region `}elsif($file=~/bw$/){#keep the first line, cut in the drawBdg_slim if overhang
     @data_region=`bigWigToBedGraph -chrom=$chr -start=$s -end=$e $file stdout `
   }
   if(scalar @data_region <=2){print STDERR "too few data at $file\n"; return (\@data_region,$outfile)}
   open OUT, ">$outfile" or die "$!";
   print OUT "track type=wiggle_0 color=255,200,0\nvariableStep chrom=$chr span=1\n";
   my ($chr0,$start0,$end0,$v0)=split/[ \t]+/,$data_region[0];
   if($start0 != 1){print OUT "1\t0\n"}#in case trunct display in act 
   foreach(@data_region){
      chomp;
      die "empty line\n" if($_ eq "" || !defined $_);
      my ($chr,$start,$end,$v)=split/[ \t]+/,$_;
      next if($v == 0);
      my $start_local=$start-$s+1;
      if($start_local <= 0){$start_local = 1}
      print OUT "$start_local\t$v\n";
   }
   close OUT;
   return (\@data_region,$outfile);

}

sub drawBdg_slim(){
   #drawBdg_slim version 0.1, Nov. 13, 2018; with barplot (bold/thin) and curve two modes
   my ($data,$peaks,$y,$region,$mark,$color,$ppx,$range,$stroke,$mode) = @_;
   $mode ||= "barplot";
   $stroke||="bold";
   my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
   my $len = $e_r-$s_r+1;
   my ($min,$mid,$max)=split/-/,$range;
   $min||=0;
   $max||=15;
   my $f=$boxH/$max;# $boxH defined at the top head of this script
   my $mid_height=$boxH*$mid/$max;
   $svg->line("x1",$space,"y1",$y+$boxH,"x2",$space+$len*$ppx,"y2",$y+$boxH,"stroke","grey","stroke-width",0.5);
   $svg->text("x",$space-18,"y",$y+$boxH-2,"-cdata",$mark,"text-anchor","end","stroke","black");

  #y axis ticks and text
   $svg->line("x1",$space,"y1",$y+$boxH+5,"x2",$space,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
   for(my $i=0;$i<=$boxH;$i+=$boxH/2){
     $svg->line("x1",$space,"y1",$y+$boxH-$i,"x2",$space-2,"y2",$y+$boxH-$i,style=>{"stroke-width",1,"stroke","black"});
     if($i == 0){
       $svg->text("x",$space-5,"y",$y+$boxH-$i-2,"-cdata",int $min,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
     }elsif($i==$boxH){
       $svg->text("x",$space-5,"y",$y+$boxH-$i,"-cdata",int $max,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
     }
   }

   #draw data with barplot or filled curve
   if($mode eq "barplot"){
     foreach my$line(@$data){
       #$svg->rect("x",$x,"y",$y+15-$value*$scale{"all"}*$scale{$mark},"width",1,"height",$value*$scale{"all"}*$scale{$mark},"fill",$color{$mark});
       my ($chr_in,$s,$e,$value)=split/\t/,$line;
       next if($chr_r ne $chr_in || $e < $s_r || $s > $e_r);
       $s-=$s_r; if ($s < 0){ $s = 0}
       if($e > $e_r){ $e = $e_r}
       $e-=$s_r; if ($e < 0){ $e = 0}
       $value-=$min;
       if($value<0){$value=0}
       if($value>$max){$value=$max}
       if($value <= $mid){
        #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f,"fill","grey","stroke","none")
        $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$mid_height,"width",($e-$s+1)*$ppx,"height",$value*$f,"fill","grey","stroke","none")
        #$svg->rect("x",$x,"y",$y+$box_height-$mid_height,"width",1,"height",$value*$f,"fill","grey","stroke","none")
       }else{
           if($stroke eq "bold"){
             $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$value*$f,"width",($e-$s+1)*$ppx,"height",$value*$f-$mid_height,"fill",$color,"stroke",$color{$mark},"stroke-width",0.2) #use this for short region
           }elsif($stroke eq "thin"){
              $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$value*$f,"width",($e-$s+1)*$ppx,"height",$value*$f-$mid_height,"fill",$color,"stroke","none","fill-opacity",0.8,"stroke-opacity",0.8) #use this for big region
            }
          }
     }
   }elsif($mode eq "curve"){
      my (@x,@y);
      my (@x_axis,@y_axis);
      my (@x_sinked,@y_sinked);
      push @x, $space;
      push @y, $y+$boxH+0.5;
      foreach my $i(0..$#$data-1){
      #foreach my$line(@$data_unpack){ #don't unpack, polygon too much points, adobe illustrator will fail to display
        my ($chr_in,$s,$e,$value)=split/\t/,$data->[$i];
        my ($chr_in_next,$s_next,$e_next,$value_next)=split/\t/,$data->[$i+1];
        die "unsorted region $data->[$i] and $data->[$i+1] when draw $mark" if($e > $s_next);
        next if($chr_r ne $chr_in || $e < $s_r || $s > $e_r);
        $s-=$s_r; if ($s < 0){ $s = 0}
        if($e > $e_r){ $e = $e_r}
        $e-=$s_r; if ($e < 0){ $e = 0}

        $s_next-=$s_r; if ($s_next < 0){ $s_next = 0}
        if($e_next > $e_r){ $e_next = $e_r}
        $e_next-=$s_r; if ($e_next < 0){ $e_next = 0}

        $value-=$min;
        if($value<0){$value=0}
        if($value>$max){$value=$max}
        if($value < $mid){#the sinked part
            push @x_sinked, $space+0.5*($e+$s)*$ppx;
            push @y_sinked, $y+$boxH-$mid_height;
        }else{#the true signal part
              push @x, $space+0.5*($e+$s)*$ppx; #use midpoint to represent bed region
              push @y, $y+$boxH-$value*$f;
          }

        #fill intermedia region with value = 0
        if($e+1 == $s_next || $e == $s_next){}else{
           push @x, $space+($e+0.09*($s_next-$e))*$ppx; #use midpoint to represent bed region
           push @y, $y+$boxH;

           push @x, $space+0.5*($e+$s_next)*$ppx; #use midpoint to represent bed region
           push @y, $y+$boxH;

           push @x, $space+($e+0.98*($s_next-$e))*$ppx; #use midpoint to represent bed region
           push @y, $y+$boxH;
        }

      }#foreach line
          
      ##get coords
      push @x, $space+$len*$ppx;
      push @y, $y+$boxH;
      my $points=$svg->get_path(
                                 x=>\@x,
                                 y=>\@y,
                                 -type=>"polygon",
                                 #-closed=>'true',
      );
      ##draw polygon
      $svg->polygon(
              %$points,
              style => {
                           #'opacity' => 0.4, 
                           'fill-opacity' => 0.7, #transparency can only be observed when plot multi-layers
                           'fill' =>$color,
                           #'stroke' => "none",
                           'stroke' => $color,  #stroke color can add good transparency effect to plot
                           'stroke-opacity' => 1,
                           'stroke-width' => 1,
                       }
      )
       
    }else{die "unknow plotting type $mode"}

   #draw peak regions
   if($peaks ne "none"){
     foreach my $id(keys %{$peaks}){
       my($chr_in,$s,$e,$strand)=split/-|:/,$peaks->{$id};
       if($s< $s_r && $e > $s_r){$s=$s_r}
       if($s< $e_r && $e > $e_r){$e=$e_r}
       $svg->rect("x",$space+($s-$s_r)*$ppx,"y",$y+$boxH+2,"width",($e-$s+1)*$ppx,"height",4,"fill",$color);
       #$svg->rect("x",$space+($s-$s_r)*$ppx,"y",$y+$boxH+2,"width",($e-$s+1)*$ppx,"height",4,"fill","black");
     }
   }
 
   return $y+$boxH+4+3;
}#sub end



sub bedExtract(){
   my ($file,$region,$outdir,$prefix,$ftype)=@_;
   die "file $file not exists\n" if(! -f $file);
   my ($chr,$s,$e)=split/:|-/,$region;
   my %bed_region;
   my $outfile="$outdir/${prefix}_${ftype}_${chr}_${s}_${e}.gff";
   open BED, $file or die "$!";
   open GFFO, ">$outfile" or die "$!";
   while(<BED>){
      chomp;
      next if($_ eq "" || $_=~/^#/);
      my ($chr_in,$start,$end,$id,undef,$strand)=split/[\t ]+/,$_;
      if($strand eq "."){$strand = "+"};
      next if($chr_in ne $chr);
      ($start,$end)=($start<$end)?($start,$end):($end,$start);
      next if( $start >= $e || $end <= $s); #not overlap at all
      if($start < $s ){$start = $s}; # cut pending ends
      if($end > $e){$end = $e};
      if(!exists $bed_region{$id}){ $bed_region{$id}="$chr_in:$start-$end:$strand" }else{die "dup $id in bed file $file\n"}
      my $start_real=$start-$s+1;
      my $end_real=$end-$s+1;
      print GFFO "$chr_in\tbedExtract\t$ftype\t$start_real\t$end_real\t.\t$strand\t.\tID=$id\n";
   }
   close BED;
   close GFFO;
   return (\%bed_region,$outfile);

}


sub drawBed(){  #for annotation tracks only, not for peaks
  my($table,$region,$y,$ftype,$col)=@_;
  my ($chr,$start_seg,$end_seg)=split/:|-/,$region;
  my $x_left=$space;
  my $x_right=$space+($end_seg-$start_seg+1)*$ppx;
  my $y_te=$y+20;
  my $barH=4;
  #draw one line anno
  
  #$svg->line("x1",$x_left,"y1",$y_te,"x2",$x_right,"y2",$y_te,"stroke","grey","stroke-width",0.5);#te border line
  $svg->text("x",$x_left-5,"y",$y_te,"font-family","Arial-Narrow","font-size",8,"font-weight","bold","text-anchor","end","-cdata",$ftype);
  ##$svg->text("x",$x_left-5,"y",$y_te-2,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","+");
  #$svg->text("x",$x_left-15,"y",$y_te-2,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","->");
  #$svg->text("x",$x_left-15,"y",$y_te+5,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","<-");
  ##$svg->text("x",$x_left-5,"y",$y_te+5,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata"," -");

  foreach my $id(keys %{$table}){
     my($chr,$start,$end,$strand)=split/-|:/,$table->{$id};
     #if($strand eq "+"){
     #   $svg->rect("x",$x_left+($start-$start_seg)*$ppx,"y",$y_te-6,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$ftype},"stroke",$color{$ftype}   );
     #}elsif($strand eq "-"){
     #   $svg->rect("x",$x_left+($start-$start_seg)*$ppx,"y",$y_te+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$ftype},"stroke",$color{$ftype}   );
     # }
     $svg->rect("x",$x_left+($start-$start_seg)*$ppx,"y",$y_te-6,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$col,"stroke",$col);

  }
  return $y_te-10;
}


sub getBlastm8(){
   my ($file,$file_next,$region,$region_next,$outdir,$anchorId,$anchorId_next,$method,$evalue)=@_;
   $evalue||=1e-5; #too small evalue like 1e-10 will lost links
   die "file $file not exists\n" if(! -f $file);
   die "fai file ${file}.fai not exists\n" if(! -f ${file}.".fai");
   die "file $file_next not exists\n" if(! -f $file_next);
   die "fai file ${file_next}.fai not exists\n" if(! -f ${file_next}.".fai");
   my ($chr,$s,$e)=split/:|-/,$region;
   my ($chr_next,$s_next,$e_next)=split/:|-/,$region_next;
   #my $outfile="$outdir/${prefix}_${ftype}_${chr}_${s}_${e}.gff";
   #extract fa
   my $prefix="${anchorId}_${chr}_${s}_${e}";
   my $prefix_next="${anchorId_next}_${chr_next}_${s_next}_${e_next}";
   readpipe("samtools faidx $file $region > $outdir/${prefix}.fa");
   readpipe("samtools faidx $file_next $region_next > $outdir/${prefix_next}.fa");

   #blast first to next, get blastm8 
   my $outfile="$outdir/${prefix}_VS_${prefix_next}.blast";
   if($method eq "blast+"){
         readpipe("cat $outdir/${prefix_next}.fa |makeblastdb -dbtype nucl -out $outdir/${prefix_next}.fa -title ${prefix_next}.fa");
         readpipe("blastn -task blastn -db $outdir/${prefix_next}.fa -out $outfile -outfmt 6 -num_threads 2 -evalue $evalue -query $outdir/${prefix}.fa");

   }elsif($method eq "blast"){
        readpipe("formatdb -i $outdir/${prefix_next}.fa -p F");   #formatdb write out files as the same location as input file
        readpipe("rm formatdb.log");
        readpipe("blastall -p blastn -U -F F -i $outdir/${prefix}.fa -d $outdir/${prefix_next}.fa -o $outfile -m 8 -e $evalue");
   }else{die "unknow method $method"}

   #filter m8 ?
     #readpipe("cat $outfile|perl bin/m8filter.pl ...");

   #parse blastm8  into hash table
   my $m8_tab=&readBlastm8($outfile); 

   return ($m8_tab, $outfile, "$outdir/${prefix}.fa", "$outdir/${prefix_next}.fa");
}


sub readBlastm8(){
   #blast m8 (12 columns):
   ## qid sid identity aligned_length mismatch gap qstart qend sstart send e-value bit-score
   ##   the columns are separated by tab
   #
   ##MSPcrunch(ACT, 8 columns):
   ## score identity qstart qend  qid sstart send sid
   ##
   ##  the columns are separated by spaces
   my $file=shift;
   open IN, $file or die"$!";   
   my $cnt=0;
   my %tab;
   while(<IN>){
     chomp;
     next if($_ eq "" || $_=~/^#/);
     $cnt++;
     my ($qid,$sid,$identity,$aligned_length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$e_value,$bit_score)=split/[\t| ]+/,$_;
     $tab{$cnt}="$bit_score\t$identity\t$qstart\t$qend\t$qid\t$sstart\t$send\t$sid";
   }
   close IN;
   return \%tab;
}


sub drawBlastm8(){
 my ($tab,$y,$col, $region,$region_next)=@_;
 my $h=50;
 #my($queChr,$queStart,$queEnd)=split/:|-/,$region;
 #my($refChr,$refStart,$refEnd)=split/:|-/,$region_next;
 foreach my $id(keys %{$tab}){
   my ($bit_score,$identity,$qstart,$qend,$qid,$sstart,$send,$sid)=split/\t/,$tab->{$id};
=pod
   $svg->rect("x",$x_left+($refStart-$left_up)*$ppx,"y",$y-5,"height",$barH,"width",($refEnd-$refStart+1)*$ppx,"fill","grey");
   if($refEnd-$refStart+1>300){
      for($refStart+1..$refEnd){
         if(($_-$refStart)%300==0){
            my($x1,$y1,$x2,$y2,$x3,$y3)=($_,$y,$_-200,$y-5,$_-200,$y+5);
            $svg->line("x1",$x_left+($x1-$left_up)*$ppx,"y1",$y1,"x2",$x_left+($x2-$left_up)*$ppx,"y2",$y2,"stroke","white","stroke-dasharray","1,1");
            $svg->line("x1",$x_left+($x1-$left_up)*$ppx,"y1",$y1,"x2",$x_left+($x3-$left_up)*$ppx,"y2",$y3,"stroke","white","stroke-dasharray","1,1");
         }
      }
   }
   $svg->rect("x",$x_left+($queStart-$left_down)*$ppx,"y",$y+$space-5,"height",$barH,"width",($queEnd-$queStart+1)*$ppx,"fill","grey");
   if($queEnd-$queStart+1>300){
      for($queStart+1..$queEnd){
         if(($_-$queStart)%300==0){
            my($x1,$y1,$x2,$y2,$x3,$y3)=($_,$y+60,$_-200,$y+60-5,$_-200,$y+60+5);
            $svg->line("x1",$x_left+($x1-$left_down)*$ppx,"y1",$y1,"x2",$x_left+($x2-$left_down)*$ppx,"y2",$y2,"stroke","white","stroke-dasharray","1,1");
            $svg->line("x1",$x_left+($x1-$left_down)*$ppx,"y1",$y1,"x2",$x_left+($x3-$left_down)*$ppx,"y2",$y3,"stroke","white","stroke-dasharray","1,1");
         }
     }
  }
=cut
  my @x=($space+$qstart*$ppx, $space+$qend*$ppx, $space+$send*$ppx, $space+$sstart*$ppx);
  my @y=($y+5,$y+5,$y+$h-5,$y+$h-5);
  my $points=$svg->get_path(
                     x=>\@x,
                     y=>\@y,
                     '-type'=>"polygon",
  );
  $svg->polygon(
                       %{$points},
                       'style'=>{
                        #"fill"=>"green",
                        "fill"=>$col,
                        #"fill"=>"#694d9f",
                        "fill-opacity"=>1,
                        #"fill-opacity"=>0.6,
                        #"fill-opacity"=>0.15,
                       }
               );
  
 }#foreach m8 tab
 return $y+$h+2;
}


sub svgWrite(){
  my ($svg,$outf)=@_;
  print STDERR "generating svg and output..\n";
  die "check outfile $outf" if(!defined $outf || $outf eq "");
  open OUT, ">$outf" or die "$!";
  print  OUT $svg->xmlify();  #don't use stdout if _exit(0) has been called, use FH instead or output file will be truncted. (don't know why)
  close OUT;
  print STDERR "write file done\n";
  return 1;
}



