use strict;
use warnings;
use Data::Dumper;
use SVG;
#usage: perl thisfile.pl genome.info orthTab
#Note:  chrs order and ref_species and ref_chr will be guessed  from orthTab file,
#       genome.info file is species_order insensitive

#v1.1 modification: use orthData instead of orthTab+FPKM




my ($info,$order)=&readInfo($ARGV[0]);
#my ($tab,$archi)=&readOrthTab($ARGV[1]);
my ($tab,$archi)=&readOrthData($ARGV[1]);
#my ($fpkm1,$fpkm2)=&readFPKM($ARGV[2]);
my ($refSpecies,$refChr)=split/:/,$archi->[0]; #orthTab first group as ref
#print STDERR Dumper $archi; exit;
#print Dumper $info,$tab,$archi; exit;
#print STDERR Dumper $fpkm1,$fpkm2;

my $space=100;
my $barH=10;
my ($height,$width)=(250,800);
#my ($height,$width)=(250,1600);
my $num=scalar @{$archi};
my $dist_chr=$height/$num;  #y   
my $len=&getLenAll($info,$archi);
#print STDERR "%{$info}\n@{$archi}\n$num\n@{$len}\n";
my $maxLen=&max($len);
my $ppx=$width/$maxLen; # x 

print STDERR "Guessed and calculated from orthoData file $ARGV[1]: $num species, maxLen=$maxLen, with ppx=$ppx; for @{$archi}\n ";
my $svg=SVG->new('height'=>$height+$space,'width'=>$width+2*$space);
#1, draw overall architechture
  my $y_next=$space;
  foreach(@{$archi}){
     my($species,$chr)=split/:/,$_;
     my $len_chr=shift @{$len};
     $svg->rect("x",$space,"y",$y_next,"height",$barH,"width",$len_chr*$ppx,"rx",8,"ry",8,"fill","white","stroke","none","stroke-width",2,"fill-opacity",0.2);
     $svg->text("x",$space/4,"y",$y_next+5,"height",$barH,"width",$len_chr*$ppx,"-cdata",$species."|".$refChr,);
     $y_next+=$dist_chr;
  } 

#2, draw anchors & links and fpkm
  my $maxFPKM=100;
  my $boxH=30;
  foreach(@{$tab}){
    $y_next=$space;
    my $y_fpkm=$space-20;
    my (@mids_x,@mids_y);
    my @anchors=split/\t/,$_; 
    my $flag;
    foreach(@anchors){
      $flag++;
      my ($species,$chr,$s,$e,$id,$strand,$v)=split/\|/,$_; #print STDERR "$species,$chr,$s,$e,$id,$strand\n";      
      if($species eq "-" || $chr eq "NONE" ){$y_next+=$dist_chr;next}
      #if($species eq "-" || $chr eq "NONE" || $chr ne $refChr){$y_next+=$dist_chr;next}
      $svg->rect("x",$space+$s*$ppx,"y",$y_next,"height",$barH,"width",($e-$s)*$ppx,"fill","black","stroke","none");

      my $FPKM=$v;
      #if(exists $fpkm1->{$id}){
      #  $FPKM=$fpkm1->{$id}
      #}elsif(exists $fpkm2->{$id}){
      #    $FPKM=$fpkm2->{$id}
      #  }
      my $height_fpkm=($FPKM>$maxFPKM)?($boxH):($FPKM*$boxH/$maxFPKM);
      $svg->rect("x",$space+($s+($e-$s)*0.5)*$ppx,"y",$y_fpkm-$height_fpkm,"height",$height_fpkm,"width",1,"fill","green","stroke","none","stroke-width",2);


      push @mids_x,$space+($s+($e-$s)/2)*$ppx;
      if($flag % 2!=0){push @mids_y,$y_next+$barH}else{push @mids_y,$y_next}
      $y_next=$y_next+$dist_chr;
      $y_fpkm+=$dist_chr+80;
    }

    my $points=$svg->get_path(
       "x" => \@mids_x,
       "y" => \@mids_y,
       "-type" => "polyline",
       "-closed" => "false",
    );
    $svg->polyline(
        %$points,       
        "style" => {
                     "stroke" => "#777777",
                     "stroke-opacity" => "0.5",
                     #"stroke" => "none",
                     "fill" => "none"
                     #"fill" => "#777777"
                   }
    );

  }

#3,draw feature elements

#=pod
## tigr6, bracRS2
#heterochromatin annotation and centremere positons
my $y_feature=$space-5;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   my $bed_corehet=&readBed($info->{$species}->{"coreHet"})->{$chr};
   my $bed_cent=&readBed($info->{$species}->{"cent"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed_corehet}){
    my($chr_in,$s,$e,$id)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","black","stroke","none","stoke-width",0.1);
       #$svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","black","stroke","black","stoke-width",0.1,"stroke-opacity",0.1);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   foreach(@{$bed_cent}){
     my($chr_in,$s,$e,$id)=split/\t/,$_;
     #my ($cents,$cente)=($s,$e);
     #my $centm=$cents+0.5*($cente-$cents);
     #$cents=$space+$cents*$ppx-5;
     #$cente=$space+$cente*$ppx+5;
     #$centm=$space+$centm*$ppx;
     #my $centy1=$y_feature+0.8*$barH+8;
     #my $centy2=$centy1-8;
     #$svg->ellipse("cx",$centm,"cy",$y_feature+0.65*$boxH,"rx",3,"ry",5,"fill","orange");
     #die "cent chr diffs" if($chr ne $chr_in);

     if($chr eq $chr_in){
        #print STDERR "$chr_in,$s,$e,$id?\n";
        $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","orange","stroke","none","stoke-width",0.1);
        #$svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","black","stroke","black","stoke-width",0.1,"stroke-opacity",0.1);
     }else{die "chr not eq: $chr with $chr_in\n"}


   }

   $y_feature+=$dist_chr+15;
   #last;
}
#=cut

=pod
##anchor blocks tigr6, brachy1.5
my $y_feature=$space-5;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"ancBlock"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$id,$v1,$v2)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","blue","stroke","none","stoke-width",0.1);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+15;
   #last;
}
=cut


##LTR tigr6, brachy1.5
 $y_feature=$space-10;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"teBed_LTR"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$id)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","blue","stroke","none","stoke-width",0);
       #$svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","blue","stroke","blue","stoke-width",0.5);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+25;
   #last;



}

##k92r1 tigr6, brachy1.5
$y_feature=$space-15;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"k92r1"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$id)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","red","stroke","none");
       #$svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","red","stroke","red","stroke-width",0.5);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+35;
   #last;
}

##k92r2 tigr6, brachy1.5(tilling stage)
$y_feature=$space-20;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"k92r2"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$id)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","red","stroke","none","stroke-width",0);
       #$svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","red","stroke","red","stroke-width",0.5);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+45;
   #last;
}


=pod
##chromHMM segs A/B tigr6 brac1.5
$y_feature=$space-25;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"chromHMM"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$type)=split/\t/,$_;
    if($chr eq $chr_in){
       my $color;
       #if($type eq "A"){$color="green"}elsif($type eq "B"){$color = "red"}else{$color="grey"}
       if($species eq "japo"){if($type eq "B"){$color="red"}else{$color="none"}}elsif($species eq "brac"){if($type eq "B" || $type eq "A"){$color="red"}else{$color="none"}}

       #($type eq "A")?($color="green"):($color="red");
       #print STDERR "$chr_in,$s,$e,$id?\n";
       next if($type eq "A");
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill",$color,"stroke",'none',"stroke-width",0.1);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+55;
   #last;
}


##LHRs tigr6
$y_feature=$space-30;
foreach my $str(@{$archi}){
   my ($species,$chr)=split/:/,$str;
   #print STDERR "\n$species,$chr,",$info->{$species}->{"teBed"},"\n";
   my $bed=&readBed($info->{$species}->{"LHRs"})->{$chr};
   #print STDERR Dumper $bed;last;
   foreach(@{$bed}){
    my($chr_in,$s,$e,$id)=split/\t/,$_;
    if($chr eq $chr_in){
       #print STDERR "$chr_in,$s,$e,$id?\n";
       $svg->rect("x",$space+$s*$ppx,"y",$y_feature,"height",$barH*0.4,"width",($e-$s)*$ppx,"fill","navy","stroke","navy","stroke-width",0.5);
     }else{die "chr not eq: $chr with $chr_in\n"}
   }
   $y_feature+=$dist_chr+65;
   last;
}
=cut



#=pod
#4, draw block shadows
#format :  japo|H1_japo|1230100|1310000	brac|H1_brac|748297|852229	brac|H1_brac|548297|652229
#          japo|H1_japo|1430100|1560000	brac|H1_brac|548997|656229	brac|H1_brac|848297|852229

#my @region=("japo|chr04|1700000|4500000\tbrac|chr04|944841|1815604","japo|chr04|28601378|30902378\tbrac|chr04|16056638|18341819");
#my @region=("japo|chr04|1700000|4500000\tbrac|chr04|944841|1815604","japo|chr04|30701378|30902378\tbrac|chr04|18056638|18241819");
my @region=("japo|chr04|1700000|4500000\tbrac|Chr04|1021209|1930025","japo|chr04|30701378|30902378\tbrac|Chr04|18140667|18321175");

#print STDERR "@test\n";
my ($ref_spec,$ref_chr)=split/:/,$archi->[0];
my ($que_spec,$que_chr)=split/:/,$archi->[1];
foreach my $ctrl(@region){
  my (@x,@y);
  my $y_next=$space;
  my @species=split/\t/,$ctrl;
  my $flag = 1;    
  foreach my$ctrl_in(@species){
    my($species,$chr,$s,$e)=split/\|/,$ctrl_in;
    if($chr ne $ref_chr && $chr ne $que_chr){$flag = 0;last};

    unshift @x,$space+$s*$ppx;
    push @x, $space+$e*$ppx;
    unshift @y, $y_next;
    push @y, $y_next;
    
    unshift @x,$space+$s*$ppx;
    push @x, $space+$e*$ppx;
    unshift @y, $y_next+$barH;
    push @y, $y_next+$barH;

    $y_next+=$dist_chr;
  }

  if($flag == 1){
    print STDERR "draw box $ctrl\n";
    my $points=$svg->get_path(
         "x" => \@x,
         "y" => \@y,
         "-type" => "polyline",
         "-closed" => "true",
    );
    $svg->polyline(
          %$points,
          "style" => {
                       "stroke" => "none",
                       "fill" => "navy",
                       "fill-opacity" => 0.5
                     }
    );
  }


}
#=cut

print $svg -> xmlify();



###sub###
sub readOrthTab(){  #standard 6*n column orthologus table format, '- or NONE' if not exists
                #filter anchors if chr ne chr_ref to NONE
                # guess and return species_num, species_order&species_chrs
  my $file=shift;
  my (@anchors,%info,@archi);
  open TAB, $file or die "$!";
  while(<TAB>){
    chomp;
    next if($_ eq "" || $_ =~/^#/);
    my @box=split/\t/,$_;    
    die "err at line $_\n" if(scalar @box %6 != 0);
    my $cnt;
    my $str1="";
    my $str2="";
    foreach(@box){
      $cnt++;
      if($cnt%6==0){$str1.=$_."\t"}else{$str1.=$_."|"}
      if($cnt%6==1){$str2.=$_}
      if($cnt%6==2){$str2.=":".$_."\t"}
    }
    push @anchors,$str1;
    $info{$str2}++;  #include all possible architechture
  }
  close TAB;
  my @temp=sort{ $info{$b} <=> $info{$a}} keys %info;  #guess the most frequency one
  @archi=split/\t/,$temp[0];
  return \@anchors,\@archi;
}


sub readOrthData(){  #standard 7*n (6+1) column orthologus data format, '- or NONE' if not exists 
                #filter anchors if chr ne chr_ref to NONE
                # guess and return species_num, species_order&species_chrs
  my $file=shift;
  my (@anchors,%archi,@archi,%value);
  open TAB, $file or die "$!";
  while(<TAB>){
    chomp;
    next if($_ eq "" || $_ =~/^#/);
    my @box=split/[\t ]+/,$_;    
    die "err at line $_\n" if(scalar @box %7 != 0);
    my $cnt;
    my $str1="";
    my $str2="";
    my $str3="";
    foreach(@box){
      $cnt++;
      if($cnt%7==0 ){if($cnt != scalar @box){$str1.=$_."\t"}else{$str1.=$_}}else{$str1.=$_."|"}
      if($cnt%7==1){$str2.=$_}
      if($cnt%7==2){$str2.=":".$_."\t"}
    }
    push @anchors,$str1;
    $archi{$str2}++;  #include all possible architechture
  }
  close TAB;
  my @temp=sort{ $archi{$b} <=> $archi{$a}} keys %archi;  #guess the most frequency one
  @archi=split/\t/,$temp[0];
  return \@anchors,\@archi;
}


sub readInfo_old(){ #genome info file
  my %info;
  #my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n/,$_;
   my $species=shift @box;
   #push @order,$species;
   foreach(@box){my ($type,$file)=split/[\t ]+/,$_;$info{$species}->{$type}=$file}
  }
  close INF;
  $/="\n";
  return \%info;
} 

sub readInfo(){ #genome info file
  my %info;
  my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n+/,$_;
   my $species=shift @box;
   push @order,$species;
   foreach(@box){
     my ($type,$file)=split/[\t ]+/,$_;
     next if($type =~/^"/);
     if(!exists $info{$species}->{$type}){ $info{$species}->{$type}=$file }else{die "dup $type in species $species\n"}
   }
  }
  close INF;
  $/="\n";
  return (\%info,\@order);
}



sub readLen(){
  my %len;
  my $file=shift @_;
  open LEN, $file or die "$!";
  while(<LEN>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,$len)=split/\t/,$_;
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
   my $sizef=$info->{$species}->{"chrSize"};
   #if ( -e $sizef ){}else{print "\n$species\t$sizef\n"; die "not exists $sizef\n"}
   if($chr eq "NONE"){$chr = $refChr}
   my $len=&readLen($sizef)->{$chr} ;  # here $_ will be unpredictable when call another sub in a sub
   push @len, $len;
  }
  return \@len;
}


sub max(){
  my $list=shift;
  my $max=0;
  foreach(@{$list}){
    if($_ > $max){$max = $_}
  }
  return $max;
}

sub readFPKM(){
   my $file=shift;
   my %fpkm1;
   my %fpkm2;
   open IN, $file or die "$!";
   while(<IN>){
    chomp;
    #print STDERR "$_\n";
    next if($_ eq "" || $_=~/^#/);
    my ($anch1,$fpkm1,$anch2,$fpkm2)=split/\t/,$_;
    #print STDERR "$anch1,$fpkm1,$anch2,$fpkm2\n";
    $fpkm1{$anch1}=$fpkm1;
    $fpkm2{$anch2}=$fpkm2;
    #if(!exists $fpkm1{$anch1}){$fpkm1{$anch1}=$fpkm1}else{die "dup $anch1\n"}
    #if(!exists $fpkm2{$anch2}){$fpkm2{$anch2}=$fpkm2}else{die "dup $anch2\n"}
   }
   close IN;
   return \%fpkm1,\%fpkm2;

}

sub readBed(){
   my $file = shift;
   my %bed;
   open BED, $file or die"$! $file";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/); 
     my($chr,$s,$e,$id)=split/[\t ]+/,$_;
     push @{$bed{$chr}},"$chr\t$s\t$e\t$id";

   }
   close BED;
   return \%bed;
}



sub readBed_cntid(){
   my $file = shift;
   my %bed;
   open BED, $file or die"$! $file";
   my $cnt;
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/); 
     $cnt++;
     my($chr,$s,$e)=split/[\t ]+/,$_;
     push @{$bed{$chr}},"$chr\t$s\t$e\tid$cnt";

   }
   close BED;
   return \%bed;
}
