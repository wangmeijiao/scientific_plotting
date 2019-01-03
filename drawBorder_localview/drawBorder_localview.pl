use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use SVG;
use POSIX qw[ _exit ];

#Dec 15, 2016
#tailored for border view

#Oct 31 2016
# For short region blast detail view (like act + igv)
# read whole data file and extract data when svg drawing, region : "japo|H1_japo|xxxx-xxxx,zs97|H1_zs97|xxxx-xxxx"

my ($info_file,$blast,$species,$region,$line,$outfile); #  use $species to choose species and order for drawing
                                
GetOptions("info=s",\$info_file,"blast=s",\$blast,"species=s",\$species,"line=s",\$line,"region=s",\$region,"out=s",\$outfile);
my @species=split/-/,$species;
my $n_spec=scalar @species;


my %region;
my @temp = split/,/,$region;
my $maxLen;
my $n_spec_r;
foreach my $ctrl(@temp){
   my ($spec,$chr,$range) = split/\|/,$ctrl;
   if(!exists $region{$spec}){$region{$spec}="$chr:$range"}else{die"dup spec $ctrl"}
   my ($s,$e) = split/-/,$range;
   if(!defined $maxLen){$maxLen = abs($e-$s+1)}else{if($maxLen < abs($e-$s+1)){ $maxLen = abs($e-$s+1)}}
   $n_spec_r++;
}
die "species and region number diffs" if($n_spec != $n_spec_r);

my %line;
my @temp1 = split/,/,$line;
my $n_spec_r1;
foreach my $ctrl(@temp1){
   my ($spec,$chr,$line_x) = split/\|/,$ctrl;
   if(!exists $line{$spec}){$line{$spec}=$line_x}else{die"dup spec $ctrl"}
   $n_spec_r1++;
}
die "species and region number diffs" if($n_spec != $n_spec_r1);
#print STDERR Dumper %line;

my ($info,$order)=&readInfo($info_file);
#print Dumper $info,$order; exit;


##0.initialization and global plot parameters
my $space=100;
my ($width,$height)=(800,300*$n_spec);
my $svg=SVG->new("width",$width+2*$space,"height",$height+2*$space);
my %color=(
            "H3K9me2"=>"#840228",
            "k92"=>"#840228",
            "H3K4me3"=>"green",
            "k43"=>"green",

             "mappability" => "grey",
             "NPS" => "blue",             
             "DHsite" => "grey",
 
            "gene"=>"green",
            "mrna"=>"green",
            "mRNA"=>"green",
            "siRNA"=>"purple",

            "DNATE"=>"orange",
            "MITE"=>"orange",
            "LINE"=>"#228fbd",
            "SINE"=>"#008792",
            "RC"=>"#ea66a6",
            "LTR"=>"navy",
            "otherTE" => "#896a45",
            #"otherTE" => "purple",
            "gap" => "black",

            #"link" => "#6950a1",
            "link" => "#694d9f",
            #"link" => "#444693",
            #"link" => "black",
            #"link" => "red",
            #"link_rev" => "#6950a1",
            "link_rev" => "#694d9f",
            #"link_rev" => "#444693",
            #"link_rev" => "#585eaa",
            #"link_rev" => "blue",

             "CHH" => "#FF1FDF",
             #"CHH" => "pink",
             "CHG" => "#ef4136",
             "CG" => "#121a2a",
             #"CG" => "#ed1941",

             "1" => "navy",
             "2" => "navy",
             "3" => "navy",
             "4" => "navy",
             "5" => "navy",
             "6" => "navy",
             "7" => "navy",
             "8" => "navy",
             "9" => "navy",


          );

my %range=(

           "japo" => {"k43","0-0-5","k92","2-2-35", "mrna", "0-0-25" ,"siRNA", "0-0-15","CHH","0-0-100","CHG","0-0-100","CG","0-0-100","mappability","0-0-1","nps","0-0-25","DHsite","0-0-35"},
           "9311" => {"k43","0-0-15","k92","2-2-35", "mrna", "0-0-25" ,"siRNA", "0-0-15","CHH","0-0-100","CHG","0-0-100","CG","0-0-100","mappability","0-0-1","nps","0-0-25","DHsite","0-0-35"},
           #"japo" => {"k43","1-1-45","k92","2-2-55", "mrna", "0-0-25" ,"siRNA", "0-0-25","CHH","0-0-100","CHG","0-0-100","CG","0-0-100","mappability","0-0-1","nps","0-0-25","DHsite","0-0-15"},
           "glab" => {"k43","0-0-5","k92","1-1-25", "mrna", "0-0-25" ,"siRNA", "0-0-15","CHH","0-0-100","CHG","0-0-100","CG","0-0-100"},
           "punc" => {"k43","5-5-45","k92","2-2-55", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "brac" => {"k43", "5-5-45","k92","2-2-20", "mrna", "0-0-25" ,"siRNA", "0-0-10","CHH","0-0-100","CHG","0-0-100","CG","0-0-100"},
           "lper" => {"k43","5-5-45","k92","5-5-35", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "atau"   => {"k43","2-2-45","k92","2-2-35", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "bd"   => {"k43","2-2-45","k92","2-2-35", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "zea" => {"k43","1-1-15","k92","0-0-10", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "sorg" => {"k43","2-2-30","k92","2-2-15", "mrna", "0-0-25" ,"siRNA", "0-0-25"},
           "orop" => {"k43","5-5-45","k92","2-2-35", "mrna", "0-0-25" ,"siRNA", "0-0-25"},



         );


my %blast_score=(  #to filter blast links by score( act read blast m8 alignment length as score instead of the last column)
                   "japo" => 100,
                   #"japo" => 300, #for japo-brac
                   #"japo" => 100, #for japo-brac
                   "9311" => 100,
                   "zs97" => 242,
                   "glab" => 100,
                   #"brac" => 100,
                   


                );



my $ppx=$width/$maxLen;
my $y_next=$space;

my $dist_spec=100;



my $cnt=-1;
foreach my $spec(@species){
   print STDERR "\nplotting for $spec\n";
   $cnt++;
   my $len=&readLen($info->{$spec}->{"faSize"});
   my ($chr_r,$start_r,$end_r) = split/:|-/,$region{$spec};
   my $len_r = $end_r-$start_r+1;

   #1 data tracks
   #
   my $y_start_hilight = $y_next;

#=pod
   ## k43 and peak
   if(&isValid($info->{$spec}->{"k43_densi"})){
     print STDERR "#H3k4me3 track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"k43_densi"}),"none",$region{$spec},$y_next,"H3K4me3",$range{$spec}->{"k43"},"barplot")
   }
   ## mrna
   if(&isValid($info->{$spec}->{"mrna_densi"})){
     print STDERR "#mRNA track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"mrna_densi"}),"none",$region{$spec},$y_next,"mRNA",$range{$spec}->{"mrna"},"barplot" )
   }

   ## siRNA
   if(&isValid($info->{$spec}->{"siRNA_densi"})){
     print STDERR "#siRNA track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"siRNA_densi"}),"none",$region{$spec},$y_next,"siRNA",$range{$spec}->{"siRNA"},"barplot")
   }

#   if(&isValid($info->{$spec}->{"siRNA_densi_blakeleafM1"})){
#     print STDERR "#siRNA track\n";
#     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"siRNA_densi_blakeleafM1"}),"none",$region{$spec},$y_next,"siRNA",$range{$spec}->{"siRNA"},"barplot")
#   }

#   if(&isValid($info->{$spec}->{"siRNA_densi_blakeleafm1"})){
#     print STDERR "#siRNA track\n";
#     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"siRNA_densi_blakeleafm1"}),"none",$region{$spec},$y_next,"siRNA",$range{$spec}->{"siRNA"},"barplot")
#   }

#   if(&isValid($info->{$spec}->{"siRNA_densi_blakepanicleM1"})){
#     print STDERR "#siRNA track\n";
#     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"siRNA_densi_blakepanicleM1"}),"none",$region{$spec},$y_next,"siRNA",$range{$spec}->{"siRNA"},"barplot")
#   }




   ## DNA meth
   if(&isValid($info->{$spec}->{"CHH"})){
     print STDERR "#CHH track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"CHH"}),"none",$region{$spec},$y_next,"CHH",$range{$spec}->{"CHH"},"barplot")
   }

   if(&isValid($info->{$spec}->{"CHG"})){
     print STDERR "#CHG track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"CHG"}),"none",$region{$spec},$y_next,"CHG",$range{$spec}->{"CHG"},"barplot")
   }

   if(&isValid($info->{$spec}->{"CG"})){
     print STDERR "#CG track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"CG"}),"none",$region{$spec},$y_next,"CG",$range{$spec}->{"CG"},"barplot")
   }


   ## k92 and peak
   if(&isValid($info->{$spec}->{"k92_densi"})){
     print STDERR "#H3k9me2 track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"k92_densi"}),"none",$region{$spec},$y_next,"H3K9me2",$range{$spec}->{"k92"},"barplot")
   }

   #nps
   if(&isValid($info->{$spec}->{"nps"})){
     print STDERR "#nps track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"nps"}),"none",$region{$spec},$y_next,"NPS",$range{$spec}->{"nps"},"barplot")
   }

    #DHsite
   if(&isValid($info->{$spec}->{"DHsite"})){
     print STDERR "#DHsite track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"DHsite"}),"none",$region{$spec},$y_next,"DHsite",$range{$spec}->{"DHsite"},"barplot")
   }
 

   #mappability
   if(&isValid($info->{$spec}->{"mappability"})){
     print STDERR "#mappability track\n";
     $y_next=&drawBdg_slim(&readBdg($info->{$spec}->{"mappability"}),"none",$region{$spec},$y_next,"mappability",$range{$spec}->{"mappability"},"barplot")
   }
 
  
#=cut


   #2 gene LTR DNATE embl annotation tracks and coordinates
   #chromsome and coordinates
   my $y_plus = $y_next;
   $y_next=&chrPaint_modi(undef,$y_next,$spec,$len_r,1);
   $y_next=&drawCoords($y_next,$region{$spec});
   #$y_next=&drawCoords($y_next,$len_r);
   my $y_minus = $y_next;
   $y_next=&chrPaint_modi(undef,$y_next,$spec,$len_r,0);

   #gene embl
   if(&isValid($info->{$spec}->{"geneEmbl"})){
      print STDERR "#gene_embl ";
      $y_next=&drawEmbl_act_gene(&readEmbl_gene($info->{$spec}->{"geneEmbl"},"gene"),$region{$spec},$y_plus,$y_minus,"gene")
   }

   #LTR embl
   if(&isValid($info->{$spec}->{"LTR_embl"})){
     print STDERR "LTR_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_LTR($info->{$spec}->{"LTR_embl"},"note"),$region{$spec},$y_plus,$y_minus,"LTR")
   }

   #DNATE embl
   if(&isValid($info->{$spec}->{"DNATE_embl"})){
     print STDERR "DNATE_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"DNATE_embl"},"note"),$region{$spec},$y_plus,$y_minus,"DNATE")
   }
   #MITE embl
   if(&isValid($info->{$spec}->{"MITE_embl"})){
     print STDERR "MITE_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"MITE_embl"},"ID"),$region{$spec},$y_plus,$y_minus,"MITE")
   }


   #LINE embl
   if(&isValid($info->{$spec}->{"LINE_embl"})){
     print STDERR "LINE_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"LINE_embl"},"ID"),$region{$spec},$y_plus,$y_minus,"LINE")
   }

   #SINE embl
   if(&isValid($info->{$spec}->{"SINE_embl"})){
     print STDERR "SINE_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"SINE_embl"},"ID"),$region{$spec},$y_plus,$y_minus,"SINE")
   }

   #RC embl
   if(&isValid($info->{$spec}->{"RC_embl"})){
     print STDERR "RC_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"RC_embl"},"ID"),$region{$spec},$y_plus,$y_minus,"RC")
   }



   # other TE embl
   if(&isValid($info->{$spec}->{"otherTE_embl"})){
     print STDERR "otherTE_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"otherTE_embl"},"ID"),$region{$spec},$y_plus,$y_minus,"otherTE")
   }

   # gap
   if(&isValid($info->{$spec}->{"gap_embl"})){
     print STDERR "gap_embl ";
     $y_next=&drawEmbl_act_repeat(&readEmbl_DNATE($info->{$spec}->{"gap_embl"},"note"),$region{$spec},$y_plus,$y_minus,"gap")
   }



   print STDERR "\n";

   my $y_end_hilight = $y_next;


   #3 border line
   if(exists $line{$spec}){
      #my $regionBed = &readBed($info->{$spec}->{"regionBar"});
      #&drawHiLightBox($regionBed,$y_start_hilight,$y_end_hilight);
      #$y_next=&chrPaint($regionBed,$y_next,$spec,$len->{"H1_$spec"})
      my ($chr_r,$s_r,$e_r) = split/:|-/,$region{$spec};
      my $pos = $line{$spec};
      die "border line out of range:$pos $s_r $e_r at $spec" if($pos <= $s_r || $pos >= $e_r);
      $pos-=$s_r;
      $svg->line("x1",$space+$pos*$ppx,"y1",$y_start_hilight,"x2",$space+$pos*$ppx,"y2",$y_end_hilight,style=>{"stroke-width",2,"stroke","black","stroke-dasharray","5,5"});

   }


#=pod
   #4 draw blast m8
   if($spec ne $species[-1]){
    my $spec_next=$species[$cnt+1];
    
    my $prefix = $info->{$spec}->{"fa"};
    $prefix=~s/\.fa$//;
    my $prefix_next = $info->{$spec_next}->{"fa"};
    $prefix_next=~s/\.fa$//;
    my $m8_file=$prefix."_VS_".$prefix_next.".blast";

#    my $chr = "H1_$spec";
#    my $chr_next = "H1_$spec_next";
#    my $m8_file=$chr."_VS_".$chr_next.".blast";   #H1_japo-VS-H1_niva.blast

    $m8_file=$blast."/".$m8_file;
    print STDERR "$m8_file\n\n";
    $y_next=&drawBlastm8(&readBlastm8($m8_file),$region{$spec},$region{$spec_next},$y_next,$dist_spec,$blast_score{$spec})
   }
#=cut
   #$y_next+=2*$dist_spec;

}#while end



&svgWrite($svg,$outfile);
POSIX::_exit(0); #to speed up exit by ommitting perl_destruct for big hash

#print STDERR "generating svg..\n";
#print $svg->xmlify();
#print STDERR "done\n";



##sub
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



sub drawBed(){  #for annotation tracks only, not for peaks
  my($table,$y,$ftype)=@_;
  my $x_left=$space;
  my $barH=4;
  #draw one line anno

 #$svg->line("x1",$x_left,"y1",$y_te,"x2",$x_right,"y2",$y_te,"stroke","grey","stroke-width",0.5);#te border line
  #$svg->text("x",$x_left-$space*0.4,"y",$y+6,"font-family","Arial-Narrow","font-size",8,"font-weight","bold","text-anchor","start","-cdata",$ftype);
  $svg->text("x",$space-15,"y",$y+6,"width",10,"height",5,"stroke","black","-cdata",$ftype,"text-anchor","end","stroke","black");
  #$svg->text("x",$x_left-5,"y",$y_te-2,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","+");
  #$svg->text("x",$x_left-15,"y",$y_te-2,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","->");
  #$svg->text("x",$x_left-15,"y",$y_te+5,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata","<-");
  #$svg->text("x",$x_left-5,"y",$y_te+5,"font-family","Arial-Narrow","font-size",6,"font-weight","bold","text-anchor","start","-cdata"," -");

  foreach my $ctrl(@{$table}){
     my($chr,$start,$end,$name,undef,$strand)=split/\|/,$ctrl;
     #if($strand eq "+"){
     #   $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_te-6,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$feature_line[0]},"stroke",$color{$feature_line[0]}   );
     #}elsif($strand eq "-"){
     #   $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_te+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$feature_line[0]},"stroke",$color{$feature_line[0]}   );
     # }
     $svg->rect("x",$x_left+$start*$ppx,"y",$y+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$ftype},"stroke",$color{$ftype} );

  }
  return $y+6+$barH*2;
}



sub drawEmbl_act_gene(){
  my ($embl,$region,$y_plus,$y_minus,$type)=@_;
  #print STDERR Dumper $embl;
  my $y;
  my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
  my $h_track=20;
  #$svg->text("x",$space-15,"y",$y_plus+5,"width",10,"height",5,"stroke","black","-cdata","gene","text-anchor","end","stroke","black");
  #$svg->text("x",$space-15,"y",$y_minus+5,"width",10,"height",5,"stroke","black","-cdata","gene","text-anchor","end","stroke","black");
  foreach my $id(keys %{$embl}){
    if(exists $embl->{$id}->{"mRNA"}){
      my $strand=$embl->{$id}->{"mRNA"}->{"strand"};
      if($strand eq "+"){$y = $y_plus}elsif($strand eq "-"){$y = $y_minus}else{die "unknow strand $strand at $id in $type drawing"}
      my $color;
      if(exists $embl->{$id}->{"mRNA"}->{"colour"}){ $color = $color{$embl->{$id}->{"mRNA"}->{"colour"}}}else{$color = $color{$type}}
      #if(exists $embl->{$id}->{"mRNA"}->{"colour"}){ $color = $color{$embl->{$id}->{"mRNA"}->{"colour"}}}else{$color = $color{$type}}
      my $coord=$embl->{$id}->{"mRNA"}->{"coord"};
      my ($left,$right)=&getRange($coord); #print STDERR "$left,$right,@$coord\n";
      my $coord_sort = &segSort($coord);

      if(scalar @$coord_sort == 0){ die "empty segs at $id"}elsif(scalar @$coord_sort == 1){
         my ($s,$e)=split/\.\./,$coord_sort->[0];
         next if($e < $s_r || $s > $e_r);
         $s-=$s_r; if ($s < 0){ $s = 0}
         $e-=$s_r; if ($e < 0){ $e = 0}
         $svg->rect("x",$space+$s*$ppx,"y",$y, "width",($e-$s+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5});
         $svg->text("x",$space+(0.5*($s+$e))*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle"});
      }elsif(scalar @$coord_sort >= 2){
        foreach (my $i=0;$i<=$#$coord_sort-1;$i++){
          my ($s1,$e1)=split/\.\./,$coord_sort->[$i];
          next if($e1 < $s_r || $s1 > $e_r);
          $s1-=$s_r; if ($s1 < 0){ $s1 = 0}
          $e1-=$s_r; if ($e1 < 0){ $e1 = 0}
          $svg->rect("x",$space+$s1*$ppx,"y",$y, "width",($e1-$s1+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5});
          my ($s2,$e2)=split/\.\./,$coord_sort->[$i+1];
          next if($e2 < $s_r || $s2 > $e_r);
          $s2-=$s_r; if ($s2 < 0){ $s2 = 0}
          $e2-=$s_r; if ($e2 < 0){ $e2 = 0}
          if($e1 > $s2){die "e1 > s2 (segments overlap): $e1 > $s2  at $id in $type drawing"}
          if($s1 > $s2){die "s1 > s2 (segments wrong sorted): $s1 > $s2  at $id in $type drawing"}
          #print STDERR "$id: $s1-$e1 -> $s2-$e2\n";
          $svg->rect("x",$space+$s2*$ppx,"y",$y, "width",($e2-$s2+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5});

          ###transcript links
          my $mid = 0.5*($s2+$e1)+1;
          $svg->line("x1",$space+$e1*$ppx,"y1",$y+0.5*$h_track,"x2",$space+$mid*$ppx,"y2",$y,style=>{"stroke-width",0.8,"stroke",$color});       
          $svg->line("x1",$space+$mid*$ppx,"y1",$y,"x2",$space+$s2*$ppx,"y2",$y+0.5*$h_track,style=>{"stroke-width",1,"stroke",$color});       
          ### gene id  only draw feature length > 10 ppx 
          next if($i>0);
          $svg->text("x",$space+(0.5*($left+$right)-$s_r)*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle"});
          #if(abs($right-$left)*$ppx >= 4 ){$svg->text("x",$space+(0.5*($left+$right)-$s_r)*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle","stroke","black"})}
       }#foreach segs end
     }# segs >= 2

     ###transcript arrow
      next if($left < $s_r || $right > $e_r);
      $left-=$s_r; if ($left < 0){ $left = 0}
      $right-=$s_r; if ($right < 0){ $right = 0}
     if($strand eq "+"){
        #print STDERR "$id: $strand  $left,$right  @$coord_sort\n";
        &drawTriAngle_fill($y,$right,5,$h_track,"grey",0);
     }elsif($strand eq "-"){
         #print STDERR "$id: $strand  $left,$right   @$coord_sort\n";
         &drawTriAngle_fill($y,$left,5,$h_track,"grey",1);
        }else{die "unknow strand: $strand at $id"}

    }else{die "not found mRNA record at $id"}
  }#foreach id
  return $y_minus+$h_track;
}#sub end



sub drawEmbl_act_repeat(){
  my ($embl,$region,$y_plus,$y_minus,$type)=@_;
  #if($type eq "DNATE"){print STDERR Dumper $embl};
  my $y;
  my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
  my $h_track=20;
  #$svg->text("x",$space-15,"y",$y_plus+5,"width",10,"height",5,"stroke","black","-cdata","gene","text-anchor","end","stroke","black");
  #$svg->text("x",$space-15,"y",$y_minus+5,"width",10,"height",5,"stroke","black","-cdata","gene","text-anchor","end","stroke","black");
  foreach my $id(keys %{$embl}){
      my $strand=$embl->{$id}->{"strand"};
      if($strand eq "+"){$y = $y_plus}elsif($strand eq "-"){$y = $y_minus}else{die "unknow strand $strand at $id in $type drawing"}
      my $color;
      if(exists $embl->{$id}->{"colour"}){ $color = $color{$embl->{$id}->{"colour"}}}else{$color = $color{$type}}
      if(exists $embl->{$id}->{"color"}){ $color = $color{$embl->{$id}->{"color"}}}else{$color = $color{$type}}
      #if(exists $embl->{$id}->{"mRNA"}->{"colour"}){ $color = $color{$embl->{$id}->{"mRNA"}->{"colour"}}}else{$color = $color{$type}}
      my $coord=$embl->{$id}->{"coord"};
      my ($left,$right)=&getRange($coord); #print STDERR "$left,$right,@$coord\n";
      my $coord_sort = &segSort($coord);

      if(scalar @$coord_sort == 0){ die "empty segs at $id"}elsif(scalar @$coord_sort == 1){
         my ($s,$e)=split/\.\./,$coord_sort->[0];
         next if($e < $s_r || $s > $e_r);
         $s-=$s_r; if ($s < 0){ $s = 0}
         $e-=$s_r; if ($e < 0){ $e = 0}
         $svg->rect("x",$space+$s*$ppx,"y",$y, "width",($e-$s+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5,"fill-opacity",0.8});
         if(abs($e-$s)*$ppx >= 2 ){$svg->text("x",$space+(0.5*($s+$e))*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle"})};
      }elsif(scalar @$coord_sort >= 2){
        foreach (my $i=0;$i<=$#$coord_sort-1;$i++){
          my ($s1,$e1)=split/\.\./,$coord_sort->[$i];
          next if($e1 < $s_r || $s1 > $e_r);
          $s1-=$s_r; if ($s1 < 0){ $s1 = 0}
          $e1-=$s_r; if ($e1 < 0){ $e1 = 0}
          $svg->rect("x",$space+$s1*$ppx,"y",$y, "width",($e1-$s1+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5,"fill-opacity",0.8});
          my ($s2,$e2)=split/\.\./,$coord_sort->[$i+1];
          next if($e2 < $s_r || $s2 > $e_r);
          $s2-=$s_r; if ($s2 < 0){ $s2 = 0}
          $e2-=$s_r; if ($e2 < 0){ $e2 = 0}
          if($e1 > $s2){die "e1 > s2 (segments overlap): $e1 > $s2  at $id in $type drawing"}
          if($s1 > $s2){die "s1 > s2 (segments wrong sorted): $s1 > $s2  at $id in $type drawing"}
          #print STDERR "$id: $s1-$e1 -> $s2-$e2\n";
          $svg->rect("x",$space+$s2*$ppx,"y",$y, "width",($e2-$s2+1)*$ppx,"height",$h_track, style=>{"fill",$color,"stroke","black","stroke-width",0.5,"fill-opacity",0.8});

          ###transcript links
          my $mid = 0.5*($s2+$e1)+1;
          $svg->line("x1",$space+$e1*$ppx,"y1",$y+0.5*$h_track,"x2",$space+$mid*$ppx,"y2",$y,style=>{"stroke-width",0.8,"stroke",$color});       
          $svg->line("x1",$space+$mid*$ppx,"y1",$y,"x2",$space+$s2*$ppx,"y2",$y+0.5*$h_track,style=>{"stroke-width",1,"stroke",$color});       
          ### gene id  only draw feature length > 10 ppx 
          #if(abs($right-$left)*$ppx >= 4 ){$svg->text("x",$space+(0.5*($left+$right)-$s_r)*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle","stroke","black"})}
          next if($i>0);
          $svg->text("x",$space+(0.5*($left+$right)-$s_r)*$ppx,"y",$y+0.25*$h_track,"-cdata",$id,style=>{"text-anchor","middle"});

       }#foreach segs end
     }# segs >= 2

     ###transcript arrow
      next if($left < $s_r || $right > $e_r);
      $left-=$s_r; if ($left < 0){ $left = 0}
      $right-=$s_r; if ($right < 0){ $right = 0}
     if($strand eq "+"){
        #print STDERR "$id: $strand  $left,$right  @$coord_sort\n";
        &drawTriAngle_fill($y,$right,5,$h_track,"grey",0);
     }elsif($strand eq "-"){
         #print STDERR "$id: $strand  $left,$right   @$coord_sort\n";
         &drawTriAngle_fill($y,$left,5,$h_track,"grey",1);
        }else{die "unknow strand: $strand at $id"}

  }#foreach id
  return $y_minus+$h_track;
}#sub end



sub drawTriAngle_fill(){
   my ($y,$s,$len,$h, $c,$flip)=@_;
   $c||="black";
   $len||=5;
   my (@x,@y);
   if($flip == 0){
     @x=($space+$s*$ppx,$space+$s*$ppx,$space+$s*$ppx+$len);
     @y=($y,$y+$h,$y+0.5*$h)
   }elsif($flip == 1){
       @x=($space+$s*$ppx,$space+$s*$ppx,$space+$s*$ppx-$len);
       @y=($y,$y+$h,$y+0.5*$h)
     }
   #print STDERR "$s|@x|@y\n";
   ##get coords
   my $points=$svg->get_path(
                              x=>\@x,
                              y=>\@y,
                              -type=>"polygon",
                              #-closed=>'true',
   );
   #print STDERR Dumper $points;
   ##draw polygon
   $svg->polygon(
           %$points,
           style => {
                        #'opacity' => 0.4, 
                        'fill-opacity' => 0.8,
                        'fill' =>$c,
                        'stroke' => $c,
                        #'stroke' => $c,
                        #'stroke-opacity' => 1,
                    }
   );
   #draw id
   #if($e-$s > 20000){$svg->text("x",$space+$mid*$ppx,"y",$y+0.2*$h,"-cdata",$text,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","middle"})}

   return 1;

}



sub drawEmbl(){
  my ($file,$y)=@_;
  my $embl=&readEmbl($file,"gene");
  #print Dumper \$embl;
  my $h_track=10;
  $svg->text("x",$space-15,"y",$y+5,"width",10,"height",5,"stroke","black","-cdata","gene","text-anchor","end","stroke","black");
  foreach my $id(keys %{$embl}){
    if(exists $embl->{$id}->{"mRNA"}){
      my $strand=$embl->{$id}->{"mRNA"}->{"strand"};
      #my $name=$embl->{$id}->{"mRNA"}->{"name"};
      my $coord=$embl->{$id}->{"mRNA"}->{"coord"};
      my ($left,$right)=&getRange_1($coord);
      #$svg->line("x1",$space+$left*$ppx,"y1",$y+0.5*$h_track,"x2",$space+$right*$ppx,"y2",$y+0.5*$h_track,style=>{"stroke-width",1,"stroke","grey"});
      $svg->rect("x",$space+$left*$ppx,"y",$y+0.25*$h_track, "width",($right-$left+1)*$ppx,"height",0.5*$h_track, style=>{"fill","grey","stroke","grey"});
      foreach my$region(@$coord){
        my ($s,$e)=split/\.\./,$region;
        die "s >= e: $s >= $e at $id mRNA" if($s >= $e);
        $svg->rect("x",$space+$s*$ppx,"y",$y+0.25*$h_track, "width",($e-$s+1)*$ppx,"height",0.5*$h_track, style=>{"fill",$color{"gene"},"stroke",$color{"gene"}});
        ###transcript arrow
        if($strand eq "+"){
          $svg->line("x1",$space+$left*$ppx,"y1",$y+0.25*$h_track,"x2",$space+$left*$ppx,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
          $svg->line("x1",$space+$left*$ppx,"y1",$y,"x2",$space+$left*$ppx+5,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
          $svg->line("x1",$space+$left*$ppx+5,"y1",$y,"x2",$space+$left*$ppx+2,"y2",$y-2,style=>{"stroke-width",1,"stroke","black"});
          $svg->line("x1",$space+$left*$ppx+5,"y1",$y,"x2",$space+$left*$ppx+2,"y2",$y+2,style=>{"stroke-width",1,"stroke","black"});
        }elsif($strand eq "-"){
             $svg->line("x1",$space+$right*$ppx,"y1",$y+0.25*$h_track,"x2",$space+$right*$ppx,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
             $svg->line("x1",$space+$right*$ppx,"y1",$y,"x2",$space+$right*$ppx-5,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
             $svg->line("x1",$space+$right*$ppx-5,"y1",$y,"x2",$space+$right*$ppx-2,"y2",$y-2,style=>{"stroke-width",1,"stroke","black"});
             $svg->line("x1",$space+$right*$ppx-5,"y1",$y,"x2",$space+$right*$ppx-2,"y2",$y+2,style=>{"stroke-width",1,"stroke","black"});

           }else{die "unknow strand: $strand at $id"}
        ### gene id   
        #$svg->text("x",$space+0.5*($left+$right)*$ppx,"y",$y+1.25*$h_track,"-cdata",$id,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","middle"})
      }
    }else{die "not found mRNA record at $id"}
  }#foreach id
  return $y+$h_track+0.5*$h_track;
}#sub end



sub readEmbl_gene(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^mRNA/ || $_=~/^CDS/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{$feature{$k}.="-".$v}#die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions; 
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if($feature_name eq "mRNA"){
            $gene{$geneID}->{"mRNA"}=\%feature
          }elsif($feature_name eq "CDS"){
              $gene{$geneID}->{"CDS"}=\%feature
            }else{die "unknow feature $feature_name"}
          
          #inistialize for next chunk
          @block=();
          push @block, $_;          
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  return \%gene;

}


sub readEmbl_LTR(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^LTR/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          $cnt++;
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{ $feature{$k}.="-".$v}#$k.="_".&generate_random_string(4);$feature{$k}=$v}
            #if(!exists $feature{$k}){$feature{$k}=$v}else{die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions; 
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if(exists $gene{$geneID}){
             #die "dup LTR ID $geneID";
             $geneID.="_".&generate_random_string(4);
             $gene{$geneID}=\%feature;
          }else{
            $gene{$geneID}=\%feature;
          }
          #if($feature_name eq "mRNA"){
          #  $gene{$geneID}->{"mRNA"}=\%feature
          #}elsif($feature_name eq "CDS"){
          #    $gene{$geneID}->{"CDS"}=\%feature
          #  }else{die "unknow feature $feature_name"}
          
          #inistialize for next chunk
          @block=();
          push @block, $_;          
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  return \%gene;

}


sub readEmbl_DNATE(){
  my ($file,$keyword)=@_;
  my %gene;
  my $cnt;
  my @block;
  open EMBL, $file or die "$!";
  while(<EMBL>){
      chomp;
      next if ($_ eq "" || $_=~/^#/ || $_=~/^\s+$/ || !($_=~/^FT/)); #lines must start with FT
      $_=~s/^FT +//;
      if( $_=~/^repeat_region/ || $_=~/^gap/ || eof (EMBL) ){
        if(scalar @block == 0){push @block, $_}else{
          if(eof (EMBL)){push @block, $_}; #add final line
          #start to fill gene hash table
          my @box=split/[\t ]+/,(shift @block);
          my $feature_name=$box[0];
          my $feature_coord=$box[1];
          my ($strand,$start,$end,@regions)=&parseCoord($feature_coord);
          my $geneID;
          my %feature;
          foreach my $ctrl(@block){
            $ctrl=~s/^\///;
            my ($k,$v)=split/=/,$ctrl;
            $v=~s/\"//g;
            my @temp=split/:/,$v;
            $v=$temp[0];
            if($k eq $keyword){$geneID=$v;next} #mush has /ID=xxx or /gene=xxx
            if(!exists $feature{$k}){$feature{$k}=$v}else{ $feature{$k}.="-".$v}#die "dup key $k at $ctrl"}
          }
          $feature{"strand"}=$strand;
          $feature{"coord"}=\@regions; 
          die "can't find ID line at $_" if(!defined $geneID || $geneID eq "");
          if(exists $gene{$geneID}){die "dup LTR ID $geneID"}else{
           $gene{$geneID}=\%feature;
          }
          #if($feature_name eq "mRNA"){
          #  $gene{$geneID}->{"mRNA"}=\%feature
          #}elsif($feature_name eq "CDS"){
          #    $gene{$geneID}->{"CDS"}=\%feature
          #  }else{die "unknow feature $feature_name"}
          
          #inistialize for next chunk
          @block=();
          push @block, $_;          
        }
      }elsif($_=~/^\//){push @block, $_}else{die "unknown line $_"}

  }#while end
  close EMBL;
  #print Dumper \%gene;
  return \%gene;

}

sub generate_random_string{
    #http://stackoverflow.com/questions/13687643/generate-unique-random-strings
    my $length_of_randomstring = shift; # the length of 
                                        # the random string to generate

    my @chars=('a'..'z','A'..'Z','0'..'9','_');
    my $random_string;
    for(1..$length_of_randomstring){
        # rand @chars will generate a random 
        # number between 0 and scalar @chars
        $random_string.=$chars[rand @chars];
    }

    return $random_string;
}


sub parseCoord(){
#general code for parsing embl feature structure line
 my $string=shift;
 if($string=~/^join\((.+)\)$/){
        my @temp=split/\.\./,$1;
        return "+",$temp[0],$temp[$#temp],(split/,/,$1);
    }elsif($string=~/^[0-9]+\.\.[0-9]+$/){
         my @temp=split/\.\./,$string;
         return "+",$temp[0],$temp[$#temp],$string;
       }elsif($string=~/^complement\(join\((.+)\)\)$/){
             my @temp=split/\.\./,$1;
             return "-",$temp[0],$temp[$#temp],split/,/,$1;
          }elsif($string=~/^complement\(([0-9]+\.\.[0-9]+)\)$/){
             my @temp=split/\.\./,$1;
             return "-",$temp[0],$temp[$#temp],$1;
            }else{die  "err when parse coordination $string\n"}
}#sub end here

sub getRange(){
   #used in drawEmbl
   my $coord=shift;
   my ($left,$right);
   foreach my $ctrl(@$coord){
      my ($s,$e)=split/\.\./,$ctrl;
      die "start >= end: $s >= $e at $coord" if($s >= $e);
      if(!defined $left){$left=$s}else{if($left > $s){$left=$s}};
      if(!defined $right){$right=$e}else{if($right < $e){$right=$e}};
   }
   return ($left,$right);
}

sub segSort(){
  #sort by start, non overlap by default
  my $seg = shift;
  my @seg_sort = sort {
     my ($s1, $e1) = split/\.\./,$a;
     die "start >= end: $s1 >= $e1 " if($s1 >= $e1);
     my ($s2, $e2) = split/\.\./,$b;
     die "start >= end: $s2 >= $e2 " if($s2 >= $e2);
     if($s1 < $s2){return -1}elsif($s1 > $s2){return 1}else{return 0}

  } @$seg;

  return \@seg_sort;

}


sub readBdg(){
   my $file=shift;
   my @bdg;
   open BDG, $file or die "$!";
   while(<BDG>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_ =~/^\s+$/);
     my($chr,$s,$e,$v)=split/[\t ]+/,$_;
     die "$s >= $e at $_" if($s >= $e);
     push @bdg,$_;
   }
   close BDG;
   return \@bdg;
}



sub drawBdg(){
   my ($table,$region,$y,$mark)=@_;
   my ($chr,$start_seg,$end_seg)=split/:|-/,$region;
   my $x_left=$space;
   my $x_right=$space+($end_seg-$start_seg+1)*$ppx;
   my $boxH=50;
   foreach (@{$table}){
       chomp;
       my ($chr,$start,$end,$value)=split/[ \t]+/,$_;
       my $x=$x_left+($start-$start_seg+1)*$ppx;
       my $width=($end-$start+1)*$ppx;
       #$svg->rect("x",$x,"y",$y+15-$value*$scale{"all"}*$scale{$mark},"width",1,"height",$value*$scale{"all"}*$scale{$mark},"fill",$color{$mark});
       if($value>$boxH){$value=$boxH}
       $svg->rect("x",$x,"y",$y+$boxH-$value,"width",$width,"height",$value,"fill",$color{$mark});
       $x+=$ppx;
   }
   return $y+$boxH+2;
}



sub drawBdg_slim(){
   my ($data,$peaks,$region,$y,$mark,$range,$mode)=@_;
   my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
   my $len = $e_r-$s_r+1;
   my ($min,$mid,$max)=split/-/,$range;
   $min||=0;
   $max||=15;
   my $boxH=15;
   my $f=15/$max;
   my $mid_height=$boxH*$mid/$max;
   #my $x_left=$space;
   #my $x_right=$width;
   #my $x=$x_left;
   #my $switch=shift;
   $svg->line("x1",$space,"y1",$y+$boxH,"x2",$space+$len*$ppx,"y2",$y+$boxH,"stroke","grey","stroke-width",0.5);
   $svg->text("x",$space-15,"y",$y+$boxH-2,"width",10,"height",5,"-cdata",$mark,"text-anchor","end","stroke","black");

  #y axis ticks and text
   $svg->line("x1",$space,"y1",$y+$boxH+5,"x2",$space,"y2",$y,style=>{"stroke-width",1,"stroke","black"});
   for(my $i=0;$i<=$boxH;$i+=$boxH/3){
     $svg->line("x1",$space,"y1",$y+$boxH-$i,"x2",$space-2,"y2",$y+$boxH-$i,style=>{"stroke-width",1,"stroke","black"});
     if($i == 0){
       $svg->text("x",$space-5,"y",$y+$boxH-$i-2,"-cdata",int $min,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
     }elsif($i==$boxH){
       $svg->text("x",$space-5,"y",$y+$boxH-$i-2,"-cdata",int $max,style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end"})
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
           #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f,"fill",$color{$mark},"stroke","none")
           $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH-$value*$f,"width",($e-$s+1)*$ppx,"height",$value*$f-$mid_height,"fill",$color{$mark},"stroke",$color{$mark},"stroke-width",0.2)
           #$svg->rect("x",$x,"y",$y+$box_height-$value*$f,"width",1,"height",$value*$f-$mid_height,"fill",$color{$mark},"stroke","none")
          }
     }
   }elsif($mode eq "curve"){
      



    }else{die "unknow plotting type $mode"}

   #draw peak regions
   if($peaks ne "none"){
     foreach my $id(keys %{$peaks}){
       my($chr_in,$s,$e,$strand)=split/-|:/,$peaks->{$id};
       #if($s<$start && $e >$start){$s=$start}
       #if($s<$end && $e > $end){$e=$end}
       $svg->rect("x",$space+$s*$ppx,"y",$y+$boxH+2,"width",($e-$s+1)*$ppx,"height",4,"fill","black");
     }
   }

   return $y+$boxH+0.5*$boxH; #spaced to another track block
}#sub end



sub readBed(){
   my $file=shift;
   my @bed;
   open BED,$file or die"$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my($chr,$s,$e,$id,$col,$strand)=split/[\t ]+/,$_;
     die "s >= e: $s , $e" if($s >= $e);
     push @bed,join("|",($chr,$s,$e,$id,$col,$strand));
   }
   close BED;
   return \@bed;
}

sub readLen(){
  my %len;
  my $file=shift;
  open LEN, $file or die"$!";
  while(<LEN>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my($chr,$len)=split/[\t ]+/,$_;
    if(!exists $len{$chr}){$len{$chr}=$len}else{die"dup $chr\n"}
  }
  close LEN;
  return \%len;
}



sub chrPaint(){
  my ($bed,$y,$spec,$len)=@_;
  my $chrH=20;

  #draw name
  $svg->text("x",0.6*$space,"y",$y+0.6*$chrH,"-cdata",$spec,"font-family","Arial","font-size",12,"stroke","black","fill","black","text-anchor","left");

  #draw chr
  $svg->rect(x=>$space,y=>$y,height=>$chrH,width=>$len*$ppx,"fill","#694d9f","fill-opacity",0.1,"stroke","black");

  #paint chr with bed regions
  foreach (@{$bed}){
    my($chr,$s,$e,$id,$col,$strand)=split/\|/,$_;
    $svg->rect("x",$space+$s*$ppx,"y",$y,"height",$chrH,"width",($e-$s+1)*$ppx,"fill",$col,"stroke","none");
    $svg->text("x",$space+0.5*($s+$e)*$ppx,"y",$y+0.6*$chrH,"-cdata",$id,style=>{"font-family","Arial-Narrow","font-size",10,"text-anchor","middle"},"stroke","white")
  }


  return $y+$chrH+0.5*$chrH;

}

sub chrPaint_modi(){
  my ($bed,$y,$spec,$len,$flag)=@_;
  my $chrH=20;

  #draw name
  if($flag == 1){$svg->text("x",0.6*$space,"y",$y+0.6*$chrH,"-cdata",$spec,"font-family","Arial","font-size",12,"stroke","black","fill","black","text-anchor","left")};

  #draw chr
  $svg->rect(x=>$space,y=>$y,height=>$chrH,width=>$len*$ppx,"fill","#694d9f","fill-opacity",0.1,"stroke","black");

  #paint chr with bed regions
  foreach (@{$bed}){
    my($chr,$s,$e,$id,$col,$strand)=split/\|/,$_;
    $svg->rect("x",$space+$s*$ppx,"y",$y,"height",$chrH,"width",($e-$s+1)*$ppx,"fill",$col,"stroke","none");
    $svg->text("x",$space+0.5*($s+$e)*$ppx,"y",$y+0.6*$chrH,"-cdata",$id,style=>{"font-family","Arial-Narrow","font-size",10,"text-anchor","middle"},"stroke","white")
  }


  return $y+$chrH+0.5*$chrH;

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
 my ($tab,$range1,$range2,$y,$h,$score)=@_;
 $h||=100;

 #print Dumper $tab;exit; 
 my ($chr1,$s1,$e1) = split/:|-/,$range1;
 my ($chr2,$s2,$e2) = split/:|-/,$range2;
 
 #my($queChr,$queStart,$queEnd)=split/:|-/,$region;
 #my($refChr,$refStart,$refEnd)=split/:|-/,$region_next;
 foreach my $id(keys %{$tab}){
   my ($bit_score,$identity,$qstart,$qend,$qid,$sstart,$send,$sid)=split/\t/,$tab->{$id};
   next if($bit_score < $score); #for filter
   my $color;
   if($qstart > $qend || $sstart > $send){$color = $color{"link_rev"}}else{$color = $color{"link"}}

   my @temp1 = split/:/,$sid;
   $sid = $temp1[0];
   my @temp2 = split/:/,$qid;
   $qid = $temp2[0];

   next if($chr1 ne $qid || $chr2 ne $sid);
   next if($qend < $s1 || $qstart > $e1);
   next if($send < $s2 || $sstart > $e2);
   $qstart-=$s1; if ($qstart < 0){ $qstart = 0}
   $qend-=$s1; if ($qend < 0){ $qend = 0}

   $sstart-=$s2; if ($sstart < 0){ $sstart = 0}
   $send-=$s2; if ($send < 0){ $send = 0}
   
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
                        "fill"=>$color,
                        #"fill"=>"#694d9f",
                        #"fill-opacity"=>1,
                        #"fill-opacity"=>0.15,
                        #"fill-opacity"=>0.5,
                        "fill-opacity"=>0.8,
                        "stroke-width"=>6,
                       }
               );

 }#foreach m8 tab
 return $y+$h+2;
}




sub boxFilter(){
  #filter empty or other illegal elements from array, keep the order.
  my $idx=shift;
  my @filtered;
  foreach(@{$idx}){
    if(!defined $_ || $_ eq "" || $_=~/^\s+$/){}else{push @filtered,$_}
  }
  return \@filtered;
}

sub drawCoords(){
  my ($y,$region)=@_;
  my $x_left=$space;
  my ($chr_r,$s_r,$e_r) = split/:|-/,$region;
  my $len = $e_r-$s_r+1;
  my ($start,$end)=($s_r,$e_r);
  my ($unit,$tip);
  if($end-$start<=10000){$unit=1000}elsif($end-$start>10000 && $end-$start<=100000){$unit=5000}elsif($end-$start>100000){$unit=50000} #the step
  if($end-$start<= 5000){$unit = 200}
  $svg->line("x1",$x_left,"y1",$y,"x2",$x_left+$len*$ppx,"y2",$y,style=>{"stroke-width",1,"stroke","grey","stroke-width",0.5,"stroke-dasharray","2,2"}); #the backbone line
  #$svg->text("x",$x_left-10,"y",$y,);
  for(my $j=$start;$j<=$end;$j+=$unit){
     $svg->line("x1",$x_left+($j-$start)*$ppx,"y1",$y-2,"x2",$x_left+($j-$start)*$ppx,"y2",$y+2,style=>{"stroke-width",1,"stroke","black"}); #the tick
     my ($text_x,$text_y)=($x_left+($j-$start)*$ppx,$y-2);
     $svg->text("x",$text_x,"y",$text_y,"-cdata",int($j/1000)."k",style=>{"font-family","Arial-Narrow","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark
  }#for end here
 return $y+5;
}

sub drawHiLightBox(){
  my ($regions,$y_up,$y_down) = @_;
  die "y_up >= y_down: $y_up >= $y_down" if($y_up >= $y_down);
  my $h = $y_down - $y_up +1;
  foreach (@{$regions}){
    my($chr,$s,$e,$id,$col,$strand)=split/\|/,$_;
    next if($id !~/GI/);
    $s+=1500; #too close with left y tick
    #if(!defined $col){$col = "#1d953f"}
    if(!defined $col || $col=~/^\s+$/){$col = "white"}
    if(!defined $strand || $strand=~/^\s+$/){$strand = "+"}
    $svg->rect("x",$space+$s*$ppx,"y",$y_up,"height",$h,"width",($e-$s+1)*$ppx,"fill","none","fill-opacity",0.5,"stroke","black","stroke-width",2,"stroke-dasharray","2,1");
   # if($label eq "yes"){$svg->text("x",$space+0.5*($s+$e)*$ppx,"y",$y+0.6*$chrH,"-cdata",$id,style=>{"font-family","Arial-Narrow","font-size",10,"text-anchor","middle"},"stroke","white")}
  }

  return 1;

}


sub isValid(){
  #check and make sure that hash key is not space, is defined and not empty
  my $v=shift;
  my $flag=1;
  if(!defined $v || $v eq "" || $v=~/^\s+$/){$flag = 0}
  return $flag;
}



sub drawLegend(){
   my ($x,$y,$w,$h)=@_;
   my $x_edge = 0.25*$w;
   my $y_edge = 0.25*$h;
   my $y_dist=$h/5;
   my $barH = 10;
   #main frame
   $svg->rect("x",$x,"y",$y,"width",$w,"height",$h,"fill","none","stroke","black","stroke-width",0.6);
   
   


   return 1;
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



