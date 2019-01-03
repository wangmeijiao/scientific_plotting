#!/bin/perl
use strict;
use warnings;
use Getopt::Long;
#use Data::Dumper;
use SVG;
use POSIX qw[ _exit ];


 my ($refLen,$queLen,$refDensi,$queDensi,$axt);
 my($data,$header,$dimXY,$zlim,$outfile,$gene,$ltr,$GI,$GI_chr);

 GetOptions( "refLen=s"=>\$refLen,
#             "queLen=s"=>\$queLen,
             "refDensi=s"=>\$refDensi,
#             "queDensi=s"=>\$queDensi,
#             "axt=s"=>\$axt,
              "data=s"=>\$data,"header=s",\$header,"dimXY=s"=>\$dimXY,"zlim=s"=>\$zlim,
              "outfile=s"=>\$outfile,
              "gene=s"=>\$gene,"ltr=s"=>\$ltr,"GI=s",\$GI, "GI_chr=s",\$GI_chr
           );
my($refBin,$queBin);
($refLen,$refBin)=split/,/,$refLen;
#($queLen,$queBin)=split/,/,$queLen;
#drawing parameters
my ($width,$height)=split/,/,$dimXY;
#my $width=800;
#my $height=1000;  
my $edge=100;
my $space = $edge;
my $x_left=$edge;
my $x_right=$edge+$width;
my $ppx=$width/$refLen; #x pix/bp
my $shift_que=100;
my $track_h=15;
#my $track_h=30;
my %color=(
           "bg_color"=>"rgb(240,240,240)",
           "fill_color"=>"rgb(114,207,142)",
           ##marks
           "H3K4me2"=>"#7fb80e",
           "k42"=>"#7fb80e",
           "k42m1"=>"#7fb80e",
           "H3K4me3"=>"#1d953f",
           "k43"=>"#1d953f",
           "k43m1"=>"#1d953f",
           "H3K36me3"=>"#375830",
           "k363"=>"#375830",
           "k363m1"=>"#375830",
           "H3K9ac"=>"#375830",
           #"H3K9ac"=>"#78a355",
           "k9ac"=>"#78a355",
           "k9acm1"=>"#78a355",
           "H3K56ac"=>"#2b6447",
           "k56ac"=>"#2b6447",
           "H4K12ac"=>"#40835e",
           "h4k12ac"=>"#40835e",
           "k12ac"=>"#40835e",
           "H3K9me2r1"=>"#d71345",
           "H3K9me2"=>"#d71345",
           "k92r1"=>"#d71345",
           "k92"=>"#d71345",
           "k92M1"=>"#d71345",
           "H3K9me2r2"=>"#ed1941",
           "k92r2"=>"#ed1941",
           "H3K27me1"=>"#892f1b",
           "k271"=>"#892f1b",
           "H3K27me2"=>"#f58f98",
           "k272"=>"#f58f98",
           "k272M1"=>"#f58f98",
           "H3K27me3"=>"#a03939",
           "k273"=>"#a03939",
           #function
           "RNASEQ"=>"green",
           "DNAmeth"=>"grey",
           "siRNA"=>"black",
           #anno
           "gene"=>"green",
           "LTR"=>"blue",
           "DNA-TEs"=>"yellow",
           #comparaOrth        
           "orthBlock"=>"grey",
           "Eigen"=>"black",
           "GI"=>"purple",
);


my %range = (
             "H3K9me2r1" => "0.3001563,3.3964547",
             "H3K9me2r2" => "0.3300492,2.6666756",
             "H3K27me1" => "0.2350127,3.6071347",
             "H3K4me3" => "0.5,3.0587623",
             #"H3K4me3" => "0.0000000,3.0587623",
             "H3K9ac" => "0.5,3.0036923",
             #"H3K9ac" => "0.0000000,3.0036923",
             "H3K36me3" => "0.0000000,2.9531159",
             "Eigen" => "-0.044409475,0.036124556",

);


my %scale_up=(
           "all"=>1,
           ##marks
           "H3K4me2"=>30,
           "k42"=>6,
           "k42m1"=>1,
           "H3K4me3"=>11,
           "k43"=>8,
           "k43m1"=>1,
           "H3K36me3"=>20,
           "k363"=>10,
           "k363m1"=>1,
           "H3K9ac"=>17,
           "k9ac"=>6,
           "k9acm1"=>1,
           "H3K56ac"=>80,
           "k56ac"=>4,
           "H4K12ac"=>25,
           "h4k12ac"=>7,
           "k12ac"=>1,
           "H3K9me2r1"=>26,
           "k92r1"=>5,
           "H3K9me2r2"=>34,
           "k92r2"=>4,
           "k92"=>1,
           "k92M1"=>1,
           "H3K27me1"=>44,
           "k271"=>4,
           "H3K27me2"=>1,
           "k272"=>1,
           "k272M1"=>1,
           "H3K27me3"=>11,
           "k273"=>14,
           #function
           "RNASEQ"=>1,
           "DNAmeth"=>1,
           "siRNA"=>1,
           "eigen_10k"=>1,
           #"eigen_10k"=>-0.01,
);

my %scale_down=(
           "all"=>1,
           ##marks
           "H3K4me2"=>114,
           "k42"=>5.5,
           "k42m1"=>1,
           "H3K4me3"=>153,
           "k43"=>5,
           "k43m1"=>1,
           "H3K36me3"=>100,
           "k363"=>4.5,
           "k363m1"=>1,
           "H3K9ac"=>110,
           "k9ac"=>4,
           "k9acm1"=>1,
           "H3K56ac"=>1,
           "k56ac"=>1,
           "H4K12ac"=>1,
           "h4k12ac"=>1,
           "k12ac"=>1,
           "H3K9me2r1"=>1,
           "k92r1"=>1,
           "H3K9me2r2"=>1,
           "k92r2"=>1,
           "k92"=>9,
           "H3K9me2"=>26,
           "k92M1"=>1,
           "H3K27me1"=>1,
           "k271"=>1,
           "H3K27me2"=>54,
           "k272"=>5,
           "k272M1"=>1,
           "H3K27me3"=>1,
           "k273"=>1,
           #function
           "RNASEQ"=>1,
           "DNAmeth"=>1,
           "siRNA"=>1,
);

#prepare drawing board
my $svg=SVG->new(height=>$height,width=>$width+2*$edge);
my $y_next=30;
my $y1=$y_next;
#draw up densi
$y_next=&drawDensiBin($refDensi,$y_next,"up",$refLen*$ppx/$refBin);
print STDERR "drawing up density done\n";
#my $y2=$y_next+5;

=pod
#draw axt backbone
$y_next=&drawAxt($axt,$y_next+5);
print STDERR "drawing axt done\n";
my $y3=$y_next;
#draw down densi
$y_next=&drawDensiBin($queDensi,$y_next,"down",$queLen*$ppx/$queBin);
print STDERR "drawing down density done\n";
my $y4=$y_next;
=cut



#geneBed
#print STDERR "#gene track\n";
#$y_next=&drawBed(&readBed($gene),$y_next+10,"gene");

#LTRBed
#print STDERR "#LTR track\n";
#$y_next=&drawBed(&readBed($ltr),$y_next+2,"LTR");


#GI of chr04
print STDERR "#GI track chromosome\n";
$y_next=&drawBed(&readBed($GI),$y_next,"GI");


#GI of H1 
print STDERR "#GI track H1\n";
$y_next=&drawBed(&readBed($GI_chr),$y_next,"GI");


my $y2=$y_next-2;
$svg->rect("x",$x_left+1700000*$ppx,"y",$y1,"height",$y2-$y1,"width",(4500000-1700000+1)*$ppx,"fill","#694d9f","fill-opacity",0.1,"stroke-dasharray","2,2","stroke","black");
#$svg->line("x1",$x_left+1700000*$ppx,"y1",$y_ref+2,"x2",$x_left+4500000*$ppx,"y2",$y_ref+2,"stroke-width",5,"stroke","red");
#$svg->ellipse("cx",$x_left+9809525*$ppx,"cy",$y2-2,"rx",4,"ry",2,"fill","red"); #9809525
$svg->ellipse("cx",$x_left+10285431*$ppx,"cy",$y2-4,"rx",4,"ry",2,"fill","red"); #9809525



#coordination line
$y_next = &drawCoords($y_next+5,$refLen);


#draw pearson matrix
##$y_next = &drawPearsonMat($data,$header,$dimXY,$zlim,$y_next);

#color bar
##&drawColorBar(&hotCold510(),$space+$width+5,$y_next-$ppx*5000000,$ppx*5000000,"v",$zlim);
#&drawColorBar(&hotCold510(),$space+$width+5,$y_next-50,50,"v",$zlim);



&svgWrite($svg,$outfile);
POSIX::_exit(0); #to speed up exit by ommitting perl_destruct for big hash

#print $svg->xmlify();
#print STDERR "all done.\n";






####sub###
sub drawDensiBin(){
  # mark1       mark2   mark3   DNAmeth mRNA    siRNA
  # 2           3        2       2      5         1
  # 2           3        2       2      5         1
  my $file=shift;
  my $y_start=shift;
  my $side=shift;
  my $step=shift;
  my $y=$y_start;
  open IN, $file or die "$!";
  my $head=<IN>;
  chomp $head;
  my @mark_order=split/[\t ]+/,$head;print STDERR "\ndata order:@mark_order\n";
  my $num_tracks = scalar @mark_order;
  die "empty densi data" if($num_tracks == 0);
  my $x=$x_left;
 #read a line and draw
  while(<IN>){
     chomp;
     next if($_ eq "" ||$_=~/^#/);
     my @data=split/[\t ]+/,$_;     
     ##draw
     my $cnt=0;
     foreach my$value(@data){
       #my $N;
       #if($side eq "up"){$N=$scale_up{$mark_order[$cnt]}}elsif($side eq "down"){$N=$scale_down{$mark_order[$cnt]}}
       next if($value eq "NA" || $value eq "N/A");
       my ($min,$max) = split/,/,$range{$mark_order[$cnt]};
       if($cnt == $num_tracks - 1){ #for eigen data that has negative points, here cnt = ncol - 1, egien is the last column
         my $y_mid = $y+0.5*$track_h;
         #my ($min,$max) = (-0.044409475,0.036124556);
         if($value <= 0){
           if($value<$min){$value=$min}
           $value=($value/$min)*(0.5*$track_h);
           $svg->rect("x",$x,"y",$y_mid,"width",$step,"height",$value,"fill",$color{$mark_order[$cnt]},"stroke",$color{$mark_order[$cnt]},"stroke-width",0.1);
         }else{ 
               if($value>$max){$value=$max}
               $value=($value/$max)*(0.5*$track_h);
               $svg->rect("x",$x,"y",$y_mid-$value,"width",$step,"height",$value,"fill",$color{$mark_order[$cnt]},"stroke",$color{$mark_order[$cnt]},"stroke-width",0.1);
              } 

       }else{
         die "value is < 0 at @data" if($value < 0);
         if($value<$min){$value=$min}
         if($value>$max){$value=$max}
         $value=(($value-$min)/($max-$min))*$track_h;
         #$svg->rect("x",$x,"y",$y+$track_h-$value,"width",$step,"height",$value,"fill",$color{$mark_order[$cnt]},"stroke","none");
         $svg->rect("x",$x,"y",$y+$track_h-$value,"width",$step,"height",$value,"fill",$color{$mark_order[$cnt]},"stroke",$color{$mark_order[$cnt]},"stroke-width",0.2);
       }
       $y+=$track_h+1;
       $cnt++;
     }
     $y=$y_start;
     $x+=$step; 
  }#while end here
  close IN;
  my $track;
  foreach my $mark(@mark_order){
   $track++;
   $svg->text("x",$x_left-60,"y",$y_start+$track*$track_h,"font-family","Arial","font-size",10,"stroke",$color{$mark},"-cdata",$mark,"text-anchor","start")
  }
  return $y_start+($track_h+1)*scalar (@mark_order); 
}#sub drawDensiMedian


sub drawAxt{
 my $file=shift;
 my $y_ref=shift;
 my $y_que=$y_ref+60;
 $svg->rect(x=>$x_left,y=>$y_ref,height=>5,width=>$refLen*$ppx);
 $svg->rect(x=>$x_left,y=>$y_que,height=>5,width=>$queLen*$ppx); 
 $/="\n\n";
 open AXT, $file or die "$!";
 while(<AXT>){
    chomp;
    my @box=split/\n/,$_;
    my $align=shift @box;
    my ($blockID,$refID,$refStart,$refEnd,$queID,$queStart,$queEnd,$strand,$score)=split/ / ,$align;
    #draw collinear 
     $svg->line(x1=>$x_left+$refStart*$ppx,y1=>$y_ref+5,x2=>$x_left+$queStart*$ppx,y2=>$y_que,stroke=>"grey","stroke-width"=>0.1);
     $svg->line(x1=>$x_left+$refEnd*$ppx,y1=>$y_ref+5,x2=>$x_left+$queEnd*$ppx,y2=>$y_que,stroke=>"grey","stroke-width"=>0.1);
 }#while end here
 close AXT;
 $/="\n";
 
 #draw knob1 region
  # 1700000-4500000:944841-1815604 
  $svg->line("x1",$x_left+1700000*$ppx,"y1",$y_ref+2,"x2",$x_left+4500000*$ppx,"y2",$y_ref+2,"stroke-width",5,"stroke","red");
  $svg->line("x1",$x_left+944841*$ppx,"y1",$y_que+2,"x2",$x_left+1815604*$ppx,"y2",$y_que+2,"stroke-width",5,"stroke","red");
 
  #draw centremere region
   #9465841-10285431: mid 9875636
   #4648900-4657405:  mid 4657405
   $svg->ellipse("cx",$x_left+10285431*$ppx,"cy",$y_ref+2,"rx",4,"ry",2,"fill","white");
   $svg->ellipse("cx",$x_left+4657405*$ppx,"cy",$y_que+2,"rx",4,"ry",2,"fill","white");
  return $y_que;
}#end of drawAxt



sub drawHiLightBox(){
  my($x1_up,$x2_up,$y1,$y2)=@_;
  $svg->rect("x",$x_left+$x1_up*$ppx,"y",$y1,"height",$y2-$y1,"width",($x2_up-$x1_up+1)*$ppx,"fill","#694d9f","fill-opacity",0.1,"stroke-dasharray","2,2","stroke","black");
  #$svg->rect("x",$x_left+$x1_down*$ppx,"y",$y3,"height",$y4-$y3,"width",($x2_down-$x1_down+1)*$ppx,"fill","#694d9f","fill-opacity",0.1,"stroke-dasharray","2,2","stroke","black");
  #my @x=($x_left+$x1_up*$ppx,$x_left+$x2_up*$ppx,$x_left+$x2_down*$ppx,$x_left+$x1_down*$ppx);
  #my @y=($y2,$y2,$y3,$y3);
  #my $points=$svg->get_path(
  #     x=>\@x,
  #     y=>\@y,
  #     "-type"=>"polygon",
  #);
  #print STDERR Dumper $points;
  #$svg->polygon(
  #    %{$points},
  #    "style"=>{
  #       "fill"=>"#694d9f",
  #       "fill-opacity"=>0.1,
  #       "stroke"=>"black",
  #       "stroke-dasharray"=>"2,2",
  #    }
  #);
  return 1;
}


sub drawPearsonMat(){
  
  #data matrix format:
  #  a    b     c
  #  0.2  0.12  0.9
  #  0.89 7.1   1.9
  
  # row order must be clustered or ordered beforehand (draw as is)
  # header=T/F
  # check and make sure that not have NA and empty values, all datapoint must be numeric (rownames col are not allowed)
  # try to guess matrix nrow and ncol after readin the data and assign the best width x height
  
  my($data,$header,$dimXY,$zlim,$y_start) = @_;
  my ($width,$height)=split/,/,$dimXY; #not include x/y margin
  $height = 150;
  my ($zmin,$zmax)=split/,/,$zlim;
  
  my ($data_idx,$header_idx,$nrow,$ncol,$min,$max)=&readTable($data,$header);
  #print Dumper $header_idx,$nrow,$ncol,$min,$max,$data_idx;
  #exit;
  
  print STDERR "start to draw pearson mat with $nrow,$ncol,$min,$max at $y_start\n";
  
  #($zmin,$zmax) = ($min,$max);
  
  my ($box_h,$box_w);
  #my $y_start=50;
  #my $edge=50;
  $box_w=$width/$ncol;
  $box_h= $box_w;
  #$box_h=$height/$nrow;
  my $x_left=$edge;
  my $x_right=$edge+$width;
  #my @color=@{&greyMono256()};
  my @color=@{&hotCold510()};
  #my @color=@{&redMono256()};
  #my @color=@{&greenMono256()};
  my $ncolor=$#color;
  
  #my $svg=SVG->new("width",$width+2*$edge+10*3,"height",$height+$y_start*2);
  
  my $cnt_line;
  my $y=$y_start;
  foreach (@{$data_idx}){
   $cnt_line++;
   #if($cnt_line >2000){last}
   if($cnt_line%1000==0){print STDERR "#"}
   my @box=split/[\t ]+/,$_;
   my $x=0;
   foreach my $data(@box){   
     my $colN; #color index
    # $data*=-1;
     if($data<=$zmin){$colN=0}elsif($data>$zmin && $data<$zmax){$colN=&round ($ncolor*($data-$zmin)/($zmax-$zmin))}else{$colN=$ncolor} #cut and map datapoint to color index
     #$svg->rect("x",$x_left+$x,"y",$y,"width",$box_w,"height",$box_h,"fill",$color[$colN],"stroke","none","stroke-opacity",0.8,"stroke-width",0);
     $svg->rect("x",$x_left+$x,"y",$y,"width",$box_w,"height",$box_h,"fill",$color[$colN],"stroke",$color[$colN],"stroke-width",1);#stroke-width, the bigger, the shapper when deal with huge data. But not too big
     $x+=$box_w;
     #$x+=$box_w+10;
   }
   $y+=$box_h; 
  }
  print STDERR "\n";
  
  #draw legends
  #my $text_x=$x_left+0.5*$box_w;
  #my $text_y=$y_start-10;
  #foreach(@{$header_idx}){
  #   $svg->text("x",$text_x,"y",$text_y,"width",$box_w,"height",10,"font-family","Arial", "font-weight","bold","text-anchor","middle","font-size",10, "-cdata", $_);
  #   $text_x+=$box_w+10;  
  #}
  
  return $y; 

}

sub readTable(){
   my $file=shift;
   my $head=shift;
   my @header;
   my @data;
   my $cnt;
   my $ncol;
   my ($min,$max);
   open DATA, $file or die "$!";
   while(<DATA>){
     chomp;
     next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
     $cnt++;
     if($cnt==1){
       if($head eq "T"){
           @header=split/[\t ]+/,$_;
           $ncol=scalar @header;
           next;
       }elsif($head eq "F"){
           my @box=split/[\t ]+/,$_;
           $ncol=scalar @box;
           if(&notNA(\@box) == 0){die "has NA at $_"}
           if(&isNumeric(\@box) == 0){die "not numeric at $_"}
           push @data,$_;
           ($min,$max)=($box[0],$box[0]);
         }else{die "unknow header swith $head"}
     }
     my @box_line=split/[\t ]+/,$_;
     my $ncol_line=scalar @box_line;
     if(!defined $ncol || $ncol_line != $ncol){die "ncol err $ncol != $ncol_line at $_"}
     if(&notNA(\@box_line) == 0){die "has NA at $_"}
     if(&isNumeric(\@box_line) == 0){die "not numeric at $_"}
     push @data,$_;
     if(! defined $min && !defined $max){
       ($min,$max)=($box_line[0],$box_line[0])
     }elsif(defined $min && defined $max){
         foreach my $ctrl(@box_line){ 
           if($ctrl < $min){$min = $ctrl}
           if($ctrl > $max){$max = $ctrl}
         }
       }else{die "err min max $min , $max"}    
   }#while end
   close DATA;
   if($head eq "T"){$cnt-=1}
   return (\@data,\@header,$cnt,$ncol,$min,$max);

}

sub isNumeric(){
  my $idx=shift;
  my $flag=1;
  foreach my $ctrl (@{$idx}){
    #signed    
    $ctrl=~s/^-//;
    #decimals
     if($ctrl =~/[^0-9\.E-]/){$flag = 0}

    #scientific notation

  }
  return $flag;
}

sub notNA(){
    my $idx=shift;
    my $flag=1;
    foreach my $ctrl (@{$idx}){if($ctrl eq "NA" || $ctrl eq "NULL" || $ctrl eq "Na" || $ctrl eq "n/a" || $ctrl eq "N/A"){$flag=0}}
    return $flag;
}

sub redMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb(255,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub redMono9(){
   my @color=("#FFFFFF", "#FFDFDF", "#FFBFBF", "#FF9F9F" ,"#FF7F7F" ,"#FF5F5F", "#FF3F3F" ,"#FF1F1F", "#FF0000");
   return \@color;  
}

sub greenMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,255,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub greyMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub hotCold510(){
  my @color;
  for(my $i=0;$i<=255;$i++){
    push @color,"rgb($i,$i,255)";
  }
  for(my $i=255;$i>=0;$i--){
    push @color,"rgb(255,$i,$i)";
  }
  return \@color;

}

sub round(){
   my $n=shift;
   my $i=int $n;
   my $d=$n-$i;
   if($d>0.6){return $i+1}else{return $i}
}

sub minMax(){
   my @data=@{shift @_};
   my ($min,$max)=(0,0);
   foreach(@data){
     if($min>$_){$min=$_}
     if($max<$_){$max=$_}
   }
   return($min,$max);
}


sub readBed(){
   my $file=shift;
   my @bed;
   open BED,$file or die"$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my($chr,$s,$e,$id)=split/[\t ]+/,$_;
     #my($chr,$s,$e,$id,$col,$strand)=split/[\t ]+/,$_;
     die "s >= e: $s , $e" if($s >= $e);
     push @bed,join("|",($chr,$s,$e,$id));
     #push @bed,join("|",($chr,$s,$e,$id,$col,$strand));
   }
   close BED;
   return \@bed;
}

sub drawBed(){  #for annotation tracks only, not for peaks
  my($table,$y,$ftype)=@_;
  my $space = $edge;
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
     my($chr,$start,$end,$name)=split/\|/,$ctrl;
     #my($chr,$start,$end,$name,undef,$strand)=split/\|/,$ctrl;
     #if($strand eq "+"){
     #   $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_te-6,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$feature_line[0]},"stroke",$color{$feature_line[0]}   );
     #}elsif($strand eq "-"){
     #   $svg->rect("x",$edge+($start-$start_seg)*$ppx,"y",$y_te+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$feature_line[0]},"stroke",$color{$feature_line[0]}   );
     # }
     $svg->rect("x",$x_left+$start*$ppx,"y",$y+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$ftype},"stroke","none" );
     #$svg->rect("x",$x_left+$start*$ppx,"y",$y+2,"width",($end-$start+1)*$ppx,"height",$barH,"fill",$color{$ftype},"stroke",$color{$ftype} );

  }
  return $y+2*$barH+2;
}



sub drawCoords(){
  my ($y,$len)=@_;
  my $space = $edge;
  my $x_left=$space;
  my ($start,$end)=(1,$len);
  my ($unit,$tip);
  if($end-$start<=10000){$unit=1000}elsif($end-$start>10000 && $end-$start<=100000){$unit=5000}elsif($end-$start>100000){$unit=1000000} #the step
  $svg->line("x1",$x_left,"y1",$y,"x2",$x_left+$len*$ppx,"y2",$y,style=>{"stroke-width",1,"stroke","grey","stroke-width",0.5,"stroke-dasharray","2,2"}); #the backbone line
  #$svg->text("x",$x_left-10,"y",$y,);
  for(my $j=$start;$j<=$end;$j+=5*$unit){
     $svg->line("x1",$x_left+($j-$start)*$ppx,"y1",$y-2,"x2",$x_left+($j-$start)*$ppx,"y2",$y,style=>{"stroke-width",1,"stroke","black"}); #the tick
     my ($text_x,$text_y)=($x_left+($j-$start)*$ppx,$y-2);
     $svg->text("x",$text_x,"y",$text_y,"-cdata",int($j/$unit)."M",style=>{"font-family","Arial","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark
  }#for end here
  
  $svg->line("x1",$x_left,"y1",$y,"x2",$x_left,"y2",$y+5000000*$ppx,style=>{"stroke-width",1,"stroke","grey","stroke-width",0.5,"stroke-dasharray","2,2"}); #the backbone line

  $svg->line("x1",$x_left-1,"y1",$y,"x2",$x_left,"y2",$y,style=>{"stroke-width",1,"stroke","black"}); #the tick
  my ($text_x,$text_y)=($x_left-1,$y);
  #$svg->text("x",$text_x,"y",$text_y,"-cdata","0M",style=>{"font-family","Arial","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark


  $svg->line("x1",$x_left-3,"y1",$y+1700000*$ppx,"x2",$x_left,"y2",$y+1700000*$ppx,style=>{"stroke-width",1,"stroke","black"}); #the tick
  ($text_x,$text_y)=($x_left-3,$y+1700000*$ppx);
  $svg->text("x",$text_x,"y",$text_y,"-cdata","1.7M",style=>{"font-family","Arial","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark

  $svg->line("x1",$x_left-3,"y1",$y+4500000*$ppx,"x2",$x_left,"y2",$y+4500000*$ppx,style=>{"stroke-width",1,"stroke","black"}); #the tick
  ($text_x,$text_y)=($x_left-3,$y+4500000*$ppx);
  $svg->text("x",$text_x,"y",$text_y,"-cdata","4.5M",style=>{"font-family","Arial","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark

  $svg->line("x1",$x_left-1,"y1",$y+5000000*$ppx,"x2",$x_left,"y2",$y+5000000*$ppx,style=>{"stroke-width",1,"stroke","black"}); #the tick
  ($text_x,$text_y)=($x_left-1,$y+5000000*$ppx);
  $svg->text("x",$text_x,"y",$text_y,"-cdata","5M",style=>{"font-family","Arial","font-size",6,"text-anchor","end","font-weight","bold"}); #the tickmark

  

 return $y;
}



sub drawColorBar(){
   my $color= shift @_;
   my $numCol=scalar @$color;
   my ($x,$y,$len,$type,$zlim)=@_;
   my ($min,$max) = split/,/,$zlim;
   my $ppx_local=$len/$numCol;
   #vertical or horizonal
   if($type eq "v"){
      foreach(reverse @$color){
        $svg->rect("x",$x,"y",$y,"height",$ppx_local,"width",15,"fill",$_,"stroke","none");#stroke "none" unwork
        $y+=$ppx_local;
      }
    }elsif($type eq "h"){
        foreach(reverse @$color){
          $svg->rect("x",$x,"y",$y,"height",15,"width",$ppx_local,"fill",$_,"stroke",$_);
          $x+=$ppx_local;
         }
     }
    #draw text
    $min = sprintf("%.1f",$min);
    $max = sprintf("%.1f",$max);
    my ($x_text,$y_text)=($x+15,$y-$len);
    $svg->text("x",$x_text,"y",$y_text,"-cdata","$max","stroke","black","font-size",8);
    #$svg->text("x",$x_text,"y",$y_text,"-cdata","$min","transform","rotate(-90,$x_text,$y_text)","stroke","black","font-size",6);
    ($x_text,$y_text)=($x+15,$y);
    $svg->text("x",$x_text,"y",$y_text,"-cdata","$min","stroke","black","font-size",8);
    #$svg->text("x",$x_text,"y",$y_text,"-cdata","$max","transform","rotate(-90,$x_text,$y_text)","stroke","black");
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


