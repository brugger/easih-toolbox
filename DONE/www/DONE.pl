#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (18 Jul 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::HTML;
use EASIH::DONE;


#EASIH::DONE::Connect('done_dev');



$0 =~ s/.*\///;

print EASIH::HTML::start('EASIH-DONE', 'EASIH-DONE.css');


if ( $EASIH::HTML::parameters{ 'QC' } ) {

  my $fid = $EASIH::HTML::parameters{ 'fid' };
  
  print_header();
  my $boxid1 = "b1_boxplot";
  print_boxplot($fid, $boxid1);
  my $base1 = "basesplits1";
  print_basesplits($fid, $base1);
  my $qvshist1 = "qvshist1";
  print_qv_histogram($fid, $qvshist1);
  my $gc1 = "gc1";
  print_gc($fid, $gc1);
  my $ad1 = "ad1";
  print_adaptor_histogram($fid, $ad1);
  my $dh1 = "dh1";
  print_dups_histogram($fid, $dh1);
  
  my ( $name, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC, $sid, $rid ) = EASIH::DONE::fetch_file_info( $fid );
  
  $name =~ s/.*\///;
  

  my @traffic_light;
  push @traffic_light,["NOT ANALYSED YET..."] if (!$Q30bases && !$duplicates && !$partial_adaptors && !$Avg_AC );
  push @traffic_light,[{value => 'Bases >=Q30:'}, {value=>'<70%', bgcolor=>'red'},{value=>'70-90%', bgcolor=>'yellow'},{value=>'>90%', bgcolor=>'#00DD00'}] if ( $Q30bases );
  push @traffic_light,[{value => 'Duplicates sequences:'}, {value=>'>10%', bgcolor=>'red'},{value=>'1-10%', bgcolor=>'yellow'},{value=>'<1%', bgcolor=>'#00DD00'}] if ( $duplicates );
  push @traffic_light, [{value => 'Partial adaptors:'}, {value=>'>10%', bgcolor=>'red'},{value=>'1-10%', bgcolor=>'yellow'},{value=>'<1%', bgcolor=>'#00DD00'}] if ( $partial_adaptors );
  push @traffic_light, [{value => 'Mean AC:'}, {value=>'<40% or > 60%', bgcolor=>'red'},{value=>'40-45% or 55-60%', bgcolor=>'yellow'},{value=>'45-55%', bgcolor=>'#00DD00'}] if ( $Avg_AC );
  

  easih_top();
  
  print " <H1>QC report for: $name [$fid]</H1>";
  

  my @data = ([{value=>'Filename:', bgcolor => 'grey'}, {value=>"$name", bgcolor => 'grey'}],
	      [{value=>'Total reads:' }, {value=>"$total_reads" }],
	      [{value=>'Read length:', bgcolor => 'grey'}, {value=>"$read_length bp", bgcolor => 'grey'}],
	      ['Sample size:', "$sample_size"]);
  
  my $run_name = EASIH::DONE::fetch_run_name( $rid );
  
  push @data, [{value=>'Run:', bgcolor => 'grey' }, {value=>"<a href=$0?rid=$rid>$run_name $rid $sid</A>", bgcolor => 'grey' }] if ($run_name);
  
  my $colour = 'red';
  $colour = 'yellow' if ($Q30bases >= 70 );
  $colour = '#00DD00'  if ($Q30bases >= 90 );
  
  push @data, [{value=>'Q30 bases:', bgcolor => $colour}, {value=>"$Q30bases %", bgcolor => $colour}] if ($Q30bases);
  
  $colour = 'red';
  $colour = 'yellow' if ($duplicates <10 );
  $colour = '#00DD00'  if ($duplicates < 1 );
  push @data, [{value=>'Duplicates:', bgcolor => $colour}, {value=>"$duplicates %", bgcolor => $colour}] if ($duplicates);
  
  $colour = 'red';
  $colour = 'yellow' if ($partial_adaptors < 10  );
  $colour = '#00DD00'  if ($partial_adaptors < 1 );
  push @data, [{value=>'Partial adaptors:', bgcolor => $colour}, {value=>"$partial_adaptors %", bgcolor => $colour}] if ($partial_adaptors);
  
  $colour = 'red';
  $colour = 'yellow' if ($Avg_AC > 40 && $Avg_AC < 60);
  $colour = '#00DD00'  if ($Avg_AC > 45 && $Avg_AC < 55);
  push @data, [{value=>'Mean AC:', bgcolor => $colour}, {value=>"$Avg_AC %", bgcolor => $colour}] if ($Avg_AC);
  
#, $partial_adaptors, $Avg_AC


  print EASIH::HTML::advanced_table(\@data, 1, 0, 0, 'lightgrey', 0, '700px') . "<br>\n";
  
  print EASIH::HTML::advanced_table(\@traffic_light, 1, 0, 0, 'lightgrey', 0, '700px');

  print "
    <div id='$boxid1' style='width: 900px; height: 300px;'></div>
    <div id='$qvshist1' style='width: 900px; height: 300px;'></div>
    <div id='$dh1' style='width: 900px; height: 300px;'></div>
    <div id='$ad1' style='width: 900px; height: 300px;'></div>
    <div id='$base1' style='width: 900px; height: 300px;'></div>     
    <div id='$gc1' style='width: 900px; height: 300px;'></div>
";


  my @mappings = (['database', 'Hits', 'non-unique hits', 'duplicates']);
  
  my @mapping_stats = EASIH::DONE::fetch_mapping_stats($fid);
  my ($max_ref, $max) = (undef, -1);
  
  foreach my $mapping ( @mapping_stats ) {
    if ( $max < $$mapping[1] ) { 
      $max = $$mapping[1]; 
      $max_ref = $$mapping[0];
    }
  }
  
  foreach my $mapping (@mapping_stats) {

  
    my ($reference, $unique_hits, $non_unique_hits, $duplicates) = @$mapping;
    
    $unique_hits     ||= 0;
    $non_unique_hits ||= 0;
    $duplicates      ||= 0;

    next if (! $unique_hits);
    if ( $reference eq "PhiX") {
      push @mappings, [{value => $reference, bgcolor=>'coral'}, 
		       {value => $unique_hits, bgcolor=>'coral'},
		       {value => $non_unique_hits, bgcolor=>'coral'}, 
		       {value => $duplicates, bgcolor=>'coral'}];
    }
    elsif ( $reference eq $max_ref) {
      push @mappings, [{value => $reference, bgcolor=>'#00d00'}, 
		       {value => $unique_hits, bgcolor=>'#00d00'},
		       {value => $non_unique_hits, bgcolor=>'#00d00'}, 
		       {value => $duplicates, bgcolor=>'#00d00'}];

    }
    else {
      push @mappings, [$reference, $unique_hits, $non_unique_hits, $duplicates];
    }      
  }


  print EASIH::HTML::advanced_table(\@mappings, 1, 1, 0, undef, undef, '700px');
  

  print "
  </div>
  </body>
</html>
";
}
elsif ( $EASIH::HTML::parameters{ 'rid' } ) {

  my $rid = $EASIH::HTML::parameters{ 'rid' };

  my $run_name = EASIH::DONE::fetch_run_name( $rid );
  easih_top();
  print "<h1> Stats for run: '$run_name'</h1>";

  my @lanes  = EASIH::DONE::fetch_illumina_lane_stats_by_rid($rid);
  my @mplexs = EASIH::DONE::fetch_illumina_sample_stats_by_rid($rid);

  print STDERR Dumper( @lanes );

  my @rundata;
  foreach my $lane ( sort { $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2] } @lanes ) {
    my @lane_data;

    foreach my $lane_data (@$lane ) {
      if ( $lane_data > 100 && $lane_data < 1000) {
	$lane_data = sprintf("%.2fK", $lane_data/1000);
      }
      elsif ( $lane_data > 1000  && $lane_data < 10000) {
	$lane_data = sprintf("%.1fK", $lane_data/1000);
      }
      elsif ( $lane_data > 10000  && $lane_data < 100000) {
	$lane_data = sprintf("%dK", $lane_data/1000);
      }
    }

    push @lane_data, @$lane[2..6], 
                     "$$lane[7]&plusmn;$$lane[8]", "$$lane[9]&plusmn;$$lane[10]",  "$$lane[11] / $$lane[12]", 
                     "$$lane[13]&plusmn;$$lane[14]","$$lane[15]&plusmn;$$lane[16]","$$lane[17]&plusmn;$$lane[18]",
                     "$$lane[19]&plusmn;$$lane[20]","$$lane[21]&plusmn;$$lane[22]";
    push @rundata, \@lane_data;
  }
  print EASIH::HTML::advanced_table(\@rundata, 1, 1, 1, undef, undef, 600);


  my @data = (["Lane", "Sample", "total", "Pass filter", "% PF", "%QV30 bases"]);
  foreach my $lane ( sort { $$a[2] <=> $$b[2] } @lanes) {

#    print STDERR Dumper( $lane );
    my $QV30 = "NA";
    $QV30 = sprintf("'%s', '%s'", $$lane[6], $$lane[5]);
    $QV30 = sprintf("%.2f %%", $$lane[6]*100/$$lane[5]) if ( $$lane[6] && $$lane[5] );

      my @lvs = ("Lane $$lane[2]", " ", $$lane[3], $$lane[4], sprintf("%.2f %%", $$lane[4]*100/$$lane[3]), $QV30);

      if ($$lane[4] < 30_000_000 || $$lane[4]*100/$$lane[3] < 75 ) {
	my @clvs;
	foreach my $lv ( @lvs ) {
	  push @clvs, {value=> $lv, bgcolor=>'#CC0000'};
	}
	@lvs = @clvs;
      }
      else {
	my @clvs;
	foreach my $lv ( @lvs ) {
	  push @clvs, {value=> $lv, bgcolor=>'#00CC00'};
	}
	@lvs = @clvs;
      }
      
      push @data, \@lvs;


    my %sample_stats;
    foreach my $mplex (@mplexs ) {

      next if ( $$lane[2] != $$mplex[2] );
      
      $sample_stats{$$mplex[4]}{ perc  } = $$mplex[6];
      $sample_stats{$$mplex[4]}{ reads } = $$mplex[7];
      $sample_stats{$$mplex[4]}{ bcode } = $$mplex[5];
      push @{$sample_stats{$$mplex[4]}{fids}}, $$mplex[1];

      $sample_stats{$$mplex[4]}{ filename } = EASIH::DONE::fetch_filename( $$mplex[1] );
    }


    
    foreach my $sample (sort keys %sample_stats ) {

      $sample_stats{$sample}{fids} = join(",", @{$sample_stats{$sample}{fids}});
      
      $sample_stats{$sample}{ filename } =~ s/^.*\///;
      $sample_stats{$sample}{ filename } =~ s/\.gz//;
      $sample_stats{$sample}{ filename } =~ s/\.fq//;
      $sample_stats{$sample}{ filename } =~ s/\.[12]\z//;

      push @data, ["<a href=$0?QC=1&fid=$sample_stats{$sample}{fids}>$sample</a>", $sample_stats{$sample}{ filename }, $sample_stats{$sample}{ bcode }, $sample_stats{$sample}{ reads }, "$sample_stats{$sample}{ perc } %"] 

    }


  }
  print EASIH::HTML::advanced_table(\@data, 1, 1, 1, undef, undef, 600);

  easih_foot();


}
elsif( $EASIH::HTML::parameters{ 'project' } ) {

  easih_top();
  print "<h1> Projects in the database </h1>";

  my @data = (['Project name']);

  foreach my $project ( sort {$$b[1] cmp $$a[1] } EASIH::DONE::fetch_projects()) {
    push @data, ["<a href='$0?pid=$$project[0]'> $$project[1]</a>"];
  }

  
#  print "<CENTER>";
  print EASIH::HTML::advanced_table(\@data, 1, 5, 1, undef,undef, 600);
#  print "</CENTER>";
  easih_foot();

  print EASIH::HTML::end();
}
elsif( $EASIH::HTML::parameters{ 'pid' } ) {

  easih_top();
  my $name = EASIH::DONE::fetch_project_name( $EASIH::HTML::parameters{ 'pid' } );
  print "<h1> Files related to project $name in the database </h1>";

  my @data = (['Sample name', 'File name']);

  foreach my $project ( sort {$$b[1] cmp $$a[1] } EASIH::DONE::fetch_files_from_project( $EASIH::HTML::parameters{ 'pid' } )) {
    $$project[4] =~ s/.*\///;
    push @data, ["<a href='$0?sid=$$project[0]'> $$project[2]</a>", "<a href='$0?QC=1&fid=$$project[3]'> $$project[4]</a>"];
  }
			
  
  print EASIH::HTML::advanced_table(\@data, 1, 5, 1, undef,undef, 600);
  easih_foot();

  print EASIH::HTML::end();
}
elsif( $EASIH::HTML::parameters{ 'runs' } ) {

  easih_top();
  print "<h1> Runs in database </h1>";

  my @runs = (['Run', 'Platform', 'Samples']) ;
  foreach my $run ( sort {$$b[1] cmp $$a[1] } EASIH::DONE::fetch_runs()) {

    my $files = int( EASIH::DONE::fetch_samples_from_run($$run[1]));
    push @runs, ["<a href='$0?rid=$$run[0]'> $$run[1]</a>", $$run[2], $files];
  }

  
#  print "<CENTER>";
  print EASIH::HTML::advanced_table(\@runs, 1, 5, 1, undef,undef, 600);
#  print "</CENTER>";
  easih_foot();

  print EASIH::HTML::end();
}
elsif( $EASIH::HTML::parameters{ 'sid' } && $EASIH::HTML::parameters{ 'sid' } eq "all" ) {

  easih_top();
  print "<h1> Samples in database </h1>";

  my @runs = (['Name']) ;
  foreach my $run ( sort {$$b[2] cmp $$a[2] } EASIH::DONE::fetch_samples()) {
    push @runs, ["<a href='$0?sid=$$run[0]'> $$run[2]</a>"];
  }

  
#  print "<CENTER>";
  print EASIH::HTML::advanced_table(\@runs, 1, 5, 1, undef,undef, 600);
#  print "</CENTER>";
  easih_foot();

  print EASIH::HTML::end();
}
elsif( $EASIH::HTML::parameters{ 'sid' } ) {

  my $name = EASIH::DONE::fetch_sample_name(  $EASIH::HTML::parameters{ 'sid' } );
  easih_top();
  print "<h1> Files belonging to sample $name in database </h1>";

  my @data = (['File name', 'Run name']);

  foreach my $project ( sort {$$b[1] cmp $$a[1] } EASIH::DONE::fetch_files_from_sample( $EASIH::HTML::parameters{ 'sid' } )) {
    $$project[1] =~ s/.*\///;
    push @data, ["<a href='$0?QC=1&fid=$$project[0]'> $$project[1]</a>", "<a href='$0?rid=$$project[2]'> $$project[3]</a>"];
  }
			
  
  print EASIH::HTML::advanced_table(\@data, 1, 5, 1, undef,undef, 600);
  easih_foot();

  print EASIH::HTML::end();




}
elsif( $EASIH::HTML::parameters{ 'fid' } ) {

  if ( $EASIH::HTML::parameters{ 'fid' } eq "all") {

    easih_top();
    print "<h1> FASTQ files in the database: </h1>";


    my @data = (['Name', 'run','Total reads', 'Q30 bases',  'Duplicates', 'Partial Adaptors', 'Average AC']) ;
    my @files = EASIH::DONE::fetch_files($EASIH::HTML::parameters{ 'pid' }) ;
    map { $$_[2] =~ s/.*\///;} @files;
    
    my %run_cache;
    foreach my $file ( sort {$$b[2] cmp $$a[2] } @files) {
      
      my $run_name = $run_cache{ $$file[4] } if ($run_cache{ $$file[4]});
      if (! $run_name ) {
	$run_name = EASIH::DONE::fetch_run_name( $$file[4] );
	$run_cache{ $$file[4]} = $run_name if ($run_name );
	$run_name ||= "UNKNOWN";
      }
      
      my @line = ( "<a href=$0?QC=1&fid=$$file[0]>$$file[2]</a>", "<a href=$0?rid=$$file[4]>$run_name</a>", $$file[5]);

      my $colour = "red";
      $colour = 'yellow'      if( $$file[8] >= 70);
      $colour = '#00DD00'  if ( $$file[8] >= 90 );
      push @line, {value => "$$file[8]%", bgcolor=>$colour};

      $colour = 'red';
      $colour = 'yellow'   if ($$file[9] <10 );
      $colour = '#00DD00'  if ($$file[9] < 1 );
      push @line, {value=>"$$file[9]%", bgcolor => $colour};
      

      $colour = 'red';
      $colour = 'yellow'  if ($$file[10] < 10  );
      $colour = '#00DD00' if ($$file[10] < 1 );
      push @line, {value=>"$$file[10]%", bgcolor => $colour};

      $colour = 'red';
      $colour = 'yellow'  if ($$file[11] > 40 && $$file[11] < 60);
      $colour = '#00DD00' if ($$file[11] > 45 && $$file[11] < 55);
      push @line, {value=>"$$file[11]%", bgcolor => $colour};

      push @data, \@line;
    }
    
  
#  print "<CENTER>";
    print EASIH::HTML::advanced_table(\@data, 1, 5, 1, undef, undef, 850);
#  print "</CENTER>";
    easih_foot();

    print EASIH::HTML::end();
  }
}
else {

  easih_top();
  print "<h1> EASIH DONE (Data Offloading 'N' Evaluation) </h1>";

  my @runs     = EASIH::DONE::fetch_runs();
  my @projects = EASIH::DONE::fetch_projects();
  my @samples  = EASIH::DONE::fetch_samples();
  my @files    = EASIH::DONE::fetch_files();

  


  my $total_reads = 0;
  map { $total_reads += $$_[5] } @files;

  if ( $total_reads > 1_000_000_000_000_000) {
    $total_reads = sprintf("%.2f Peta-reads", $total_reads/ 1_000_000_000_000_000);
  }
  elsif ( $total_reads > 1_000_000_000_000) {
    $total_reads = sprintf("%.2f Tera-reads", $total_reads/ 1_000_000_000_000);
  }
  elsif ( $total_reads > 1_000_000_000) {
    $total_reads = sprintf("%.2f Giga-reads", $total_reads/ 1_000_000_000);
  }
  elsif ( $total_reads > 1_000_000) {
    $total_reads = sprintf("%.2f Mega-reads", $total_reads/ 1_000_000);
  }


  my $total_bases = 0;
  map { $total_bases += $$_[5]*$$_[6] } @files;

  if ( $total_bases > 1_000_000_000_000_000) {
    $total_bases = sprintf("$total_bases bp (%.2f Peta-bases)", $total_bases/ 1_000_000_000_000_000);
  }
  elsif ( $total_bases > 1_000_000_000_000) {
    $total_bases = sprintf("$total_bases bp (%.2f Tera-bases) ", $total_bases/ 1_000_000_000_000);
  }
  elsif ( $total_bases > 1_000_000_000) {
    $total_bases = sprintf("$total_bases bp (%.2f Giga-bases) ", $total_bases/ 1_000_000_000);
  }
  elsif ( $total_bases > 1_000_000) {
    $total_bases = sprintf("$total_bases bp (%.2f Mega-bases) ", $total_bases/ 1_000_000);
  }

  my @data = (['Runs', 'Projects', 'Samples', 'Files', 'Sequenced reads', 'Sequenced bases'],
	      [int(@runs), int(@projects), int(@samples), int(@files), $total_reads, $total_bases]);


  print "<p>\n";
  print "Current database statistics : <BR><BR>";
  print EASIH::HTML::advanced_table(\@data, 1, 1, 1, undef, undef, 800);
  print "<BR><BR>\n";
  print "Browse the data using the top menu bar, your choises are either by <a href='$0?runs=all'>runs</a>, <a href='$0?instr=all'>instrument</a>, <a href='$0?project=all'>project</a>, <a href='$0?sid=all'>sample</a> or <a href='$0?fid=all'>files</a>";
  print "</p>\n";
  
  easih_foot();
  print EASIH::HTML::end();
}





# 
# 
# 
# Kim Brugger (19 Jul 2011)
sub easih_top {

  print qq{
  <div id='topcontainer'>
     <div id='topcontainer'> <h1> <a id='topbanner' href='$0'></a></h1>
     </div>
  </div>

  <div id='navcontainer'>
    <div id='navarea'>
      <ul>
        <li ><a href=$0?runs=all> Runs </a></li>
        <li ><a href=$0?instr=all> Instrument </a></li>
        <li ><a href=$0?project=all> Project </a></li>
        <li ><a href=$0?sid=all> Sample </a></li>
        <li ><a href=$0?fid=all> File </a></li>
      </ul>
    </div>
  </div>

  <div id='contentcontainer'>
   <div id='contentarea'>
};
  
}



# 
# 
# 
# Kim Brugger (19 Jul 2011)
sub easih_foot {

  print qq{

      <br><br>
    </div>
  </div>
  
  <div id='footercontainer'>
    <div id='footer'>
     <div class='leftcol'> EASIH QC and run tracking database &copy; 2011-. All rights reserved.</div>
    </div>
  </div>
};
  
}


# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_basesplits {
  my ($fid, $domid) = @_;

  print '<script type="text/javascript">
      var queryString = "";
      var dataUrl = "";

      function onLoadCallback() {
        if (dataUrl.length > 0) {
          var query = new google.visualization.Query(dataUrl);
          query.setQuery(queryString);
          query.send(handleQueryResponse);
        } else {
          var dataTable = new google.visualization.DataTable();

         dataTable.addColumn("number");
         dataTable.addColumn("number");
         dataTable.addColumn("number");
         dataTable.addColumn("number");
         dataTable.addColumn("number");
';

  my @data;
  my @splits = EASIH::DONE::fetch_base_distribution( $fid );
  print"          dataTable.addRows(".(@splits  ).");";

  foreach my $line ( @splits ) {
    for (my $i=0;$i< 5;$i++) {
      print "          dataTable.setValue($$line[0], $i, ". $$line[$i+1] .");\n";
    }
  }

  print '
  draw(dataTable);
        }
      }

      function draw(dataTable) {
        var vis = new google.visualization.ImageChart(document.getElementById("'.$domid .'"));
        var options = {
          chxl: "",
          chxp: "",
          chxr: "",
          chxs: "",
          chxtc: "",
          chxt: "y,x",
          cht: "bvs",
          chco: "00CC00,0000CC,000000,CC0000,999999",
          chd: "s:UNYOYV,TZUTUS,999999",
          chdl: "A|C|G|T|N",
          chp: "0,0.017",
          chm: "r,FF0000,0,0,0",
          chtt: "Base splits pr cycle"
        };
        vis.draw(dataTable, options);
      }

      function handleQueryResponse(response) {
        if (response.isError()) {
          alert("Error in query: " + response.getMessage() + " " + response.getDetailedMessage());
          return;
        }
        draw(response.getDataTable());
      }

      google.load("visualization", "1", {packages:["imagechart"]});
      google.setOnLoadCallback(onLoadCallback);

</script> 
';  
}





# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_adaptor_histogram {
  my ($fid, $domid) = @_;

  my @data;
  my @adaptors = EASIH::DONE::fetch_adaptors( $fid );

  my %splits;
  foreach my $split ( @adaptors ) {
    $splits{$$split[0]} = $$split[1];
  }

  return if ( !@adaptors);

  print '
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn("string", "Position");
        data.addColumn("number", "Counts");
';

  print"          data.addRows(".(keys %splits ).");\n";


  for (my $i=0; $i < $adaptors[-1][0]; $i++){
    
    $splits{ $i } ||= 0;

    print "          data.setValue($i, 0, '$i');\n";
    print "          data.setValue($i, 1, $splits{$i});\n";
  }

  print '
  var chart = new google.visualization.ColumnChart(document.getElementById("'.$domid.'"));
        chart.draw(data, {title: "Partial adaptor mapping",
                          hAxis: {title: "Count"},
                          vAxis: {title: "Adaptor"}
                         });
      }
    </script>
';

}


# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_qv_histogram {
  my ($fid, $domid) = @_;


  my @data;
  my @splits = EASIH::DONE::fetch_qvs_histogram( $fid );

  my %splits;
  foreach my $split ( @splits ) {
    $splits{$$split[0]} = $$split[1];
  }

  return if ( !@splits);

  print '
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn("string", "QV Score");
        data.addColumn("number", "Counts");
';

  print"          data.addRows(".($splits[-1][0]  ).");\n";

  for (my $i=0;$i<@splits; $i++){
    
    $splits{ $i } ||= 0;

    print "          data.setValue($i, 0, '$i');\n";
    print "          data.setValue($i, 1, $splits{$i});\n";
  }

  print '
  var chart = new google.visualization.ColumnChart(document.getElementById("'.$domid.'"));
        chart.draw(data, {title: "QV split",
                          hAxis: {title: "Count"},
                          vAxis: {title: "QVs"}
                         });
      }
    </script>
';

}


# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_dups_histogram {
  my ($fid, $domid) = @_;


  my @dups = EASIH::DONE::fetch_duplicates( $fid );

  return if (!@dups);

  print '
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn("string", "QV Score");
        data.addColumn("number", "Counts");
';

  my %splits;
  foreach my $dup ( @dups ) {
    $splits{$$dup[0]} = $$dup[1];
  }

  print"          data.addRows(".($dups[-1][0] +1 ).");\n";

  for (my $i=2;$i<= $dups[-1][0]; $i++){
    
    $splits{ $i } ||= 0;

    print "          data.setValue($i, 0, '$i');\n";
    print "          data.setValue($i, 1, $splits{$i});\n";
  }

  print '
  var chart = new google.visualization.ColumnChart(document.getElementById("'.$domid.'"));
        chart.draw(data, {title: "Duplicates",
                          hAxis: {title: "Duplicates"},
                          vAxis: {title: "Sequences"}
                         });
      }
    </script>
';

}



# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_gc {
  my ($fid, $domid) = @_;


print '
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn("string", "GC");
        data.addColumn("number", "Percentage");
';

  my @data;
  my @gcs = EASIH::DONE::fetch_gc_distribution( $fid );
  
  my %splits;
  foreach my $gc ( @gcs ) {
    $splits{$$gc[0]} = $$gc[1];
  }
  
  print"          data.addRows(21);\n";
  
  for (my $i=0;$i < 21; $i++){
    my $bin = $i*5;
    $splits{ $bin } ||= 0;
    
    print "          data.setValue($i, 0, '$bin');\n";
    print "          data.setValue($i, 1, $splits{$bin});\n";
  }

print '
        var chart = new google.visualization.LineChart(document.getElementById("'.$domid.'"));
        chart.draw(data, {title: "%GC"});
      }
    </script>
';
}
# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_boxplot {
  my ($fid, $domid) = @_;

  print '
     <script type="text/javascript">
      function drawVisualization() {
        // Populate the data table.
        var dataTable = google.visualization.arrayToDataTable([';


  my @data;
  foreach my $line (EASIH::DONE::fetch_qvs_boxplot( $fid )) {
    next if ( $$line[1] < 0);
    push @data,  "[". join(",", "'".($$line[0])."'", $$line[2],$$line[3],$$line[5], $$line[6]) ."]";
    
  }

  print join(",\n", @data) . '], true)';

  print "
        // Draw the chart.
        var chart = new google.visualization.CandlestickChart(document.getElementById(\'$domid\'));
        chart.draw(dataTable, {legend:\'none\'});
      }

      google.setOnLoadCallback(drawVisualization);
    </script>
";
  
}



# 
# 
# 
# Kim Brugger (06 Jul 2011)
sub print_header {

print q{
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <script type="text/javascript">
      google.load("visualization", "1", {packages: ["corechart"]});
    </script>};
  
}
