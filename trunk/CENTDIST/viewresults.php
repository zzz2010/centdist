<?php
echo "TESTING";
$userdir='./userdata/';
if($_GET['uid']) 
 {	
	$ip=base64_decode($_GET['uid']);
 }
else if ($HTTP_SERVER_VARS["HTTP_X_FORWARDED_FOR"])
{
$ip = $HTTP_SERVER_VARS["HTTP_X_FORWARDED_FOR"];
}
elseif ($HTTP_SERVER_VARS["HTTP_CLIENT_IP"])
{
$ip = $HTTP_SERVER_VARS["HTTP_CLIENT_IP"];
}
elseif ($HTTP_SERVER_VARS["REMOTE_ADDR"])
{
$ip = $HTTP_SERVER_VARS["REMOTE_ADDR"];
}
elseif (getenv("HTTP_X_FORWARDED_FOR"))
{
$ip = getenv("HTTP_X_FORWARDED_FOR");
}
elseif (getenv("HTTP_CLIENT_IP"))
{
$ip = getenv("HTTP_CLIENT_IP");
}
elseif (getenv("REMOTE_ADDR"))
{
$ip = getenv("REMOTE_ADDR");
}
else
{
$ip = "Unknown";
}
      $userdir=$userdir.$ip;

        //  header( 'Location: http://compbio.ddns.comp.nus.edu.sg/~chipseq/Pomoda/'.$userdir.'/'.$_GET['runid'] ) ;
        $comp=explode("__", $_GET['runid']);
        $hashvalue=$comp[0];
        $resultLink='http://compbio.ddns.comp.nus.edu.sg/~chipseq/Pomoda/'.$userdir.'/'.$_GET['runid'];
        $fastafilelink= 'http://compbio.ddns.comp.nus.edu.sg/~chipseq/tmpfasta/'.$hashvalue.'.fa';
?>
<P><a href=<? echo $fastafilelink ?> >download sequence file(fasta format)</a></P>
<P><a href=<? echo $resultLink ?> >download raw result files</a></P>
<P><a href=<? echo "../Pomoda/showresult.php?runid=".$_GET['runid']?>>view de novo motif results</a></P>
<?php
$rundir="../Pomoda/$userdir/".$_GET['runid']."/scan";
//echo "TESTING $rundir<BR>";
function getcol($data,$colname){
	$out=array();
	foreach($data as $d){
		array_push($out,$d[$colname]);
	}
	return $out;
}
function image($link,$height){
	return "<img src=$link height=$height>";
}
function ahref($href,$s){
	return "<a href=$href>$s</a>";
}
function ahrefimage($link,$height){
	return ahref($link,image($link.".tn.jpg",$height));
}
function readtable($FileName){
$FileHandle = fopen($FileName,"r");
$FileContent = fread ($FileHandle,filesize ($FileName));
fclose($FileHandle);
// You can replace the \t with whichever delimiting character you are using
$lines = explode("\n", $FileContent);
$headerrow=1;
$header=array();
$data=array();
foreach($lines as $CurrValue)
{
	$cells=explode("\t",$CurrValue);
	if ($headerrow==1){
		$header=$cells;
		$headerrow=0;
		continue;
	}
	$i=0;
	$x=array();
	if (count($cells)<2) continue;
	foreach($cells as $cell){
		$x[$header[$i]]=$cell;
		$i++;
	}
	array_push($data,$x);
}
return $data;
}
$table=readtable("$rundir/results.txt");
array_multisort(getcol($table,"SCORE"),SORT_NUMERIC,SORT_DESC,$table);
print ("<TABLE border=1>\n");
print("<TR>");
print "<TD>RANK</TD>";
print "<TD>NAME</TD>";
print "<TD>SCORE<br></TD>";
print "<TD>Z0SCORE</TD>";
print "<TD>FULLWINDOW</TD>";
print "<TD>DIST HISTOGRAM</TD>";
print ("</TR>");
$rank=1;
foreach($table as $d){
	print ("<TR>");
	print "<TD>".$rank."</TD>";
	print "<TD>".$d["NAME"]."</TD>";
	print "<TD>".($d["SCORE"])."</TD>";
	print "<TD>".$d["Z0SCORE"]."</TD>";
	print "<TD>".$d["W"]."</TD>";
	print "<TD>".ahrefimage($rundir."/".$d["NAME"].".pos_disthist.jpg",'100%')."</TD>";
	print ("</TR>");
	print ("\n");
	$rank++;
}
print ("</TABLE>");
?>
