<?php

$fIn1 = fopen("genotypes.txt", "r");
$fOut1 = fopen("genotypes.2.txt","w");
$fIn2 = fopen("outcomevars.txt", "r");
$fOut2 = fopen("outcomevars.2.txt","w");

for($i=0;$i<=2000;$i++) {

	$sLn1 = fgets($fIn1);
	$sLn2 = fgets($fIn2);
	if ($i<=500 || $i>=1500) {
		if ($sLn1 !=false ) {
			fwrite($fOut1 , $sLn1 );
		}
		
		if ($sLn2!=false ) {
			fwrite($fOut2 , $sLn2);
		}
	
	}
}

?>