#!/usr/bin/env php

<?php

$uhfftffmr_path = 'src/uhfftfmmr.c';
$uhfftffmc_path = 'src/uhfftfmmc.c';


//get the content of the files as a string
$fft_r_contents = file_get_contents(__DIR__.'/'.$uhfftffmr_path);
$fft_c_contents = file_get_contents(__DIR__.'/'.$uhfftffmc_path);

//output file name
$f_name_r = "uhfft_r";
$f_name_c = "uhfft_c";

//default translation between C and C++
$default_rep = array(
	'REAL' => 'FFT_precision',         
	'Complex *x' => 'FFT_precision *Re, FFT_precision *Im',
	'"uhfftfmmr.h"' => '"'.$f_name_r.'.h"',
	'"uhfftfmmc.h"' => '"'.$f_name_c.'.h"',
);

//translate
translate($fft_r_contents, $default_rep, $f_name_r);
translate($fft_c_contents, $default_rep, $f_name_c);

$includes = '#include "bhfmm/fft/FFT_Settings.h"'."\n";
$function_param = '('.$default_rep['Complex *x'].');';
$r_h = $includes;
$c_h = $includes;
for($ord=6;$ord<17;$ord++) {
	$r_h .= 'void MFFTR' .$ord*2 .$function_param."\n";
	$r_h .= 'void MIFFTR'.$ord*2 .$function_param."\n";
	$c_h .= 'void MFFTC' .$ord*2 .$function_param."\n";
	$c_h .= 'void MIFFTC'.$ord*2 .$function_param."\n";
}
file_put_contents(__DIR__.'/output/'.$f_name_r.'.h',$r_h);
file_put_contents(__DIR__.'/output/'.$f_name_c.'.h',$c_h);

/**
 * Translate a function from C to C++
 * 
 * @param string $file_string        file_get_contents of the file containing the function
 * @param array  $default_replace    Basic translation between C and matlab (see the $default_rep array below)
 */
function translate($file_string, $default_replace, $output_name) {
	
	//do the default replacement
	$string1 = strtr($file_string,$default_replace);
	
	$translated = preg_replace('#(Re|Im)\(x\[(\d{1,3})\]\)#Us','$1[$2]',$string1);
	
	//write the result to a file in the output dir
	file_put_contents(__DIR__.'/output/'.$output_name.'.cpp',$translated);
}



?>
