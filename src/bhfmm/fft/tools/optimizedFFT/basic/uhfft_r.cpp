/*
*	UHFFT Library of FFT modules for FMM
*/


#include "uhfft_r.h"


/*
*	Number of additions = 64
*	Number of multiplications = 14
*	Number of sign changes = 13
*	Number of assigns = 90
*	Total number of operations = 181
*/
void	MFFTR12(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, /*tmp23r,*/ tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, /*tmp32r,*/ tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, /*tmp38r,*/ tmp39r, /*tmp40r,*/ tmp41r, /*tmp42r,*/
		/*tmp43r,*/ /*tmp44r,*/ /*tmp45r,*/ /*tmp46r,*/ /*tmp47r,*/ tmp48r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, /*tmp23i,*/ tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, /*tmp32i,*/ tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, /*tmp38i,*/ tmp39i, /*tmp40i,*/ tmp41i, /*tmp42i,*/
		/*tmp43i,*/ /*tmp44i,*/ /*tmp45i,*/ /*tmp46i,*/ /*tmp47i,*/ tmp48i;

	const FFT_precision	C3 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =                    2;	/* INTEGERCONST	*/

	tmp0i = C2*Im[3];
	tmp0r = -C2*Re[3];
	tmp1i = Im[4]-Im[2];
	tmp1r = Re[4]+Re[2];
	tmp2i = Im[4]+Im[2];
	tmp2r = Re[4]-Re[2];
	tmp3i = Im[1]+Im[5];
	tmp3r = Re[1]-Re[5];
	tmp4i = Im[1]-Im[5];
	tmp4r = Re[1]+Re[5];
	tmp5i = tmp1i+tmp3i;
	tmp5r = tmp1r+tmp3r;
	tmp6i = tmp1i-tmp3i;
	tmp6r = tmp1r-tmp3r;
	tmp7i = tmp2i-tmp4r;
	tmp7r = tmp2r+tmp4i;
	tmp8i = tmp2i+tmp4r;
	tmp8r = tmp2r-tmp4i;
	tmp9i = -(Im[4]-Im[2]);
	tmp9r = Re[4]+Re[2];
	tmp10i = -(Im[4]+Im[2]);
	tmp10r = Re[4]-Re[2];
	tmp11i = Im[5]+Im[1];
	tmp11r = Re[5]-Re[1];
	tmp12i = Im[5]-Im[1];
	tmp12r = Re[5]+Re[1];
	tmp13i = tmp9i+tmp11i;
	tmp13r = tmp9r+tmp11r;
	tmp14i = tmp9i-tmp11i;
	tmp14r = tmp9r-tmp11r;
	tmp15i = tmp10i-tmp12r;
	tmp15r = tmp10r+tmp12i;
	tmp16i = tmp10i+tmp12r;
	tmp16r = tmp10r-tmp12i;
	tmp17i = tmp5i+tmp13i;
	tmp17r = tmp5r+tmp13r;
	tmp18i = tmp5i-tmp13i;
	tmp18r = tmp5r-tmp13r;
	tmp19i = -C3*tmp17i;
	tmp19r = -C3*tmp17r;
	tmp20i = -C6*tmp18r;
	tmp20r = C6*tmp18i;
	tmp21i = tmp19i+tmp0i;
	tmp21r = tmp19r+Re[0];
	tmp22i = tmp21i+tmp20i;
	tmp22r = tmp21r+tmp20r;
	tmp24i = tmp0i+tmp17i;
	tmp24r = Re[0]+tmp17r;
	tmp25i = tmp8i+tmp16i;
	tmp25r = tmp8r+tmp16r;
	tmp26i = tmp8i-tmp16i;
	tmp26r = tmp8r-tmp16r;
	tmp27i = -C3*tmp25i;
	tmp27r = -C3*tmp25r;
	tmp28i = -C6*tmp26r;
	tmp28r = C6*tmp26i;
	tmp29i = tmp27i+tmp0r;
	tmp29r = tmp27r+Re[0];
	tmp30i = tmp29i+tmp28i;
	tmp30r = tmp29r+tmp28r;
	tmp31i = tmp29i-tmp28i;
	tmp31r = tmp29r-tmp28r;
	tmp33i = tmp6i+tmp14i;
	tmp33r = tmp6r+tmp14r;
	tmp34i = tmp6i-tmp14i;
	tmp34r = tmp6r-tmp14r;
	tmp35i = -C3*tmp33i;
	tmp35r = -C3*tmp33r;
	tmp36i = -C6*tmp34r;
	tmp36r = C6*tmp34i;
	tmp37i = tmp35i-tmp0i;
	tmp37r = tmp35r+Re[0];
	tmp39i = tmp37i-tmp36i;
	tmp39r = tmp37r-tmp36r;
	tmp41i = tmp7i+tmp15i;
	tmp41r = tmp7r+tmp15r;
	tmp48i = -(tmp0r-tmp41i);
	tmp48r = Re[0]+tmp41r;
	Re[0] = tmp24r;
	Im[0] = tmp24i;
	Re[1] = tmp30r;
	Im[1] = tmp30i;
	Re[2] = tmp39r;
	Im[2] = tmp39i;
	Re[3] = tmp48r;
	Im[3] = tmp48i;
	Re[4] = tmp22r;
	Im[4] = tmp22i;
	Re[5] = tmp31r;
	Im[5] = tmp31i;
}

/*
*	Number of additions = 27
*	Number of multiplications = 20
*	Number of sign changes = 10
*	Number of assigns = 58
*	Total number of operations = 115
*/
void	MIFFTR12(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, /*tmp24r,*/ tmp25r, tmp26r, tmp27r, tmp28r,
		/*tmp29r,*/ tmp30r, /*tmp31r,*/ /*tmp32r,*/ /*tmp33r,*/ /*tmp34r,*/ /*tmp35r,*/
		/*tmp36r,*/ tmp37r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, /*tmp14i,*/
		tmp15i, tmp16i, /*tmp17i,*/ /*tmp18i,*/ /*tmp19i,*/ /*tmp20i,*/ /*tmp21i,*/
		/*tmp22i,*/ tmp23i;

	const FFT_precision	C3 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =                    2;	/* INTEGERCONST	*/

	tmp0r = C2*Re[0];
	tmp0i = C2*Im[0];
	tmp1r = C2*Re[3];
	tmp1i = -C2*Im[3];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = C2*Re[4];
	tmp2i = C2*Im[4];
	tmp5r = C2*Re[1];
	tmp3i = C2*Im[1];
	tmp6r = tmp4r+tmp5r;
	tmp7r = tmp4r-tmp5r;
	tmp8r = C2*Re[2];
	tmp4i = -C2*Im[2];
	tmp9r = C2*Re[5];
	tmp5i = C2*Im[5];
	tmp10r = tmp8r+tmp9r;
	tmp11r = tmp8r-tmp9r;
	tmp12r = tmp6r+tmp10r;
	tmp13r = tmp6r-tmp10r;
	tmp14r = -C3*tmp12r;
	tmp6i = C6*tmp13r;
	tmp15r = tmp14r+tmp2r;
	tmp16r = tmp2r+tmp12r;
	tmp7i = tmp2i+tmp4i;
	tmp17r = tmp3i+tmp5i;
	tmp8i = tmp2i-tmp4i;
	tmp18r = tmp3i-tmp5i;
	tmp9i = -C3*tmp7i;
	tmp19r = -C3*tmp17r;
	tmp10i = C6*tmp18r;
	tmp20r = -C6*tmp8i;
	tmp11i = tmp9i+tmp0i;
	tmp21r = tmp19r+tmp1i;
	tmp12i = tmp11i+tmp10i;
	tmp22r = tmp21r+tmp20r;
	tmp13i = tmp11i-tmp10i;
	tmp23r = tmp21r-tmp20r;
	tmp25r = tmp7r+tmp11r;
	tmp26r = tmp7r-tmp11r;
	tmp27r = -C3*tmp25r;
	tmp15i = C6*tmp26r;
	tmp28r = tmp27r+tmp3r;
	tmp16i = tmp2i+tmp4i;
	tmp30r = -(tmp3i+tmp5i);
	tmp23i = tmp0i+tmp16i;
	tmp37r = -(tmp1i-tmp30r);
	Re[0] = tmp16r;
	Im[0] = 0.0;
	Re[1] = tmp22r;
	Im[1] = tmp12i;
	Re[2] = tmp28r;
	Im[2] = -tmp15i;
	Re[3] = tmp37r;
	Im[3] = tmp23i;
	Re[4] = tmp15r;
	Im[4] = tmp6i;
	Re[5] = tmp23r;
	Im[5] = tmp13i;
}

/*
*	Number of additions = 44
*	Number of multiplications = 48
*	Number of sign changes = 6
*	Number of assigns = 55
*	Total number of operations = 153
*/
void	MFFTR14(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i;

	const FFT_precision	C13 =     0.22252093395631;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.43388373911756;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =     0.62348980185873;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.78183148246803;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.90096886790242;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =     0.97492791218182;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =                    2;	/* INTEGERCONST	*/

	tmp0r = C16*Re[6];
	tmp0i = -C16*Im[6];
	tmp1r = C16*Re[4];
	tmp1i = -C16*Im[4];
	tmp2r = C16*Re[2];
	tmp2i = C16*Im[2];
	tmp3r = C9*tmp0r-C11*tmp1r-C13*tmp2r+Re[0];
	tmp4r = tmp0r+tmp1r;
	tmp5r = -C11*tmp0r-C13*tmp1r+C9*tmp2r+Re[0];
	tmp6r = tmp4r+tmp2r;
	tmp7r = -C13*tmp0r+C9*tmp1r-C11*tmp2r+Re[0];
	tmp8r = C10*tmp0i+C12*tmp1i+C14*tmp2i;
	tmp9r = C12*tmp0i+C14*tmp1i-C10*tmp2i;
	tmp10r = C14*tmp0i-C10*tmp1i-C12*tmp2i;
	tmp11r = tmp3r+tmp8r;
	tmp12r = tmp3r-tmp8r;
	tmp13r = tmp5r+tmp9r;
	tmp14r = tmp5r-tmp9r;
	tmp15r = tmp7r+tmp10r;
	tmp16r = tmp7r-tmp10r;
	tmp17r = Re[0]+tmp6r;
	tmp3i = C16*Im[1];
	tmp18r = C16*Re[1];
	tmp4i = C16*Im[3];
	tmp19r = C16*Re[3];
	tmp5i = C16*Im[5];
	tmp20r = -C16*Re[5];
	tmp6i = C9*tmp3i-C11*tmp4i-C13*tmp5i;
	tmp7i = tmp3i+tmp4i;
	tmp8i = -C11*tmp3i-C13*tmp4i+C9*tmp5i;
	tmp9i = tmp7i+tmp5i;
	tmp10i = -C13*tmp3i+C9*tmp4i-C11*tmp5i;
	tmp11i = -C10*tmp18r-C12*tmp19r-C14*tmp20r;
	tmp12i = -C12*tmp18r-C14*tmp19r+C10*tmp20r;
	tmp13i = -C14*tmp18r+C10*tmp19r+C12*tmp20r;
	tmp14i = tmp6i+tmp11i;
	tmp15i = tmp6i-tmp11i;
	tmp16i = tmp8i+tmp12i;
	tmp17i = tmp8i-tmp12i;
	tmp18i = tmp10i+tmp13i;
	tmp19i = tmp10i-tmp13i;
	Re[0] = tmp17r;
	Im[0] = tmp9i;
	Re[1] = tmp14r;
	Im[1] = -tmp17i;
	Re[2] = tmp11r;
	Im[2] = tmp14i;
	Re[3] = tmp16r;
	Im[3] = -tmp19i;
	Re[4] = tmp15r;
	Im[4] = tmp18i;
	Re[5] = tmp12r;
	Im[5] = -tmp15i;
	Re[6] = tmp13r;
	Im[6] = tmp16i;
}

/*
*	Number of additions = 127
*	Number of multiplications = 72
*	Number of sign changes = 6
*	Number of assigns = 105
*	Total number of operations = 310
*/
void	MIFFTR14(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		/*tmp43r,*/ /*tmp44r,*/ tmp45r, tmp46r, /*tmp47r,*/ /*tmp48r,*/ tmp49r,
		tmp50r, /*tmp51r,*/ /*tmp52r,*/ tmp53r, tmp54r /*tmp55r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, /*tmp7i,*/
		tmp8i, /*tmp9i,*/ tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, /*tmp20i,*/ tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, /*tmp28i,*/
		tmp29i, /*tmp30i,*/ tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, /*tmp41i,*/ /*tmp42i,*/
		/*tmp43i,*/ /*tmp44i,*/ tmp45i, tmp46i, /*tmp47i,*/ /*tmp48i,*/ tmp49i,
		tmp50i, /*tmp51i,*/ /*tmp52i,*/ tmp53i, tmp54i /*tmp55i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C13 =     0.22252093395631;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.43388373911756;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =     0.62348980185873;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.78183148246803;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.90096886790242;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =     0.97492791218182;	/* FFT_precisionCONST	*/

	tmp0i = -(Im[1]-Im[6]);
	tmp0r = Re[1]+Re[6];
	tmp1i = -(Im[1]+Im[6]);
	tmp1r = Re[1]-Re[6];
	tmp2i = -(Im[3]-Im[4]);
	tmp2r = Re[3]+Re[4];
	tmp3i = -(Im[3]+Im[4]);
	tmp3r = Re[3]-Re[4];
	tmp4i = Im[2]-Im[5];
	tmp4r = Re[2]+Re[5];
	tmp5i = Im[2]+Im[5];
	tmp5r = Re[2]-Re[5];
	tmp6i = C9*tmp0i-C11*tmp2i-C13*tmp4i+Im[0];
	tmp6r = C9*tmp0r-C11*tmp2r-C13*tmp4r+Re[0];
	tmp7r = tmp0r+tmp2r;
	tmp8i = -C11*tmp0i-C13*tmp2i+C9*tmp4i+Im[0];
	tmp8r = -C11*tmp0r-C13*tmp2r+C9*tmp4r+Re[0];
	tmp9r = tmp7r+tmp4r;
	tmp10i = -C13*tmp0i+C9*tmp2i-C11*tmp4i+Im[0];
	tmp10r = -C13*tmp0r+C9*tmp2r-C11*tmp4r+Re[0];
	tmp11i = C10*tmp1r+C12*tmp3r+C14*tmp5r;
	tmp11r = -C10*tmp1i-C12*tmp3i-C14*tmp5i;
	tmp12i = C12*tmp1r+C14*tmp3r-C10*tmp5r;
	tmp12r = -C12*tmp1i-C14*tmp3i+C10*tmp5i;
	tmp13i = C14*tmp1r-C10*tmp3r-C12*tmp5r;
	tmp13r = -C14*tmp1i+C10*tmp3i+C12*tmp5i;
	tmp14i = tmp6i+tmp11i;
	tmp14r = tmp6r+tmp11r;
	tmp15i = tmp6i-tmp11i;
	tmp15r = tmp6r-tmp11r;
	tmp16i = tmp8i+tmp12i;
	tmp16r = tmp8r+tmp12r;
	tmp17i = tmp8i-tmp12i;
	tmp17r = tmp8r-tmp12r;
	tmp18i = tmp10i+tmp13i;
	tmp18r = tmp10r+tmp13r;
	tmp19i = tmp10i-tmp13i;
	tmp19r = tmp10r-tmp13r;
	tmp20r = Re[0]+tmp9r;
	tmp21i = Im[1]-Im[6];
	tmp21r = Re[1]+Re[6];
	tmp22i = Im[1]+Im[6];
	tmp22r = Re[1]-Re[6];
	tmp23i = Im[3]-Im[4];
	tmp23r = Re[3]+Re[4];
	tmp24i = Im[3]+Im[4];
	tmp24r = Re[3]-Re[4];
	tmp25i = -(Im[2]-Im[5]);
	tmp25r = Re[2]+Re[5];
	tmp26i = -(Im[2]+Im[5]);
	tmp26r = Re[2]-Re[5];
	tmp27i = C9*tmp21i-C11*tmp23i-C13*tmp25i-Im[0];
	tmp27r = C9*tmp21r-C11*tmp23r-C13*tmp25r+Re[0];
	tmp28r = tmp21r+tmp23r;
	tmp29i = -C11*tmp21i-C13*tmp23i+C9*tmp25i-Im[0];
	tmp29r = -C11*tmp21r-C13*tmp23r+C9*tmp25r+Re[0];
	tmp30r = tmp28r+tmp25r;
	tmp31i = -C13*tmp21i+C9*tmp23i-C11*tmp25i-Im[0];
	tmp31r = -C13*tmp21r+C9*tmp23r-C11*tmp25r+Re[0];
	tmp32i = C10*tmp22r+C12*tmp24r+C14*tmp26r;
	tmp32r = -C10*tmp22i-C12*tmp24i-C14*tmp26i;
	tmp33i = C12*tmp22r+C14*tmp24r-C10*tmp26r;
	tmp33r = -C12*tmp22i-C14*tmp24i+C10*tmp26i;
	tmp34i = C14*tmp22r-C10*tmp24r-C12*tmp26r;
	tmp34r = -C14*tmp22i+C10*tmp24i+C12*tmp26i;
	tmp35i = tmp27i+tmp32i;
	tmp35r = tmp27r+tmp32r;
	tmp36i = tmp27i-tmp32i;
	tmp36r = tmp27r-tmp32r;
	tmp37i = tmp29i+tmp33i;
	tmp37r = tmp29r+tmp33r;
	tmp38i = tmp29i-tmp33i;
	tmp38r = tmp29r-tmp33r;
	tmp39i = tmp31i+tmp34i;
	tmp39r = tmp31r+tmp34r;
	tmp40i = tmp31i-tmp34i;
	tmp40r = tmp31r-tmp34r;
	tmp41r = Re[0]+tmp30r;
	tmp42r = tmp20r+tmp41r;
	tmp45i = tmp17i-tmp38i;
	tmp45r = tmp17r-tmp38r;
	tmp46i = tmp14i+tmp35i;
	tmp46r = tmp14r+tmp35r;
	tmp49i = tmp19i-tmp40i;
	tmp49r = tmp19r-tmp40r;
	tmp50i = tmp18i+tmp39i;
	tmp50r = tmp18r+tmp39r;
	tmp53i = tmp15i-tmp36i;
	tmp53r = tmp15r-tmp36r;
	tmp54i = tmp16i+tmp37i;
	tmp54r = tmp16r+tmp37r;
	Re[0] = tmp42r;
	Im[0] = C0;
	Re[1] = tmp45r;
	Im[1] = tmp45i;
	Re[2] = tmp46r;
	Im[2] = tmp46i;
	Re[3] = tmp49r;
	Im[3] = tmp49i;
	Re[4] = tmp50r;
	Im[4] = tmp50i;
	Re[5] = tmp53r;
	Im[5] = tmp53i;
	Re[6] = tmp54r;
	Im[6] = tmp54i;
}

/*
*	Number of additions = 108
*	Number of multiplications = 26
*	Number of sign changes = 4
*	Number of assigns = 130
*	Total number of operations = 268
*/
void	MFFTR16(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r /*tmp63r,*/
		/*tmp64r,*/ /*tmp65r,*/ /*tmp66r,*/ /*tmp67r,*/ /*tmp68r,*/ /*tmp69r,*/ /*tmp70r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i /*tmp51i,*/ /*tmp52i,*/ /*tmp53i,*/ /*tmp54i,*/ /*tmp55i,*/ /*tmp56i,*/
		/*tmp57i,*/ /*tmp58i*/;

	const FFT_precision	C3 =     0.38268343236509;	/* FFT_precisionCONST	*/
	const FFT_precision	C4 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =     0.92387953251129;	/* FFT_precisionCONST	*/
	const FFT_precision	C5 =                    2;	/* INTEGERCONST	*/

	tmp0r = C5*Re[4];
	tmp0i = C5*Im[4];
	tmp1r = Re[0]+tmp0r;
	tmp2r = Re[0]-tmp0r;
	tmp3r = Re[0]+tmp0i;
	tmp4r = Re[0]-tmp0i;
	tmp1i = Im[2]-Im[6];
	tmp5r = Re[2]+Re[6];
	tmp2i = Im[2]+Im[6];
	tmp6r = Re[2]-Re[6];
	tmp3i = Im[6]-Im[2];
	tmp7r = Re[6]+Re[2];
	tmp4i = Im[6]+Im[2];
	tmp8r = Re[6]-Re[2];
	tmp5i = tmp1i+tmp3i;
	tmp9r = tmp5r+tmp7r;
	tmp6i = tmp1i-tmp3i;
	tmp10r = tmp5r-tmp7r;
	tmp7i = tmp2i-tmp8r;
	tmp11r = tmp6r+tmp4i;
	tmp8i = tmp2i+tmp8r;
	tmp12r = tmp6r-tmp4i;
	tmp9i = C4*(tmp7i-tmp11r);
	tmp13r = C4*(tmp11r+tmp7i);
	tmp10i = -C4*(tmp8i+tmp12r);
	tmp14r = -C4*(tmp12r-tmp8i);
	tmp15r = tmp1r+tmp9r;
	tmp16r = tmp1r-tmp9r;
	tmp17r = tmp3r+tmp13r;
	tmp18r = tmp3r-tmp13r;
	tmp19r = tmp2r+tmp6i;
	tmp20r = tmp2r-tmp6i;
	tmp21r = tmp4r+tmp14r;
	tmp22r = tmp4r-tmp14r;
	tmp11i = Im[1]+Im[7];
	tmp23r = Re[1]-Re[7];
	tmp12i = Im[1]-Im[7];
	tmp24r = Re[1]+Re[7];
	tmp13i = Im[5]+Im[3];
	tmp25r = Re[5]-Re[3];
	tmp14i = Im[5]-Im[3];
	tmp26r = Re[5]+Re[3];
	tmp15i = tmp11i+tmp13i;
	tmp27r = tmp23r+tmp25r;
	tmp16i = tmp11i-tmp13i;
	tmp28r = tmp23r-tmp25r;
	tmp17i = tmp12i-tmp26r;
	tmp29r = tmp24r+tmp14i;
	tmp18i = tmp12i+tmp26r;
	tmp30r = tmp24r-tmp14i;
	tmp31r = tmp27r;
	tmp19i = tmp15i;
	tmp32r = C2*tmp29r+C3*tmp17i;
	tmp20i = C2*tmp17i-C3*tmp29r;
	tmp33r = C4*(tmp28r+tmp16i);
	tmp21i = C4*(tmp16i-tmp28r);
	tmp34r = C3*tmp30r+C2*tmp18i;
	tmp22i = C3*tmp18i-C2*tmp30r;
	tmp23i = Im[3]+Im[5];
	tmp35r = Re[3]-Re[5];
	tmp24i = Im[3]-Im[5];
	tmp36r = Re[3]+Re[5];
	tmp25i = Im[7]+Im[1];
	tmp37r = Re[7]-Re[1];
	tmp26i = Im[7]-Im[1];
	tmp38r = Re[7]+Re[1];
	tmp27i = tmp23i+tmp25i;
	tmp39r = tmp35r+tmp37r;
	tmp28i = tmp23i-tmp25i;
	tmp40r = tmp35r-tmp37r;
	tmp29i = tmp24i-tmp38r;
	tmp41r = tmp36r+tmp26i;
	tmp30i = tmp24i+tmp38r;
	tmp42r = tmp36r-tmp26i;
	tmp43r = tmp39r;
	tmp31i = tmp27i;
	tmp44r = C3*tmp41r+C2*tmp29i;
	tmp32i = C3*tmp29i-C2*tmp41r;
	tmp45r = -C4*(tmp40r-tmp28i);
	tmp33i = -C4*(tmp28i+tmp40r);
	tmp46r = -C2*tmp42r-C3*tmp30i;
	tmp34i = -C2*tmp30i+C3*tmp42r;
	tmp47r = tmp31r+tmp43r;
	tmp35i = tmp19i+tmp31i;
	tmp48r = tmp32r+tmp44r;
	tmp36i = tmp20i+tmp32i;
	tmp49r = tmp33r+tmp45r;
	tmp37i = tmp21i+tmp33i;
	tmp50r = tmp34r+tmp46r;
	tmp38i = tmp22i+tmp34i;
	tmp51r = -tmp19i+tmp31i;
	tmp39i = tmp31r-tmp43r;
	tmp52r = -tmp20i+tmp32i;
	tmp40i = tmp32r-tmp44r;
	tmp53r = -tmp21i+tmp33i;
	tmp41i = tmp33r-tmp45r;
	tmp54r = -tmp22i+tmp34i;
	tmp42i = tmp34r-tmp46r;
	tmp55r = tmp15r+tmp47r;
	tmp43i = tmp5i+tmp35i;
	tmp56r = tmp17r+tmp48r;
	tmp44i = tmp9i+tmp36i;
	tmp57r = tmp19r+tmp49r;
	tmp45i = -tmp10r+tmp37i;
	tmp58r = tmp21r+tmp50r;
	tmp46i = tmp10i+tmp38i;
	tmp59r = tmp16r-tmp51r;
	tmp47i = -tmp5i-tmp39i;
	tmp60r = tmp18r-tmp52r;
	tmp48i = -tmp9i-tmp40i;
	tmp61r = tmp20r-tmp53r;
	tmp49i = tmp10r-tmp41i;
	tmp62r = tmp22r-tmp54r;
	tmp50i = -tmp10i-tmp42i;
	Re[0] = tmp55r;
	Im[0] = tmp43i;
	Re[1] = tmp56r;
	Im[1] = tmp44i;
	Re[2] = tmp57r;
	Im[2] = tmp45i;
	Re[3] = tmp58r;
	Im[3] = tmp46i;
	Re[4] = tmp59r;
	Im[4] = tmp47i;
	Re[5] = tmp60r;
	Im[5] = tmp48i;
	Re[6] = tmp61r;
	Im[6] = tmp49i;
	Re[7] = tmp62r;
	Im[7] = tmp50i;
}

/*
*	Number of additions = 58
*	Number of multiplications = 40
*	Number of sign changes = 11
*	Number of assigns = 97
*	Total number of operations = 206
*/
void	MIFFTR16(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r /*tmp48r,*/ /*tmp49r,*/
		/*tmp50r,*/ /*tmp51r,*/ /*tmp52r,*/ /*tmp53r,*/ /*tmp54r,*/ /*tmp55r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, /*tmp12i,*/ tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, /*tmp18i,*/ tmp19i, tmp20i, tmp21i,
		/*tmp22i,*/ tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, /*tmp30i,*/ tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i /*tmp38i,*/ /*tmp39i,*/ /*tmp40i,*/ /*tmp41i,*/ /*tmp42i,*/
		/*tmp43i,*/ /*tmp44i,*/ /*tmp45i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C3 =     0.38268343236509;	/* FFT_precisionCONST	*/
	const FFT_precision	C4 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =     0.92387953251129;	/* FFT_precisionCONST	*/
	const FFT_precision	C5 =                    2;	/* INTEGERCONST	*/

	tmp0r = C5*Re[0];
	tmp0i = C5*Im[0];
	tmp1r = C5*Re[4];
	tmp1i = C5*Im[4];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = C5*Re[2];
	tmp2i = C5*Im[2];
	tmp5r = C5*Re[6];
	tmp3i = C5*Im[6];
	tmp6r = tmp4r+tmp5r;
	tmp7r = tmp4r-tmp5r;
	tmp4i = C4*(tmp2i-tmp3i);
	tmp8r = -C4*(tmp3i+tmp2i);
	tmp5i = -C4*(tmp2i-tmp3i);
	tmp9r = -C4*(tmp3i+tmp2i);
	tmp10r = tmp2r+tmp6r;
	tmp11r = tmp2r-tmp6r;
	tmp6i = tmp0i+tmp4i;
	tmp12r = -(tmp1i-tmp8r);
	tmp7i = tmp0i-tmp4i;
	tmp13r = -(tmp1i+tmp8r);
	tmp8i = tmp0i+tmp5i;
	tmp14r = tmp1i+tmp9r;
	tmp9i = tmp0i-tmp5i;
	tmp15r = tmp1i-tmp9r;
	tmp16r = C5*Re[1];
	tmp10i = C5*Im[1];
	tmp17r = C5*Re[5];
	tmp11i = C5*Im[5];
	tmp18r = tmp16r+tmp17r;
	tmp19r = tmp16r-tmp17r;
	tmp20r = tmp18r;
	tmp21r = -C2*tmp11i-C3*tmp10i;
	tmp13i = C2*tmp10i-C3*tmp11i;
	tmp22r = C4*tmp19r;
	tmp14i = C4*tmp19r;
	tmp23r = C3*tmp11i-C2*tmp10i;
	tmp15i = C3*tmp10i+C2*tmp11i;
	tmp24r = C5*Re[3];
	tmp16i = C5*Im[3];
	tmp25r = C5*Re[7];
	tmp17i = C5*Im[7];
	tmp26r = tmp24r+tmp25r;
	tmp27r = tmp24r-tmp25r;
	tmp28r = tmp26r;
	tmp29r = -C3*tmp17i-C2*tmp16i;
	tmp19i = C3*tmp16i-C2*tmp17i;
	tmp30r = -C4*tmp27r;
	tmp20i = C4*tmp27r;
	tmp31r = -C2*tmp17i+C3*tmp16i;
	tmp21i = -C2*tmp16i-C3*tmp17i;
	tmp32r = tmp20r+tmp28r;
	tmp33r = tmp21r+tmp29r;
	tmp23i = tmp13i+tmp19i;
	tmp34r = tmp22r+tmp30r;
	tmp24i = tmp14i+tmp20i;
	tmp35r = tmp23r+tmp31r;
	tmp25i = tmp15i+tmp21i;
	tmp36r = C0;
	tmp26i = -(tmp20r-tmp28r);
	tmp37r = tmp13i-tmp19i;
	tmp27i = -(tmp21r-tmp29r);
	tmp38r = tmp14i-tmp20i;
	tmp28i = -(tmp22r-tmp30r);
	tmp39r = tmp15i-tmp21i;
	tmp29i = -(tmp23r-tmp31r);
	tmp40r = tmp10r+tmp32r;
	tmp41r = tmp12r+tmp33r;
	tmp31i = tmp6i+tmp23i;
	tmp42r = tmp3r+tmp34r;
	tmp32i = tmp7r+tmp24i;
	tmp43r = tmp14r+tmp35r;
	tmp33i = tmp8i+tmp25i;
	tmp44r = tmp11r-tmp36r;
	tmp34i = -tmp26i;
	tmp45r = tmp13r-tmp37r;
	tmp35i = tmp7i-tmp27i;
	tmp46r = tmp3r-tmp38r;
	tmp36i = -tmp7r-tmp28i;
	tmp47r = tmp15r-tmp39r;
	tmp37i = tmp9i-tmp29i;
	Re[0] = tmp40r;
	Im[0] = 0.0;
	Re[1] = tmp41r;
	Im[1] = tmp31i;
	Re[2] = tmp42r;
	Im[2] = tmp32i;
	Re[3] = tmp43r;
	Im[3] = tmp33i;
	Re[4] = tmp44r;
	Im[4] = tmp34i;
	Re[5] = tmp45r;
	Im[5] = tmp35i;
	Re[6] = tmp46r;
	Im[6] = tmp36i;
	Re[7] = tmp47r;
	Im[7] = tmp37i;
}

/*
*	Number of additions = 148
*	Number of multiplications = 80
*	Number of sign changes = 40
*	Number of assigns = 214
*	Total number of operations = 482
*/
void	MFFTR18(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, /*tmp91r,*/
		/*tmp92r,*/ tmp93r, tmp94r, /*tmp95r,*/ /*tmp96r,*/ tmp97r, tmp98r,
		/*tmp99r,*/ /*tmp100r,*/ tmp101r, tmp102r, /*tmp103r,*/ /*tmp104r,*/ tmp105r,
		tmp106r /*tmp107r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, tmp86i, tmp87i, tmp88i, /*tmp89i,*/ /*tmp90i,*/ tmp91i,
		tmp92i, /*tmp93i,*/ /*tmp94i,*/ tmp95i, tmp96i, /*tmp97i,*/ /*tmp98i,*/
		tmp99i, tmp100i, /*tmp101i,*/ /*tmp102i,*/ tmp103i, tmp104i /*tmp105i*/;

	const FFT_precision	C9 =     0.17364817766693;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.34202014332567;	/* FFT_precisionCONST	*/
	const FFT_precision	C1 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.64278760968654;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =     0.76604444311898;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.93969262078591;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.98480775301221;	/* FFT_precisionCONST	*/
	const FFT_precision	C5 =                    2;	/* INTEGERCONST	*/

	tmp0r = C5*Re[6];
	tmp0i = -C5*Im[6];
	tmp1r = -C1*tmp0r;
	tmp2r = C6*tmp0i;
	tmp3r = tmp1r+Re[0];
	tmp4r = tmp3r+tmp2r;
	tmp5r = tmp3r-tmp2r;
	tmp6r = Re[0]+tmp0r;
	tmp1i = Im[4]-Im[2];
	tmp7r = Re[4]+Re[2];
	tmp2i = Im[4]+Im[2];
	tmp8r = Re[4]-Re[2];
	tmp3i = -C1*tmp1i;
	tmp9r = -C1*tmp7r;
	tmp4i = -C6*tmp8r;
	tmp10r = C6*tmp2i;
	tmp5i = tmp3i-Im[8];
	tmp11r = tmp9r+Re[8];
	tmp6i = tmp5i+tmp4i;
	tmp12r = tmp11r+tmp10r;
	tmp7i = tmp5i-tmp4i;
	tmp13r = tmp11r-tmp10r;
	tmp8i = -(Im[8]-tmp1i);
	tmp14r = Re[8]+tmp7r;
	tmp9i = -(Im[4]-Im[8]);
	tmp15r = Re[4]+Re[8];
	tmp10i = -(Im[4]+Im[8]);
	tmp16r = Re[4]-Re[8];
	tmp11i = -C1*tmp9i;
	tmp17r = -C1*tmp15r;
	tmp12i = -C6*tmp16r;
	tmp18r = C6*tmp10i;
	tmp13i = tmp11i+Im[2];
	tmp19r = tmp17r+Re[2];
	tmp14i = tmp13i+tmp12i;
	tmp20r = tmp19r+tmp18r;
	tmp15i = tmp13i-tmp12i;
	tmp21r = tmp19r-tmp18r;
	tmp16i = Im[2]+tmp9i;
	tmp22r = Re[2]+tmp15r;
	tmp17i = C7*tmp6i-C8*tmp12r;
	tmp23r = C7*tmp12r+C8*tmp6i;
	tmp18i = C9*tmp7i-C10*tmp13r;
	tmp24r = C9*tmp13r+C10*tmp7i;
	tmp19i = C9*tmp14i-C10*tmp20r;
	tmp25r = C9*tmp20r+C10*tmp14i;
	tmp20i = -C11*tmp15i-C12*tmp21r;
	tmp26r = -C11*tmp21r+C12*tmp15i;
	tmp21i = tmp8i+tmp16i;
	tmp27r = tmp14r+tmp22r;
	tmp22i = tmp8i-tmp16i;
	tmp28r = tmp14r-tmp22r;
	tmp23i = -C1*tmp21i;
	tmp29r = -C1*tmp27r;
	tmp24i = -C6*tmp28r;
	tmp30r = C6*tmp22i;
	tmp31r = tmp29r+tmp6r;
	tmp25i = tmp23i+tmp24i;
	tmp32r = tmp31r+tmp30r;
	tmp26i = tmp23i-tmp24i;
	tmp33r = tmp31r-tmp30r;
	tmp34r = tmp6r+tmp27r;
	tmp27i = tmp17i+tmp19i;
	tmp35r = tmp23r+tmp25r;
	tmp28i = tmp17i-tmp19i;
	tmp36r = tmp23r-tmp25r;
	tmp29i = -C1*tmp27i;
	tmp37r = -C1*tmp35r;
	tmp30i = -C6*tmp36r;
	tmp38r = C6*tmp28i;
	tmp39r = tmp37r+tmp4r;
	tmp31i = tmp29i+tmp30i;
	tmp40r = tmp39r+tmp38r;
	tmp32i = tmp29i-tmp30i;
	tmp41r = tmp39r-tmp38r;
	tmp42r = tmp4r+tmp35r;
	tmp33i = tmp18i+tmp20i;
	tmp43r = tmp24r+tmp26r;
	tmp34i = tmp18i-tmp20i;
	tmp44r = tmp24r-tmp26r;
	tmp35i = -C1*tmp33i;
	tmp45r = -C1*tmp43r;
	tmp36i = -C6*tmp44r;
	tmp46r = C6*tmp34i;
	tmp47r = tmp45r+tmp5r;
	tmp37i = tmp35i+tmp36i;
	tmp48r = tmp47r+tmp46r;
	tmp38i = tmp35i-tmp36i;
	tmp49r = tmp47r-tmp46r;
	tmp50r = tmp5r+tmp43r;
	tmp39i = C5*Im[3];
	tmp51r = C5*Re[3];
	tmp40i = -C1*tmp39i;
	tmp41i = -C6*tmp51r;
	tmp42i = tmp40i+tmp41i;
	tmp43i = tmp40i-tmp41i;
	tmp44i = Im[5]+Im[7];
	tmp52r = -(Re[5]-Re[7]);
	tmp45i = Im[5]-Im[7];
	tmp53r = -(Re[5]+Re[7]);
	tmp46i = -C1*tmp44i;
	tmp54r = -C1*tmp52r;
	tmp47i = -C6*tmp53r;
	tmp55r = C6*tmp45i;
	tmp48i = tmp46i+Im[1];
	tmp56r = tmp54r+Re[1];
	tmp49i = tmp48i+tmp47i;
	tmp57r = tmp56r+tmp55r;
	tmp50i = tmp48i-tmp47i;
	tmp58r = tmp56r-tmp55r;
	tmp51i = Im[1]+tmp44i;
	tmp59r = Re[1]+tmp52r;
	tmp52i = Im[5]+Im[1];
	tmp60r = Re[5]-Re[1];
	tmp53i = Im[5]-Im[1];
	tmp61r = Re[5]+Re[1];
	tmp54i = -C1*tmp52i;
	tmp62r = -C1*tmp60r;
	tmp55i = -C6*tmp61r;
	tmp63r = C6*tmp53i;
	tmp56i = tmp54i+Im[7];
	tmp64r = tmp62r-Re[7];
	tmp57i = tmp56i+tmp55i;
	tmp65r = tmp64r+tmp63r;
	tmp58i = tmp56i-tmp55i;
	tmp66r = tmp64r-tmp63r;
	tmp59i = Im[7]+tmp52i;
	tmp67r = -(Re[7]-tmp60r);
	tmp60i = C7*tmp49i-C8*tmp57r;
	tmp68r = C7*tmp57r+C8*tmp49i;
	tmp61i = C9*tmp50i-C10*tmp58r;
	tmp69r = C9*tmp58r+C10*tmp50i;
	tmp62i = C9*tmp57i-C10*tmp65r;
	tmp70r = C9*tmp65r+C10*tmp57i;
	tmp63i = -C11*tmp58i-C12*tmp66r;
	tmp71r = -C11*tmp66r+C12*tmp58i;
	tmp64i = tmp51i+tmp59i;
	tmp72r = tmp59r+tmp67r;
	tmp65i = tmp51i-tmp59i;
	tmp73r = tmp59r-tmp67r;
	tmp66i = -C1*tmp64i;
	tmp74r = -C1*tmp72r;
	tmp67i = -C6*tmp73r;
	tmp75r = C6*tmp65i;
	tmp68i = tmp66i+tmp39i;
	tmp69i = tmp68i+tmp67i;
	tmp76r = tmp74r+tmp75r;
	tmp70i = tmp68i-tmp67i;
	tmp77r = tmp74r-tmp75r;
	tmp71i = tmp39i+tmp64i;
	tmp72i = tmp60i+tmp62i;
	tmp78r = tmp68r+tmp70r;
	tmp73i = tmp60i-tmp62i;
	tmp79r = tmp68r-tmp70r;
	tmp74i = -C1*tmp72i;
	tmp80r = -C1*tmp78r;
	tmp75i = -C6*tmp79r;
	tmp81r = C6*tmp73i;
	tmp76i = tmp74i+tmp42i;
	tmp77i = tmp76i+tmp75i;
	tmp82r = tmp80r+tmp81r;
	tmp78i = tmp76i-tmp75i;
	tmp83r = tmp80r-tmp81r;
	tmp79i = tmp42i+tmp72i;
	tmp80i = tmp61i+tmp63i;
	tmp84r = tmp69r+tmp71r;
	tmp81i = tmp61i-tmp63i;
	tmp85r = tmp69r-tmp71r;
	tmp82i = -C1*tmp80i;
	tmp86r = -C1*tmp84r;
	tmp83i = -C6*tmp85r;
	tmp87r = C6*tmp81i;
	tmp84i = tmp82i+tmp43i;
	tmp85i = tmp84i+tmp83i;
	tmp88r = tmp86r+tmp87r;
	tmp86i = tmp84i-tmp83i;
	tmp89r = tmp86r-tmp87r;
	tmp87i = tmp43i+tmp80i;
	tmp88i = tmp21i+tmp71i;
	tmp90r = tmp34r+tmp72r;
	tmp91i = tmp37i-tmp85i;
	tmp93r = tmp48r-tmp88r;
	tmp92i = tmp27i+tmp79i;
	tmp94r = tmp42r+tmp78r;
	tmp95i = tmp26i-tmp70i;
	tmp97r = tmp33r-tmp77r;
	tmp96i = tmp33i+tmp87i;
	tmp98r = tmp50r+tmp84r;
	tmp99i = tmp32i-tmp78i;
	tmp101r = tmp41r-tmp83r;
	tmp100i = tmp25i+tmp69i;
	tmp102r = tmp32r+tmp76r;
	tmp103i = tmp38i-tmp86i;
	tmp105r = tmp49r-tmp89r;
	tmp104i = tmp31i+tmp77i;
	tmp106r = tmp40r+tmp82r;
	Re[0] = tmp90r;
	Im[0] = tmp88i;
	Re[1] = tmp93r;
	Im[1] = tmp91i;
	Re[2] = tmp94r;
	Im[2] = tmp92i;
	Re[3] = tmp97r;
	Im[3] = tmp95i;
	Re[4] = tmp98r;
	Im[4] = tmp96i;
	Re[5] = tmp101r;
	Im[5] = tmp99i;
	Re[6] = tmp102r;
	Im[6] = tmp100i;
	Re[7] = tmp105r;
	Im[7] = tmp103i;
	Re[8] = tmp106r;
	Im[8] = tmp104i;
}

/*
*	Number of additions = 175
*	Number of multiplications = 80
*	Number of sign changes = 45
*	Number of assigns = 241
*	Total number of operations = 541
*/
void	MIFFTR18(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, /*tmp105r,*/
		/*tmp106r,*/ tmp107r, tmp108r, /*tmp109r,*/ /*tmp110r,*/ tmp111r, tmp112r,
		/*tmp113r,*/ /*tmp114r,*/ tmp115r, tmp116r, /*tmp117r,*/ /*tmp118r,*/ tmp119r,
		tmp120r /*tmp121r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, /*tmp35i,*/
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, tmp86i, /*tmp87i,*/ tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i, /*tmp104i,*/ /*tmp105i,*/
		/*tmp106i,*/ tmp107i, tmp108i, /*tmp109i,*/ /*tmp110i,*/ tmp111i, tmp112i,
		/*tmp113i,*/ /*tmp114i,*/ tmp115i, tmp116i, /*tmp117i,*/ /*tmp118i,*/ tmp119i,
		tmp120i /*tmp121i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C9 =     0.17364817766693;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.34202014332567;	/* FFT_precisionCONST	*/
	const FFT_precision	C1 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.64278760968654;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =     0.76604444311898;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.93969262078591;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.98480775301221;	/* FFT_precisionCONST	*/

	tmp0i = -(Im[3]-Im[6]);
	tmp0r = Re[3]+Re[6];
	tmp1i = -(Im[3]+Im[6]);
	tmp1r = Re[3]-Re[6];
	tmp2i = -C1*tmp0i;
	tmp2r = -C1*tmp0r;
	tmp3i = C6*tmp1r;
	tmp3r = -C6*tmp1i;
	tmp4i = tmp2i+Im[0];
	tmp4r = tmp2r+Re[0];
	tmp5i = tmp4i+tmp3i;
	tmp5r = tmp4r+tmp3r;
	tmp6i = tmp4i-tmp3i;
	tmp6r = tmp4r-tmp3r;
	tmp7i = Im[0]+tmp0i;
	tmp7r = Re[0]+tmp0r;
	tmp8i = Im[4]-Im[7];
	tmp8r = Re[4]+Re[7];
	tmp9i = Im[4]+Im[7];
	tmp9r = Re[4]-Re[7];
	tmp10i = -C1*tmp8i;
	tmp10r = -C1*tmp8r;
	tmp11i = C6*tmp9r;
	tmp11r = -C6*tmp9i;
	tmp12i = tmp10i-Im[1];
	tmp12r = tmp10r+Re[1];
	tmp13i = tmp12i+tmp11i;
	tmp13r = tmp12r+tmp11r;
	tmp14i = tmp12i-tmp11i;
	tmp14r = tmp12r-tmp11r;
	tmp15i = -(Im[1]-tmp8i);
	tmp15r = Re[1]+tmp8r;
	tmp16i = -(Im[5]-Im[8]);
	tmp16r = Re[5]+Re[8];
	tmp17i = -(Im[5]+Im[8]);
	tmp17r = Re[5]-Re[8];
	tmp18i = -C1*tmp16i;
	tmp18r = -C1*tmp16r;
	tmp19i = C6*tmp17r;
	tmp19r = -C6*tmp17i;
	tmp20i = tmp18i+Im[2];
	tmp20r = tmp18r+Re[2];
	tmp21i = tmp20i+tmp19i;
	tmp21r = tmp20r+tmp19r;
	tmp22i = tmp20i-tmp19i;
	tmp22r = tmp20r-tmp19r;
	tmp23i = Im[2]+tmp16i;
	tmp23r = Re[2]+tmp16r;
	tmp24i = C7*tmp13i+C8*tmp13r;
	tmp24r = C7*tmp13r-C8*tmp13i;
	tmp25i = C9*tmp14i+C10*tmp14r;
	tmp25r = C9*tmp14r-C10*tmp14i;
	tmp26i = C9*tmp21i+C10*tmp21r;
	tmp26r = C9*tmp21r-C10*tmp21i;
	tmp27i = -C11*tmp22i+C12*tmp22r;
	tmp27r = -C11*tmp22r-C12*tmp22i;
	tmp28i = tmp15i+tmp23i;
	tmp28r = tmp15r+tmp23r;
	tmp29i = tmp15i-tmp23i;
	tmp29r = tmp15r-tmp23r;
	tmp30i = -C1*tmp28i;
	tmp30r = -C1*tmp28r;
	tmp31i = C6*tmp29r;
	tmp31r = -C6*tmp29i;
	tmp32i = tmp30i+tmp7i;
	tmp32r = tmp30r+tmp7r;
	tmp33i = tmp32i+tmp31i;
	tmp33r = tmp32r+tmp31r;
	tmp34i = tmp32i-tmp31i;
	tmp34r = tmp32r-tmp31r;
	tmp35r = tmp7r+tmp28r;
	tmp36i = tmp24i+tmp26i;
	tmp36r = tmp24r+tmp26r;
	tmp37i = tmp24i-tmp26i;
	tmp37r = tmp24r-tmp26r;
	tmp38i = -C1*tmp36i;
	tmp38r = -C1*tmp36r;
	tmp39i = C6*tmp37r;
	tmp39r = -C6*tmp37i;
	tmp40i = tmp38i+tmp5i;
	tmp40r = tmp38r+tmp5r;
	tmp41i = tmp40i+tmp39i;
	tmp41r = tmp40r+tmp39r;
	tmp42i = tmp40i-tmp39i;
	tmp42r = tmp40r-tmp39r;
	tmp43i = tmp5i+tmp36i;
	tmp43r = tmp5r+tmp36r;
	tmp44i = tmp25i+tmp27i;
	tmp44r = tmp25r+tmp27r;
	tmp45i = tmp25i-tmp27i;
	tmp45r = tmp25r-tmp27r;
	tmp46i = -C1*tmp44i;
	tmp46r = -C1*tmp44r;
	tmp47i = C6*tmp45r;
	tmp47r = -C6*tmp45i;
	tmp48i = tmp46i+tmp6i;
	tmp48r = tmp46r+tmp6r;
	tmp49i = tmp48i+tmp47i;
	tmp49r = tmp48r+tmp47r;
	tmp50i = tmp48i-tmp47i;
	tmp50r = tmp48r-tmp47r;
	tmp51i = tmp6i+tmp44i;
	tmp51r = tmp6r+tmp44r;
	tmp52i = Im[3]-Im[6];
	tmp52r = Re[3]+Re[6];
	tmp53i = Im[3]+Im[6];
	tmp53r = Re[3]-Re[6];
	tmp54i = -C1*tmp52i;
	tmp54r = -C1*tmp52r;
	tmp55i = C6*tmp53r;
	tmp55r = -C6*tmp53i;
	tmp56i = tmp54i-Im[0];
	tmp56r = tmp54r+Re[0];
	tmp57i = tmp56i+tmp55i;
	tmp57r = tmp56r+tmp55r;
	tmp58i = tmp56i-tmp55i;
	tmp58r = tmp56r-tmp55r;
	tmp59i = -(Im[0]-tmp52i);
	tmp59r = Re[0]+tmp52r;
	tmp60i = -(Im[4]-Im[7]);
	tmp60r = Re[4]+Re[7];
	tmp61i = -(Im[4]+Im[7]);
	tmp61r = Re[4]-Re[7];
	tmp62i = -C1*tmp60i;
	tmp62r = -C1*tmp60r;
	tmp63i = C6*tmp61r;
	tmp63r = -C6*tmp61i;
	tmp64i = tmp62i+Im[1];
	tmp64r = tmp62r+Re[1];
	tmp65i = tmp64i+tmp63i;
	tmp65r = tmp64r+tmp63r;
	tmp66i = tmp64i-tmp63i;
	tmp66r = tmp64r-tmp63r;
	tmp67i = Im[1]+tmp60i;
	tmp67r = Re[1]+tmp60r;
	tmp68i = Im[5]-Im[8];
	tmp68r = Re[5]+Re[8];
	tmp69i = Im[5]+Im[8];
	tmp69r = Re[5]-Re[8];
	tmp70i = -C1*tmp68i;
	tmp70r = -C1*tmp68r;
	tmp71i = C6*tmp69r;
	tmp71r = -C6*tmp69i;
	tmp72i = tmp70i-Im[2];
	tmp72r = tmp70r+Re[2];
	tmp73i = tmp72i+tmp71i;
	tmp73r = tmp72r+tmp71r;
	tmp74i = tmp72i-tmp71i;
	tmp74r = tmp72r-tmp71r;
	tmp75i = -(Im[2]-tmp68i);
	tmp75r = Re[2]+tmp68r;
	tmp76i = C7*tmp65i+C8*tmp65r;
	tmp76r = C7*tmp65r-C8*tmp65i;
	tmp77i = C9*tmp66i+C10*tmp66r;
	tmp77r = C9*tmp66r-C10*tmp66i;
	tmp78i = C9*tmp73i+C10*tmp73r;
	tmp78r = C9*tmp73r-C10*tmp73i;
	tmp79i = -C11*tmp74i+C12*tmp74r;
	tmp79r = -C11*tmp74r-C12*tmp74i;
	tmp80i = tmp67i+tmp75i;
	tmp80r = tmp67r+tmp75r;
	tmp81i = tmp67i-tmp75i;
	tmp81r = tmp67r-tmp75r;
	tmp82i = -C1*tmp80i;
	tmp82r = -C1*tmp80r;
	tmp83i = C6*tmp81r;
	tmp83r = -C6*tmp81i;
	tmp84i = tmp82i+tmp59i;
	tmp84r = tmp82r+tmp59r;
	tmp85i = tmp84i+tmp83i;
	tmp85r = tmp84r+tmp83r;
	tmp86i = tmp84i-tmp83i;
	tmp86r = tmp84r-tmp83r;
	tmp87r = tmp59r+tmp80r;
	tmp88i = tmp76i+tmp78i;
	tmp88r = tmp76r+tmp78r;
	tmp89i = tmp76i-tmp78i;
	tmp89r = tmp76r-tmp78r;
	tmp90i = -C1*tmp88i;
	tmp90r = -C1*tmp88r;
	tmp91i = C6*tmp89r;
	tmp91r = -C6*tmp89i;
	tmp92i = tmp90i+tmp57i;
	tmp92r = tmp90r+tmp57r;
	tmp93i = tmp92i+tmp91i;
	tmp93r = tmp92r+tmp91r;
	tmp94i = tmp92i-tmp91i;
	tmp94r = tmp92r-tmp91r;
	tmp95i = tmp57i+tmp88i;
	tmp95r = tmp57r+tmp88r;
	tmp96i = tmp77i+tmp79i;
	tmp96r = tmp77r+tmp79r;
	tmp97i = tmp77i-tmp79i;
	tmp97r = tmp77r-tmp79r;
	tmp98i = -C1*tmp96i;
	tmp98r = -C1*tmp96r;
	tmp99i = C6*tmp97r;
	tmp99r = -C6*tmp97i;
	tmp100i = tmp98i+tmp58i;
	tmp100r = tmp98r+tmp58r;
	tmp101i = tmp100i+tmp99i;
	tmp101r = tmp100r+tmp99r;
	tmp102i = tmp100i-tmp99i;
	tmp102r = tmp100r-tmp99r;
	tmp103i = tmp58i+tmp96i;
	tmp103r = tmp58r+tmp96r;
	tmp104r = tmp35r+tmp87r;
	tmp107i = tmp49i-tmp101i;
	tmp107r = tmp49r-tmp101r;
	tmp108i = tmp43i+tmp95i;
	tmp108r = tmp43r+tmp95r;
	tmp111i = tmp34i-tmp86i;
	tmp111r = tmp34r-tmp86r;
	tmp112i = tmp51i+tmp103i;
	tmp112r = tmp51r+tmp103r;
	tmp115i = tmp42i-tmp94i;
	tmp115r = tmp42r-tmp94r;
	tmp116i = tmp33i+tmp85i;
	tmp116r = tmp33r+tmp85r;
	tmp119i = tmp50i-tmp102i;
	tmp119r = tmp50r-tmp102r;
	tmp120i = tmp41i+tmp93i;
	tmp120r = tmp41r+tmp93r;
	Re[0] = tmp104r;
	Im[0] = C0;
	Re[1] = tmp107r;
	Im[1] = tmp107i;
	Re[2] = tmp108r;
	Im[2] = tmp108i;
	Re[3] = tmp111r;
	Im[3] = tmp111i;
	Re[4] = tmp112r;
	Im[4] = tmp112i;
	Re[5] = tmp115r;
	Im[5] = tmp115i;
	Re[6] = tmp116r;
	Im[6] = tmp116i;
	Re[7] = tmp119r;
	Im[7] = tmp119i;
	Re[8] = tmp120r;
	Im[8] = tmp120i;
}

/*
*	Number of additions = 126
*	Number of multiplications = 44
*	Number of sign changes = 18
*	Number of assigns = 166
*	Total number of operations = 354
*/
void	MFFTR20(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		/*tmp71r,*/ tmp72r, /*tmp73r,*/ tmp74r, tmp75r, tmp76r, tmp77r,
		/*tmp78r,*/ tmp79r, tmp80r, /*tmp81r,*/ tmp82r, tmp83r, tmp84r,
		tmp85r, /*tmp86r,*/ tmp87r, /*tmp88r,*/ tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, /*tmp95r,*/ /*tmp96r,*/ tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, /*tmp103r,*/ tmp104r /*tmp105r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i;

	const FFT_precision	C1 =                 0.25;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.55901699437495;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.58778525229247;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.95105651629515;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =                    2;	/* INTEGERCONST	*/

	tmp0r = C12*Re[4];
	tmp0i = -C12*Im[4];
	tmp1r = C12*Re[8];
	tmp1i = -C12*Im[8];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = -C1*tmp2r;
	tmp5r = C15*tmp3r;
	tmp6r = tmp4r+Re[0];
	tmp7r = tmp6r+tmp5r;
	tmp8r = tmp6r-tmp5r;
	tmp9r = C8*tmp0i+C10*tmp1i;
	tmp10r = C10*tmp0i-C8*tmp1i;
	tmp11r = tmp7r+tmp9r;
	tmp12r = tmp7r-tmp9r;
	tmp13r = tmp8r+tmp10r;
	tmp14r = tmp8r-tmp10r;
	tmp15r = Re[0]+tmp2r;
	tmp2i = Im[1]+Im[9];
	tmp16r = Re[1]+Re[9];
	tmp3i = Im[1]-Im[9];
	tmp17r = Re[1]-Re[9];
	tmp4i = Im[3]+Im[7];
	tmp18r = -(Re[3]+Re[7]);
	tmp5i = Im[3]-Im[7];
	tmp19r = -(Re[3]-Re[7]);
	tmp6i = tmp2i+tmp4i;
	tmp20r = tmp16r+tmp18r;
	tmp7i = tmp2i-tmp4i;
	tmp21r = tmp16r-tmp18r;
	tmp8i = -C1*tmp6i;
	tmp22r = -C1*tmp20r;
	tmp9i = C15*tmp7i;
	tmp23r = C15*tmp21r;
	tmp10i = tmp8i+Im[5];
	tmp24r = tmp22r+Re[5];
	tmp11i = tmp10i+tmp9i;
	tmp25r = tmp24r+tmp23r;
	tmp12i = tmp10i-tmp9i;
	tmp26r = tmp24r-tmp23r;
	tmp13i = -C8*tmp17r-C10*tmp19r;
	tmp27r = C8*tmp3i+C10*tmp5i;
	tmp14i = -C10*tmp17r+C8*tmp19r;
	tmp28r = C10*tmp3i-C8*tmp5i;
	tmp15i = tmp11i+tmp13i;
	tmp29r = tmp25r+tmp27r;
	tmp16i = tmp11i-tmp13i;
	tmp30r = tmp25r-tmp27r;
	tmp17i = tmp12i+tmp14i;
	tmp31r = tmp26r+tmp28r;
	tmp18i = tmp12i-tmp14i;
	tmp32r = tmp26r-tmp28r;
	tmp19i = Im[5]+tmp6i;
	tmp33r = Re[5]+tmp20r;
	tmp34r = C12*Re[6];
	tmp20i = C12*Im[6];
	tmp35r = C12*Re[2];
	tmp21i = C12*Im[2];
	tmp36r = tmp34r+tmp35r;
	tmp37r = tmp34r-tmp35r;
	tmp38r = -C1*tmp36r;
	tmp39r = C15*tmp37r;
	tmp40r = tmp38r+tmp39r;
	tmp41r = tmp38r-tmp39r;
	tmp42r = C8*tmp20i+C10*tmp21i;
	tmp43r = C10*tmp20i-C8*tmp21i;
	tmp44r = tmp40r+tmp42r;
	tmp45r = tmp40r-tmp42r;
	tmp46r = tmp41r+tmp43r;
	tmp47r = tmp41r-tmp43r;
	tmp22i = Im[9]+Im[1];
	tmp48r = -(Re[9]+Re[1]);
	tmp23i = Im[9]-Im[1];
	tmp49r = -(Re[9]-Re[1]);
	tmp24i = Im[7]+Im[3];
	tmp50r = Re[7]+Re[3];
	tmp25i = Im[7]-Im[3];
	tmp51r = Re[7]-Re[3];
	tmp26i = tmp22i+tmp24i;
	tmp52r = tmp48r+tmp50r;
	tmp27i = tmp22i-tmp24i;
	tmp53r = tmp48r-tmp50r;
	tmp28i = -C1*tmp26i;
	tmp54r = -C1*tmp52r;
	tmp29i = C15*tmp27i;
	tmp55r = C15*tmp53r;
	tmp30i = tmp28i+Im[5];
	tmp56r = tmp54r-Re[5];
	tmp31i = tmp30i+tmp29i;
	tmp57r = tmp56r+tmp55r;
	tmp32i = tmp30i-tmp29i;
	tmp58r = tmp56r-tmp55r;
	tmp33i = -C8*tmp49r-C10*tmp51r;
	tmp59r = C8*tmp23i+C10*tmp25i;
	tmp34i = -C10*tmp49r+C8*tmp51r;
	tmp60r = C10*tmp23i-C8*tmp25i;
	tmp35i = tmp31i+tmp33i;
	tmp61r = tmp57r+tmp59r;
	tmp36i = tmp31i-tmp33i;
	tmp62r = tmp57r-tmp59r;
	tmp37i = tmp32i+tmp34i;
	tmp63r = tmp58r+tmp60r;
	tmp38i = tmp32i-tmp34i;
	tmp64r = tmp58r-tmp60r;
	tmp39i = Im[5]+tmp26i;
	tmp65r = -(Re[5]-tmp52r);
	tmp66r = tmp15r+tmp36r;
	tmp67r = tmp15r-tmp36r;
	tmp40i = tmp19i+tmp39i;
	tmp68r = tmp33r+tmp65r;
	tmp41i = tmp19i-tmp39i;
	tmp69r = tmp33r-tmp65r;
	tmp70r = tmp66r+tmp68r;
	tmp72r = tmp67r+tmp41i;
	tmp74r = tmp12r+tmp45r;
	tmp75r = tmp12r-tmp45r;
	tmp42i = tmp16i+tmp36i;
	tmp76r = tmp30r+tmp62r;
	tmp43i = tmp16i-tmp36i;
	tmp77r = tmp30r-tmp62r;
	tmp79r = tmp74r-tmp76r;
	tmp80r = tmp75r+tmp43i;
	tmp82r = tmp14r+tmp47r;
	tmp83r = tmp14r-tmp47r;
	tmp44i = tmp18i+tmp38i;
	tmp84r = tmp32r+tmp64r;
	tmp45i = tmp18i-tmp38i;
	tmp85r = tmp32r-tmp64r;
	tmp87r = tmp82r-tmp84r;
	tmp89r = tmp83r-tmp45i;
	tmp90r = tmp13r+tmp46r;
	tmp91r = tmp13r-tmp46r;
	tmp46i = tmp17i+tmp37i;
	tmp92r = tmp31r+tmp63r;
	tmp47i = tmp17i-tmp37i;
	tmp93r = tmp31r-tmp63r;
	tmp94r = tmp90r+tmp92r;
	tmp97r = tmp91r-tmp47i;
	tmp98r = tmp11r+tmp44r;
	tmp99r = tmp11r-tmp44r;
	tmp48i = tmp15i+tmp35i;
	tmp100r = tmp29r+tmp61r;
	tmp49i = tmp15i-tmp35i;
	tmp101r = tmp29r-tmp61r;
	tmp102r = tmp98r+tmp100r;
	tmp104r = tmp99r+tmp49i;
	Re[0] = tmp70r;
	Im[0] = tmp40i;
	Re[1] = tmp80r;
	Im[1] = -tmp77r;
	Re[2] = tmp87r;
	Im[2] = -tmp44i;
	Re[3] = tmp97r;
	Im[3] = tmp93r;
	Re[4] = tmp102r;
	Im[4] = tmp48i;
	Re[5] = tmp72r;
	Im[5] = -tmp69r;
	Re[6] = tmp79r;
	Im[6] = -tmp42i;
	Re[7] = tmp89r;
	Im[7] = tmp85r;
	Re[8] = tmp94r;
	Im[8] = tmp46i;
	Re[9] = tmp104r;
	Im[9] = -tmp101r;
}

/*
*	Number of additions = 185
*	Number of multiplications = 48
*	Number of sign changes = 18
*	Number of assigns = 221
*	Total number of operations = 472
*/
void	MIFFTR20(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, /*tmp77r,*/
		tmp78r, /*tmp79r,*/ tmp80r, tmp81r, tmp82r, tmp83r, /*tmp84r,*/
		tmp85r, tmp86r, /*tmp87r,*/ tmp88r, tmp89r, tmp90r, tmp91r,
		/*tmp92r,*/ tmp93r, /*tmp94r,*/ tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, /*tmp101r,*/ /*tmp102r,*/ tmp103r, tmp104r, tmp105r,
		tmp106r, tmp107r, tmp108r, /*tmp109r,*/ tmp110r /*tmp111r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, /*tmp72i,*/ tmp73i, /*tmp74i,*/ tmp75i, /*tmp76i,*/ /*tmp77i,*/
		tmp78i, /*tmp79i,*/ tmp80i, tmp81i, tmp82i, tmp83i, /*tmp84i,*/
		tmp85i, tmp86i, /*tmp87i,*/ tmp88i, tmp89i, tmp90i, tmp91i,
		/*tmp92i,*/ tmp93i, /*tmp94i,*/ tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, tmp100i, /*tmp101i,*/ /*tmp102i,*/ tmp103i, tmp104i, tmp105i,
		tmp106i, tmp107i, tmp108i, /*tmp109i,*/ tmp110i /*tmp111i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C1 =                 0.25;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.55901699437495;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.58778525229247;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.95105651629515;	/* FFT_precisionCONST	*/

	tmp0i = -(Im[6]-Im[4]);
	tmp0r = Re[6]+Re[4];
	tmp1i = -(Im[6]+Im[4]);
	tmp1r = Re[6]-Re[4];
	tmp2i = -(Im[2]-Im[8]);
	tmp2r = Re[2]+Re[8];
	tmp3i = -(Im[2]+Im[8]);
	tmp3r = Re[2]-Re[8];
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = -C1*tmp4i;
	tmp6r = -C1*tmp4r;
	tmp7i = C15*tmp5i;
	tmp7r = C15*tmp5r;
	tmp8i = tmp6i+Im[0];
	tmp8r = tmp6r+Re[0];
	tmp9i = tmp8i+tmp7i;
	tmp9r = tmp8r+tmp7r;
	tmp10i = tmp8i-tmp7i;
	tmp10r = tmp8r-tmp7r;
	tmp11i = C8*tmp1r+C10*tmp3r;
	tmp11r = -C8*tmp1i-C10*tmp3i;
	tmp12i = C10*tmp1r-C8*tmp3r;
	tmp12r = -C10*tmp1i+C8*tmp3i;
	tmp13i = tmp9i+tmp11i;
	tmp13r = tmp9r+tmp11r;
	tmp14i = tmp9i-tmp11i;
	tmp14r = tmp9r-tmp11r;
	tmp15i = tmp10i+tmp12i;
	tmp15r = tmp10r+tmp12r;
	tmp16i = tmp10i-tmp12i;
	tmp16r = tmp10r-tmp12r;
	tmp17i = Im[0]+tmp4i;
	tmp17r = Re[0]+tmp4r;
	tmp18i = Im[1]+Im[9];
	tmp18r = Re[1]+Re[9];
	tmp19i = Im[1]-Im[9];
	tmp19r = Re[1]-Re[9];
	tmp20i = -(Im[7]+Im[3]);
	tmp20r = Re[7]+Re[3];
	tmp21i = -(Im[7]-Im[3]);
	tmp21r = Re[7]-Re[3];
	tmp22i = tmp18i+tmp20i;
	tmp22r = tmp18r+tmp20r;
	tmp23i = tmp18i-tmp20i;
	tmp23r = tmp18r-tmp20r;
	tmp24i = -C1*tmp22i;
	tmp24r = -C1*tmp22r;
	tmp25i = C15*tmp23i;
	tmp25r = C15*tmp23r;
	tmp26i = tmp24i+Im[5];
	tmp26r = tmp24r+Re[5];
	tmp27i = tmp26i+tmp25i;
	tmp27r = tmp26r+tmp25r;
	tmp28i = tmp26i-tmp25i;
	tmp28r = tmp26r-tmp25r;
	tmp29i = C8*tmp19r+C10*tmp21r;
	tmp29r = -C8*tmp19i-C10*tmp21i;
	tmp30i = C10*tmp19r-C8*tmp21r;
	tmp30r = -C10*tmp19i+C8*tmp21i;
	tmp31i = tmp27i+tmp29i;
	tmp31r = tmp27r+tmp29r;
	tmp32i = tmp27i-tmp29i;
	tmp32r = tmp27r-tmp29r;
	tmp33i = tmp28i+tmp30i;
	tmp33r = tmp28r+tmp30r;
	tmp34i = tmp28i-tmp30i;
	tmp34r = tmp28r-tmp30r;
	tmp35i = Im[5]+tmp22i;
	tmp35r = Re[5]+tmp22r;
	tmp36i = Im[6]-Im[4];
	tmp36r = Re[6]+Re[4];
	tmp37i = Im[6]+Im[4];
	tmp37r = Re[6]-Re[4];
	tmp38i = Im[2]-Im[8];
	tmp38r = Re[2]+Re[8];
	tmp39i = Im[2]+Im[8];
	tmp39r = Re[2]-Re[8];
	tmp40i = tmp36i+tmp38i;
	tmp40r = tmp36r+tmp38r;
	tmp41i = tmp36i-tmp38i;
	tmp41r = tmp36r-tmp38r;
	tmp42i = -C1*tmp40i;
	tmp42r = -C1*tmp40r;
	tmp43i = C15*tmp41i;
	tmp43r = C15*tmp41r;
	tmp44i = tmp42i-Im[0];
	tmp44r = tmp42r+Re[0];
	tmp45i = tmp44i+tmp43i;
	tmp45r = tmp44r+tmp43r;
	tmp46i = tmp44i-tmp43i;
	tmp46r = tmp44r-tmp43r;
	tmp47i = C8*tmp37r+C10*tmp39r;
	tmp47r = -C8*tmp37i-C10*tmp39i;
	tmp48i = C10*tmp37r-C8*tmp39r;
	tmp48r = -C10*tmp37i+C8*tmp39i;
	tmp49i = tmp45i+tmp47i;
	tmp49r = tmp45r+tmp47r;
	tmp50i = tmp45i-tmp47i;
	tmp50r = tmp45r-tmp47r;
	tmp51i = tmp46i+tmp48i;
	tmp51r = tmp46r+tmp48r;
	tmp52i = tmp46i-tmp48i;
	tmp52r = tmp46r-tmp48r;
	tmp53i = -(Im[0]-tmp40i);
	tmp53r = Re[0]+tmp40r;
	tmp54i = -(Im[1]+Im[9]);
	tmp54r = Re[1]+Re[9];
	tmp55i = -(Im[1]-Im[9]);
	tmp55r = Re[1]-Re[9];
	tmp56i = Im[7]+Im[3];
	tmp56r = Re[7]+Re[3];
	tmp57i = Im[7]-Im[3];
	tmp57r = Re[7]-Re[3];
	tmp58i = tmp54i+tmp56i;
	tmp58r = tmp54r+tmp56r;
	tmp59i = tmp54i-tmp56i;
	tmp59r = tmp54r-tmp56r;
	tmp60i = -C1*tmp58i;
	tmp60r = -C1*tmp58r;
	tmp61i = C15*tmp59i;
	tmp61r = C15*tmp59r;
	tmp62i = tmp60i-Im[5];
	tmp62r = tmp60r+Re[5];
	tmp63i = tmp62i+tmp61i;
	tmp63r = tmp62r+tmp61r;
	tmp64i = tmp62i-tmp61i;
	tmp64r = tmp62r-tmp61r;
	tmp65i = C8*tmp55r+C10*tmp57r;
	tmp65r = -C8*tmp55i-C10*tmp57i;
	tmp66i = C10*tmp55r-C8*tmp57r;
	tmp66r = -C10*tmp55i+C8*tmp57i;
	tmp67i = tmp63i+tmp65i;
	tmp67r = tmp63r+tmp65r;
	tmp68i = tmp63i-tmp65i;
	tmp68r = tmp63r-tmp65r;
	tmp69i = tmp64i+tmp66i;
	tmp69r = tmp64r+tmp66r;
	tmp70i = tmp64i-tmp66i;
	tmp70r = tmp64r-tmp66r;
	tmp71i = -(Im[5]-tmp58i);
	tmp71r = Re[5]+tmp58r;
	tmp72r = tmp17r+tmp53r;
	tmp73i = tmp17i-tmp53i;
	tmp73r = tmp17r-tmp53r;
	tmp74r = tmp35r+tmp71r;
	tmp75i = tmp35i-tmp71i;
	tmp75r = tmp35r-tmp71r;
	tmp76r = tmp72r+tmp74r;
	tmp78i = tmp73i+tmp75r;
	tmp78r = tmp73r-tmp75i;
	tmp80i = tmp14i+tmp50i;
	tmp80r = tmp14r+tmp50r;
	tmp81i = tmp14i-tmp50i;
	tmp81r = tmp14r-tmp50r;
	tmp82i = tmp32i+tmp68i;
	tmp82r = tmp32r+tmp68r;
	tmp83i = tmp32i-tmp68i;
	tmp83r = tmp32r-tmp68r;
	tmp85i = tmp80i-tmp82i;
	tmp85r = tmp80r-tmp82r;
	tmp86i = tmp81i+tmp83r;
	tmp86r = tmp81r-tmp83i;
	tmp88i = tmp16i+tmp52i;
	tmp88r = tmp16r+tmp52r;
	tmp89i = tmp16i-tmp52i;
	tmp89r = tmp16r-tmp52r;
	tmp90i = tmp34i+tmp70i;
	tmp90r = tmp34r+tmp70r;
	tmp91i = tmp34i-tmp70i;
	tmp91r = tmp34r-tmp70r;
	tmp93i = tmp88i-tmp90i;
	tmp93r = tmp88r-tmp90r;
	tmp95i = tmp89i-tmp91r;
	tmp95r = tmp89r+tmp91i;
	tmp96i = tmp15i+tmp51i;
	tmp96r = tmp15r+tmp51r;
	tmp97i = tmp15i-tmp51i;
	tmp97r = tmp15r-tmp51r;
	tmp98i = tmp33i+tmp69i;
	tmp98r = tmp33r+tmp69r;
	tmp99i = tmp33i-tmp69i;
	tmp99r = tmp33r-tmp69r;
	tmp100i = tmp96i+tmp98i;
	tmp100r = tmp96r+tmp98r;
	tmp103i = tmp97i-tmp99r;
	tmp103r = tmp97r+tmp99i;
	tmp104i = tmp13i+tmp49i;
	tmp104r = tmp13r+tmp49r;
	tmp105i = tmp13i-tmp49i;
	tmp105r = tmp13r-tmp49r;
	tmp106i = tmp31i+tmp67i;
	tmp106r = tmp31r+tmp67r;
	tmp107i = tmp31i-tmp67i;
	tmp107r = tmp31r-tmp67r;
	tmp108i = tmp104i+tmp106i;
	tmp108r = tmp104r+tmp106r;
	tmp110i = tmp105i+tmp107r;
	tmp110r = tmp105r-tmp107i;
	Re[0] = tmp76r;
	Im[0] = C0;
	Re[1] = tmp86r;
	Im[1] = tmp86i;
	Re[2] = tmp93r;
	Im[2] = tmp93i;
	Re[3] = tmp103r;
	Im[3] = tmp103i;
	Re[4] = tmp108r;
	Im[4] = tmp108i;
	Re[5] = tmp78r;
	Im[5] = tmp78i;
	Re[6] = tmp85r;
	Im[6] = tmp85i;
	Re[7] = tmp95r;
	Im[7] = tmp95i;
	Re[8] = tmp100r;
	Im[8] = tmp100i;
	Re[9] = tmp110r;
	Im[9] = tmp110i;
}

/*
*	Number of additions = 114
*	Number of multiplications = 120
*	Number of sign changes = 10
*	Number of assigns = 91
*	Total number of operations = 335
*/
void	MFFTR22(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i;

	const FFT_precision	C19 =     0.14231483827329;	/* FFT_precisionCONST	*/
	const FFT_precision	C22 =     0.28173255684143;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.41541501300189;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =      0.5406408174556;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.65486073394529;	/* FFT_precisionCONST	*/
	const FFT_precision	C18 =     0.75574957435426;	/* FFT_precisionCONST	*/
	const FFT_precision	C13 =     0.84125353283118;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =     0.90963199535452;	/* FFT_precisionCONST	*/
	const FFT_precision	C21 =      0.9594929736145;	/* FFT_precisionCONST	*/
	const FFT_precision	C20 =     0.98982144188093;	/* FFT_precisionCONST	*/
	const FFT_precision	C24 =                    2;	/* INTEGERCONST	*/

	tmp0r = C24*Re[10];
	tmp0i = -C24*Im[10];
	tmp1r = C24*Re[2];
	tmp1i = C24*Im[2];
	tmp2r = C24*Re[4];
	tmp2i = C24*Im[4];
	tmp3r = C24*Re[8];
	tmp3i = C24*Im[8];
	tmp4r = C24*Re[6];
	tmp4i = -C24*Im[6];
	tmp5r = C13*tmp0r+C15*tmp1r-C17*tmp2r-C19*tmp3r-C21*tmp4r+Re[0];
	tmp6r = tmp0r+tmp1r;
	tmp7r = C15*tmp0r-C17*tmp1r-C19*tmp2r-C21*tmp3r+C13*tmp4r+Re[0];
	tmp8r = tmp6r+tmp2r;
	tmp9r = -C17*tmp0r-C19*tmp1r-C21*tmp2r+C13*tmp3r+C15*tmp4r+Re[0];
	tmp10r = tmp8r+tmp3r;
	tmp11r = -C19*tmp0r-C21*tmp1r+C13*tmp2r+C15*tmp3r-C17*tmp4r+Re[0];
	tmp12r = tmp10r+tmp4r;
	tmp13r = -C21*tmp0r+C13*tmp1r+C15*tmp2r-C17*tmp3r-C19*tmp4r+Re[0];
	tmp14r = C14*tmp0i+C16*tmp1i+C18*tmp2i-C20*tmp3i+C22*tmp4i;
	tmp15r = C16*tmp0i+C18*tmp1i-C20*tmp2i+C22*tmp3i-C14*tmp4i;
	tmp16r = C18*tmp0i-C20*tmp1i+C22*tmp2i-C14*tmp3i-C16*tmp4i;
	tmp17r = -C20*tmp0i+C22*tmp1i-C14*tmp2i-C16*tmp3i-C18*tmp4i;
	tmp18r = C22*tmp0i-C14*tmp1i-C16*tmp2i-C18*tmp3i+C20*tmp4i;
	tmp19r = tmp5r+tmp14r;
	tmp20r = tmp5r-tmp14r;
	tmp21r = tmp7r+tmp15r;
	tmp22r = tmp7r-tmp15r;
	tmp23r = tmp9r+tmp16r;
	tmp24r = tmp9r-tmp16r;
	tmp25r = tmp11r+tmp17r;
	tmp26r = tmp11r-tmp17r;
	tmp27r = tmp13r+tmp18r;
	tmp28r = tmp13r-tmp18r;
	tmp29r = Re[0]+tmp12r;
	tmp5i = C24*Im[1];
	tmp30r = C24*Re[1];
	tmp6i = C24*Im[9];
	tmp31r = -C24*Re[9];
	tmp7i = C24*Im[7];
	tmp32r = -C24*Re[7];
	tmp8i = C24*Im[3];
	tmp33r = -C24*Re[3];
	tmp9i = C24*Im[5];
	tmp34r = C24*Re[5];
	tmp10i = C13*tmp5i+C15*tmp6i-C17*tmp7i-C19*tmp8i-C21*tmp9i;
	tmp11i = tmp5i+tmp6i;
	tmp12i = C15*tmp5i-C17*tmp6i-C19*tmp7i-C21*tmp8i+C13*tmp9i;
	tmp13i = tmp11i+tmp7i;
	tmp14i = -C17*tmp5i-C19*tmp6i-C21*tmp7i+C13*tmp8i+C15*tmp9i;
	tmp15i = tmp13i+tmp8i;
	tmp16i = -C19*tmp5i-C21*tmp6i+C13*tmp7i+C15*tmp8i-C17*tmp9i;
	tmp17i = tmp15i+tmp9i;
	tmp18i = -C21*tmp5i+C13*tmp6i+C15*tmp7i-C17*tmp8i-C19*tmp9i;
	tmp19i = -C14*tmp30r-C16*tmp31r-C18*tmp32r+C20*tmp33r-C22*tmp34r;
	tmp20i = -C16*tmp30r-C18*tmp31r+C20*tmp32r-C22*tmp33r+C14*tmp34r;
	tmp21i = -C18*tmp30r+C20*tmp31r-C22*tmp32r+C14*tmp33r+C16*tmp34r;
	tmp22i = C20*tmp30r-C22*tmp31r+C14*tmp32r+C16*tmp33r+C18*tmp34r;
	tmp23i = -C22*tmp30r+C14*tmp31r+C16*tmp32r+C18*tmp33r-C20*tmp34r;
	tmp24i = tmp10i+tmp19i;
	tmp25i = tmp10i-tmp19i;
	tmp26i = tmp12i+tmp20i;
	tmp27i = tmp12i-tmp20i;
	tmp28i = tmp14i+tmp21i;
	tmp29i = tmp14i-tmp21i;
	tmp30i = tmp16i+tmp22i;
	tmp31i = tmp16i-tmp22i;
	tmp32i = tmp18i+tmp23i;
	tmp33i = tmp18i-tmp23i;
	Re[0] = tmp29r;
	Im[0] = tmp17i;
	Re[1] = tmp28r;
	Im[1] = -tmp33i;
	Re[2] = tmp19r;
	Im[2] = tmp24i;
	Re[3] = tmp24r;
	Im[3] = -tmp29i;
	Re[4] = tmp21r;
	Im[4] = tmp26i;
	Re[5] = tmp25r;
	Im[5] = -tmp30i;
	Re[6] = tmp26r;
	Im[6] = tmp31i;
	Re[7] = tmp22r;
	Im[7] = -tmp27i;
	Re[8] = tmp23r;
	Im[8] = tmp28i;
	Re[9] = tmp20r;
	Im[9] = -tmp25i;
	Re[10] = tmp27r;
	Im[10] = tmp32i;
}

/*
*	Number of additions = 291
*	Number of multiplications = 200
*	Number of sign changes = 10
*	Number of assigns = 173
*	Total number of operations = 674
*/
void	MIFFTR22(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		/*tmp71r,*/ /*tmp72r,*/ tmp73r, tmp74r, /*tmp75r,*/ /*tmp76r,*/ tmp77r,
		tmp78r, /*tmp79r,*/ /*tmp80r,*/ tmp81r, tmp82r, /*tmp83r,*/ /*tmp84r,*/
		tmp85r, tmp86r, /*tmp87r,*/ /*tmp88r,*/ tmp89r, tmp90r /*tmp91r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, /*tmp11i,*/ tmp12i, /*tmp13i,*/ tmp14i,
		/*tmp15i,*/ tmp16i, /*tmp17i,*/ tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, /*tmp34i,*/ tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, /*tmp46i,*/ tmp47i, /*tmp48i,*/ tmp49i,
		/*tmp50i,*/ tmp51i, /*tmp52i,*/ tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, /*tmp69i,*/ /*tmp70i,*/
		/*tmp71i,*/ /*tmp72i,*/ tmp73i, tmp74i, /*tmp75i,*/ /*tmp76i,*/ tmp77i,
		tmp78i, /*tmp79i,*/ /*tmp80i,*/ tmp81i, tmp82i, /*tmp83i,*/ /*tmp84i,*/
		tmp85i, tmp86i, /*tmp87i,*/ /*tmp88i,*/ tmp89i, tmp90i /*tmp91i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C19 =     0.14231483827329;	/* FFT_precisionCONST	*/
	const FFT_precision	C22 =     0.28173255684143;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.41541501300189;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =      0.5406408174556;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.65486073394529;	/* FFT_precisionCONST	*/
	const FFT_precision	C18 =     0.75574957435426;	/* FFT_precisionCONST	*/
	const FFT_precision	C13 =     0.84125353283118;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =     0.90963199535452;	/* FFT_precisionCONST	*/
	const FFT_precision	C21 =      0.9594929736145;	/* FFT_precisionCONST	*/
	const FFT_precision	C20 =     0.98982144188093;	/* FFT_precisionCONST	*/

	tmp0i = -(Im[1]-Im[10]);
	tmp0r = Re[1]+Re[10];
	tmp1i = -(Im[1]+Im[10]);
	tmp1r = Re[1]-Re[10];
	tmp2i = Im[2]-Im[9];
	tmp2r = Re[2]+Re[9];
	tmp3i = Im[2]+Im[9];
	tmp3r = Re[2]-Re[9];
	tmp4i = Im[4]-Im[7];
	tmp4r = Re[4]+Re[7];
	tmp5i = Im[4]+Im[7];
	tmp5r = Re[4]-Re[7];
	tmp6i = Im[8]-Im[3];
	tmp6r = Re[8]+Re[3];
	tmp7i = Im[8]+Im[3];
	tmp7r = Re[8]-Re[3];
	tmp8i = -(Im[5]-Im[6]);
	tmp8r = Re[5]+Re[6];
	tmp9i = -(Im[5]+Im[6]);
	tmp9r = Re[5]-Re[6];
	tmp10i = C13*tmp0i+C15*tmp2i-C17*tmp4i-C19*tmp6i-C21*tmp8i+Im[0];
	tmp10r = C13*tmp0r+C15*tmp2r-C17*tmp4r-C19*tmp6r-C21*tmp8r+Re[0];
	tmp11r = tmp0r+tmp2r;
	tmp12i = C15*tmp0i-C17*tmp2i-C19*tmp4i-C21*tmp6i+C13*tmp8i+Im[0];
	tmp12r = C15*tmp0r-C17*tmp2r-C19*tmp4r-C21*tmp6r+C13*tmp8r+Re[0];
	tmp13r = tmp11r+tmp4r;
	tmp14i = -C17*tmp0i-C19*tmp2i-C21*tmp4i+C13*tmp6i+C15*tmp8i+Im[0];
	tmp14r = -C17*tmp0r-C19*tmp2r-C21*tmp4r+C13*tmp6r+C15*tmp8r+Re[0];
	tmp15r = tmp13r+tmp6r;
	tmp16i = -C19*tmp0i-C21*tmp2i+C13*tmp4i+C15*tmp6i-C17*tmp8i+Im[0];
	tmp16r = -C19*tmp0r-C21*tmp2r+C13*tmp4r+C15*tmp6r-C17*tmp8r+Re[0];
	tmp17r = tmp15r+tmp8r;
	tmp18i = -C21*tmp0i+C13*tmp2i+C15*tmp4i-C17*tmp6i-C19*tmp8i+Im[0];
	tmp18r = -C21*tmp0r+C13*tmp2r+C15*tmp4r-C17*tmp6r-C19*tmp8r+Re[0];
	tmp19i = C14*tmp1r+C16*tmp3r+C18*tmp5r-C20*tmp7r+C22*tmp9r;
	tmp19r = -C14*tmp1i-C16*tmp3i-C18*tmp5i+C20*tmp7i-C22*tmp9i;
	tmp20i = C16*tmp1r+C18*tmp3r-C20*tmp5r+C22*tmp7r-C14*tmp9r;
	tmp20r = -C16*tmp1i-C18*tmp3i+C20*tmp5i-C22*tmp7i+C14*tmp9i;
	tmp21i = C18*tmp1r-C20*tmp3r+C22*tmp5r-C14*tmp7r-C16*tmp9r;
	tmp21r = -C18*tmp1i+C20*tmp3i-C22*tmp5i+C14*tmp7i+C16*tmp9i;
	tmp22i = -C20*tmp1r+C22*tmp3r-C14*tmp5r-C16*tmp7r-C18*tmp9r;
	tmp22r = C20*tmp1i-C22*tmp3i+C14*tmp5i+C16*tmp7i+C18*tmp9i;
	tmp23i = C22*tmp1r-C14*tmp3r-C16*tmp5r-C18*tmp7r+C20*tmp9r;
	tmp23r = -C22*tmp1i+C14*tmp3i+C16*tmp5i+C18*tmp7i-C20*tmp9i;
	tmp24i = tmp10i+tmp19i;
	tmp24r = tmp10r+tmp19r;
	tmp25i = tmp10i-tmp19i;
	tmp25r = tmp10r-tmp19r;
	tmp26i = tmp12i+tmp20i;
	tmp26r = tmp12r+tmp20r;
	tmp27i = tmp12i-tmp20i;
	tmp27r = tmp12r-tmp20r;
	tmp28i = tmp14i+tmp21i;
	tmp28r = tmp14r+tmp21r;
	tmp29i = tmp14i-tmp21i;
	tmp29r = tmp14r-tmp21r;
	tmp30i = tmp16i+tmp22i;
	tmp30r = tmp16r+tmp22r;
	tmp31i = tmp16i-tmp22i;
	tmp31r = tmp16r-tmp22r;
	tmp32i = tmp18i+tmp23i;
	tmp32r = tmp18r+tmp23r;
	tmp33i = tmp18i-tmp23i;
	tmp33r = tmp18r-tmp23r;
	tmp34r = Re[0]+tmp17r;
	tmp35i = Im[1]-Im[10];
	tmp35r = Re[1]+Re[10];
	tmp36i = Im[1]+Im[10];
	tmp36r = Re[1]-Re[10];
	tmp37i = -(Im[2]-Im[9]);
	tmp37r = Re[2]+Re[9];
	tmp38i = -(Im[2]+Im[9]);
	tmp38r = Re[2]-Re[9];
	tmp39i = -(Im[4]-Im[7]);
	tmp39r = Re[4]+Re[7];
	tmp40i = -(Im[4]+Im[7]);
	tmp40r = Re[4]-Re[7];
	tmp41i = -(Im[8]-Im[3]);
	tmp41r = Re[8]+Re[3];
	tmp42i = -(Im[8]+Im[3]);
	tmp42r = Re[8]-Re[3];
	tmp43i = Im[5]-Im[6];
	tmp43r = Re[5]+Re[6];
	tmp44i = Im[5]+Im[6];
	tmp44r = Re[5]-Re[6];
	tmp45i = C13*tmp35i+C15*tmp37i-C17*tmp39i-C19*tmp41i-C21*tmp43i-Im[0];
	tmp45r = C13*tmp35r+C15*tmp37r-C17*tmp39r-C19*tmp41r-C21*tmp43r+Re[0];
	tmp46r = tmp35r+tmp37r;
	tmp47i = C15*tmp35i-C17*tmp37i-C19*tmp39i-C21*tmp41i+C13*tmp43i-Im[0];
	tmp47r = C15*tmp35r-C17*tmp37r-C19*tmp39r-C21*tmp41r+C13*tmp43r+Re[0];
	tmp48r = tmp46r+tmp39r;
	tmp49i = -C17*tmp35i-C19*tmp37i-C21*tmp39i+C13*tmp41i+C15*tmp43i-Im[0];
	tmp49r = -C17*tmp35r-C19*tmp37r-C21*tmp39r+C13*tmp41r+C15*tmp43r+Re[0];
	tmp50r = tmp48r+tmp41r;
	tmp51i = -C19*tmp35i-C21*tmp37i+C13*tmp39i+C15*tmp41i-C17*tmp43i-Im[0];
	tmp51r = -C19*tmp35r-C21*tmp37r+C13*tmp39r+C15*tmp41r-C17*tmp43r+Re[0];
	tmp52r = tmp50r+tmp43r;
	tmp53i = -C21*tmp35i+C13*tmp37i+C15*tmp39i-C17*tmp41i-C19*tmp43i-Im[0];
	tmp53r = -C21*tmp35r+C13*tmp37r+C15*tmp39r-C17*tmp41r-C19*tmp43r+Re[0];
	tmp54i = C14*tmp36r+C16*tmp38r+C18*tmp40r-C20*tmp42r+C22*tmp44r;
	tmp54r = -C14*tmp36i-C16*tmp38i-C18*tmp40i+C20*tmp42i-C22*tmp44i;
	tmp55i = C16*tmp36r+C18*tmp38r-C20*tmp40r+C22*tmp42r-C14*tmp44r;
	tmp55r = -C16*tmp36i-C18*tmp38i+C20*tmp40i-C22*tmp42i+C14*tmp44i;
	tmp56i = C18*tmp36r-C20*tmp38r+C22*tmp40r-C14*tmp42r-C16*tmp44r;
	tmp56r = -C18*tmp36i+C20*tmp38i-C22*tmp40i+C14*tmp42i+C16*tmp44i;
	tmp57i = -C20*tmp36r+C22*tmp38r-C14*tmp40r-C16*tmp42r-C18*tmp44r;
	tmp57r = C20*tmp36i-C22*tmp38i+C14*tmp40i+C16*tmp42i+C18*tmp44i;
	tmp58i = C22*tmp36r-C14*tmp38r-C16*tmp40r-C18*tmp42r+C20*tmp44r;
	tmp58r = -C22*tmp36i+C14*tmp38i+C16*tmp40i+C18*tmp42i-C20*tmp44i;
	tmp59i = tmp45i+tmp54i;
	tmp59r = tmp45r+tmp54r;
	tmp60i = tmp45i-tmp54i;
	tmp60r = tmp45r-tmp54r;
	tmp61i = tmp47i+tmp55i;
	tmp61r = tmp47r+tmp55r;
	tmp62i = tmp47i-tmp55i;
	tmp62r = tmp47r-tmp55r;
	tmp63i = tmp49i+tmp56i;
	tmp63r = tmp49r+tmp56r;
	tmp64i = tmp49i-tmp56i;
	tmp64r = tmp49r-tmp56r;
	tmp65i = tmp51i+tmp57i;
	tmp65r = tmp51r+tmp57r;
	tmp66i = tmp51i-tmp57i;
	tmp66r = tmp51r-tmp57r;
	tmp67i = tmp53i+tmp58i;
	tmp67r = tmp53r+tmp58r;
	tmp68i = tmp53i-tmp58i;
	tmp68r = tmp53r-tmp58r;
	tmp69r = Re[0]+tmp52r;
	tmp70r = tmp34r+tmp69r;
	tmp73i = tmp33i-tmp68i;
	tmp73r = tmp33r-tmp68r;
	tmp74i = tmp24i+tmp59i;
	tmp74r = tmp24r+tmp59r;
	tmp77i = tmp29i-tmp64i;
	tmp77r = tmp29r-tmp64r;
	tmp78i = tmp26i+tmp61i;
	tmp78r = tmp26r+tmp61r;
	tmp81i = tmp30i-tmp65i;
	tmp81r = tmp30r-tmp65r;
	tmp82i = tmp31i+tmp66i;
	tmp82r = tmp31r+tmp66r;
	tmp85i = tmp27i-tmp62i;
	tmp85r = tmp27r-tmp62r;
	tmp86i = tmp28i+tmp63i;
	tmp86r = tmp28r+tmp63r;
	tmp89i = tmp25i-tmp60i;
	tmp89r = tmp25r-tmp60r;
	tmp90i = tmp32i+tmp67i;
	tmp90r = tmp32r+tmp67r;
	Re[0] = tmp70r;
	Im[0] = C0;
	Re[1] = tmp73r;
	Im[1] = tmp73i;
	Re[2] = tmp74r;
	Im[2] = tmp74i;
	Re[3] = tmp77r;
	Im[3] = tmp77i;
	Re[4] = tmp78r;
	Im[4] = tmp78i;
	Re[5] = tmp81r;
	Im[5] = tmp81i;
	Re[6] = tmp82r;
	Im[6] = tmp82i;
	Re[7] = tmp85r;
	Im[7] = tmp85i;
	Re[8] = tmp86r;
	Im[8] = tmp86i;
	Re[9] = tmp89r;
	Im[9] = tmp89i;
	Re[10] = tmp90r;
	Im[10] = tmp90i;
}

/*
*	Number of additions = 204
*	Number of multiplications = 42
*	Number of sign changes = 35
*	Number of assigns = 258
*	Total number of operations = 539
*/
void	MFFTR24(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, /*tmp81r,*/ tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, /*tmp88r,*/ tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, /*tmp98r,*/
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, /*tmp105r,*/
		tmp106r, tmp107r, tmp108r, tmp109r, tmp110r, tmp111r, /*tmp112r,*/
		tmp113r, /*tmp114r,*/ tmp115r, tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, /*tmp121r,*/ /*tmp122r,*/ tmp123r, /*tmp124r,*/ /*tmp125r,*/ /*tmp126r,*/
		/*tmp127r,*/ /*tmp128r,*/ /*tmp129r,*/ tmp130r, tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r, /*tmp136r,*/ tmp137r /*tmp138r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, /*tmp69i,*/ tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, /*tmp76i,*/ tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, /*tmp86i,*/ tmp87i, tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, /*tmp93i,*/ tmp94i, tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, /*tmp100i,*/ tmp101i, /*tmp102i,*/ tmp103i, tmp104i, tmp105i,
		tmp106i, tmp107i, tmp108i, /*tmp109i,*/ /*tmp110i,*/ tmp111i, /*tmp112i,*/
		/*tmp113i,*/ /*tmp114i,*/ /*tmp115i,*/ /*tmp116i,*/ /*tmp117i,*/ tmp118i, tmp119i,
		tmp120i, tmp121i, tmp122i, tmp123i, /*tmp124i,*/ tmp125i /*tmp126i*/;

	const FFT_precision	C4 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C3 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =                    2;	/* INTEGERCONST	*/

	tmp0r = C2*Re[6];
	tmp0i = -C2*Im[6];
	tmp1r = Re[0]+tmp0r;
	tmp2r = Re[0]-tmp0r;
	tmp3r = Re[0]+tmp0i;
	tmp4r = Re[0]-tmp0i;
	tmp1i = Im[9]+Im[3];
	tmp5r = Re[9]-Re[3];
	tmp2i = Im[9]-Im[3];
	tmp6r = Re[9]+Re[3];
	tmp3i = Im[3]+Im[9];
	tmp7r = Re[3]-Re[9];
	tmp4i = Im[3]-Im[9];
	tmp8r = Re[3]+Re[9];
	tmp5i = tmp1i+tmp3i;
	tmp9r = tmp5r+tmp7r;
	tmp6i = tmp1i-tmp3i;
	tmp10r = tmp5r-tmp7r;
	tmp7i = tmp2i-tmp8r;
	tmp11r = tmp6r+tmp4i;
	tmp8i = tmp2i+tmp8r;
	tmp12r = tmp6r-tmp4i;
	tmp9i = C3*(tmp7i-tmp11r);
	tmp13r = C3*(tmp11r+tmp7i);
	tmp10i = -C3*(tmp8i+tmp12r);
	tmp14r = -C3*(tmp12r-tmp8i);
	tmp15r = tmp1r+tmp9r;
	tmp16r = tmp1r-tmp9r;
	tmp17r = tmp3r+tmp13r;
	tmp18r = tmp3r-tmp13r;
	tmp19r = tmp2r+tmp6i;
	tmp20r = tmp2r-tmp6i;
	tmp21r = tmp4r+tmp14r;
	tmp22r = tmp4r-tmp14r;
	tmp11i = -(Im[8]-Im[4]);
	tmp23r = Re[8]+Re[4];
	tmp12i = -(Im[8]+Im[4]);
	tmp24r = Re[8]-Re[4];
	tmp13i = Im[10]-Im[2];
	tmp25r = Re[10]+Re[2];
	tmp14i = Im[10]+Im[2];
	tmp26r = Re[10]-Re[2];
	tmp15i = tmp11i+tmp13i;
	tmp27r = tmp23r+tmp25r;
	tmp16i = tmp11i-tmp13i;
	tmp28r = tmp23r-tmp25r;
	tmp17i = tmp12i-tmp26r;
	tmp29r = tmp24r+tmp14i;
	tmp18i = tmp12i+tmp26r;
	tmp30r = tmp24r-tmp14i;
	tmp19i = Im[1]+Im[11];
	tmp31r = Re[1]-Re[11];
	tmp20i = Im[1]-Im[11];
	tmp32r = Re[1]+Re[11];
	tmp21i = Im[5]+Im[7];
	tmp33r = -(Re[5]-Re[7]);
	tmp22i = Im[5]-Im[7];
	tmp34r = -(Re[5]+Re[7]);
	tmp23i = tmp19i+tmp21i;
	tmp35r = tmp31r+tmp33r;
	tmp24i = tmp19i-tmp21i;
	tmp36r = tmp31r-tmp33r;
	tmp25i = tmp20i-tmp34r;
	tmp37r = tmp32r+tmp22i;
	tmp26i = tmp20i+tmp34r;
	tmp38r = tmp32r-tmp22i;
	tmp27i = C3*(tmp25i-tmp37r);
	tmp39r = C3*(tmp37r+tmp25i);
	tmp28i = -C3*(tmp26i+tmp38r);
	tmp40r = -C3*(tmp38r-tmp26i);
	tmp29i = tmp15i+tmp23i;
	tmp41r = tmp27r+tmp35r;
	tmp30i = tmp15i-tmp23i;
	tmp42r = tmp27r-tmp35r;
	tmp31i = tmp17i+tmp27i;
	tmp43r = tmp29r+tmp39r;
	tmp32i = tmp17i-tmp27i;
	tmp44r = tmp29r-tmp39r;
	tmp33i = tmp16i-tmp36r;
	tmp45r = tmp28r+tmp24i;
	tmp34i = tmp16i+tmp36r;
	tmp46r = tmp28r-tmp24i;
	tmp35i = tmp18i+tmp28i;
	tmp47r = tmp30r+tmp40r;
	tmp36i = tmp18i-tmp28i;
	tmp48r = tmp30r-tmp40r;
	tmp37i = Im[8]-Im[4];
	tmp49r = Re[8]+Re[4];
	tmp38i = Im[8]+Im[4];
	tmp50r = Re[8]-Re[4];
	tmp39i = Im[2]-Im[10];
	tmp51r = Re[2]+Re[10];
	tmp40i = Im[2]+Im[10];
	tmp52r = Re[2]-Re[10];
	tmp41i = tmp37i+tmp39i;
	tmp53r = tmp49r+tmp51r;
	tmp42i = tmp37i-tmp39i;
	tmp54r = tmp49r-tmp51r;
	tmp43i = tmp38i-tmp52r;
	tmp55r = tmp50r+tmp40i;
	tmp44i = tmp38i+tmp52r;
	tmp56r = tmp50r-tmp40i;
	tmp45i = Im[7]+Im[5];
	tmp57r = -(Re[7]-Re[5]);
	tmp46i = Im[7]-Im[5];
	tmp58r = -(Re[7]+Re[5]);
	tmp47i = Im[11]+Im[1];
	tmp59r = Re[11]-Re[1];
	tmp48i = Im[11]-Im[1];
	tmp60r = Re[11]+Re[1];
	tmp49i = tmp45i+tmp47i;
	tmp61r = tmp57r+tmp59r;
	tmp50i = tmp45i-tmp47i;
	tmp62r = tmp57r-tmp59r;
	tmp51i = tmp46i-tmp60r;
	tmp63r = tmp58r+tmp48i;
	tmp52i = tmp46i+tmp60r;
	tmp64r = tmp58r-tmp48i;
	tmp53i = C3*(tmp51i-tmp63r);
	tmp65r = C3*(tmp63r+tmp51i);
	tmp54i = -C3*(tmp52i+tmp64r);
	tmp66r = -C3*(tmp64r-tmp52i);
	tmp55i = tmp41i+tmp49i;
	tmp67r = tmp53r+tmp61r;
	tmp56i = tmp41i-tmp49i;
	tmp68r = tmp53r-tmp61r;
	tmp57i = tmp43i+tmp53i;
	tmp69r = tmp55r+tmp65r;
	tmp58i = tmp43i-tmp53i;
	tmp70r = tmp55r-tmp65r;
	tmp59i = tmp42i-tmp62r;
	tmp71r = tmp54r+tmp50i;
	tmp60i = tmp42i+tmp62r;
	tmp72r = tmp54r-tmp50i;
	tmp61i = tmp44i+tmp54i;
	tmp73r = tmp56r+tmp66r;
	tmp62i = tmp44i-tmp54i;
	tmp74r = tmp56r-tmp66r;
	tmp63i = tmp29i+tmp55i;
	tmp75r = tmp41r+tmp67r;
	tmp64i = tmp29i-tmp55i;
	tmp76r = tmp41r-tmp67r;
	tmp65i = -C4*tmp63i;
	tmp77r = -C4*tmp75r;
	tmp66i = -C7*tmp76r;
	tmp78r = C7*tmp64i;
	tmp67i = tmp65i+tmp5i;
	tmp79r = tmp77r+tmp15r;
	tmp68i = tmp67i+tmp66i;
	tmp80r = tmp79r+tmp78r;
	tmp70i = tmp5i+tmp63i;
	tmp82r = tmp15r+tmp75r;
	tmp71i = tmp35i+tmp61i;
	tmp83r = tmp47r+tmp73r;
	tmp72i = tmp35i-tmp61i;
	tmp84r = tmp47r-tmp73r;
	tmp73i = -C4*tmp71i;
	tmp85r = -C4*tmp83r;
	tmp74i = -C7*tmp84r;
	tmp86r = C7*tmp72i;
	tmp75i = tmp73i+tmp10i;
	tmp87r = tmp85r+tmp21r;
	tmp77i = tmp75i-tmp74i;
	tmp89r = tmp87r-tmp86r;
	tmp78i = tmp10i+tmp71i;
	tmp90r = tmp21r+tmp83r;
	tmp79i = tmp34i+tmp60i;
	tmp91r = tmp46r+tmp72r;
	tmp80i = tmp34i-tmp60i;
	tmp92r = tmp46r-tmp72r;
	tmp81i = -C4*tmp79i;
	tmp93r = -C4*tmp91r;
	tmp82i = -C7*tmp92r;
	tmp94r = C7*tmp80i;
	tmp83i = tmp81i+tmp10r;
	tmp95r = tmp93r+tmp20r;
	tmp84i = tmp83i+tmp82i;
	tmp96r = tmp95r+tmp94r;
	tmp85i = tmp83i-tmp82i;
	tmp97r = tmp95r-tmp94r;
	tmp87i = tmp31i+tmp57i;
	tmp99r = tmp43r+tmp69r;
	tmp88i = tmp31i-tmp57i;
	tmp100r = tmp43r-tmp69r;
	tmp89i = -C4*tmp87i;
	tmp101r = -C4*tmp99r;
	tmp90i = -C7*tmp100r;
	tmp102r = C7*tmp88i;
	tmp91i = tmp89i+tmp9i;
	tmp103r = tmp101r+tmp17r;
	tmp92i = tmp91i+tmp90i;
	tmp104r = tmp103r+tmp102r;
	tmp94i = tmp9i+tmp87i;
	tmp106r = tmp17r+tmp99r;
	tmp95i = tmp30i+tmp56i;
	tmp107r = tmp42r+tmp68r;
	tmp96i = tmp30i-tmp56i;
	tmp108r = tmp42r-tmp68r;
	tmp97i = -C4*tmp95i;
	tmp109r = -C4*tmp107r;
	tmp98i = -C7*tmp108r;
	tmp110r = C7*tmp96i;
	tmp99i = tmp97i-tmp5i;
	tmp111r = tmp109r+tmp16r;
	tmp101i = tmp99i-tmp98i;
	tmp113r = tmp111r-tmp110r;
	tmp103i = tmp36i+tmp62i;
	tmp115r = tmp48r+tmp74r;
	tmp104i = tmp36i-tmp62i;
	tmp116r = tmp48r-tmp74r;
	tmp105i = -C4*tmp103i;
	tmp117r = -C4*tmp115r;
	tmp106i = -C7*tmp116r;
	tmp118r = C7*tmp104i;
	tmp107i = tmp105i-tmp10i;
	tmp119r = tmp117r+tmp22r;
	tmp108i = tmp107i+tmp106i;
	tmp120r = tmp119r+tmp118r;
	tmp111i = tmp33i+tmp59i;
	tmp123r = tmp45r+tmp71r;
	tmp118i = -(tmp10r-tmp111i);
	tmp130r = tmp19r+tmp123r;
	tmp119i = tmp32i+tmp58i;
	tmp131r = tmp44r+tmp70r;
	tmp120i = tmp32i-tmp58i;
	tmp132r = tmp44r-tmp70r;
	tmp121i = -C4*tmp119i;
	tmp133r = -C4*tmp131r;
	tmp122i = -C7*tmp132r;
	tmp134r = C7*tmp120i;
	tmp123i = tmp121i-tmp9i;
	tmp135r = tmp133r+tmp18r;
	tmp125i = tmp123i-tmp122i;
	tmp137r = tmp135r-tmp134r;
	Re[0] = tmp82r;
	Im[0] = tmp70i;
	Re[1] = tmp89r;
	Im[1] = tmp77i;
	Re[2] = tmp96r;
	Im[2] = tmp84i;
	Re[3] = tmp106r;
	Im[3] = tmp94i;
	Re[4] = tmp113r;
	Im[4] = tmp101i;
	Re[5] = tmp120r;
	Im[5] = tmp108i;
	Re[6] = tmp130r;
	Im[6] = tmp118i;
	Re[7] = tmp137r;
	Im[7] = tmp125i;
	Re[8] = tmp80r;
	Im[8] = tmp68i;
	Re[9] = tmp90r;
	Im[9] = tmp78i;
	Re[10] = tmp97r;
	Im[10] = tmp85i;
	Re[11] = tmp104r;
	Im[11] = tmp92i;
}

/*
*	Number of additions = 111
*	Number of multiplications = 60
*	Number of sign changes = 39
*	Number of assigns = 182
*	Total number of operations = 392
*/
void	MIFFTR24(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, /*tmp58r,*/ tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, /*tmp68r,*/ tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, /*tmp75r,*/ tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, /*tmp81r,*/ tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, /*tmp88r,*/ /*tmp89r,*/ tmp90r, /*tmp91r,*/
		/*tmp92r,*/ /*tmp93r,*/ /*tmp94r,*/ /*tmp95r,*/ /*tmp96r,*/ tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, /*tmp103r,*/ tmp104r /*tmp105r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		/*tmp36i,*/ tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, /*tmp46i,*/ tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, /*tmp53i,*/ tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, /*tmp62i,*/ /*tmp63i,*/
		tmp64i, /*tmp65i,*/ /*tmp66i,*/ /*tmp67i,*/ /*tmp68i,*/ /*tmp69i,*/ /*tmp70i,*/
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, /*tmp77i,*/
		tmp78i /*tmp79i*/;

	const FFT_precision	C4 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C3 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =                    2;	/* INTEGERCONST	*/

	tmp0r = C2*Re[0];
	tmp0i = C2*Im[0];
	tmp1r = C2*Re[6];
	tmp1i = -C2*Im[6];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = C2*Re[9];
	tmp2i = C2*Im[9];
	tmp5r = C2*Re[3];
	tmp3i = C2*Im[3];
	tmp6r = tmp4r+tmp5r;
	tmp7r = tmp4r-tmp5r;
	tmp4i = C3*(tmp2i-tmp3i);
	tmp8r = -C3*(tmp3i+tmp2i);
	tmp5i = -C3*(tmp2i-tmp3i);
	tmp9r = -C3*(tmp3i+tmp2i);
	tmp10r = tmp2r+tmp6r;
	tmp11r = tmp2r-tmp6r;
	tmp6i = tmp0i+tmp4i;
	tmp12r = -(tmp1i-tmp8r);
	tmp7i = tmp0i-tmp4i;
	tmp13r = -(tmp1i+tmp8r);
	tmp8i = tmp0i+tmp5i;
	tmp14r = tmp1i+tmp9r;
	tmp9i = tmp0i-tmp5i;
	tmp15r = tmp1i-tmp9r;
	tmp16r = C2*Re[4];
	tmp10i = -C2*Im[4];
	tmp17r = C2*Re[10];
	tmp11i = C2*Im[10];
	tmp18r = tmp16r+tmp17r;
	tmp19r = tmp16r-tmp17r;
	tmp20r = C2*Re[1];
	tmp12i = C2*Im[1];
	tmp21r = C2*Re[7];
	tmp13i = -C2*Im[7];
	tmp22r = tmp20r+tmp21r;
	tmp23r = tmp20r-tmp21r;
	tmp14i = C3*(tmp12i-tmp13i);
	tmp24r = -C3*(tmp13i+tmp12i);
	tmp15i = -C3*(tmp12i-tmp13i);
	tmp25r = -C3*(tmp13i+tmp12i);
	tmp26r = tmp18r+tmp22r;
	tmp27r = tmp18r-tmp22r;
	tmp16i = tmp10i+tmp14i;
	tmp28r = -(tmp11i-tmp24r);
	tmp17i = tmp10i-tmp14i;
	tmp29r = -(tmp11i+tmp24r);
	tmp18i = tmp10i+tmp15i;
	tmp30r = tmp11i+tmp25r;
	tmp19i = tmp10i-tmp15i;
	tmp31r = tmp11i-tmp25r;
	tmp32r = C2*Re[8];
	tmp20i = C2*Im[8];
	tmp33r = C2*Re[2];
	tmp21i = C2*Im[2];
	tmp34r = tmp32r+tmp33r;
	tmp35r = tmp32r-tmp33r;
	tmp36r = C2*Re[5];
	tmp22i = -C2*Im[5];
	tmp37r = C2*Re[11];
	tmp23i = C2*Im[11];
	tmp38r = tmp36r+tmp37r;
	tmp39r = tmp36r-tmp37r;
	tmp24i = C3*(tmp22i-tmp23i);
	tmp40r = -C3*(tmp23i+tmp22i);
	tmp25i = -C3*(tmp22i-tmp23i);
	tmp41r = -C3*(tmp23i+tmp22i);
	tmp42r = tmp34r+tmp38r;
	tmp43r = tmp34r-tmp38r;
	tmp26i = tmp20i+tmp24i;
	tmp44r = -(tmp21i-tmp40r);
	tmp27i = tmp20i-tmp24i;
	tmp45r = -(tmp21i+tmp40r);
	tmp28i = tmp20i+tmp25i;
	tmp46r = tmp21i+tmp41r;
	tmp29i = tmp20i-tmp25i;
	tmp47r = tmp21i-tmp41r;
	tmp48r = tmp26r+tmp42r;
	tmp49r = tmp26r-tmp42r;
	tmp50r = -C4*tmp48r;
	tmp30i = C7*tmp49r;
	tmp51r = tmp50r+tmp10r;
	tmp52r = tmp10r+tmp48r;
	tmp31i = tmp18i+tmp28i;
	tmp53r = tmp30r+tmp46r;
	tmp32i = tmp18i-tmp28i;
	tmp54r = tmp30r-tmp46r;
	tmp33i = -C4*tmp31i;
	tmp55r = -C4*tmp53r;
	tmp34i = C7*tmp54r;
	tmp56r = -C7*tmp32i;
	tmp35i = tmp33i+tmp8i;
	tmp57r = tmp55r+tmp14r;
	tmp37i = tmp35i-tmp34i;
	tmp59r = tmp57r-tmp56r;
	tmp38i = tmp8i+tmp31i;
	tmp60r = tmp14r+tmp53r;
	tmp39i = -(tmp23r+tmp39r);
	tmp61r = tmp19r+tmp35r;
	tmp40i = -(tmp23r-tmp39r);
	tmp62r = tmp19r-tmp35r;
	tmp41i = -C4*tmp39i;
	tmp63r = -C4*tmp61r;
	tmp42i = C7*tmp62r;
	tmp64r = -C7*tmp40i;
	tmp43i = tmp41i-tmp7r;
	tmp65r = tmp63r+tmp3r;
	tmp44i = tmp43i+tmp42i;
	tmp66r = tmp65r+tmp64r;
	tmp45i = tmp43i-tmp42i;
	tmp67r = tmp65r-tmp64r;
	tmp47i = tmp16i+tmp26i;
	tmp69r = tmp28r+tmp44r;
	tmp48i = tmp16i-tmp26i;
	tmp70r = tmp28r-tmp44r;
	tmp49i = -C4*tmp47i;
	tmp71r = -C4*tmp69r;
	tmp50i = C7*tmp70r;
	tmp72r = -C7*tmp48i;
	tmp51i = tmp49i+tmp6i;
	tmp73r = tmp71r+tmp12r;
	tmp52i = tmp51i+tmp50i;
	tmp74r = tmp73r+tmp72r;
	tmp54i = tmp6i+tmp47i;
	tmp76r = tmp12r+tmp69r;
	tmp77r = tmp27r+tmp43r;
	tmp78r = tmp27r-tmp43r;
	tmp79r = -C4*tmp77r;
	tmp55i = C7*tmp78r;
	tmp80r = tmp79r+tmp11r;
	tmp56i = tmp19i+tmp29i;
	tmp82r = tmp31r+tmp47r;
	tmp57i = tmp19i-tmp29i;
	tmp83r = tmp31r-tmp47r;
	tmp58i = -C4*tmp56i;
	tmp84r = -C4*tmp82r;
	tmp59i = C7*tmp83r;
	tmp85r = -C7*tmp57i;
	tmp60i = tmp58i+tmp9i;
	tmp86r = tmp84r+tmp15r;
	tmp61i = tmp60i+tmp59i;
	tmp87r = tmp86r+tmp85r;
	tmp64i = tmp23r+tmp39r;
	tmp90r = tmp19r+tmp35r;
	tmp71i = tmp7r+tmp64i;
	tmp97r = tmp3r+tmp90r;
	tmp72i = tmp17i+tmp27i;
	tmp98r = tmp29r+tmp45r;
	tmp73i = tmp17i-tmp27i;
	tmp99r = tmp29r-tmp45r;
	tmp74i = -C4*tmp72i;
	tmp100r = -C4*tmp98r;
	tmp75i = C7*tmp99r;
	tmp101r = -C7*tmp73i;
	tmp76i = tmp74i+tmp7i;
	tmp102r = tmp100r+tmp13r;
	tmp78i = tmp76i-tmp75i;
	tmp104r = tmp102r-tmp101r;
	Re[0] = tmp52r;
	Im[0] = 0.0;
	Re[1] = tmp59r;
	Im[1] = tmp37i;
	Re[2] = tmp66r;
	Im[2] = tmp44i;
	Re[3] = tmp76r;
	Im[3] = tmp54i;
	Re[4] = tmp80r;
	Im[4] = -tmp55i;
	Re[5] = tmp87r;
	Im[5] = tmp61i;
	Re[6] = tmp97r;
	Im[6] = tmp71i;
	Re[7] = tmp104r;
	Im[7] = tmp78i;
	Re[8] = tmp51r;
	Im[8] = tmp30i;
	Re[9] = tmp60r;
	Im[9] = tmp38i;
	Re[10] = tmp67r;
	Im[10] = tmp45i;
	Re[11] = tmp74r;
	Im[11] = tmp52i;
}

/*
*	Number of additions = 161
*	Number of multiplications = 168
*	Number of sign changes = 12
*	Number of assigns = 109
*	Total number of operations = 450
*/
void	MFFTR26(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i;

	const FFT_precision	C23 =     0.12053668025532;	/* FFT_precisionCONST	*/
	const FFT_precision	C26 =     0.23931566428756;	/* FFT_precisionCONST	*/
	const FFT_precision	C19 =     0.35460488704254;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =     0.46472317204377;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.56806474673116;	/* FFT_precisionCONST	*/
	const FFT_precision	C22 =      0.6631226582408;	/* FFT_precisionCONST	*/
	const FFT_precision	C21 =      0.7485107481711;	/* FFT_precisionCONST	*/
	const FFT_precision	C18 =     0.82298386589366;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.88545602565321;	/* FFT_precisionCONST	*/
	const FFT_precision	C20 =     0.93501624268541;	/* FFT_precisionCONST	*/
	const FFT_precision	C25 =     0.97094181742605;	/* FFT_precisionCONST	*/
	const FFT_precision	C24 =     0.99270887409805;	/* FFT_precisionCONST	*/
	const FFT_precision	C28 =                    2;	/* INTEGERCONST	*/

	tmp0r = C28*Re[12];
	tmp0i = -C28*Im[12];
	tmp1r = C28*Re[2];
	tmp1i = C28*Im[2];
	tmp2r = C28*Re[4];
	tmp2i = C28*Im[4];
	tmp3r = C28*Re[8];
	tmp3i = C28*Im[8];
	tmp4r = C28*Re[10];
	tmp4i = -C28*Im[10];
	tmp5r = C28*Re[6];
	tmp5i = C28*Im[6];
	tmp6r = C15*tmp0r+C17*tmp1r-C19*tmp2r-C21*tmp3r+C23*tmp4r-C25*tmp5r+Re[0];
	tmp7r = tmp0r+tmp1r;
	tmp8r = C17*tmp0r-C19*tmp1r-C21*tmp2r+C23*tmp3r-C25*tmp4r+C15*tmp5r+Re[0];
	tmp9r = tmp7r+tmp2r;
	tmp10r = -C19*tmp0r-C21*tmp1r+C23*tmp2r-C25*tmp3r+C15*tmp4r+C17*tmp5r+Re[0];
	tmp11r = tmp9r+tmp3r;
	tmp12r = -C21*tmp0r+C23*tmp1r-C25*tmp2r+C15*tmp3r+C17*tmp4r-C19*tmp5r+Re[0];
	tmp13r = tmp11r+tmp4r;
	tmp14r = C23*tmp0r-C25*tmp1r+C15*tmp2r+C17*tmp3r-C19*tmp4r-C21*tmp5r+Re[0];
	tmp15r = tmp13r+tmp5r;
	tmp16r = -C25*tmp0r+C15*tmp1r+C17*tmp2r-C19*tmp3r-C21*tmp4r+C23*tmp5r+Re[0];
	tmp17r = C16*tmp0i+C18*tmp1i+C20*tmp2i-C22*tmp3i+C24*tmp4i+C26*tmp5i;
	tmp18r = C18*tmp0i+C20*tmp1i-C22*tmp2i+C24*tmp3i+C26*tmp4i-C16*tmp5i;
	tmp19r = C20*tmp0i-C22*tmp1i+C24*tmp2i+C26*tmp3i-C16*tmp4i-C18*tmp5i;
	tmp20r = -C22*tmp0i+C24*tmp1i+C26*tmp2i-C16*tmp3i-C18*tmp4i-C20*tmp5i;
	tmp21r = C24*tmp0i+C26*tmp1i-C16*tmp2i-C18*tmp3i-C20*tmp4i+C22*tmp5i;
	tmp22r = C26*tmp0i-C16*tmp1i-C18*tmp2i-C20*tmp3i+C22*tmp4i-C24*tmp5i;
	tmp23r = tmp6r+tmp17r;
	tmp24r = tmp6r-tmp17r;
	tmp25r = tmp8r+tmp18r;
	tmp26r = tmp8r-tmp18r;
	tmp27r = tmp10r+tmp19r;
	tmp28r = tmp10r-tmp19r;
	tmp29r = tmp12r+tmp20r;
	tmp30r = tmp12r-tmp20r;
	tmp31r = tmp14r+tmp21r;
	tmp32r = tmp14r-tmp21r;
	tmp33r = tmp16r+tmp22r;
	tmp34r = tmp16r-tmp22r;
	tmp35r = Re[0]+tmp15r;
	tmp6i = C28*Im[1];
	tmp36r = C28*Re[1];
	tmp7i = C28*Im[11];
	tmp37r = -C28*Re[11];
	tmp8i = C28*Im[9];
	tmp38r = -C28*Re[9];
	tmp9i = C28*Im[5];
	tmp39r = -C28*Re[5];
	tmp10i = C28*Im[3];
	tmp40r = C28*Re[3];
	tmp11i = C28*Im[7];
	tmp41r = -C28*Re[7];
	tmp12i = C15*tmp6i+C17*tmp7i-C19*tmp8i-C21*tmp9i+C23*tmp10i-C25*tmp11i;
	tmp13i = tmp6i+tmp7i;
	tmp14i = C17*tmp6i-C19*tmp7i-C21*tmp8i+C23*tmp9i-C25*tmp10i+C15*tmp11i;
	tmp15i = tmp13i+tmp8i;
	tmp16i = -C19*tmp6i-C21*tmp7i+C23*tmp8i-C25*tmp9i+C15*tmp10i+C17*tmp11i;
	tmp17i = tmp15i+tmp9i;
	tmp18i = -C21*tmp6i+C23*tmp7i-C25*tmp8i+C15*tmp9i+C17*tmp10i-C19*tmp11i;
	tmp19i = tmp17i+tmp10i;
	tmp20i = C23*tmp6i-C25*tmp7i+C15*tmp8i+C17*tmp9i-C19*tmp10i-C21*tmp11i;
	tmp21i = tmp19i+tmp11i;
	tmp22i = -C25*tmp6i+C15*tmp7i+C17*tmp8i-C19*tmp9i-C21*tmp10i+C23*tmp11i;
	tmp23i = -C16*tmp36r-C18*tmp37r-C20*tmp38r+C22*tmp39r-C24*tmp40r-C26*tmp41r;
	tmp24i = -C18*tmp36r-C20*tmp37r+C22*tmp38r-C24*tmp39r-C26*tmp40r+C16*tmp41r;
	tmp25i = -C20*tmp36r+C22*tmp37r-C24*tmp38r-C26*tmp39r+C16*tmp40r+C18*tmp41r;
	tmp26i = C22*tmp36r-C24*tmp37r-C26*tmp38r+C16*tmp39r+C18*tmp40r+C20*tmp41r;
	tmp27i = -C24*tmp36r-C26*tmp37r+C16*tmp38r+C18*tmp39r+C20*tmp40r-C22*tmp41r;
	tmp28i = -C26*tmp36r+C16*tmp37r+C18*tmp38r+C20*tmp39r-C22*tmp40r+C24*tmp41r;
	tmp29i = tmp12i+tmp23i;
	tmp30i = tmp12i-tmp23i;
	tmp31i = tmp14i+tmp24i;
	tmp32i = tmp14i-tmp24i;
	tmp33i = tmp16i+tmp25i;
	tmp34i = tmp16i-tmp25i;
	tmp35i = tmp18i+tmp26i;
	tmp36i = tmp18i-tmp26i;
	tmp37i = tmp20i+tmp27i;
	tmp38i = tmp20i-tmp27i;
	tmp39i = tmp22i+tmp28i;
	tmp40i = tmp22i-tmp28i;
	Re[0] = tmp35r;
	Im[0] = tmp21i;
	Re[1] = tmp34r;
	Im[1] = -tmp40i;
	Re[2] = tmp23r;
	Im[2] = tmp29i;
	Re[3] = tmp29r;
	Im[3] = -tmp35i;
	Re[4] = tmp25r;
	Im[4] = tmp31i;
	Re[5] = tmp28r;
	Im[5] = -tmp34i;
	Re[6] = tmp31r;
	Im[6] = tmp37i;
	Re[7] = tmp32r;
	Im[7] = -tmp38i;
	Re[8] = tmp27r;
	Im[8] = tmp33i;
	Re[9] = tmp26r;
	Im[9] = -tmp32i;
	Re[10] = tmp30r;
	Im[10] = tmp36i;
	Re[11] = tmp24r;
	Im[11] = -tmp30i;
	Re[12] = tmp33r;
	Im[12] = tmp39i;
}

/*
*	Number of additions = 397
*	Number of multiplications = 288
*	Number of sign changes = 12
*	Number of assigns = 207
*	Total number of operations = 904
*/
void	MIFFTR26(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		/*tmp85r,*/ /*tmp86r,*/ tmp87r, tmp88r, /*tmp89r,*/ /*tmp90r,*/ tmp91r,
		tmp92r, /*tmp93r,*/ /*tmp94r,*/ tmp95r, tmp96r, /*tmp97r,*/ /*tmp98r,*/
		tmp99r, tmp100r, /*tmp101r,*/ /*tmp102r,*/ tmp103r, tmp104r, /*tmp105r,*/
		/*tmp106r,*/ tmp107r, tmp108r /*tmp109r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, /*tmp13i,*/ tmp14i,
		/*tmp15i,*/ tmp16i, /*tmp17i,*/ tmp18i, /*tmp19i,*/ tmp20i, /*tmp21i,*/
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, /*tmp41i,*/ tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, /*tmp55i,*/ tmp56i,
		/*tmp57i,*/ tmp58i, /*tmp59i,*/ tmp60i, /*tmp61i,*/ tmp62i, /*tmp63i,*/
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, /*tmp83i,*/ /*tmp84i,*/
		/*tmp85i,*/ /*tmp86i,*/ tmp87i, tmp88i, /*tmp89i,*/ /*tmp90i,*/ tmp91i,
		tmp92i, /*tmp93i,*/ /*tmp94i,*/ tmp95i, tmp96i, /*tmp97i,*/ /*tmp98i,*/
		tmp99i, tmp100i, /*tmp101i,*/ /*tmp102i,*/ tmp103i, tmp104i, /*tmp105i,*/
		/*tmp106i,*/ tmp107i, tmp108i /*tmp109i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C23 =     0.12053668025532;	/* FFT_precisionCONST	*/
	const FFT_precision	C26 =     0.23931566428756;	/* FFT_precisionCONST	*/
	const FFT_precision	C19 =     0.35460488704254;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =     0.46472317204377;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.56806474673116;	/* FFT_precisionCONST	*/
	const FFT_precision	C22 =      0.6631226582408;	/* FFT_precisionCONST	*/
	const FFT_precision	C21 =      0.7485107481711;	/* FFT_precisionCONST	*/
	const FFT_precision	C18 =     0.82298386589366;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.88545602565321;	/* FFT_precisionCONST	*/
	const FFT_precision	C20 =     0.93501624268541;	/* FFT_precisionCONST	*/
	const FFT_precision	C25 =     0.97094181742605;	/* FFT_precisionCONST	*/
	const FFT_precision	C24 =     0.99270887409805;	/* FFT_precisionCONST	*/

	tmp0i = -(Im[1]-Im[12]);
	tmp0r = Re[1]+Re[12];
	tmp1i = -(Im[1]+Im[12]);
	tmp1r = Re[1]-Re[12];
	tmp2i = Im[2]-Im[11];
	tmp2r = Re[2]+Re[11];
	tmp3i = Im[2]+Im[11];
	tmp3r = Re[2]-Re[11];
	tmp4i = Im[4]-Im[9];
	tmp4r = Re[4]+Re[9];
	tmp5i = Im[4]+Im[9];
	tmp5r = Re[4]-Re[9];
	tmp6i = Im[8]-Im[5];
	tmp6r = Re[8]+Re[5];
	tmp7i = Im[8]+Im[5];
	tmp7r = Re[8]-Re[5];
	tmp8i = -(Im[3]-Im[10]);
	tmp8r = Re[3]+Re[10];
	tmp9i = -(Im[3]+Im[10]);
	tmp9r = Re[3]-Re[10];
	tmp10i = Im[6]-Im[7];
	tmp10r = Re[6]+Re[7];
	tmp11i = Im[6]+Im[7];
	tmp11r = Re[6]-Re[7];
	tmp12i = C15*tmp0i+C17*tmp2i-C19*tmp4i-C21*tmp6i+C23*tmp8i-C25*tmp10i+Im[0];
	tmp12r = C15*tmp0r+C17*tmp2r-C19*tmp4r-C21*tmp6r+C23*tmp8r-C25*tmp10r+Re[0];
	tmp13r = tmp0r+tmp2r;
	tmp14i = C17*tmp0i-C19*tmp2i-C21*tmp4i+C23*tmp6i-C25*tmp8i+C15*tmp10i+Im[0];
	tmp14r = C17*tmp0r-C19*tmp2r-C21*tmp4r+C23*tmp6r-C25*tmp8r+C15*tmp10r+Re[0];
	tmp15r = tmp13r+tmp4r;
	tmp16i = -C19*tmp0i-C21*tmp2i+C23*tmp4i-C25*tmp6i+C15*tmp8i+C17*tmp10i+Im[0];
	tmp16r = -C19*tmp0r-C21*tmp2r+C23*tmp4r-C25*tmp6r+C15*tmp8r+C17*tmp10r+Re[0];
	tmp17r = tmp15r+tmp6r;
	tmp18i = -C21*tmp0i+C23*tmp2i-C25*tmp4i+C15*tmp6i+C17*tmp8i-C19*tmp10i+Im[0];
	tmp18r = -C21*tmp0r+C23*tmp2r-C25*tmp4r+C15*tmp6r+C17*tmp8r-C19*tmp10r+Re[0];
	tmp19r = tmp17r+tmp8r;
	tmp20i = C23*tmp0i-C25*tmp2i+C15*tmp4i+C17*tmp6i-C19*tmp8i-C21*tmp10i+Im[0];
	tmp20r = C23*tmp0r-C25*tmp2r+C15*tmp4r+C17*tmp6r-C19*tmp8r-C21*tmp10r+Re[0];
	tmp21r = tmp19r+tmp10r;
	tmp22i = -C25*tmp0i+C15*tmp2i+C17*tmp4i-C19*tmp6i-C21*tmp8i+C23*tmp10i+Im[0];
	tmp22r = -C25*tmp0r+C15*tmp2r+C17*tmp4r-C19*tmp6r-C21*tmp8r+C23*tmp10r+Re[0];
	tmp23i = C16*tmp1r+C18*tmp3r+C20*tmp5r-C22*tmp7r+C24*tmp9r+C26*tmp11r;
	tmp23r = -C16*tmp1i-C18*tmp3i-C20*tmp5i+C22*tmp7i-C24*tmp9i-C26*tmp11i;
	tmp24i = C18*tmp1r+C20*tmp3r-C22*tmp5r+C24*tmp7r+C26*tmp9r-C16*tmp11r;
	tmp24r = -C18*tmp1i-C20*tmp3i+C22*tmp5i-C24*tmp7i-C26*tmp9i+C16*tmp11i;
	tmp25i = C20*tmp1r-C22*tmp3r+C24*tmp5r+C26*tmp7r-C16*tmp9r-C18*tmp11r;
	tmp25r = -C20*tmp1i+C22*tmp3i-C24*tmp5i-C26*tmp7i+C16*tmp9i+C18*tmp11i;
	tmp26i = -C22*tmp1r+C24*tmp3r+C26*tmp5r-C16*tmp7r-C18*tmp9r-C20*tmp11r;
	tmp26r = C22*tmp1i-C24*tmp3i-C26*tmp5i+C16*tmp7i+C18*tmp9i+C20*tmp11i;
	tmp27i = C24*tmp1r+C26*tmp3r-C16*tmp5r-C18*tmp7r-C20*tmp9r+C22*tmp11r;
	tmp27r = -C24*tmp1i-C26*tmp3i+C16*tmp5i+C18*tmp7i+C20*tmp9i-C22*tmp11i;
	tmp28i = C26*tmp1r-C16*tmp3r-C18*tmp5r-C20*tmp7r+C22*tmp9r-C24*tmp11r;
	tmp28r = -C26*tmp1i+C16*tmp3i+C18*tmp5i+C20*tmp7i-C22*tmp9i+C24*tmp11i;
	tmp29i = tmp12i+tmp23i;
	tmp29r = tmp12r+tmp23r;
	tmp30i = tmp12i-tmp23i;
	tmp30r = tmp12r-tmp23r;
	tmp31i = tmp14i+tmp24i;
	tmp31r = tmp14r+tmp24r;
	tmp32i = tmp14i-tmp24i;
	tmp32r = tmp14r-tmp24r;
	tmp33i = tmp16i+tmp25i;
	tmp33r = tmp16r+tmp25r;
	tmp34i = tmp16i-tmp25i;
	tmp34r = tmp16r-tmp25r;
	tmp35i = tmp18i+tmp26i;
	tmp35r = tmp18r+tmp26r;
	tmp36i = tmp18i-tmp26i;
	tmp36r = tmp18r-tmp26r;
	tmp37i = tmp20i+tmp27i;
	tmp37r = tmp20r+tmp27r;
	tmp38i = tmp20i-tmp27i;
	tmp38r = tmp20r-tmp27r;
	tmp39i = tmp22i+tmp28i;
	tmp39r = tmp22r+tmp28r;
	tmp40i = tmp22i-tmp28i;
	tmp40r = tmp22r-tmp28r;
	tmp41r = Re[0]+tmp21r;
	tmp42i = Im[1]-Im[12];
	tmp42r = Re[1]+Re[12];
	tmp43i = Im[1]+Im[12];
	tmp43r = Re[1]-Re[12];
	tmp44i = -(Im[2]-Im[11]);
	tmp44r = Re[2]+Re[11];
	tmp45i = -(Im[2]+Im[11]);
	tmp45r = Re[2]-Re[11];
	tmp46i = -(Im[4]-Im[9]);
	tmp46r = Re[4]+Re[9];
	tmp47i = -(Im[4]+Im[9]);
	tmp47r = Re[4]-Re[9];
	tmp48i = -(Im[8]-Im[5]);
	tmp48r = Re[8]+Re[5];
	tmp49i = -(Im[8]+Im[5]);
	tmp49r = Re[8]-Re[5];
	tmp50i = Im[3]-Im[10];
	tmp50r = Re[3]+Re[10];
	tmp51i = Im[3]+Im[10];
	tmp51r = Re[3]-Re[10];
	tmp52i = -(Im[6]-Im[7]);
	tmp52r = Re[6]+Re[7];
	tmp53i = -(Im[6]+Im[7]);
	tmp53r = Re[6]-Re[7];
	tmp54i = C15*tmp42i+C17*tmp44i-C19*tmp46i-C21*tmp48i+C23*tmp50i-C25*tmp52i-Im[0];
	tmp54r = C15*tmp42r+C17*tmp44r-C19*tmp46r-C21*tmp48r+C23*tmp50r-C25*tmp52r+Re[0];
	tmp55r = tmp42r+tmp44r;
	tmp56i = C17*tmp42i-C19*tmp44i-C21*tmp46i+C23*tmp48i-C25*tmp50i+C15*tmp52i-Im[0];
	tmp56r = C17*tmp42r-C19*tmp44r-C21*tmp46r+C23*tmp48r-C25*tmp50r+C15*tmp52r+Re[0];
	tmp57r = tmp55r+tmp46r;
	tmp58i = -C19*tmp42i-C21*tmp44i+C23*tmp46i-C25*tmp48i+C15*tmp50i+C17*tmp52i-Im[0];
	tmp58r = -C19*tmp42r-C21*tmp44r+C23*tmp46r-C25*tmp48r+C15*tmp50r+C17*tmp52r+Re[0];
	tmp59r = tmp57r+tmp48r;
	tmp60i = -C21*tmp42i+C23*tmp44i-C25*tmp46i+C15*tmp48i+C17*tmp50i-C19*tmp52i-Im[0];
	tmp60r = -C21*tmp42r+C23*tmp44r-C25*tmp46r+C15*tmp48r+C17*tmp50r-C19*tmp52r+Re[0];
	tmp61r = tmp59r+tmp50r;
	tmp62i = C23*tmp42i-C25*tmp44i+C15*tmp46i+C17*tmp48i-C19*tmp50i-C21*tmp52i-Im[0];
	tmp62r = C23*tmp42r-C25*tmp44r+C15*tmp46r+C17*tmp48r-C19*tmp50r-C21*tmp52r+Re[0];
	tmp63r = tmp61r+tmp52r;
	tmp64i = -C25*tmp42i+C15*tmp44i+C17*tmp46i-C19*tmp48i-C21*tmp50i+C23*tmp52i-Im[0];
	tmp64r = -C25*tmp42r+C15*tmp44r+C17*tmp46r-C19*tmp48r-C21*tmp50r+C23*tmp52r+Re[0];
	tmp65i = C16*tmp43r+C18*tmp45r+C20*tmp47r-C22*tmp49r+C24*tmp51r+C26*tmp53r;
	tmp65r = -C16*tmp43i-C18*tmp45i-C20*tmp47i+C22*tmp49i-C24*tmp51i-C26*tmp53i;
	tmp66i = C18*tmp43r+C20*tmp45r-C22*tmp47r+C24*tmp49r+C26*tmp51r-C16*tmp53r;
	tmp66r = -C18*tmp43i-C20*tmp45i+C22*tmp47i-C24*tmp49i-C26*tmp51i+C16*tmp53i;
	tmp67i = C20*tmp43r-C22*tmp45r+C24*tmp47r+C26*tmp49r-C16*tmp51r-C18*tmp53r;
	tmp67r = -C20*tmp43i+C22*tmp45i-C24*tmp47i-C26*tmp49i+C16*tmp51i+C18*tmp53i;
	tmp68i = -C22*tmp43r+C24*tmp45r+C26*tmp47r-C16*tmp49r-C18*tmp51r-C20*tmp53r;
	tmp68r = C22*tmp43i-C24*tmp45i-C26*tmp47i+C16*tmp49i+C18*tmp51i+C20*tmp53i;
	tmp69i = C24*tmp43r+C26*tmp45r-C16*tmp47r-C18*tmp49r-C20*tmp51r+C22*tmp53r;
	tmp69r = -C24*tmp43i-C26*tmp45i+C16*tmp47i+C18*tmp49i+C20*tmp51i-C22*tmp53i;
	tmp70i = C26*tmp43r-C16*tmp45r-C18*tmp47r-C20*tmp49r+C22*tmp51r-C24*tmp53r;
	tmp70r = -C26*tmp43i+C16*tmp45i+C18*tmp47i+C20*tmp49i-C22*tmp51i+C24*tmp53i;
	tmp71i = tmp54i+tmp65i;
	tmp71r = tmp54r+tmp65r;
	tmp72i = tmp54i-tmp65i;
	tmp72r = tmp54r-tmp65r;
	tmp73i = tmp56i+tmp66i;
	tmp73r = tmp56r+tmp66r;
	tmp74i = tmp56i-tmp66i;
	tmp74r = tmp56r-tmp66r;
	tmp75i = tmp58i+tmp67i;
	tmp75r = tmp58r+tmp67r;
	tmp76i = tmp58i-tmp67i;
	tmp76r = tmp58r-tmp67r;
	tmp77i = tmp60i+tmp68i;
	tmp77r = tmp60r+tmp68r;
	tmp78i = tmp60i-tmp68i;
	tmp78r = tmp60r-tmp68r;
	tmp79i = tmp62i+tmp69i;
	tmp79r = tmp62r+tmp69r;
	tmp80i = tmp62i-tmp69i;
	tmp80r = tmp62r-tmp69r;
	tmp81i = tmp64i+tmp70i;
	tmp81r = tmp64r+tmp70r;
	tmp82i = tmp64i-tmp70i;
	tmp82r = tmp64r-tmp70r;
	tmp83r = Re[0]+tmp63r;
	tmp84r = tmp41r+tmp83r;
	tmp87i = tmp40i-tmp82i;
	tmp87r = tmp40r-tmp82r;
	tmp88i = tmp29i+tmp71i;
	tmp88r = tmp29r+tmp71r;
	tmp91i = tmp35i-tmp77i;
	tmp91r = tmp35r-tmp77r;
	tmp92i = tmp31i+tmp73i;
	tmp92r = tmp31r+tmp73r;
	tmp95i = tmp34i-tmp76i;
	tmp95r = tmp34r-tmp76r;
	tmp96i = tmp37i+tmp79i;
	tmp96r = tmp37r+tmp79r;
	tmp99i = tmp38i-tmp80i;
	tmp99r = tmp38r-tmp80r;
	tmp100i = tmp33i+tmp75i;
	tmp100r = tmp33r+tmp75r;
	tmp103i = tmp32i-tmp74i;
	tmp103r = tmp32r-tmp74r;
	tmp104i = tmp36i+tmp78i;
	tmp104r = tmp36r+tmp78r;
	tmp107i = tmp30i-tmp72i;
	tmp107r = tmp30r-tmp72r;
	tmp108i = tmp39i+tmp81i;
	tmp108r = tmp39r+tmp81r;
	Re[0] = tmp84r;
	Im[0] = C0;
	Re[1] = tmp87r;
	Im[1] = tmp87i;
	Re[2] = tmp88r;
	Im[2] = tmp88i;
	Re[3] = tmp91r;
	Im[3] = tmp91i;
	Re[4] = tmp92r;
	Im[4] = tmp92i;
	Re[5] = tmp95r;
	Im[5] = tmp95i;
	Re[6] = tmp96r;
	Im[6] = tmp96i;
	Re[7] = tmp99r;
	Im[7] = tmp99i;
	Re[8] = tmp100r;
	Im[8] = tmp100i;
	Re[9] = tmp103r;
	Im[9] = tmp103i;
	Re[10] = tmp104r;
	Im[10] = tmp104i;
	Re[11] = tmp107r;
	Im[11] = tmp107i;
	Re[12] = tmp108r;
	Im[12] = tmp108i;
}

/*
*	Number of additions = 220
*	Number of multiplications = 120
*	Number of sign changes = 16
*	Number of assigns = 209
*	Total number of operations = 565
*/
void	MFFTR28(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, /*tmp82r,*/ tmp83r, /*tmp84r,*/
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, /*tmp90r,*/ /*tmp91r,*/
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, /*tmp97r,*/ tmp98r,
		/*tmp99r,*/ tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, /*tmp105r,*/
		tmp106r, tmp107r, /*tmp108r,*/ tmp109r, tmp110r, tmp111r, tmp112r,
		tmp113r, /*tmp114r,*/ tmp115r, /*tmp116r,*/ tmp117r, tmp118r, tmp119r,
		tmp120r, tmp121r, /*tmp122r,*/ /*tmp123r,*/ tmp124r, tmp125r, tmp126r,
		tmp127r, tmp128r, /*tmp129r,*/ tmp130r, /*tmp131r,*/ tmp132r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i;

	const FFT_precision	C13 =     0.22252093395631;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.43388373911756;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =     0.62348980185873;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.78183148246803;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.90096886790242;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =     0.97492791218182;	/* FFT_precisionCONST	*/
	const FFT_precision	C16 =                    2;	/* INTEGERCONST	*/

	tmp0r = C16*Re[8];
	tmp0i = C16*Im[8];
	tmp1r = C16*Re[4];
	tmp1i = -C16*Im[4];
	tmp2r = C16*Re[12];
	tmp2i = -C16*Im[12];
	tmp3r = C9*tmp0r-C11*tmp1r-C13*tmp2r+Re[0];
	tmp4r = tmp0r+tmp1r;
	tmp5r = -C11*tmp0r-C13*tmp1r+C9*tmp2r+Re[0];
	tmp6r = tmp4r+tmp2r;
	tmp7r = -C13*tmp0r+C9*tmp1r-C11*tmp2r+Re[0];
	tmp8r = C10*tmp0i+C12*tmp1i+C14*tmp2i;
	tmp9r = C12*tmp0i+C14*tmp1i-C10*tmp2i;
	tmp10r = C14*tmp0i-C10*tmp1i-C12*tmp2i;
	tmp11r = tmp3r+tmp8r;
	tmp12r = tmp3r-tmp8r;
	tmp13r = tmp5r+tmp9r;
	tmp14r = tmp5r-tmp9r;
	tmp15r = tmp7r+tmp10r;
	tmp16r = tmp7r-tmp10r;
	tmp17r = Re[0]+tmp6r;
	tmp3i = Im[1]+Im[13];
	tmp18r = Re[1]+Re[13];
	tmp4i = Im[1]-Im[13];
	tmp19r = Re[1]-Re[13];
	tmp5i = Im[11]+Im[3];
	tmp20r = -(Re[11]+Re[3]);
	tmp6i = Im[11]-Im[3];
	tmp21r = -(Re[11]-Re[3]);
	tmp7i = Im[9]+Im[5];
	tmp22r = Re[9]+Re[5];
	tmp8i = Im[9]-Im[5];
	tmp23r = Re[9]-Re[5];
	tmp9i = C9*tmp3i-C11*tmp5i-C13*tmp7i+Im[7];
	tmp24r = C9*tmp18r-C11*tmp20r-C13*tmp22r-Re[7];
	tmp10i = tmp3i+tmp5i;
	tmp25r = tmp18r+tmp20r;
	tmp11i = -C11*tmp3i-C13*tmp5i+C9*tmp7i+Im[7];
	tmp26r = -C11*tmp18r-C13*tmp20r+C9*tmp22r-Re[7];
	tmp12i = tmp10i+tmp7i;
	tmp27r = tmp25r+tmp22r;
	tmp13i = -C13*tmp3i+C9*tmp5i-C11*tmp7i+Im[7];
	tmp28r = -C13*tmp18r+C9*tmp20r-C11*tmp22r-Re[7];
	tmp14i = -C10*tmp19r-C12*tmp21r-C14*tmp23r;
	tmp29r = C10*tmp4i+C12*tmp6i+C14*tmp8i;
	tmp15i = -C12*tmp19r-C14*tmp21r+C10*tmp23r;
	tmp30r = C12*tmp4i+C14*tmp6i-C10*tmp8i;
	tmp16i = -C14*tmp19r+C10*tmp21r+C12*tmp23r;
	tmp31r = C14*tmp4i-C10*tmp6i-C12*tmp8i;
	tmp17i = tmp9i+tmp14i;
	tmp32r = tmp24r+tmp29r;
	tmp18i = tmp9i-tmp14i;
	tmp33r = tmp24r-tmp29r;
	tmp19i = tmp11i+tmp15i;
	tmp34r = tmp26r+tmp30r;
	tmp20i = tmp11i-tmp15i;
	tmp35r = tmp26r-tmp30r;
	tmp21i = tmp13i+tmp16i;
	tmp36r = tmp28r+tmp31r;
	tmp22i = tmp13i-tmp16i;
	tmp37r = tmp28r-tmp31r;
	tmp23i = Im[7]+tmp12i;
	tmp38r = -(Re[7]-tmp27r);
	tmp39r = C16*Re[6];
	tmp24i = -C16*Im[6];
	tmp40r = C16*Re[10];
	tmp25i = C16*Im[10];
	tmp41r = C16*Re[2];
	tmp26i = C16*Im[2];
	tmp42r = C9*tmp39r-C11*tmp40r-C13*tmp41r;
	tmp43r = tmp39r+tmp40r;
	tmp44r = -C11*tmp39r-C13*tmp40r+C9*tmp41r;
	tmp45r = tmp43r+tmp41r;
	tmp46r = -C13*tmp39r+C9*tmp40r-C11*tmp41r;
	tmp47r = C10*tmp24i+C12*tmp25i+C14*tmp26i;
	tmp48r = C12*tmp24i+C14*tmp25i-C10*tmp26i;
	tmp49r = C14*tmp24i-C10*tmp25i-C12*tmp26i;
	tmp50r = tmp42r+tmp47r;
	tmp51r = tmp42r-tmp47r;
	tmp52r = tmp44r+tmp48r;
	tmp53r = tmp44r-tmp48r;
	tmp54r = tmp46r+tmp49r;
	tmp55r = tmp46r-tmp49r;
	tmp27i = Im[13]+Im[1];
	tmp56r = -(Re[13]+Re[1]);
	tmp28i = Im[13]-Im[1];
	tmp57r = -(Re[13]-Re[1]);
	tmp29i = Im[3]+Im[11];
	tmp58r = Re[3]+Re[11];
	tmp30i = Im[3]-Im[11];
	tmp59r = Re[3]-Re[11];
	tmp31i = Im[5]+Im[9];
	tmp60r = -(Re[5]+Re[9]);
	tmp32i = Im[5]-Im[9];
	tmp61r = -(Re[5]-Re[9]);
	tmp33i = C9*tmp27i-C11*tmp29i-C13*tmp31i+Im[7];
	tmp62r = C9*tmp56r-C11*tmp58r-C13*tmp60r+Re[7];
	tmp34i = tmp27i+tmp29i;
	tmp63r = tmp56r+tmp58r;
	tmp35i = -C11*tmp27i-C13*tmp29i+C9*tmp31i+Im[7];
	tmp64r = -C11*tmp56r-C13*tmp58r+C9*tmp60r+Re[7];
	tmp36i = tmp34i+tmp31i;
	tmp65r = tmp63r+tmp60r;
	tmp37i = -C13*tmp27i+C9*tmp29i-C11*tmp31i+Im[7];
	tmp66r = -C13*tmp56r+C9*tmp58r-C11*tmp60r+Re[7];
	tmp38i = -C10*tmp57r-C12*tmp59r-C14*tmp61r;
	tmp67r = C10*tmp28i+C12*tmp30i+C14*tmp32i;
	tmp39i = -C12*tmp57r-C14*tmp59r+C10*tmp61r;
	tmp68r = C12*tmp28i+C14*tmp30i-C10*tmp32i;
	tmp40i = -C14*tmp57r+C10*tmp59r+C12*tmp61r;
	tmp69r = C14*tmp28i-C10*tmp30i-C12*tmp32i;
	tmp41i = tmp33i+tmp38i;
	tmp70r = tmp62r+tmp67r;
	tmp42i = tmp33i-tmp38i;
	tmp71r = tmp62r-tmp67r;
	tmp43i = tmp35i+tmp39i;
	tmp72r = tmp64r+tmp68r;
	tmp44i = tmp35i-tmp39i;
	tmp73r = tmp64r-tmp68r;
	tmp45i = tmp37i+tmp40i;
	tmp74r = tmp66r+tmp69r;
	tmp46i = tmp37i-tmp40i;
	tmp75r = tmp66r-tmp69r;
	tmp47i = Im[7]+tmp36i;
	tmp76r = Re[7]+tmp65r;
	tmp77r = tmp17r+tmp45r;
	tmp78r = tmp17r-tmp45r;
	tmp48i = tmp23i+tmp47i;
	tmp79r = tmp38r+tmp76r;
	tmp49i = tmp23i-tmp47i;
	tmp80r = tmp38r-tmp76r;
	tmp81r = tmp77r+tmp79r;
	tmp83r = tmp78r+tmp49i;
	tmp85r = tmp15r+tmp54r;
	tmp86r = tmp15r-tmp54r;
	tmp50i = tmp21i+tmp45i;
	tmp87r = tmp36r+tmp74r;
	tmp51i = tmp21i-tmp45i;
	tmp88r = tmp36r-tmp74r;
	tmp89r = tmp85r+tmp87r;
	tmp92r = tmp86r-tmp51i;
	tmp93r = tmp14r+tmp53r;
	tmp94r = tmp14r-tmp53r;
	tmp52i = tmp20i+tmp44i;
	tmp95r = tmp35r+tmp73r;
	tmp53i = tmp20i-tmp44i;
	tmp96r = tmp35r-tmp73r;
	tmp98r = tmp93r-tmp95r;
	tmp100r = tmp94r-tmp53i;
	tmp101r = tmp12r+tmp51r;
	tmp102r = tmp12r-tmp51r;
	tmp54i = tmp18i+tmp42i;
	tmp103r = tmp33r+tmp71r;
	tmp55i = tmp18i-tmp42i;
	tmp104r = tmp33r-tmp71r;
	tmp106r = tmp101r-tmp103r;
	tmp107r = tmp102r+tmp55i;
	tmp109r = tmp11r+tmp50r;
	tmp110r = tmp11r-tmp50r;
	tmp56i = tmp17i+tmp41i;
	tmp111r = tmp32r+tmp70r;
	tmp57i = tmp17i-tmp41i;
	tmp112r = tmp32r-tmp70r;
	tmp113r = tmp109r+tmp111r;
	tmp115r = tmp110r+tmp57i;
	tmp117r = tmp13r+tmp52r;
	tmp118r = tmp13r-tmp52r;
	tmp58i = tmp19i+tmp43i;
	tmp119r = tmp34r+tmp72r;
	tmp59i = tmp19i-tmp43i;
	tmp120r = tmp34r-tmp72r;
	tmp121r = tmp117r+tmp119r;
	tmp124r = tmp118r-tmp59i;
	tmp125r = tmp16r+tmp55r;
	tmp126r = tmp16r-tmp55r;
	tmp60i = tmp22i+tmp46i;
	tmp127r = tmp37r+tmp75r;
	tmp61i = tmp22i-tmp46i;
	tmp128r = tmp37r-tmp75r;
	tmp130r = tmp125r-tmp127r;
	tmp132r = tmp126r-tmp61i;
	Re[0] = tmp81r;
	Im[0] = tmp48i;
	Re[1] = tmp92r;
	Im[1] = tmp88r;
	Re[2] = tmp98r;
	Im[2] = -tmp52i;
	Re[3] = tmp107r;
	Im[3] = -tmp104r;
	Re[4] = tmp113r;
	Im[4] = tmp56i;
	Re[5] = tmp124r;
	Im[5] = tmp120r;
	Re[6] = tmp130r;
	Im[6] = -tmp60i;
	Re[7] = tmp83r;
	Im[7] = -tmp80r;
	Re[8] = tmp89r;
	Im[8] = tmp50i;
	Re[9] = tmp100r;
	Im[9] = tmp96r;
	Re[10] = tmp106r;
	Im[10] = -tmp54i;
	Re[11] = tmp115r;
	Im[11] = -tmp112r;
	Re[12] = tmp121r;
	Im[12] = tmp58i;
	Re[13] = tmp132r;
	Im[13] = tmp128r;
}

/*
*	Number of additions = 321
*	Number of multiplications = 144
*	Number of sign changes = 14
*	Number of assigns = 277
*	Total number of operations = 756
*/
void	MIFFTR28(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, /*tmp89r,*/ tmp90r, /*tmp91r,*/
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, /*tmp97r,*/ /*tmp98r,*/
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, /*tmp104r,*/ tmp105r,
		/*tmp106r,*/ tmp107r, tmp108r, tmp109r, tmp110r, tmp111r, /*tmp112r,*/
		tmp113r, tmp114r, /*tmp115r,*/ tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, /*tmp121r,*/ tmp122r, /*tmp123r,*/ tmp124r, tmp125r, tmp126r,
		tmp127r, tmp128r, /*tmp129r,*/ /*tmp130r,*/ tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r, /*tmp136r,*/ tmp137r, /*tmp138r,*/ tmp139r;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, /*tmp84i,*/
		tmp85i, /*tmp86i,*/ tmp87i, /*tmp88i,*/ /*tmp89i,*/ tmp90i, /*tmp91i,*/
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, /*tmp97i,*/ /*tmp98i,*/
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i, /*tmp104i,*/ tmp105i,
		/*tmp106i,*/ tmp107i, tmp108i, tmp109i, tmp110i, tmp111i, /*tmp112i,*/
		tmp113i, tmp114i, /*tmp115i,*/ tmp116i, tmp117i, tmp118i, tmp119i,
		tmp120i, /*tmp121i,*/ tmp122i, /*tmp123i,*/ tmp124i, tmp125i, tmp126i,
		tmp127i, tmp128i, /*tmp129i,*/ /*tmp130i,*/ tmp131i, tmp132i, tmp133i,
		tmp134i, tmp135i, /*tmp136i,*/ tmp137i, /*tmp138i,*/ tmp139i;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C13 =     0.22252093395631;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =     0.43388373911756;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =     0.62348980185873;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.78183148246803;	/* FFT_precisionCONST	*/
	const FFT_precision	C11 =     0.90096886790242;	/* FFT_precisionCONST	*/
	const FFT_precision	C14 =     0.97492791218182;	/* FFT_precisionCONST	*/

	tmp0i = Im[8]-Im[6];
	tmp0r = Re[8]+Re[6];
	tmp1i = Im[8]+Im[6];
	tmp1r = Re[8]-Re[6];
	tmp2i = -(Im[10]-Im[4]);
	tmp2r = Re[10]+Re[4];
	tmp3i = -(Im[10]+Im[4]);
	tmp3r = Re[10]-Re[4];
	tmp4i = -(Im[2]-Im[12]);
	tmp4r = Re[2]+Re[12];
	tmp5i = -(Im[2]+Im[12]);
	tmp5r = Re[2]-Re[12];
	tmp6i = C9*tmp0i-C11*tmp2i-C13*tmp4i+Im[0];
	tmp6r = C9*tmp0r-C11*tmp2r-C13*tmp4r+Re[0];
	tmp7i = tmp0i+tmp2i;
	tmp7r = tmp0r+tmp2r;
	tmp8i = -C11*tmp0i-C13*tmp2i+C9*tmp4i+Im[0];
	tmp8r = -C11*tmp0r-C13*tmp2r+C9*tmp4r+Re[0];
	tmp9i = tmp7i+tmp4i;
	tmp9r = tmp7r+tmp4r;
	tmp10i = -C13*tmp0i+C9*tmp2i-C11*tmp4i+Im[0];
	tmp10r = -C13*tmp0r+C9*tmp2r-C11*tmp4r+Re[0];
	tmp11i = C10*tmp1r+C12*tmp3r+C14*tmp5r;
	tmp11r = -C10*tmp1i-C12*tmp3i-C14*tmp5i;
	tmp12i = C12*tmp1r+C14*tmp3r-C10*tmp5r;
	tmp12r = -C12*tmp1i-C14*tmp3i+C10*tmp5i;
	tmp13i = C14*tmp1r-C10*tmp3r-C12*tmp5r;
	tmp13r = -C14*tmp1i+C10*tmp3i+C12*tmp5i;
	tmp14i = tmp6i+tmp11i;
	tmp14r = tmp6r+tmp11r;
	tmp15i = tmp6i-tmp11i;
	tmp15r = tmp6r-tmp11r;
	tmp16i = tmp8i+tmp12i;
	tmp16r = tmp8r+tmp12r;
	tmp17i = tmp8i-tmp12i;
	tmp17r = tmp8r-tmp12r;
	tmp18i = tmp10i+tmp13i;
	tmp18r = tmp10r+tmp13r;
	tmp19i = tmp10i-tmp13i;
	tmp19r = tmp10r-tmp13r;
	tmp20i = Im[0]+tmp9i;
	tmp20r = Re[0]+tmp9r;
	tmp21i = Im[1]+Im[13];
	tmp21r = Re[1]+Re[13];
	tmp22i = Im[1]-Im[13];
	tmp22r = Re[1]-Re[13];
	tmp23i = -(Im[3]+Im[11]);
	tmp23r = Re[3]+Re[11];
	tmp24i = -(Im[3]-Im[11]);
	tmp24r = Re[3]-Re[11];
	tmp25i = Im[9]+Im[5];
	tmp25r = Re[9]+Re[5];
	tmp26i = Im[9]-Im[5];
	tmp26r = Re[9]-Re[5];
	tmp27i = C9*tmp21i-C11*tmp23i-C13*tmp25i-Im[7];
	tmp27r = C9*tmp21r-C11*tmp23r-C13*tmp25r+Re[7];
	tmp28i = tmp21i+tmp23i;
	tmp28r = tmp21r+tmp23r;
	tmp29i = -C11*tmp21i-C13*tmp23i+C9*tmp25i-Im[7];
	tmp29r = -C11*tmp21r-C13*tmp23r+C9*tmp25r+Re[7];
	tmp30i = tmp28i+tmp25i;
	tmp30r = tmp28r+tmp25r;
	tmp31i = -C13*tmp21i+C9*tmp23i-C11*tmp25i-Im[7];
	tmp31r = -C13*tmp21r+C9*tmp23r-C11*tmp25r+Re[7];
	tmp32i = C10*tmp22r+C12*tmp24r+C14*tmp26r;
	tmp32r = -C10*tmp22i-C12*tmp24i-C14*tmp26i;
	tmp33i = C12*tmp22r+C14*tmp24r-C10*tmp26r;
	tmp33r = -C12*tmp22i-C14*tmp24i+C10*tmp26i;
	tmp34i = C14*tmp22r-C10*tmp24r-C12*tmp26r;
	tmp34r = -C14*tmp22i+C10*tmp24i+C12*tmp26i;
	tmp35i = tmp27i+tmp32i;
	tmp35r = tmp27r+tmp32r;
	tmp36i = tmp27i-tmp32i;
	tmp36r = tmp27r-tmp32r;
	tmp37i = tmp29i+tmp33i;
	tmp37r = tmp29r+tmp33r;
	tmp38i = tmp29i-tmp33i;
	tmp38r = tmp29r-tmp33r;
	tmp39i = tmp31i+tmp34i;
	tmp39r = tmp31r+tmp34r;
	tmp40i = tmp31i-tmp34i;
	tmp40r = tmp31r-tmp34r;
	tmp41i = -(Im[7]-tmp30i);
	tmp41r = Re[7]+tmp30r;
	tmp42i = -(Im[8]-Im[6]);
	tmp42r = Re[8]+Re[6];
	tmp43i = -(Im[8]+Im[6]);
	tmp43r = Re[8]-Re[6];
	tmp44i = Im[10]-Im[4];
	tmp44r = Re[10]+Re[4];
	tmp45i = Im[10]+Im[4];
	tmp45r = Re[10]-Re[4];
	tmp46i = Im[2]-Im[12];
	tmp46r = Re[2]+Re[12];
	tmp47i = Im[2]+Im[12];
	tmp47r = Re[2]-Re[12];
	tmp48i = C9*tmp42i-C11*tmp44i-C13*tmp46i-Im[0];
	tmp48r = C9*tmp42r-C11*tmp44r-C13*tmp46r+Re[0];
	tmp49i = tmp42i+tmp44i;
	tmp49r = tmp42r+tmp44r;
	tmp50i = -C11*tmp42i-C13*tmp44i+C9*tmp46i-Im[0];
	tmp50r = -C11*tmp42r-C13*tmp44r+C9*tmp46r+Re[0];
	tmp51i = tmp49i+tmp46i;
	tmp51r = tmp49r+tmp46r;
	tmp52i = -C13*tmp42i+C9*tmp44i-C11*tmp46i-Im[0];
	tmp52r = -C13*tmp42r+C9*tmp44r-C11*tmp46r+Re[0];
	tmp53i = C10*tmp43r+C12*tmp45r+C14*tmp47r;
	tmp53r = -C10*tmp43i-C12*tmp45i-C14*tmp47i;
	tmp54i = C12*tmp43r+C14*tmp45r-C10*tmp47r;
	tmp54r = -C12*tmp43i-C14*tmp45i+C10*tmp47i;
	tmp55i = C14*tmp43r-C10*tmp45r-C12*tmp47r;
	tmp55r = -C14*tmp43i+C10*tmp45i+C12*tmp47i;
	tmp56i = tmp48i+tmp53i;
	tmp56r = tmp48r+tmp53r;
	tmp57i = tmp48i-tmp53i;
	tmp57r = tmp48r-tmp53r;
	tmp58i = tmp50i+tmp54i;
	tmp58r = tmp50r+tmp54r;
	tmp59i = tmp50i-tmp54i;
	tmp59r = tmp50r-tmp54r;
	tmp60i = tmp52i+tmp55i;
	tmp60r = tmp52r+tmp55r;
	tmp61i = tmp52i-tmp55i;
	tmp61r = tmp52r-tmp55r;
	tmp62i = -(Im[0]-tmp51i);
	tmp62r = Re[0]+tmp51r;
	tmp63i = -(Im[1]+Im[13]);
	tmp63r = Re[1]+Re[13];
	tmp64i = -(Im[1]-Im[13]);
	tmp64r = Re[1]-Re[13];
	tmp65i = Im[3]+Im[11];
	tmp65r = Re[3]+Re[11];
	tmp66i = Im[3]-Im[11];
	tmp66r = Re[3]-Re[11];
	tmp67i = -(Im[9]+Im[5]);
	tmp67r = Re[9]+Re[5];
	tmp68i = -(Im[9]-Im[5]);
	tmp68r = Re[9]-Re[5];
	tmp69i = C9*tmp63i-C11*tmp65i-C13*tmp67i+Im[7];
	tmp69r = C9*tmp63r-C11*tmp65r-C13*tmp67r+Re[7];
	tmp70i = tmp63i+tmp65i;
	tmp70r = tmp63r+tmp65r;
	tmp71i = -C11*tmp63i-C13*tmp65i+C9*tmp67i+Im[7];
	tmp71r = -C11*tmp63r-C13*tmp65r+C9*tmp67r+Re[7];
	tmp72i = tmp70i+tmp67i;
	tmp72r = tmp70r+tmp67r;
	tmp73i = -C13*tmp63i+C9*tmp65i-C11*tmp67i+Im[7];
	tmp73r = -C13*tmp63r+C9*tmp65r-C11*tmp67r+Re[7];
	tmp74i = C10*tmp64r+C12*tmp66r+C14*tmp68r;
	tmp74r = -C10*tmp64i-C12*tmp66i-C14*tmp68i;
	tmp75i = C12*tmp64r+C14*tmp66r-C10*tmp68r;
	tmp75r = -C12*tmp64i-C14*tmp66i+C10*tmp68i;
	tmp76i = C14*tmp64r-C10*tmp66r-C12*tmp68r;
	tmp76r = -C14*tmp64i+C10*tmp66i+C12*tmp68i;
	tmp77i = tmp69i+tmp74i;
	tmp77r = tmp69r+tmp74r;
	tmp78i = tmp69i-tmp74i;
	tmp78r = tmp69r-tmp74r;
	tmp79i = tmp71i+tmp75i;
	tmp79r = tmp71r+tmp75r;
	tmp80i = tmp71i-tmp75i;
	tmp80r = tmp71r-tmp75r;
	tmp81i = tmp73i+tmp76i;
	tmp81r = tmp73r+tmp76r;
	tmp82i = tmp73i-tmp76i;
	tmp82r = tmp73r-tmp76r;
	tmp83i = Im[7]+tmp72i;
	tmp83r = Re[7]+tmp72r;
	tmp84r = tmp20r+tmp62r;
	tmp85i = tmp20i-tmp62i;
	tmp85r = tmp20r-tmp62r;
	tmp86r = tmp41r+tmp83r;
	tmp87i = tmp41i-tmp83i;
	tmp87r = tmp41r-tmp83r;
	tmp88r = tmp84r+tmp86r;
	tmp90i = tmp85i+tmp87r;
	tmp90r = tmp85r-tmp87i;
	tmp92i = tmp18i+tmp60i;
	tmp92r = tmp18r+tmp60r;
	tmp93i = tmp18i-tmp60i;
	tmp93r = tmp18r-tmp60r;
	tmp94i = tmp39i+tmp81i;
	tmp94r = tmp39r+tmp81r;
	tmp95i = tmp39i-tmp81i;
	tmp95r = tmp39r-tmp81r;
	tmp96i = tmp92i+tmp94i;
	tmp96r = tmp92r+tmp94r;
	tmp99i = tmp93i-tmp95r;
	tmp99r = tmp93r+tmp95i;
	tmp100i = tmp17i+tmp59i;
	tmp100r = tmp17r+tmp59r;
	tmp101i = tmp17i-tmp59i;
	tmp101r = tmp17r-tmp59r;
	tmp102i = tmp38i+tmp80i;
	tmp102r = tmp38r+tmp80r;
	tmp103i = tmp38i-tmp80i;
	tmp103r = tmp38r-tmp80r;
	tmp105i = tmp100i-tmp102i;
	tmp105r = tmp100r-tmp102r;
	tmp107i = tmp101i-tmp103r;
	tmp107r = tmp101r+tmp103i;
	tmp108i = tmp15i+tmp57i;
	tmp108r = tmp15r+tmp57r;
	tmp109i = tmp15i-tmp57i;
	tmp109r = tmp15r-tmp57r;
	tmp110i = tmp36i+tmp78i;
	tmp110r = tmp36r+tmp78r;
	tmp111i = tmp36i-tmp78i;
	tmp111r = tmp36r-tmp78r;
	tmp113i = tmp108i-tmp110i;
	tmp113r = tmp108r-tmp110r;
	tmp114i = tmp109i+tmp111r;
	tmp114r = tmp109r-tmp111i;
	tmp116i = tmp14i+tmp56i;
	tmp116r = tmp14r+tmp56r;
	tmp117i = tmp14i-tmp56i;
	tmp117r = tmp14r-tmp56r;
	tmp118i = tmp35i+tmp77i;
	tmp118r = tmp35r+tmp77r;
	tmp119i = tmp35i-tmp77i;
	tmp119r = tmp35r-tmp77r;
	tmp120i = tmp116i+tmp118i;
	tmp120r = tmp116r+tmp118r;
	tmp122i = tmp117i+tmp119r;
	tmp122r = tmp117r-tmp119i;
	tmp124i = tmp16i+tmp58i;
	tmp124r = tmp16r+tmp58r;
	tmp125i = tmp16i-tmp58i;
	tmp125r = tmp16r-tmp58r;
	tmp126i = tmp37i+tmp79i;
	tmp126r = tmp37r+tmp79r;
	tmp127i = tmp37i-tmp79i;
	tmp127r = tmp37r-tmp79r;
	tmp128i = tmp124i+tmp126i;
	tmp128r = tmp124r+tmp126r;
	tmp131i = tmp125i-tmp127r;
	tmp131r = tmp125r+tmp127i;
	tmp132i = tmp19i+tmp61i;
	tmp132r = tmp19r+tmp61r;
	tmp133i = tmp19i-tmp61i;
	tmp133r = tmp19r-tmp61r;
	tmp134i = tmp40i+tmp82i;
	tmp134r = tmp40r+tmp82r;
	tmp135i = tmp40i-tmp82i;
	tmp135r = tmp40r-tmp82r;
	tmp137i = tmp132i-tmp134i;
	tmp137r = tmp132r-tmp134r;
	tmp139i = tmp133i-tmp135r;
	tmp139r = tmp133r+tmp135i;
	Re[0] = tmp88r;
	Im[0] = C0;
	Re[1] = tmp99r;
	Im[1] = tmp99i;
	Re[2] = tmp105r;
	Im[2] = tmp105i;
	Re[3] = tmp114r;
	Im[3] = tmp114i;
	Re[4] = tmp120r;
	Im[4] = tmp120i;
	Re[5] = tmp131r;
	Im[5] = tmp131i;
	Re[6] = tmp137r;
	Im[6] = tmp137i;
	Re[7] = tmp90r;
	Im[7] = tmp90i;
	Re[8] = tmp96r;
	Im[8] = tmp96i;
	Re[9] = tmp107r;
	Im[9] = tmp107i;
	Re[10] = tmp113r;
	Im[10] = tmp113i;
	Re[11] = tmp122r;
	Im[11] = tmp122i;
	Re[12] = tmp128r;
	Im[12] = tmp128i;
	Re[13] = tmp139r;
	Im[13] = tmp139i;
}

/*
*	Number of additions = 280
*	Number of multiplications = 108
*	Number of sign changes = 52
*	Number of assigns = 378
*	Total number of operations = 818
*/
void	MFFTR30(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, /*tmp105r,*/
		/*tmp106r,*/ tmp107r, tmp108r, /*tmp109r,*/ tmp110r, tmp111r, tmp112r,
		tmp113r, tmp114r, tmp115r, tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, tmp121r, tmp122r, tmp123r, tmp124r, /*tmp125r,*/ /*tmp126r,*/
		tmp127r, /*tmp128r,*/ tmp129r, tmp130r, tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r, tmp136r, tmp137r, tmp138r, tmp139r, tmp140r,
		tmp141r, tmp142r, tmp143r, tmp144r, /*tmp145r,*/ tmp146r, /*tmp147r,*/
		/*tmp148r,*/ tmp149r, tmp150r, tmp151r, tmp152r, tmp153r, tmp154r,
		tmp155r, tmp156r, tmp157r, tmp158r, tmp159r, tmp160r, tmp161r,
		tmp162r, tmp163r, /*tmp164r,*/ tmp165r, tmp166r, /*tmp167r,*/ /*tmp168r,*/
		tmp169r, tmp170r, tmp171r, tmp172r, tmp173r, tmp174r, tmp175r,
		tmp176r, tmp177r, tmp178r, tmp179r, tmp180r, tmp181r, tmp182r,
		tmp183r, /*tmp184r,*/ tmp185r, tmp186r, /*tmp187r,*/ tmp188r /*tmp189r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, tmp100i, tmp101i, tmp102i, /*tmp103i,*/ /*tmp104i,*/ tmp105i,
		tmp106i, /*tmp107i,*/ tmp108i, tmp109i, tmp110i, tmp111i, tmp112i,
		tmp113i, tmp114i, tmp115i, tmp116i, tmp117i, tmp118i, tmp119i,
		tmp120i, tmp121i, tmp122i, /*tmp123i,*/ /*tmp124i,*/ tmp125i, /*tmp126i,*/
		tmp127i, tmp128i, tmp129i, tmp130i, tmp131i, tmp132i, tmp133i,
		tmp134i, tmp135i, tmp136i, tmp137i, tmp138i, tmp139i, tmp140i,
		tmp141i, tmp142i, /*tmp143i,*/ tmp144i, /*tmp145i,*/ /*tmp146i,*/ tmp147i,
		tmp148i, tmp149i, tmp150i, tmp151i, tmp152i, tmp153i, tmp154i,
		tmp155i, tmp156i, tmp157i, tmp158i, tmp159i, tmp160i, tmp161i,
		/*tmp162i,*/ tmp163i, tmp164i, /*tmp165i,*/ /*tmp166i,*/ tmp167i, tmp168i,
		tmp169i, tmp170i, tmp171i, tmp172i, tmp173i, tmp174i, tmp175i,
		tmp176i, tmp177i, tmp178i, tmp179i, tmp180i, tmp181i, /*tmp182i,*/
		tmp183i, tmp184i, /*tmp185i,*/ tmp186i /*tmp187i*/;

	const FFT_precision	C1 =                 0.25;	/* FFT_precisionCONST	*/
	const FFT_precision	C13 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.55901699437495;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.58778525229247;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.95105651629515;	/* FFT_precisionCONST	*/
	const FFT_precision	C12 =                    2;	/* INTEGERCONST	*/

	tmp0r = C12*Re[6];
	tmp0i = C12*Im[6];
	tmp1r = C12*Re[12];
	tmp1i = C12*Im[12];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = -C1*tmp2r;
	tmp5r = C15*tmp3r;
	tmp6r = tmp4r+Re[0];
	tmp7r = tmp6r+tmp5r;
	tmp8r = tmp6r-tmp5r;
	tmp9r = C8*tmp0i+C10*tmp1i;
	tmp10r = C10*tmp0i-C8*tmp1i;
	tmp11r = tmp7r+tmp9r;
	tmp12r = tmp7r-tmp9r;
	tmp13r = tmp8r+tmp10r;
	tmp14r = tmp8r-tmp10r;
	tmp15r = Re[0]+tmp2r;
	tmp2i = Im[1]+Im[11];
	tmp16r = Re[1]-Re[11];
	tmp3i = Im[1]-Im[11];
	tmp17r = Re[1]+Re[11];
	tmp4i = Im[7]+Im[13];
	tmp18r = Re[7]+Re[13];
	tmp5i = Im[7]-Im[13];
	tmp19r = Re[7]-Re[13];
	tmp6i = tmp2i+tmp4i;
	tmp20r = tmp16r+tmp18r;
	tmp7i = tmp2i-tmp4i;
	tmp21r = tmp16r-tmp18r;
	tmp8i = -C1*tmp6i;
	tmp22r = -C1*tmp20r;
	tmp9i = C15*tmp7i;
	tmp23r = C15*tmp21r;
	tmp10i = tmp8i+Im[5];
	tmp24r = tmp22r-Re[5];
	tmp11i = tmp10i+tmp9i;
	tmp25r = tmp24r+tmp23r;
	tmp12i = tmp10i-tmp9i;
	tmp26r = tmp24r-tmp23r;
	tmp13i = -C8*tmp17r-C10*tmp19r;
	tmp27r = C8*tmp3i+C10*tmp5i;
	tmp14i = -C10*tmp17r+C8*tmp19r;
	tmp28r = C10*tmp3i-C8*tmp5i;
	tmp15i = tmp11i+tmp13i;
	tmp29r = tmp25r+tmp27r;
	tmp16i = tmp11i-tmp13i;
	tmp30r = tmp25r-tmp27r;
	tmp17i = tmp12i+tmp14i;
	tmp31r = tmp26r+tmp28r;
	tmp18i = tmp12i-tmp14i;
	tmp32r = tmp26r-tmp28r;
	tmp19i = Im[5]+tmp6i;
	tmp33r = -(Re[5]-tmp20r);
	tmp20i = -(Im[4]-Im[14]);
	tmp34r = Re[4]+Re[14];
	tmp21i = -(Im[4]+Im[14]);
	tmp35r = Re[4]-Re[14];
	tmp22i = Im[2]+Im[8];
	tmp36r = Re[2]+Re[8];
	tmp23i = Im[2]-Im[8];
	tmp37r = Re[2]-Re[8];
	tmp24i = tmp20i+tmp22i;
	tmp38r = tmp34r+tmp36r;
	tmp25i = tmp20i-tmp22i;
	tmp39r = tmp34r-tmp36r;
	tmp26i = -C1*tmp24i;
	tmp40r = -C1*tmp38r;
	tmp27i = C15*tmp25i;
	tmp41r = C15*tmp39r;
	tmp28i = tmp26i-Im[10];
	tmp42r = tmp40r+Re[10];
	tmp29i = tmp28i+tmp27i;
	tmp43r = tmp42r+tmp41r;
	tmp30i = tmp28i-tmp27i;
	tmp44r = tmp42r-tmp41r;
	tmp31i = -C8*tmp35r-C10*tmp37r;
	tmp45r = C8*tmp21i+C10*tmp23i;
	tmp32i = -C10*tmp35r+C8*tmp37r;
	tmp46r = C10*tmp21i-C8*tmp23i;
	tmp33i = tmp29i+tmp31i;
	tmp47r = tmp43r+tmp45r;
	tmp34i = tmp29i-tmp31i;
	tmp48r = tmp43r-tmp45r;
	tmp35i = tmp30i+tmp32i;
	tmp49r = tmp44r+tmp46r;
	tmp36i = tmp30i-tmp32i;
	tmp50r = tmp44r-tmp46r;
	tmp37i = -(Im[10]-tmp24i);
	tmp51r = Re[10]+tmp38r;
	tmp38i = C12*Im[9];
	tmp52r = -C12*Re[9];
	tmp39i = C12*Im[3];
	tmp53r = -C12*Re[3];
	tmp40i = tmp38i+tmp39i;
	tmp41i = tmp38i-tmp39i;
	tmp42i = -C1*tmp40i;
	tmp43i = C15*tmp41i;
	tmp44i = tmp42i+tmp43i;
	tmp45i = tmp42i-tmp43i;
	tmp46i = -C8*tmp52r-C10*tmp53r;
	tmp47i = -C10*tmp52r+C8*tmp53r;
	tmp48i = tmp44i+tmp46i;
	tmp49i = tmp44i-tmp46i;
	tmp50i = tmp45i+tmp47i;
	tmp51i = tmp45i-tmp47i;
	tmp52i = -(Im[14]-Im[4]);
	tmp54r = Re[14]+Re[4];
	tmp53i = -(Im[14]+Im[4]);
	tmp55r = Re[14]-Re[4];
	tmp54i = -(Im[8]+Im[2]);
	tmp56r = Re[8]+Re[2];
	tmp55i = -(Im[8]-Im[2]);
	tmp57r = Re[8]-Re[2];
	tmp56i = tmp52i+tmp54i;
	tmp58r = tmp54r+tmp56r;
	tmp57i = tmp52i-tmp54i;
	tmp59r = tmp54r-tmp56r;
	tmp58i = -C1*tmp56i;
	tmp60r = -C1*tmp58r;
	tmp59i = C15*tmp57i;
	tmp61r = C15*tmp59r;
	tmp60i = tmp58i+Im[10];
	tmp62r = tmp60r+Re[10];
	tmp61i = tmp60i+tmp59i;
	tmp63r = tmp62r+tmp61r;
	tmp62i = tmp60i-tmp59i;
	tmp64r = tmp62r-tmp61r;
	tmp63i = -C8*tmp55r-C10*tmp57r;
	tmp65r = C8*tmp53i+C10*tmp55i;
	tmp64i = -C10*tmp55r+C8*tmp57r;
	tmp66r = C10*tmp53i-C8*tmp55i;
	tmp65i = tmp61i+tmp63i;
	tmp67r = tmp63r+tmp65r;
	tmp66i = tmp61i-tmp63i;
	tmp68r = tmp63r-tmp65r;
	tmp67i = tmp62i+tmp64i;
	tmp69r = tmp64r+tmp66r;
	tmp68i = tmp62i-tmp64i;
	tmp70r = tmp64r-tmp66r;
	tmp69i = Im[10]+tmp56i;
	tmp71r = Re[10]+tmp58r;
	tmp70i = Im[11]+Im[1];
	tmp72r = Re[11]-Re[1];
	tmp71i = Im[11]-Im[1];
	tmp73r = Re[11]+Re[1];
	tmp72i = Im[13]+Im[7];
	tmp74r = -(Re[13]+Re[7]);
	tmp73i = Im[13]-Im[7];
	tmp75r = -(Re[13]-Re[7]);
	tmp74i = tmp70i+tmp72i;
	tmp76r = tmp72r+tmp74r;
	tmp75i = tmp70i-tmp72i;
	tmp77r = tmp72r-tmp74r;
	tmp76i = -C1*tmp74i;
	tmp78r = -C1*tmp76r;
	tmp77i = C15*tmp75i;
	tmp79r = C15*tmp77r;
	tmp78i = tmp76i+Im[5];
	tmp80r = tmp78r+Re[5];
	tmp79i = tmp78i+tmp77i;
	tmp81r = tmp80r+tmp79r;
	tmp80i = tmp78i-tmp77i;
	tmp82r = tmp80r-tmp79r;
	tmp81i = -C8*tmp73r-C10*tmp75r;
	tmp83r = C8*tmp71i+C10*tmp73i;
	tmp82i = -C10*tmp73r+C8*tmp75r;
	tmp84r = C10*tmp71i-C8*tmp73i;
	tmp83i = tmp79i+tmp81i;
	tmp85r = tmp81r+tmp83r;
	tmp84i = tmp79i-tmp81i;
	tmp86r = tmp81r-tmp83r;
	tmp85i = tmp80i+tmp82i;
	tmp87r = tmp82r+tmp84r;
	tmp86i = tmp80i-tmp82i;
	tmp88r = tmp82r-tmp84r;
	tmp87i = Im[5]+tmp74i;
	tmp89r = Re[5]+tmp76r;
	tmp88i = tmp69i+tmp37i;
	tmp90r = tmp71r+tmp51r;
	tmp89i = tmp69i-tmp37i;
	tmp91r = tmp71r-tmp51r;
	tmp90i = -C13*tmp88i;
	tmp92r = -C13*tmp90r;
	tmp91i = -C17*tmp91r;
	tmp93r = C17*tmp89i;
	tmp94r = tmp92r+tmp15r;
	tmp92i = tmp90i+tmp91i;
	tmp95r = tmp94r+tmp93r;
	tmp93i = tmp90i-tmp91i;
	tmp96r = tmp94r-tmp93r;
	tmp97r = tmp15r+tmp90r;
	tmp94i = tmp19i+tmp87i;
	tmp98r = tmp33r+tmp89r;
	tmp95i = tmp19i-tmp87i;
	tmp99r = tmp33r-tmp89r;
	tmp96i = -C13*tmp94i;
	tmp100r = -C13*tmp98r;
	tmp97i = -C17*tmp99r;
	tmp101r = C17*tmp95i;
	tmp98i = tmp96i+tmp40i;
	tmp99i = tmp98i+tmp97i;
	tmp102r = tmp100r+tmp101r;
	tmp100i = tmp98i-tmp97i;
	tmp103r = tmp100r-tmp101r;
	tmp101i = tmp40i+tmp94i;
	tmp102i = tmp88i+tmp101i;
	tmp104r = tmp97r+tmp98r;
	tmp105i = tmp93i-tmp100i;
	tmp107r = tmp96r-tmp103r;
	tmp106i = tmp92i+tmp99i;
	tmp108r = tmp95r+tmp102r;
	tmp108i = tmp65i+tmp33i;
	tmp110r = tmp67r+tmp47r;
	tmp109i = tmp65i-tmp33i;
	tmp111r = tmp67r-tmp47r;
	tmp110i = -C13*tmp108i;
	tmp112r = -C13*tmp110r;
	tmp111i = -C17*tmp111r;
	tmp113r = C17*tmp109i;
	tmp114r = tmp112r+tmp11r;
	tmp112i = tmp110i+tmp111i;
	tmp115r = tmp114r+tmp113r;
	tmp113i = tmp110i-tmp111i;
	tmp116r = tmp114r-tmp113r;
	tmp117r = tmp11r+tmp110r;
	tmp114i = tmp15i+tmp83i;
	tmp118r = tmp29r+tmp85r;
	tmp115i = tmp15i-tmp83i;
	tmp119r = tmp29r-tmp85r;
	tmp116i = -C13*tmp114i;
	tmp120r = -C13*tmp118r;
	tmp117i = -C17*tmp119r;
	tmp121r = C17*tmp115i;
	tmp118i = tmp116i+tmp48i;
	tmp119i = tmp118i+tmp117i;
	tmp122r = tmp120r+tmp121r;
	tmp120i = tmp118i-tmp117i;
	tmp123r = tmp120r-tmp121r;
	tmp121i = tmp48i+tmp114i;
	tmp122i = tmp108i+tmp121i;
	tmp124r = tmp117r+tmp118r;
	tmp125i = tmp113i-tmp120i;
	tmp127r = tmp116r-tmp123r;
	tmp127i = tmp112i-tmp119i;
	tmp129r = tmp115r-tmp122r;
	tmp128i = tmp67i+tmp35i;
	tmp130r = tmp69r+tmp49r;
	tmp129i = tmp67i-tmp35i;
	tmp131r = tmp69r-tmp49r;
	tmp130i = -C13*tmp128i;
	tmp132r = -C13*tmp130r;
	tmp131i = -C17*tmp131r;
	tmp133r = C17*tmp129i;
	tmp134r = tmp132r+tmp13r;
	tmp132i = tmp130i+tmp131i;
	tmp135r = tmp134r+tmp133r;
	tmp133i = tmp130i-tmp131i;
	tmp136r = tmp134r-tmp133r;
	tmp137r = tmp13r+tmp130r;
	tmp134i = tmp17i+tmp85i;
	tmp138r = tmp31r+tmp87r;
	tmp135i = tmp17i-tmp85i;
	tmp139r = tmp31r-tmp87r;
	tmp136i = -C13*tmp134i;
	tmp140r = -C13*tmp138r;
	tmp137i = -C17*tmp139r;
	tmp141r = C17*tmp135i;
	tmp138i = tmp136i+tmp50i;
	tmp139i = tmp138i+tmp137i;
	tmp142r = tmp140r+tmp141r;
	tmp140i = tmp138i-tmp137i;
	tmp143r = tmp140r-tmp141r;
	tmp141i = tmp50i+tmp134i;
	tmp142i = tmp128i+tmp141i;
	tmp144r = tmp137r+tmp138r;
	tmp144i = tmp133i+tmp140i;
	tmp146r = tmp136r+tmp143r;
	tmp147i = tmp132i-tmp139i;
	tmp149r = tmp135r-tmp142r;
	tmp148i = tmp68i+tmp36i;
	tmp150r = tmp70r+tmp50r;
	tmp149i = tmp68i-tmp36i;
	tmp151r = tmp70r-tmp50r;
	tmp150i = -C13*tmp148i;
	tmp152r = -C13*tmp150r;
	tmp151i = -C17*tmp151r;
	tmp153r = C17*tmp149i;
	tmp154r = tmp152r+tmp14r;
	tmp152i = tmp150i+tmp151i;
	tmp155r = tmp154r+tmp153r;
	tmp153i = tmp150i-tmp151i;
	tmp156r = tmp154r-tmp153r;
	tmp157r = tmp14r+tmp150r;
	tmp154i = tmp18i+tmp86i;
	tmp158r = tmp32r+tmp88r;
	tmp155i = tmp18i-tmp86i;
	tmp159r = tmp32r-tmp88r;
	tmp156i = -C13*tmp154i;
	tmp160r = -C13*tmp158r;
	tmp157i = -C17*tmp159r;
	tmp161r = C17*tmp155i;
	tmp158i = tmp156i+tmp51i;
	tmp159i = tmp158i+tmp157i;
	tmp162r = tmp160r+tmp161r;
	tmp160i = tmp158i-tmp157i;
	tmp163r = tmp160r-tmp161r;
	tmp161i = tmp51i+tmp154i;
	tmp163i = tmp148i-tmp161i;
	tmp165r = tmp157r-tmp158r;
	tmp164i = tmp153i+tmp160i;
	tmp166r = tmp156r+tmp163r;
	tmp167i = tmp152i-tmp159i;
	tmp169r = tmp155r-tmp162r;
	tmp168i = tmp66i+tmp34i;
	tmp170r = tmp68r+tmp48r;
	tmp169i = tmp66i-tmp34i;
	tmp171r = tmp68r-tmp48r;
	tmp170i = -C13*tmp168i;
	tmp172r = -C13*tmp170r;
	tmp171i = -C17*tmp171r;
	tmp173r = C17*tmp169i;
	tmp174r = tmp172r+tmp12r;
	tmp172i = tmp170i+tmp171i;
	tmp175r = tmp174r+tmp173r;
	tmp173i = tmp170i-tmp171i;
	tmp176r = tmp174r-tmp173r;
	tmp177r = tmp12r+tmp170r;
	tmp174i = tmp16i+tmp84i;
	tmp178r = tmp30r+tmp86r;
	tmp175i = tmp16i-tmp84i;
	tmp179r = tmp30r-tmp86r;
	tmp176i = -C13*tmp174i;
	tmp180r = -C13*tmp178r;
	tmp177i = -C17*tmp179r;
	tmp181r = C17*tmp175i;
	tmp178i = tmp176i+tmp49i;
	tmp179i = tmp178i+tmp177i;
	tmp182r = tmp180r+tmp181r;
	tmp180i = tmp178i-tmp177i;
	tmp183r = tmp180r-tmp181r;
	tmp181i = tmp49i+tmp174i;
	tmp183i = tmp168i-tmp181i;
	tmp185r = tmp177r-tmp178r;
	tmp184i = tmp173i+tmp180i;
	tmp186r = tmp176r+tmp183r;
	tmp186i = tmp172i+tmp179i;
	tmp188r = tmp175r+tmp182r;
	Re[0] = tmp104r;
	Im[0] = tmp102i;
	Re[1] = tmp129r;
	Im[1] = tmp127i;
	Re[2] = tmp146r;
	Im[2] = tmp144i;
	Re[3] = tmp165r;
	Im[3] = tmp163i;
	Re[4] = tmp188r;
	Im[4] = tmp186i;
	Re[5] = tmp107r;
	Im[5] = tmp105i;
	Re[6] = tmp124r;
	Im[6] = tmp122i;
	Re[7] = tmp149r;
	Im[7] = tmp147i;
	Re[8] = tmp166r;
	Im[8] = tmp164i;
	Re[9] = tmp185r;
	Im[9] = tmp183i;
	Re[10] = tmp108r;
	Im[10] = tmp106i;
	Re[11] = tmp127r;
	Im[11] = tmp125i;
	Re[12] = tmp144r;
	Im[12] = tmp142i;
	Re[13] = tmp169r;
	Im[13] = tmp167i;
	Re[14] = tmp186r;
	Im[14] = tmp184i;
}

/*
*	Number of additions = 339
*	Number of multiplications = 112
*	Number of sign changes = 57
*	Number of assigns = 433
*	Total number of operations = 941
*/
void	MIFFTR30(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, tmp105r,
		tmp106r, tmp107r, tmp108r, tmp109r, tmp110r, tmp111r, tmp112r,
		tmp113r, tmp114r, tmp115r, tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, tmp121r, tmp122r, tmp123r, tmp124r, /*tmp125r,*/ /*tmp126r,*/
		tmp127r, tmp128r, /*tmp129r,*/ tmp130r, tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r, tmp136r, tmp137r, tmp138r, tmp139r, tmp140r,
		tmp141r, tmp142r, tmp143r, tmp144r, tmp145r, tmp146r, /*tmp147r,*/
		/*tmp148r,*/ tmp149r, /*tmp150r,*/ tmp151r, tmp152r, tmp153r, tmp154r,
		tmp155r, tmp156r, tmp157r, tmp158r, tmp159r, tmp160r, tmp161r,
		tmp162r, tmp163r, tmp164r, tmp165r, tmp166r, tmp167r, tmp168r,
		/*tmp169r,*/ tmp170r, /*tmp171r,*/ /*tmp172r,*/ tmp173r, tmp174r, tmp175r,
		tmp176r, tmp177r, tmp178r, tmp179r, tmp180r, tmp181r, tmp182r,
		tmp183r, tmp184r, tmp185r, tmp186r, tmp187r, tmp188r, tmp189r,
		/*tmp190r,*/ tmp191r, tmp192r, /*tmp193r,*/ /*tmp194r,*/ tmp195r, tmp196r,
		tmp197r, tmp198r, tmp199r, tmp200r, tmp201r, tmp202r, tmp203r,
		tmp204r, tmp205r, tmp206r, tmp207r, tmp208r, tmp209r, tmp210r,
		tmp211r, /*tmp212r,*/ tmp213r, tmp214r, /*tmp215r,*/ tmp216r /*tmp217r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i, tmp104i, tmp105i,
		tmp106i, tmp107i, tmp108i, tmp109i, tmp110i, tmp111i, tmp112i,
		tmp113i, tmp114i, /*tmp115i,*/ tmp116i, tmp117i, tmp118i, tmp119i,
		tmp120i, tmp121i, tmp122i, /*tmp123i,*/ /*tmp124i,*/ /*tmp125i,*/ /*tmp126i,*/
		tmp127i, tmp128i, /*tmp129i,*/ tmp130i, tmp131i, tmp132i, tmp133i,
		tmp134i, tmp135i, tmp136i, tmp137i, tmp138i, tmp139i, tmp140i,
		tmp141i, tmp142i, tmp143i, tmp144i, tmp145i, tmp146i, /*tmp147i,*/
		/*tmp148i,*/ tmp149i, /*tmp150i,*/ tmp151i, tmp152i, tmp153i, tmp154i,
		tmp155i, tmp156i, tmp157i, tmp158i, tmp159i, tmp160i, tmp161i,
		tmp162i, tmp163i, tmp164i, tmp165i, tmp166i, tmp167i, tmp168i,
		/*tmp169i,*/ tmp170i, /*tmp171i,*/ /*tmp172i,*/ tmp173i, tmp174i, tmp175i,
		tmp176i, tmp177i, tmp178i, tmp179i, tmp180i, tmp181i, tmp182i,
		tmp183i, tmp184i, tmp185i, tmp186i, tmp187i, tmp188i, tmp189i,
		/*tmp190i,*/ tmp191i, tmp192i, /*tmp193i,*/ /*tmp194i,*/ tmp195i, tmp196i,
		tmp197i, tmp198i, tmp199i, tmp200i, tmp201i, tmp202i, tmp203i,
		tmp204i, tmp205i, tmp206i, tmp207i, tmp208i, tmp209i, tmp210i,
		tmp211i, /*tmp212i,*/ tmp213i, tmp214i, /*tmp215i,*/ tmp216i /*tmp217i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C1 =                 0.25;	/* FFT_precisionCONST	*/
	const FFT_precision	C13 =                  0.5;	/* FFT_precisionCONST	*/
	const FFT_precision	C15 =     0.55901699437495;	/* FFT_precisionCONST	*/
	const FFT_precision	C10 =     0.58778525229247;	/* FFT_precisionCONST	*/
	const FFT_precision	C17 =     0.86602540378444;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.95105651629515;	/* FFT_precisionCONST	*/

	tmp0i = Im[6]-Im[9];
	tmp0r = Re[6]+Re[9];
	tmp1i = Im[6]+Im[9];
	tmp1r = Re[6]-Re[9];
	tmp2i = Im[12]-Im[3];
	tmp2r = Re[12]+Re[3];
	tmp3i = Im[12]+Im[3];
	tmp3r = Re[12]-Re[3];
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = -C1*tmp4i;
	tmp6r = -C1*tmp4r;
	tmp7i = C15*tmp5i;
	tmp7r = C15*tmp5r;
	tmp8i = tmp6i+Im[0];
	tmp8r = tmp6r+Re[0];
	tmp9i = tmp8i+tmp7i;
	tmp9r = tmp8r+tmp7r;
	tmp10i = tmp8i-tmp7i;
	tmp10r = tmp8r-tmp7r;
	tmp11i = C8*tmp1r+C10*tmp3r;
	tmp11r = -C8*tmp1i-C10*tmp3i;
	tmp12i = C10*tmp1r-C8*tmp3r;
	tmp12r = -C10*tmp1i+C8*tmp3i;
	tmp13i = tmp9i+tmp11i;
	tmp13r = tmp9r+tmp11r;
	tmp14i = tmp9i-tmp11i;
	tmp14r = tmp9r-tmp11r;
	tmp15i = tmp10i+tmp12i;
	tmp15r = tmp10r+tmp12r;
	tmp16i = tmp10i-tmp12i;
	tmp16r = tmp10r-tmp12r;
	tmp17i = Im[0]+tmp4i;
	tmp17r = Re[0]+tmp4r;
	tmp18i = Im[1]-Im[4];
	tmp18r = Re[1]+Re[4];
	tmp19i = Im[1]+Im[4];
	tmp19r = Re[1]-Re[4];
	tmp20i = Im[7]+Im[13];
	tmp20r = Re[7]+Re[13];
	tmp21i = Im[7]-Im[13];
	tmp21r = Re[7]-Re[13];
	tmp22i = tmp18i+tmp20i;
	tmp22r = tmp18r+tmp20r;
	tmp23i = tmp18i-tmp20i;
	tmp23r = tmp18r-tmp20r;
	tmp24i = -C1*tmp22i;
	tmp24r = -C1*tmp22r;
	tmp25i = C15*tmp23i;
	tmp25r = C15*tmp23r;
	tmp26i = tmp24i-Im[10];
	tmp26r = tmp24r+Re[10];
	tmp27i = tmp26i+tmp25i;
	tmp27r = tmp26r+tmp25r;
	tmp28i = tmp26i-tmp25i;
	tmp28r = tmp26r-tmp25r;
	tmp29i = C8*tmp19r+C10*tmp21r;
	tmp29r = -C8*tmp19i-C10*tmp21i;
	tmp30i = C10*tmp19r-C8*tmp21r;
	tmp30r = -C10*tmp19i+C8*tmp21i;
	tmp31i = tmp27i+tmp29i;
	tmp31r = tmp27r+tmp29r;
	tmp32i = tmp27i-tmp29i;
	tmp32r = tmp27r-tmp29r;
	tmp33i = tmp28i+tmp30i;
	tmp33r = tmp28r+tmp30r;
	tmp34i = tmp28i-tmp30i;
	tmp34r = tmp28r-tmp30r;
	tmp35i = -(Im[10]-tmp22i);
	tmp35r = Re[10]+tmp22r;
	tmp36i = -(Im[11]-Im[14]);
	tmp36r = Re[11]+Re[14];
	tmp37i = -(Im[11]+Im[14]);
	tmp37r = Re[11]-Re[14];
	tmp38i = Im[2]+Im[8];
	tmp38r = Re[2]+Re[8];
	tmp39i = Im[2]-Im[8];
	tmp39r = Re[2]-Re[8];
	tmp40i = tmp36i+tmp38i;
	tmp40r = tmp36r+tmp38r;
	tmp41i = tmp36i-tmp38i;
	tmp41r = tmp36r-tmp38r;
	tmp42i = -C1*tmp40i;
	tmp42r = -C1*tmp40r;
	tmp43i = C15*tmp41i;
	tmp43r = C15*tmp41r;
	tmp44i = tmp42i-Im[5];
	tmp44r = tmp42r+Re[5];
	tmp45i = tmp44i+tmp43i;
	tmp45r = tmp44r+tmp43r;
	tmp46i = tmp44i-tmp43i;
	tmp46r = tmp44r-tmp43r;
	tmp47i = C8*tmp37r+C10*tmp39r;
	tmp47r = -C8*tmp37i-C10*tmp39i;
	tmp48i = C10*tmp37r-C8*tmp39r;
	tmp48r = -C10*tmp37i+C8*tmp39i;
	tmp49i = tmp45i+tmp47i;
	tmp49r = tmp45r+tmp47r;
	tmp50i = tmp45i-tmp47i;
	tmp50r = tmp45r-tmp47r;
	tmp51i = tmp46i+tmp48i;
	tmp51r = tmp46r+tmp48r;
	tmp52i = tmp46i-tmp48i;
	tmp52r = tmp46r-tmp48r;
	tmp53i = -(Im[5]-tmp40i);
	tmp53r = Re[5]+tmp40r;
	tmp54i = -(Im[6]-Im[9]);
	tmp54r = Re[6]+Re[9];
	tmp55i = -(Im[6]+Im[9]);
	tmp55r = Re[6]-Re[9];
	tmp56i = -(Im[12]-Im[3]);
	tmp56r = Re[12]+Re[3];
	tmp57i = -(Im[12]+Im[3]);
	tmp57r = Re[12]-Re[3];
	tmp58i = tmp54i+tmp56i;
	tmp58r = tmp54r+tmp56r;
	tmp59i = tmp54i-tmp56i;
	tmp59r = tmp54r-tmp56r;
	tmp60i = -C1*tmp58i;
	tmp60r = -C1*tmp58r;
	tmp61i = C15*tmp59i;
	tmp61r = C15*tmp59r;
	tmp62i = tmp60i-Im[0];
	tmp62r = tmp60r+Re[0];
	tmp63i = tmp62i+tmp61i;
	tmp63r = tmp62r+tmp61r;
	tmp64i = tmp62i-tmp61i;
	tmp64r = tmp62r-tmp61r;
	tmp65i = C8*tmp55r+C10*tmp57r;
	tmp65r = -C8*tmp55i-C10*tmp57i;
	tmp66i = C10*tmp55r-C8*tmp57r;
	tmp66r = -C10*tmp55i+C8*tmp57i;
	tmp67i = tmp63i+tmp65i;
	tmp67r = tmp63r+tmp65r;
	tmp68i = tmp63i-tmp65i;
	tmp68r = tmp63r-tmp65r;
	tmp69i = tmp64i+tmp66i;
	tmp69r = tmp64r+tmp66r;
	tmp70i = tmp64i-tmp66i;
	tmp70r = tmp64r-tmp66r;
	tmp71i = -(Im[0]-tmp58i);
	tmp71r = Re[0]+tmp58r;
	tmp72i = -(Im[1]-Im[4]);
	tmp72r = Re[1]+Re[4];
	tmp73i = -(Im[1]+Im[4]);
	tmp73r = Re[1]-Re[4];
	tmp74i = -(Im[7]+Im[13]);
	tmp74r = Re[7]+Re[13];
	tmp75i = -(Im[7]-Im[13]);
	tmp75r = Re[7]-Re[13];
	tmp76i = tmp72i+tmp74i;
	tmp76r = tmp72r+tmp74r;
	tmp77i = tmp72i-tmp74i;
	tmp77r = tmp72r-tmp74r;
	tmp78i = -C1*tmp76i;
	tmp78r = -C1*tmp76r;
	tmp79i = C15*tmp77i;
	tmp79r = C15*tmp77r;
	tmp80i = tmp78i+Im[10];
	tmp80r = tmp78r+Re[10];
	tmp81i = tmp80i+tmp79i;
	tmp81r = tmp80r+tmp79r;
	tmp82i = tmp80i-tmp79i;
	tmp82r = tmp80r-tmp79r;
	tmp83i = C8*tmp73r+C10*tmp75r;
	tmp83r = -C8*tmp73i-C10*tmp75i;
	tmp84i = C10*tmp73r-C8*tmp75r;
	tmp84r = -C10*tmp73i+C8*tmp75i;
	tmp85i = tmp81i+tmp83i;
	tmp85r = tmp81r+tmp83r;
	tmp86i = tmp81i-tmp83i;
	tmp86r = tmp81r-tmp83r;
	tmp87i = tmp82i+tmp84i;
	tmp87r = tmp82r+tmp84r;
	tmp88i = tmp82i-tmp84i;
	tmp88r = tmp82r-tmp84r;
	tmp89i = Im[10]+tmp76i;
	tmp89r = Re[10]+tmp76r;
	tmp90i = Im[11]-Im[14];
	tmp90r = Re[11]+Re[14];
	tmp91i = Im[11]+Im[14];
	tmp91r = Re[11]-Re[14];
	tmp92i = -(Im[2]+Im[8]);
	tmp92r = Re[2]+Re[8];
	tmp93i = -(Im[2]-Im[8]);
	tmp93r = Re[2]-Re[8];
	tmp94i = tmp90i+tmp92i;
	tmp94r = tmp90r+tmp92r;
	tmp95i = tmp90i-tmp92i;
	tmp95r = tmp90r-tmp92r;
	tmp96i = -C1*tmp94i;
	tmp96r = -C1*tmp94r;
	tmp97i = C15*tmp95i;
	tmp97r = C15*tmp95r;
	tmp98i = tmp96i+Im[5];
	tmp98r = tmp96r+Re[5];
	tmp99i = tmp98i+tmp97i;
	tmp99r = tmp98r+tmp97r;
	tmp100i = tmp98i-tmp97i;
	tmp100r = tmp98r-tmp97r;
	tmp101i = C8*tmp91r+C10*tmp93r;
	tmp101r = -C8*tmp91i-C10*tmp93i;
	tmp102i = C10*tmp91r-C8*tmp93r;
	tmp102r = -C10*tmp91i+C8*tmp93i;
	tmp103i = tmp99i+tmp101i;
	tmp103r = tmp99r+tmp101r;
	tmp104i = tmp99i-tmp101i;
	tmp104r = tmp99r-tmp101r;
	tmp105i = tmp100i+tmp102i;
	tmp105r = tmp100r+tmp102r;
	tmp106i = tmp100i-tmp102i;
	tmp106r = tmp100r-tmp102r;
	tmp107i = Im[5]+tmp94i;
	tmp107r = Re[5]+tmp94r;
	tmp108i = tmp89i+tmp53i;
	tmp108r = tmp89r+tmp53r;
	tmp109i = tmp89i-tmp53i;
	tmp109r = tmp89r-tmp53r;
	tmp110i = -C13*tmp108i;
	tmp110r = -C13*tmp108r;
	tmp111i = C17*tmp109r;
	tmp111r = -C17*tmp109i;
	tmp112i = tmp110i+tmp17i;
	tmp112r = tmp110r+tmp17r;
	tmp113i = tmp112i+tmp111i;
	tmp113r = tmp112r+tmp111r;
	tmp114i = tmp112i-tmp111i;
	tmp114r = tmp112r-tmp111r;
	tmp115r = tmp17r+tmp108r;
	tmp116i = tmp35i+tmp107i;
	tmp116r = tmp35r+tmp107r;
	tmp117i = tmp35i-tmp107i;
	tmp117r = tmp35r-tmp107r;
	tmp118i = -C13*tmp116i;
	tmp118r = -C13*tmp116r;
	tmp119i = C17*tmp117r;
	tmp119r = -C17*tmp117i;
	tmp120i = tmp118i+tmp71i;
	tmp120r = tmp118r+tmp71r;
	tmp121i = tmp120i+tmp119i;
	tmp121r = tmp120r+tmp119r;
	tmp122i = tmp120i-tmp119i;
	tmp122r = tmp120r-tmp119r;
	tmp123r = tmp71r+tmp116r;
	tmp124r = tmp115r+tmp123r;
	tmp127i = tmp114i-tmp122i;
	tmp127r = tmp114r-tmp122r;
	tmp128i = tmp113i+tmp121i;
	tmp128r = tmp113r+tmp121r;
	tmp130i = tmp85i+tmp49i;
	tmp130r = tmp85r+tmp49r;
	tmp131i = tmp85i-tmp49i;
	tmp131r = tmp85r-tmp49r;
	tmp132i = -C13*tmp130i;
	tmp132r = -C13*tmp130r;
	tmp133i = C17*tmp131r;
	tmp133r = -C17*tmp131i;
	tmp134i = tmp132i+tmp13i;
	tmp134r = tmp132r+tmp13r;
	tmp135i = tmp134i+tmp133i;
	tmp135r = tmp134r+tmp133r;
	tmp136i = tmp134i-tmp133i;
	tmp136r = tmp134r-tmp133r;
	tmp137i = tmp13i+tmp130i;
	tmp137r = tmp13r+tmp130r;
	tmp138i = tmp31i+tmp103i;
	tmp138r = tmp31r+tmp103r;
	tmp139i = tmp31i-tmp103i;
	tmp139r = tmp31r-tmp103r;
	tmp140i = -C13*tmp138i;
	tmp140r = -C13*tmp138r;
	tmp141i = C17*tmp139r;
	tmp141r = -C17*tmp139i;
	tmp142i = tmp140i+tmp67i;
	tmp142r = tmp140r+tmp67r;
	tmp143i = tmp142i+tmp141i;
	tmp143r = tmp142r+tmp141r;
	tmp144i = tmp142i-tmp141i;
	tmp144r = tmp142r-tmp141r;
	tmp145i = tmp67i+tmp138i;
	tmp145r = tmp67r+tmp138r;
	tmp146i = tmp137i+tmp145i;
	tmp146r = tmp137r+tmp145r;
	tmp149i = tmp136i-tmp144i;
	tmp149r = tmp136r-tmp144r;
	tmp151i = tmp135i-tmp143i;
	tmp151r = tmp135r-tmp143r;
	tmp152i = tmp87i+tmp51i;
	tmp152r = tmp87r+tmp51r;
	tmp153i = tmp87i-tmp51i;
	tmp153r = tmp87r-tmp51r;
	tmp154i = -C13*tmp152i;
	tmp154r = -C13*tmp152r;
	tmp155i = C17*tmp153r;
	tmp155r = -C17*tmp153i;
	tmp156i = tmp154i+tmp15i;
	tmp156r = tmp154r+tmp15r;
	tmp157i = tmp156i+tmp155i;
	tmp157r = tmp156r+tmp155r;
	tmp158i = tmp156i-tmp155i;
	tmp158r = tmp156r-tmp155r;
	tmp159i = tmp15i+tmp152i;
	tmp159r = tmp15r+tmp152r;
	tmp160i = tmp33i+tmp105i;
	tmp160r = tmp33r+tmp105r;
	tmp161i = tmp33i-tmp105i;
	tmp161r = tmp33r-tmp105r;
	tmp162i = -C13*tmp160i;
	tmp162r = -C13*tmp160r;
	tmp163i = C17*tmp161r;
	tmp163r = -C17*tmp161i;
	tmp164i = tmp162i+tmp69i;
	tmp164r = tmp162r+tmp69r;
	tmp165i = tmp164i+tmp163i;
	tmp165r = tmp164r+tmp163r;
	tmp166i = tmp164i-tmp163i;
	tmp166r = tmp164r-tmp163r;
	tmp167i = tmp69i+tmp160i;
	tmp167r = tmp69r+tmp160r;
	tmp168i = tmp159i+tmp167i;
	tmp168r = tmp159r+tmp167r;
	tmp170i = tmp158i+tmp166i;
	tmp170r = tmp158r+tmp166r;
	tmp173i = tmp157i-tmp165i;
	tmp173r = tmp157r-tmp165r;
	tmp174i = tmp88i+tmp52i;
	tmp174r = tmp88r+tmp52r;
	tmp175i = tmp88i-tmp52i;
	tmp175r = tmp88r-tmp52r;
	tmp176i = -C13*tmp174i;
	tmp176r = -C13*tmp174r;
	tmp177i = C17*tmp175r;
	tmp177r = -C17*tmp175i;
	tmp178i = tmp176i+tmp16i;
	tmp178r = tmp176r+tmp16r;
	tmp179i = tmp178i+tmp177i;
	tmp179r = tmp178r+tmp177r;
	tmp180i = tmp178i-tmp177i;
	tmp180r = tmp178r-tmp177r;
	tmp181i = tmp16i+tmp174i;
	tmp181r = tmp16r+tmp174r;
	tmp182i = tmp34i+tmp106i;
	tmp182r = tmp34r+tmp106r;
	tmp183i = tmp34i-tmp106i;
	tmp183r = tmp34r-tmp106r;
	tmp184i = -C13*tmp182i;
	tmp184r = -C13*tmp182r;
	tmp185i = C17*tmp183r;
	tmp185r = -C17*tmp183i;
	tmp186i = tmp184i+tmp70i;
	tmp186r = tmp184r+tmp70r;
	tmp187i = tmp186i+tmp185i;
	tmp187r = tmp186r+tmp185r;
	tmp188i = tmp186i-tmp185i;
	tmp188r = tmp186r-tmp185r;
	tmp189i = tmp70i+tmp182i;
	tmp189r = tmp70r+tmp182r;
	tmp191i = tmp181i-tmp189i;
	tmp191r = tmp181r-tmp189r;
	tmp192i = tmp180i+tmp188i;
	tmp192r = tmp180r+tmp188r;
	tmp195i = tmp179i-tmp187i;
	tmp195r = tmp179r-tmp187r;
	tmp196i = tmp86i+tmp50i;
	tmp196r = tmp86r+tmp50r;
	tmp197i = tmp86i-tmp50i;
	tmp197r = tmp86r-tmp50r;
	tmp198i = -C13*tmp196i;
	tmp198r = -C13*tmp196r;
	tmp199i = C17*tmp197r;
	tmp199r = -C17*tmp197i;
	tmp200i = tmp198i+tmp14i;
	tmp200r = tmp198r+tmp14r;
	tmp201i = tmp200i+tmp199i;
	tmp201r = tmp200r+tmp199r;
	tmp202i = tmp200i-tmp199i;
	tmp202r = tmp200r-tmp199r;
	tmp203i = tmp14i+tmp196i;
	tmp203r = tmp14r+tmp196r;
	tmp204i = tmp32i+tmp104i;
	tmp204r = tmp32r+tmp104r;
	tmp205i = tmp32i-tmp104i;
	tmp205r = tmp32r-tmp104r;
	tmp206i = -C13*tmp204i;
	tmp206r = -C13*tmp204r;
	tmp207i = C17*tmp205r;
	tmp207r = -C17*tmp205i;
	tmp208i = tmp206i+tmp68i;
	tmp208r = tmp206r+tmp68r;
	tmp209i = tmp208i+tmp207i;
	tmp209r = tmp208r+tmp207r;
	tmp210i = tmp208i-tmp207i;
	tmp210r = tmp208r-tmp207r;
	tmp211i = tmp68i+tmp204i;
	tmp211r = tmp68r+tmp204r;
	tmp213i = tmp203i-tmp211i;
	tmp213r = tmp203r-tmp211r;
	tmp214i = tmp202i+tmp210i;
	tmp214r = tmp202r+tmp210r;
	tmp216i = tmp201i+tmp209i;
	tmp216r = tmp201r+tmp209r;
	Re[0] = tmp124r;
	Im[0] = C0;
	Re[1] = tmp151r;
	Im[1] = tmp151i;
	Re[2] = tmp170r;
	Im[2] = tmp170i;
	Re[3] = tmp191r;
	Im[3] = tmp191i;
	Re[4] = tmp216r;
	Im[4] = tmp216i;
	Re[5] = tmp127r;
	Im[5] = tmp127i;
	Re[6] = tmp146r;
	Im[6] = tmp146i;
	Re[7] = tmp173r;
	Im[7] = tmp173i;
	Re[8] = tmp192r;
	Im[8] = tmp192i;
	Re[9] = tmp213r;
	Im[9] = tmp213i;
	Re[10] = tmp128r;
	Im[10] = tmp128i;
	Re[11] = tmp149r;
	Im[11] = tmp149i;
	Re[12] = tmp168r;
	Im[12] = tmp168i;
	Re[13] = tmp195r;
	Im[13] = tmp195i;
	Re[14] = tmp214r;
	Im[14] = tmp214i;
}

/*
*	Number of additions = 320
*	Number of multiplications = 86
*	Number of sign changes = 10
*	Number of assigns = 362
*	Total number of operations = 778
*/
void	MFFTR32(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, tmp105r,
		tmp106r, tmp107r, tmp108r, tmp109r, tmp110r, tmp111r, tmp112r,
		tmp113r, tmp114r, tmp115r, tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, tmp121r, tmp122r, tmp123r, tmp124r, tmp125r, tmp126r,
		tmp127r, tmp128r, tmp129r, tmp130r, tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r, tmp136r, tmp137r, tmp138r, tmp139r, tmp140r,
		tmp141r, tmp142r, tmp143r, tmp144r, tmp145r, tmp146r, tmp147r,
		tmp148r, tmp149r, tmp150r, tmp151r, tmp152r, tmp153r, tmp154r,
		tmp155r, tmp156r, tmp157r, tmp158r, tmp159r, tmp160r, tmp161r,
		tmp162r, tmp163r, tmp164r, tmp165r, tmp166r, tmp167r, tmp168r,
		tmp169r, tmp170r /*tmp171r,*/ /*tmp172r,*/ /*tmp173r,*/ /*tmp174r,*/ /*tmp175r,*/
		/*tmp176r,*/ /*tmp177r,*/ /*tmp178r,*/ /*tmp179r,*/ /*tmp180r,*/ /*tmp181r,*/ /*tmp182r,*/
		/*tmp183r,*/ /*tmp184r,*/ /*tmp185r,*/ /*tmp186r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, tmp74i, tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i,
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, tmp97i, tmp98i,
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i, tmp104i, tmp105i,
		tmp106i, tmp107i, tmp108i, tmp109i, tmp110i, tmp111i, tmp112i,
		tmp113i, tmp114i, tmp115i, tmp116i, tmp117i, tmp118i, tmp119i,
		tmp120i, tmp121i, tmp122i, tmp123i, tmp124i, tmp125i, tmp126i,
		tmp127i, tmp128i, tmp129i, tmp130i, tmp131i, tmp132i, tmp133i,
		tmp134i, tmp135i, tmp136i, tmp137i, tmp138i, tmp139i, tmp140i,
		tmp141i, tmp142i, tmp143i, tmp144i, tmp145i, tmp146i, tmp147i,
		tmp148i, tmp149i, tmp150i, tmp151i, tmp152i, tmp153i, tmp154i,
		tmp155i, tmp156i, tmp157i, tmp158i /*tmp159i,*/ /*tmp160i,*/ /*tmp161i,*/
		/*tmp162i,*/ /*tmp163i,*/ /*tmp164i,*/ /*tmp165i,*/ /*tmp166i,*/ /*tmp167i,*/ /*tmp168i,*/
		/*tmp169i,*/ /*tmp170i,*/ /*tmp171i,*/ /*tmp172i,*/ /*tmp173i,*/ /*tmp174i*/;

	const FFT_precision	C3 =     0.19509032201613;	/* FFT_precisionCONST	*/
	const FFT_precision	C5 =     0.38268343236509;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =      0.5555702330196;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.83146961230255;	/* FFT_precisionCONST	*/
	const FFT_precision	C4 =     0.92387953251129;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =     0.98078528040323;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =                    2;	/* INTEGERCONST	*/

	tmp0r = C9*Re[8];
	tmp0i = C9*Im[8];
	tmp1r = Re[0]+tmp0r;
	tmp2r = Re[0]-tmp0r;
	tmp3r = Re[0]+tmp0i;
	tmp4r = Re[0]-tmp0i;
	tmp1i = Im[4]-Im[12];
	tmp5r = Re[4]+Re[12];
	tmp2i = Im[4]+Im[12];
	tmp6r = Re[4]-Re[12];
	tmp3i = Im[12]-Im[4];
	tmp7r = Re[12]+Re[4];
	tmp4i = Im[12]+Im[4];
	tmp8r = Re[12]-Re[4];
	tmp5i = tmp1i+tmp3i;
	tmp9r = tmp5r+tmp7r;
	tmp6i = tmp1i-tmp3i;
	tmp10r = tmp5r-tmp7r;
	tmp7i = tmp2i-tmp8r;
	tmp11r = tmp6r+tmp4i;
	tmp8i = tmp2i+tmp8r;
	tmp12r = tmp6r-tmp4i;
	tmp9i = C8*(tmp7i-tmp11r);
	tmp13r = C8*(tmp11r+tmp7i);
	tmp10i = -C8*(tmp8i+tmp12r);
	tmp14r = -C8*(tmp12r-tmp8i);
	tmp15r = tmp1r+tmp9r;
	tmp16r = tmp1r-tmp9r;
	tmp17r = tmp3r+tmp13r;
	tmp18r = tmp3r-tmp13r;
	tmp19r = tmp2r+tmp6i;
	tmp20r = tmp2r-tmp6i;
	tmp21r = tmp4r+tmp14r;
	tmp22r = tmp4r-tmp14r;
	tmp11i = Im[2]-Im[14];
	tmp23r = Re[2]+Re[14];
	tmp12i = Im[2]+Im[14];
	tmp24r = Re[2]-Re[14];
	tmp13i = Im[10]-Im[6];
	tmp25r = Re[10]+Re[6];
	tmp14i = Im[10]+Im[6];
	tmp26r = Re[10]-Re[6];
	tmp15i = tmp11i+tmp13i;
	tmp27r = tmp23r+tmp25r;
	tmp16i = tmp11i-tmp13i;
	tmp28r = tmp23r-tmp25r;
	tmp17i = tmp12i-tmp26r;
	tmp29r = tmp24r+tmp14i;
	tmp18i = tmp12i+tmp26r;
	tmp30r = tmp24r-tmp14i;
	tmp31r = tmp27r;
	tmp19i = tmp15i;
	tmp32r = C4*tmp29r+C5*tmp17i;
	tmp20i = C4*tmp17i-C5*tmp29r;
	tmp33r = C8*(tmp28r+tmp16i);
	tmp21i = C8*(tmp16i-tmp28r);
	tmp34r = C5*tmp30r+C4*tmp18i;
	tmp22i = C5*tmp18i-C4*tmp30r;
	tmp23i = Im[6]-Im[10];
	tmp35r = Re[6]+Re[10];
	tmp24i = Im[6]+Im[10];
	tmp36r = Re[6]-Re[10];
	tmp25i = Im[14]-Im[2];
	tmp37r = Re[14]+Re[2];
	tmp26i = Im[14]+Im[2];
	tmp38r = Re[14]-Re[2];
	tmp27i = tmp23i+tmp25i;
	tmp39r = tmp35r+tmp37r;
	tmp28i = tmp23i-tmp25i;
	tmp40r = tmp35r-tmp37r;
	tmp29i = tmp24i-tmp38r;
	tmp41r = tmp36r+tmp26i;
	tmp30i = tmp24i+tmp38r;
	tmp42r = tmp36r-tmp26i;
	tmp43r = tmp39r;
	tmp31i = tmp27i;
	tmp44r = C5*tmp41r+C4*tmp29i;
	tmp32i = C5*tmp29i-C4*tmp41r;
	tmp45r = -C8*(tmp40r-tmp28i);
	tmp33i = -C8*(tmp28i+tmp40r);
	tmp46r = -C4*tmp42r-C5*tmp30i;
	tmp34i = -C4*tmp30i+C5*tmp42r;
	tmp47r = tmp31r+tmp43r;
	tmp35i = tmp19i+tmp31i;
	tmp48r = tmp32r+tmp44r;
	tmp36i = tmp20i+tmp32i;
	tmp49r = tmp33r+tmp45r;
	tmp37i = tmp21i+tmp33i;
	tmp50r = tmp34r+tmp46r;
	tmp38i = tmp22i+tmp34i;
	tmp51r = -tmp19i+tmp31i;
	tmp39i = tmp31r-tmp43r;
	tmp52r = -tmp20i+tmp32i;
	tmp40i = tmp32r-tmp44r;
	tmp53r = -tmp21i+tmp33i;
	tmp41i = tmp33r-tmp45r;
	tmp54r = -tmp22i+tmp34i;
	tmp42i = tmp34r-tmp46r;
	tmp55r = tmp15r+tmp47r;
	tmp43i = tmp5i+tmp35i;
	tmp56r = tmp17r+tmp48r;
	tmp44i = tmp9i+tmp36i;
	tmp57r = tmp19r+tmp49r;
	tmp45i = -tmp10r+tmp37i;
	tmp58r = tmp21r+tmp50r;
	tmp46i = tmp10i+tmp38i;
	tmp59r = tmp16r-tmp51r;
	tmp47i = -tmp5i-tmp39i;
	tmp60r = tmp18r-tmp52r;
	tmp48i = -tmp9i-tmp40i;
	tmp61r = tmp20r-tmp53r;
	tmp49i = tmp10r-tmp41i;
	tmp62r = tmp22r-tmp54r;
	tmp50i = -tmp10i-tmp42i;
	tmp63r = tmp15r-tmp47r;
	tmp51i = tmp5i-tmp35i;
	tmp64r = tmp17r-tmp48r;
	tmp52i = tmp9i-tmp36i;
	tmp65r = tmp19r-tmp49r;
	tmp53i = -tmp10r-tmp37i;
	tmp66r = tmp21r-tmp50r;
	tmp54i = tmp10i-tmp38i;
	tmp67r = tmp16r+tmp51r;
	tmp55i = -tmp5i+tmp39i;
	tmp68r = tmp18r+tmp52r;
	tmp56i = -tmp9i+tmp40i;
	tmp69r = tmp20r+tmp53r;
	tmp57i = tmp10r+tmp41i;
	tmp70r = tmp22r+tmp54r;
	tmp58i = -tmp10i+tmp42i;
	tmp59i = Im[1]+Im[15];
	tmp71r = Re[1]-Re[15];
	tmp60i = Im[1]-Im[15];
	tmp72r = Re[1]+Re[15];
	tmp61i = Im[9]+Im[7];
	tmp73r = Re[9]-Re[7];
	tmp62i = Im[9]-Im[7];
	tmp74r = Re[9]+Re[7];
	tmp63i = tmp59i+tmp61i;
	tmp75r = tmp71r+tmp73r;
	tmp64i = tmp59i-tmp61i;
	tmp76r = tmp71r-tmp73r;
	tmp65i = tmp60i-tmp74r;
	tmp77r = tmp72r+tmp62i;
	tmp66i = tmp60i+tmp74r;
	tmp78r = tmp72r-tmp62i;
	tmp67i = Im[5]+Im[11];
	tmp79r = Re[5]-Re[11];
	tmp68i = Im[5]-Im[11];
	tmp80r = Re[5]+Re[11];
	tmp69i = Im[13]+Im[3];
	tmp81r = Re[13]-Re[3];
	tmp70i = Im[13]-Im[3];
	tmp82r = Re[13]+Re[3];
	tmp71i = tmp67i+tmp69i;
	tmp83r = tmp79r+tmp81r;
	tmp72i = tmp67i-tmp69i;
	tmp84r = tmp79r-tmp81r;
	tmp73i = tmp68i-tmp82r;
	tmp85r = tmp80r+tmp70i;
	tmp74i = tmp68i+tmp82r;
	tmp86r = tmp80r-tmp70i;
	tmp75i = C8*(tmp73i-tmp85r);
	tmp87r = C8*(tmp85r+tmp73i);
	tmp76i = -C8*(tmp74i+tmp86r);
	tmp88r = -C8*(tmp86r-tmp74i);
	tmp77i = tmp63i+tmp71i;
	tmp89r = tmp75r+tmp83r;
	tmp78i = tmp63i-tmp71i;
	tmp90r = tmp75r-tmp83r;
	tmp79i = tmp65i+tmp75i;
	tmp91r = tmp77r+tmp87r;
	tmp80i = tmp65i-tmp75i;
	tmp92r = tmp77r-tmp87r;
	tmp81i = tmp64i-tmp84r;
	tmp93r = tmp76r+tmp72i;
	tmp82i = tmp64i+tmp84r;
	tmp94r = tmp76r-tmp72i;
	tmp83i = tmp66i+tmp76i;
	tmp95r = tmp78r+tmp88r;
	tmp84i = tmp66i-tmp76i;
	tmp96r = tmp78r-tmp88r;
	tmp97r = tmp89r;
	tmp85i = tmp77i;
	tmp98r = C2*tmp91r+C3*tmp79i;
	tmp86i = C2*tmp79i-C3*tmp91r;
	tmp99r = C4*tmp93r+C5*tmp81i;
	tmp87i = C4*tmp81i-C5*tmp93r;
	tmp100r = C6*tmp95r+C7*tmp83i;
	tmp88i = C6*tmp83i-C7*tmp95r;
	tmp101r = C8*(tmp90r+tmp78i);
	tmp89i = C8*(tmp78i-tmp90r);
	tmp102r = C7*tmp92r+C6*tmp80i;
	tmp90i = C7*tmp80i-C6*tmp92r;
	tmp103r = C5*tmp94r+C4*tmp82i;
	tmp91i = C5*tmp82i-C4*tmp94r;
	tmp104r = C3*tmp96r+C2*tmp84i;
	tmp92i = C3*tmp84i-C2*tmp96r;
	tmp93i = Im[3]+Im[13];
	tmp105r = Re[3]-Re[13];
	tmp94i = Im[3]-Im[13];
	tmp106r = Re[3]+Re[13];
	tmp95i = Im[11]+Im[5];
	tmp107r = Re[11]-Re[5];
	tmp96i = Im[11]-Im[5];
	tmp108r = Re[11]+Re[5];
	tmp97i = tmp93i+tmp95i;
	tmp109r = tmp105r+tmp107r;
	tmp98i = tmp93i-tmp95i;
	tmp110r = tmp105r-tmp107r;
	tmp99i = tmp94i-tmp108r;
	tmp111r = tmp106r+tmp96i;
	tmp100i = tmp94i+tmp108r;
	tmp112r = tmp106r-tmp96i;
	tmp101i = Im[7]+Im[9];
	tmp113r = Re[7]-Re[9];
	tmp102i = Im[7]-Im[9];
	tmp114r = Re[7]+Re[9];
	tmp103i = Im[15]+Im[1];
	tmp115r = Re[15]-Re[1];
	tmp104i = Im[15]-Im[1];
	tmp116r = Re[15]+Re[1];
	tmp105i = tmp101i+tmp103i;
	tmp117r = tmp113r+tmp115r;
	tmp106i = tmp101i-tmp103i;
	tmp118r = tmp113r-tmp115r;
	tmp107i = tmp102i-tmp116r;
	tmp119r = tmp114r+tmp104i;
	tmp108i = tmp102i+tmp116r;
	tmp120r = tmp114r-tmp104i;
	tmp109i = C8*(tmp107i-tmp119r);
	tmp121r = C8*(tmp119r+tmp107i);
	tmp110i = -C8*(tmp108i+tmp120r);
	tmp122r = -C8*(tmp120r-tmp108i);
	tmp111i = tmp97i+tmp105i;
	tmp123r = tmp109r+tmp117r;
	tmp112i = tmp97i-tmp105i;
	tmp124r = tmp109r-tmp117r;
	tmp113i = tmp99i+tmp109i;
	tmp125r = tmp111r+tmp121r;
	tmp114i = tmp99i-tmp109i;
	tmp126r = tmp111r-tmp121r;
	tmp115i = tmp98i-tmp118r;
	tmp127r = tmp110r+tmp106i;
	tmp116i = tmp98i+tmp118r;
	tmp128r = tmp110r-tmp106i;
	tmp117i = tmp100i+tmp110i;
	tmp129r = tmp112r+tmp122r;
	tmp118i = tmp100i-tmp110i;
	tmp130r = tmp112r-tmp122r;
	tmp131r = tmp123r;
	tmp119i = tmp111i;
	tmp132r = C6*tmp125r+C7*tmp113i;
	tmp120i = C6*tmp113i-C7*tmp125r;
	tmp133r = C5*tmp127r+C4*tmp115i;
	tmp121i = C5*tmp115i-C4*tmp127r;
	tmp134r = -C3*tmp129r+C2*tmp117i;
	tmp122i = -C3*tmp117i-C2*tmp129r;
	tmp135r = -C8*(tmp124r-tmp112i);
	tmp123i = -C8*(tmp112i+tmp124r);
	tmp136r = -C2*tmp126r+C3*tmp114i;
	tmp124i = -C2*tmp114i-C3*tmp126r;
	tmp137r = -C4*tmp128r-C5*tmp116i;
	tmp125i = -C4*tmp116i+C5*tmp128r;
	tmp138r = -C7*tmp130r-C6*tmp118i;
	tmp126i = -C7*tmp118i+C6*tmp130r;
	tmp139r = tmp97r+tmp131r;
	tmp127i = tmp85i+tmp119i;
	tmp140r = tmp98r+tmp132r;
	tmp128i = tmp86i+tmp120i;
	tmp141r = tmp99r+tmp133r;
	tmp129i = tmp87i+tmp121i;
	tmp142r = tmp100r+tmp134r;
	tmp130i = tmp88i+tmp122i;
	tmp143r = tmp101r+tmp135r;
	tmp131i = tmp89i+tmp123i;
	tmp144r = tmp102r+tmp136r;
	tmp132i = tmp90i+tmp124i;
	tmp145r = tmp103r+tmp137r;
	tmp133i = tmp91i+tmp125i;
	tmp146r = tmp104r+tmp138r;
	tmp134i = tmp92i+tmp126i;
	tmp147r = -tmp85i+tmp119i;
	tmp135i = tmp97r-tmp131r;
	tmp148r = -tmp86i+tmp120i;
	tmp136i = tmp98r-tmp132r;
	tmp149r = -tmp87i+tmp121i;
	tmp137i = tmp99r-tmp133r;
	tmp150r = -tmp88i+tmp122i;
	tmp138i = tmp100r-tmp134r;
	tmp151r = -tmp89i+tmp123i;
	tmp139i = tmp101r-tmp135r;
	tmp152r = -tmp90i+tmp124i;
	tmp140i = tmp102r-tmp136r;
	tmp153r = -tmp91i+tmp125i;
	tmp141i = tmp103r-tmp137r;
	tmp154r = -tmp92i+tmp126i;
	tmp142i = tmp104r-tmp138r;
	tmp155r = tmp55r+tmp139r;
	tmp143i = tmp43i+tmp127i;
	tmp156r = tmp56r+tmp140r;
	tmp144i = tmp44i+tmp128i;
	tmp157r = tmp57r+tmp141r;
	tmp145i = tmp45i+tmp129i;
	tmp158r = tmp58r+tmp142r;
	tmp146i = tmp46i+tmp130i;
	tmp159r = tmp59r+tmp143r;
	tmp147i = tmp47i+tmp131i;
	tmp160r = tmp60r+tmp144r;
	tmp148i = tmp48i+tmp132i;
	tmp161r = tmp61r+tmp145r;
	tmp149i = tmp49i+tmp133i;
	tmp162r = tmp62r+tmp146r;
	tmp150i = tmp50i+tmp134i;
	tmp163r = tmp63r-tmp147r;
	tmp151i = tmp51i-tmp135i;
	tmp164r = tmp64r-tmp148r;
	tmp152i = tmp52i-tmp136i;
	tmp165r = tmp65r-tmp149r;
	tmp153i = tmp53i-tmp137i;
	tmp166r = tmp66r-tmp150r;
	tmp154i = tmp54i-tmp138i;
	tmp167r = tmp67r-tmp151r;
	tmp155i = tmp55i-tmp139i;
	tmp168r = tmp68r-tmp152r;
	tmp156i = tmp56i-tmp140i;
	tmp169r = tmp69r-tmp153r;
	tmp157i = tmp57i-tmp141i;
	tmp170r = tmp70r-tmp154r;
	tmp158i = tmp58i-tmp142i;
	Re[0] = tmp155r;
	Im[0] = tmp143i;
	Re[1] = tmp156r;
	Im[1] = tmp144i;
	Re[2] = tmp157r;
	Im[2] = tmp145i;
	Re[3] = tmp158r;
	Im[3] = tmp146i;
	Re[4] = tmp159r;
	Im[4] = tmp147i;
	Re[5] = tmp160r;
	Im[5] = tmp148i;
	Re[6] = tmp161r;
	Im[6] = tmp149i;
	Re[7] = tmp162r;
	Im[7] = tmp150i;
	Re[8] = tmp163r;
	Im[8] = tmp151i;
	Re[9] = tmp164r;
	Im[9] = tmp152i;
	Re[10] = tmp165r;
	Im[10] = tmp153i;
	Re[11] = tmp166r;
	Im[11] = tmp154i;
	Re[12] = tmp167r;
	Im[12] = tmp155i;
	Re[13] = tmp168r;
	Im[13] = tmp156i;
	Re[14] = tmp169r;
	Im[14] = tmp157i;
	Re[15] = tmp170r;
	Im[15] = tmp158i;
}

/*
*	Number of additions = 192
*	Number of multiplications = 116
*	Number of sign changes = 31
*	Number of assigns = 272
*	Total number of operations = 611
*/
void	MIFFTR32(FFT_precision *Re, FFT_precision *Im)
{
	register FFT_precision	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r, tmp74r, tmp75r, tmp76r, tmp77r,
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r,
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r,
		tmp92r, tmp93r, tmp94r, tmp95r, tmp96r, tmp97r, tmp98r,
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r, tmp104r, tmp105r,
		tmp106r, tmp107r, tmp108r, tmp109r, tmp110r, tmp111r, tmp112r,
		tmp113r, tmp114r, tmp115r, tmp116r, tmp117r, tmp118r, tmp119r,
		tmp120r, tmp121r, tmp122r, tmp123r, tmp124r, tmp125r, tmp126r,
		tmp127r, tmp128r, tmp129r, tmp130r, tmp131r, tmp132r, tmp133r,
		tmp134r, tmp135r /*tmp136r,*/ /*tmp137r,*/ /*tmp138r,*/ /*tmp139r,*/ /*tmp140r,*/
		/*tmp141r,*/ /*tmp142r,*/ /*tmp143r,*/ /*tmp144r,*/ /*tmp145r,*/ /*tmp146r,*/ /*tmp147r,*/
		/*tmp148r,*/ /*tmp149r,*/ /*tmp150r,*/ /*tmp151r*/;
	register FFT_precision	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, /*tmp12i,*/ tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, /*tmp18i,*/ tmp19i, tmp20i, tmp21i,
		/*tmp22i,*/ tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, /*tmp30i,*/ tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, /*tmp38i,*/ tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, /*tmp56i,*/
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i, /*tmp74i,*/ tmp75i, tmp76i, tmp77i,
		tmp78i, tmp79i, tmp80i, tmp81i, /*tmp82i,*/ tmp83i, tmp84i,
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i,
		tmp92i, tmp93i, tmp94i, tmp95i, tmp96i, tmp97i, /*tmp98i,*/
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i, tmp104i, tmp105i,
		tmp106i, tmp107i, tmp108i, tmp109i, tmp110i, tmp111i, tmp112i,
		tmp113i /*tmp114i,*/ /*tmp115i,*/ /*tmp116i,*/ /*tmp117i,*/ /*tmp118i,*/ /*tmp119i,*/
		/*tmp120i,*/ /*tmp121i,*/ /*tmp122i,*/ /*tmp123i,*/ /*tmp124i,*/ /*tmp125i,*/ /*tmp126i,*/
		/*tmp127i,*/ /*tmp128i,*/ /*tmp129i*/;

	const FFT_precision	C0 =                    0;	/* ZEROCONST	*/
	const FFT_precision	C3 =     0.19509032201613;	/* FFT_precisionCONST	*/
	const FFT_precision	C5 =     0.38268343236509;	/* FFT_precisionCONST	*/
	const FFT_precision	C7 =      0.5555702330196;	/* FFT_precisionCONST	*/
	const FFT_precision	C8 =     0.70710678118655;	/* FFT_precisionCONST	*/
	const FFT_precision	C6 =     0.83146961230255;	/* FFT_precisionCONST	*/
	const FFT_precision	C4 =     0.92387953251129;	/* FFT_precisionCONST	*/
	const FFT_precision	C2 =     0.98078528040323;	/* FFT_precisionCONST	*/
	const FFT_precision	C9 =                    2;	/* INTEGERCONST	*/

	tmp0r = C9*Re[0];
	tmp0i = C9*Im[0];
	tmp1r = C9*Re[8];
	tmp1i = C9*Im[8];
	tmp2r = tmp0r+tmp1r;
	tmp3r = tmp0r-tmp1r;
	tmp4r = C9*Re[4];
	tmp2i = C9*Im[4];
	tmp5r = C9*Re[12];
	tmp3i = C9*Im[12];
	tmp6r = tmp4r+tmp5r;
	tmp7r = tmp4r-tmp5r;
	tmp4i = C8*(tmp2i-tmp3i);
	tmp8r = -C8*(tmp3i+tmp2i);
	tmp5i = -C8*(tmp2i-tmp3i);
	tmp9r = -C8*(tmp3i+tmp2i);
	tmp10r = tmp2r+tmp6r;
	tmp11r = tmp2r-tmp6r;
	tmp6i = tmp0i+tmp4i;
	tmp12r = -(tmp1i-tmp8r);
	tmp7i = tmp0i-tmp4i;
	tmp13r = -(tmp1i+tmp8r);
	tmp8i = tmp0i+tmp5i;
	tmp14r = tmp1i+tmp9r;
	tmp9i = tmp0i-tmp5i;
	tmp15r = tmp1i-tmp9r;
	tmp16r = C9*Re[2];
	tmp10i = C9*Im[2];
	tmp17r = C9*Re[10];
	tmp11i = C9*Im[10];
	tmp18r = tmp16r+tmp17r;
	tmp19r = tmp16r-tmp17r;
	tmp20r = tmp18r;
	tmp21r = -C4*tmp11i-C5*tmp10i;
	tmp13i = C4*tmp10i-C5*tmp11i;
	tmp22r = C8*tmp19r;
	tmp14i = C8*tmp19r;
	tmp23r = C5*tmp11i-C4*tmp10i;
	tmp15i = C5*tmp10i+C4*tmp11i;
	tmp24r = C9*Re[6];
	tmp16i = C9*Im[6];
	tmp25r = C9*Re[14];
	tmp17i = C9*Im[14];
	tmp26r = tmp24r+tmp25r;
	tmp27r = tmp24r-tmp25r;
	tmp28r = tmp26r;
	tmp29r = -C5*tmp17i-C4*tmp16i;
	tmp19i = C5*tmp16i-C4*tmp17i;
	tmp30r = -C8*tmp27r;
	tmp20i = C8*tmp27r;
	tmp31r = -C4*tmp17i+C5*tmp16i;
	tmp21i = -C4*tmp16i-C5*tmp17i;
	tmp32r = tmp20r+tmp28r;
	tmp33r = tmp21r+tmp29r;
	tmp23i = tmp13i+tmp19i;
	tmp34r = tmp22r+tmp30r;
	tmp24i = tmp14i+tmp20i;
	tmp35r = tmp23r+tmp31r;
	tmp25i = tmp15i+tmp21i;
	tmp36r = C0;
	tmp26i = -(tmp20r-tmp28r);
	tmp37r = tmp13i-tmp19i;
	tmp27i = -(tmp21r-tmp29r);
	tmp38r = tmp14i-tmp20i;
	tmp28i = -(tmp22r-tmp30r);
	tmp39r = tmp15i-tmp21i;
	tmp29i = -(tmp23r-tmp31r);
	tmp40r = tmp10r+tmp32r;
	tmp41r = tmp12r+tmp33r;
	tmp31i = tmp6i+tmp23i;
	tmp42r = tmp3r+tmp34r;
	tmp32i = tmp7r+tmp24i;
	tmp43r = tmp14r+tmp35r;
	tmp33i = tmp8i+tmp25i;
	tmp44r = tmp11r-tmp36r;
	tmp34i = -tmp26i;
	tmp45r = tmp13r-tmp37r;
	tmp35i = tmp7i-tmp27i;
	tmp46r = tmp3r-tmp38r;
	tmp36i = -tmp7r-tmp28i;
	tmp47r = tmp15r-tmp39r;
	tmp37i = tmp9i-tmp29i;
	tmp48r = tmp10r-tmp32r;
	tmp49r = tmp12r-tmp33r;
	tmp39i = tmp6i-tmp23i;
	tmp50r = tmp3r-tmp34r;
	tmp40i = tmp7r-tmp24i;
	tmp51r = tmp14r-tmp35r;
	tmp41i = tmp8i-tmp25i;
	tmp52r = tmp11r+tmp36r;
	tmp42i = tmp26i;
	tmp53r = tmp13r+tmp37r;
	tmp43i = tmp7i+tmp27i;
	tmp54r = tmp3r+tmp38r;
	tmp44i = -tmp7r+tmp28i;
	tmp55r = tmp15r+tmp39r;
	tmp45i = tmp9i+tmp29i;
	tmp56r = C9*Re[1];
	tmp46i = C9*Im[1];
	tmp57r = C9*Re[9];
	tmp47i = C9*Im[9];
	tmp58r = tmp56r+tmp57r;
	tmp59r = tmp56r-tmp57r;
	tmp60r = C9*Re[5];
	tmp48i = C9*Im[5];
	tmp61r = C9*Re[13];
	tmp49i = C9*Im[13];
	tmp62r = tmp60r+tmp61r;
	tmp63r = tmp60r-tmp61r;
	tmp50i = C8*(tmp48i-tmp49i);
	tmp64r = -C8*(tmp49i+tmp48i);
	tmp51i = -C8*(tmp48i-tmp49i);
	tmp65r = -C8*(tmp49i+tmp48i);
	tmp66r = tmp58r+tmp62r;
	tmp67r = tmp58r-tmp62r;
	tmp52i = tmp46i+tmp50i;
	tmp68r = -(tmp47i-tmp64r);
	tmp53i = tmp46i-tmp50i;
	tmp69r = -(tmp47i+tmp64r);
	tmp54i = tmp46i+tmp51i;
	tmp70r = tmp47i+tmp65r;
	tmp55i = tmp46i-tmp51i;
	tmp71r = tmp47i-tmp65r;
	tmp72r = tmp66r;
	tmp73r = C2*tmp68r-C3*tmp52i;
	tmp57i = C2*tmp52i+C3*tmp68r;
	tmp74r = C4*tmp59r-C5*tmp63r;
	tmp58i = C4*tmp63r+C5*tmp59r;
	tmp75r = C6*tmp70r-C7*tmp54i;
	tmp59i = C6*tmp54i+C7*tmp70r;
	tmp76r = C8*tmp67r;
	tmp60i = C8*tmp67r;
	tmp77r = C7*tmp69r-C6*tmp53i;
	tmp61i = C7*tmp53i+C6*tmp69r;
	tmp78r = C5*tmp59r+C4*tmp63r;
	tmp62i = -C5*tmp63r+C4*tmp59r;
	tmp79r = C3*tmp71r-C2*tmp55i;
	tmp63i = C3*tmp55i+C2*tmp71r;
	tmp80r = C9*Re[3];
	tmp64i = C9*Im[3];
	tmp81r = C9*Re[11];
	tmp65i = C9*Im[11];
	tmp82r = tmp80r+tmp81r;
	tmp83r = tmp80r-tmp81r;
	tmp84r = C9*Re[7];
	tmp66i = C9*Im[7];
	tmp85r = C9*Re[15];
	tmp67i = C9*Im[15];
	tmp86r = tmp84r+tmp85r;
	tmp87r = tmp84r-tmp85r;
	tmp68i = C8*(tmp66i-tmp67i);
	tmp88r = -C8*(tmp67i+tmp66i);
	tmp69i = -C8*(tmp66i-tmp67i);
	tmp89r = -C8*(tmp67i+tmp66i);
	tmp90r = tmp82r+tmp86r;
	tmp91r = tmp82r-tmp86r;
	tmp70i = tmp64i+tmp68i;
	tmp92r = -(tmp65i-tmp88r);
	tmp71i = tmp64i-tmp68i;
	tmp93r = -(tmp65i+tmp88r);
	tmp72i = tmp64i+tmp69i;
	tmp94r = tmp65i+tmp89r;
	tmp73i = tmp64i-tmp69i;
	tmp95r = tmp65i-tmp89r;
	tmp96r = tmp90r;
	tmp97r = C6*tmp92r-C7*tmp70i;
	tmp75i = C6*tmp70i+C7*tmp92r;
	tmp98r = C5*tmp83r-C4*tmp87r;
	tmp76i = C5*tmp87r+C4*tmp83r;
	tmp99r = -C3*tmp94r-C2*tmp72i;
	tmp77i = -C3*tmp72i+C2*tmp94r;
	tmp100r = -C8*tmp91r;
	tmp78i = C8*tmp91r;
	tmp101r = -C2*tmp93r-C3*tmp71i;
	tmp79i = -C2*tmp71i+C3*tmp93r;
	tmp102r = -C4*tmp83r-C5*tmp87r;
	tmp80i = C4*tmp87r-C5*tmp83r;
	tmp103r = -C7*tmp95r+C6*tmp73i;
	tmp81i = -C7*tmp73i-C6*tmp95r;
	tmp104r = tmp72r+tmp96r;
	tmp105r = tmp73r+tmp97r;
	tmp83i = tmp57i+tmp75i;
	tmp106r = tmp74r+tmp98r;
	tmp84i = tmp58i+tmp76i;
	tmp107r = tmp75r+tmp99r;
	tmp85i = tmp59i+tmp77i;
	tmp108r = tmp76r+tmp100r;
	tmp86i = tmp60i+tmp78i;
	tmp109r = tmp77r+tmp101r;
	tmp87i = tmp61i+tmp79i;
	tmp110r = tmp78r+tmp102r;
	tmp88i = tmp62i+tmp80i;
	tmp111r = tmp79r+tmp103r;
	tmp89i = tmp63i+tmp81i;
	tmp112r = C0;
	tmp90i = -(tmp72r-tmp96r);
	tmp113r = tmp57i-tmp75i;
	tmp91i = -(tmp73r-tmp97r);
	tmp114r = tmp58i-tmp76i;
	tmp92i = -(tmp74r-tmp98r);
	tmp115r = tmp59i-tmp77i;
	tmp93i = -(tmp75r-tmp99r);
	tmp116r = tmp60i-tmp78i;
	tmp94i = -(tmp76r-tmp100r);
	tmp117r = tmp61i-tmp79i;
	tmp95i = -(tmp77r-tmp101r);
	tmp118r = tmp62i-tmp80i;
	tmp96i = -(tmp78r-tmp102r);
	tmp119r = tmp63i-tmp81i;
	tmp97i = -(tmp79r-tmp103r);
	tmp120r = tmp40r+tmp104r;
	tmp121r = tmp41r+tmp105r;
	tmp99i = tmp31i+tmp83i;
	tmp122r = tmp42r+tmp106r;
	tmp100i = tmp32i+tmp84i;
	tmp123r = tmp43r+tmp107r;
	tmp101i = tmp33i+tmp85i;
	tmp124r = tmp44r+tmp108r;
	tmp102i = tmp34i+tmp86i;
	tmp125r = tmp45r+tmp109r;
	tmp103i = tmp35i+tmp87i;
	tmp126r = tmp46r+tmp110r;
	tmp104i = tmp36i+tmp88i;
	tmp127r = tmp47r+tmp111r;
	tmp105i = tmp37i+tmp89i;
	tmp128r = tmp48r-tmp112r;
	tmp106i = -tmp90i;
	tmp129r = tmp49r-tmp113r;
	tmp107i = tmp39i-tmp91i;
	tmp130r = tmp50r-tmp114r;
	tmp108i = tmp40i-tmp92i;
	tmp131r = tmp51r-tmp115r;
	tmp109i = tmp41i-tmp93i;
	tmp132r = tmp52r-tmp116r;
	tmp110i = tmp42i-tmp94i;
	tmp133r = tmp53r-tmp117r;
	tmp111i = tmp43i-tmp95i;
	tmp134r = tmp54r-tmp118r;
	tmp112i = tmp44i-tmp96i;
	tmp135r = tmp55r-tmp119r;
	tmp113i = tmp45i-tmp97i;
	Re[0] = tmp120r;
	Im[0] = 0.0;
	Re[1] = tmp121r;
	Im[1] = tmp99i;
	Re[2] = tmp122r;
	Im[2] = tmp100i;
	Re[3] = tmp123r;
	Im[3] = tmp101i;
	Re[4] = tmp124r;
	Im[4] = tmp102i;
	Re[5] = tmp125r;
	Im[5] = tmp103i;
	Re[6] = tmp126r;
	Im[6] = tmp104i;
	Re[7] = tmp127r;
	Im[7] = tmp105i;
	Re[8] = tmp128r;
	Im[8] = tmp106i;
	Re[9] = tmp129r;
	Im[9] = tmp107i;
	Re[10] = tmp130r;
	Im[10] = tmp108i;
	Re[11] = tmp131r;
	Im[11] = tmp109i;
	Re[12] = tmp132r;
	Im[12] = tmp110i;
	Re[13] = tmp133r;
	Im[13] = tmp111i;
	Re[14] = tmp134r;
	Im[14] = tmp112i;
	Re[15] = tmp135r;
	Im[15] = tmp113i;
}
