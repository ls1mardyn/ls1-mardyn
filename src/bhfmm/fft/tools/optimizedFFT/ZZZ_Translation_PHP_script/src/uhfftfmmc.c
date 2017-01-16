/*
*	UHFFT Library of FFT modules for FMM
*/


#include "uhfftfmmc.h"


/*
*	Number of additions = 72
*	Number of multiplications = 16
*	Number of sign changes = 16
*	Number of assigns = 112
*	Total number of operations = 216
*/
void	MFFTC12(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i;

	const REAL	C2 =                  0.5;	/* REALCONST	*/
	const REAL	C5 =     0.86602540378444;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[18]);
	tmp0r = Re(x[0])+Re(x[18]);
	tmp1i = Im(x[0])-Im(x[18]);
	tmp1r = Re(x[0])-Re(x[18]);
	tmp2i = Im(x[0])+Re(x[18]);
	tmp2r = Re(x[0])-Im(x[18]);
	tmp3i = Im(x[0])-Re(x[18]);
	tmp3r = Re(x[0])+Im(x[18]);
	tmp4i = Im(x[24])+Im(x[6]);
	tmp4r = Re(x[24])+Re(x[6]);
	tmp5i = Im(x[24])-Im(x[6]);
	tmp5r = Re(x[24])-Re(x[6]);
	tmp6i = Im(x[24])-Re(x[6]);
	tmp6r = Re(x[24])+Im(x[6]);
	tmp7i = Im(x[24])+Re(x[6]);
	tmp7r = Re(x[24])-Im(x[6]);
	tmp8i = Im(x[12])+Im(x[30]);
	tmp8r = Re(x[12])+Re(x[30]);
	tmp9i = Im(x[12])-Im(x[30]);
	tmp9r = Re(x[12])-Re(x[30]);
	tmp10i = -(Im(x[12])+Re(x[30]));
	tmp10r = -(Re(x[12])-Im(x[30]));
	tmp11i = -(Im(x[12])-Re(x[30]));
	tmp11r = -(Re(x[12])+Im(x[30]));
	tmp12i = tmp4i+tmp8i;
	tmp12r = tmp4r+tmp8r;
	tmp13i = tmp4i-tmp8i;
	tmp13r = tmp4r-tmp8r;
	tmp14i = -C2*tmp12i;
	tmp14r = -C2*tmp12r;
	tmp15i = -C5*tmp13r;
	tmp15r = C5*tmp13i;
	tmp16i = tmp14i+tmp0i;
	tmp16r = tmp14r+tmp0r;
	tmp17i = tmp16i+tmp15i;
	tmp17r = tmp16r+tmp15r;
	tmp18i = tmp16i-tmp15i;
	tmp18r = tmp16r-tmp15r;
	tmp19i = tmp0i+tmp12i;
	tmp19r = tmp0r+tmp12r;
	tmp20i = tmp7i+tmp11i;
	tmp20r = tmp7r+tmp11r;
	tmp21i = tmp7i-tmp11i;
	tmp21r = tmp7r-tmp11r;
	tmp22i = -C2*tmp20i;
	tmp22r = -C2*tmp20r;
	tmp23i = -C5*tmp21r;
	tmp23r = C5*tmp21i;
	tmp24i = tmp22i+tmp3i;
	tmp24r = tmp22r+tmp3r;
	tmp25i = tmp24i+tmp23i;
	tmp25r = tmp24r+tmp23r;
	tmp26i = tmp24i-tmp23i;
	tmp26r = tmp24r-tmp23r;
	tmp27i = tmp3i+tmp20i;
	tmp27r = tmp3r+tmp20r;
	tmp28i = tmp5i+tmp9i;
	tmp28r = tmp5r+tmp9r;
	tmp29i = tmp5i-tmp9i;
	tmp29r = tmp5r-tmp9r;
	tmp30i = -C2*tmp28i;
	tmp30r = -C2*tmp28r;
	tmp31i = -C5*tmp29r;
	tmp31r = C5*tmp29i;
	tmp32i = tmp30i+tmp1i;
	tmp32r = tmp30r+tmp1r;
	tmp33i = tmp32i+tmp31i;
	tmp33r = tmp32r+tmp31r;
	tmp34i = tmp32i-tmp31i;
	tmp34r = tmp32r-tmp31r;
	tmp35i = tmp1i+tmp28i;
	tmp35r = tmp1r+tmp28r;
	tmp36i = tmp6i+tmp10i;
	tmp36r = tmp6r+tmp10r;
	tmp37i = tmp6i-tmp10i;
	tmp37r = tmp6r-tmp10r;
	tmp38i = -C2*tmp36i;
	tmp38r = -C2*tmp36r;
	tmp39i = -C5*tmp37r;
	tmp39r = C5*tmp37i;
	tmp40i = tmp38i+tmp2i;
	tmp40r = tmp38r+tmp2r;
	tmp41i = tmp40i+tmp39i;
	tmp41r = tmp40r+tmp39r;
	tmp42i = tmp40i-tmp39i;
	tmp42r = tmp40r-tmp39r;
	tmp43i = tmp2i+tmp36i;
	tmp43r = tmp2r+tmp36r;
	Re(x[0]) = tmp19r;
	Im(x[0]) = tmp19i;
	Re(x[6]) = tmp25r;
	Im(x[6]) = tmp25i;
	Re(x[12]) = tmp34r;
	Im(x[12]) = tmp34i;
	Re(x[18]) = tmp43r;
	Im(x[18]) = tmp43i;
	Re(x[24]) = tmp17r;
	Im(x[24]) = tmp17i;
	Re(x[30]) = tmp26r;
	Im(x[30]) = tmp26i;
	Re(x[36]) = tmp35r;
	Im(x[36]) = tmp35i;
	Re(x[42]) = tmp41r;
	Im(x[42]) = tmp41i;
	Re(x[48]) = tmp18r;
	Im(x[48]) = tmp18i;
	Re(x[54]) = tmp27r;
	Im(x[54]) = tmp27i;
	Re(x[60]) = tmp33r;
	Im(x[60]) = tmp33i;
	Re(x[66]) = tmp42r;
	Im(x[66]) = tmp42i;
}

/*
*	Number of additions = 96
*	Number of multiplications = 16
*	Number of sign changes = 12
*	Number of assigns = 136
*	Total number of operations = 260
*/
void	MIFFTC12(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i;

	const REAL	C2 =                  0.5;	/* REALCONST	*/
	const REAL	C5 =     0.86602540378444;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[36]);
	tmp0r = Re(x[0])+Re(x[36]);
	tmp1i = Im(x[0])-Im(x[36]);
	tmp1r = Re(x[0])-Re(x[36]);
	tmp2i = Im(x[54])+Im(x[18]);
	tmp2r = Re(x[54])+Re(x[18]);
	tmp3i = Im(x[54])-Im(x[18]);
	tmp3r = Re(x[54])-Re(x[18]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = tmp1i+tmp3r;
	tmp6r = tmp1r-tmp3i;
	tmp7i = tmp1i-tmp3r;
	tmp7r = tmp1r+tmp3i;
	tmp8i = Im(x[24])+Im(x[60]);
	tmp8r = Re(x[24])+Re(x[60]);
	tmp9i = Im(x[24])-Im(x[60]);
	tmp9r = Re(x[24])-Re(x[60]);
	tmp10i = Im(x[6])+Im(x[42]);
	tmp10r = Re(x[6])+Re(x[42]);
	tmp11i = Im(x[6])-Im(x[42]);
	tmp11r = Re(x[6])-Re(x[42]);
	tmp12i = tmp8i+tmp10i;
	tmp12r = tmp8r+tmp10r;
	tmp13i = tmp8i-tmp10i;
	tmp13r = tmp8r-tmp10r;
	tmp14i = tmp9i+tmp11r;
	tmp14r = tmp9r-tmp11i;
	tmp15i = tmp9i-tmp11r;
	tmp15r = tmp9r+tmp11i;
	tmp16i = Im(x[48])+Im(x[12]);
	tmp16r = Re(x[48])+Re(x[12]);
	tmp17i = Im(x[48])-Im(x[12]);
	tmp17r = Re(x[48])-Re(x[12]);
	tmp18i = Im(x[30])+Im(x[66]);
	tmp18r = Re(x[30])+Re(x[66]);
	tmp19i = Im(x[30])-Im(x[66]);
	tmp19r = Re(x[30])-Re(x[66]);
	tmp20i = tmp16i+tmp18i;
	tmp20r = tmp16r+tmp18r;
	tmp21i = tmp16i-tmp18i;
	tmp21r = tmp16r-tmp18r;
	tmp22i = tmp17i+tmp19r;
	tmp22r = tmp17r-tmp19i;
	tmp23i = tmp17i-tmp19r;
	tmp23r = tmp17r+tmp19i;
	tmp24i = tmp12i+tmp20i;
	tmp24r = tmp12r+tmp20r;
	tmp25i = tmp12i-tmp20i;
	tmp25r = tmp12r-tmp20r;
	tmp26i = -C2*tmp24i;
	tmp26r = -C2*tmp24r;
	tmp27i = C5*tmp25r;
	tmp27r = -C5*tmp25i;
	tmp28i = tmp26i+tmp4i;
	tmp28r = tmp26r+tmp4r;
	tmp29i = tmp28i+tmp27i;
	tmp29r = tmp28r+tmp27r;
	tmp30i = tmp28i-tmp27i;
	tmp30r = tmp28r-tmp27r;
	tmp31i = tmp4i+tmp24i;
	tmp31r = tmp4r+tmp24r;
	tmp32i = tmp15i+tmp23i;
	tmp32r = tmp15r+tmp23r;
	tmp33i = tmp15i-tmp23i;
	tmp33r = tmp15r-tmp23r;
	tmp34i = -C2*tmp32i;
	tmp34r = -C2*tmp32r;
	tmp35i = C5*tmp33r;
	tmp35r = -C5*tmp33i;
	tmp36i = tmp34i+tmp7i;
	tmp36r = tmp34r+tmp7r;
	tmp37i = tmp36i+tmp35i;
	tmp37r = tmp36r+tmp35r;
	tmp38i = tmp36i-tmp35i;
	tmp38r = tmp36r-tmp35r;
	tmp39i = tmp7i+tmp32i;
	tmp39r = tmp7r+tmp32r;
	tmp40i = tmp13i+tmp21i;
	tmp40r = tmp13r+tmp21r;
	tmp41i = tmp13i-tmp21i;
	tmp41r = tmp13r-tmp21r;
	tmp42i = -C2*tmp40i;
	tmp42r = -C2*tmp40r;
	tmp43i = C5*tmp41r;
	tmp43r = -C5*tmp41i;
	tmp44i = tmp42i+tmp5i;
	tmp44r = tmp42r+tmp5r;
	tmp45i = tmp44i+tmp43i;
	tmp45r = tmp44r+tmp43r;
	tmp46i = tmp44i-tmp43i;
	tmp46r = tmp44r-tmp43r;
	tmp47i = tmp5i+tmp40i;
	tmp47r = tmp5r+tmp40r;
	tmp48i = tmp14i+tmp22i;
	tmp48r = tmp14r+tmp22r;
	tmp49i = tmp14i-tmp22i;
	tmp49r = tmp14r-tmp22r;
	tmp50i = -C2*tmp48i;
	tmp50r = -C2*tmp48r;
	tmp51i = C5*tmp49r;
	tmp51r = -C5*tmp49i;
	tmp52i = tmp50i+tmp6i;
	tmp52r = tmp50r+tmp6r;
	tmp53i = tmp52i+tmp51i;
	tmp53r = tmp52r+tmp51r;
	tmp54i = tmp52i-tmp51i;
	tmp54r = tmp52r-tmp51r;
	tmp55i = tmp6i+tmp48i;
	tmp55r = tmp6r+tmp48r;
	Re(x[0]) = tmp31r;
	Im(x[0]) = tmp31i;
	Re(x[6]) = tmp37r;
	Im(x[6]) = tmp37i;
	Re(x[12]) = tmp46r;
	Im(x[12]) = tmp46i;
	Re(x[18]) = tmp55r;
	Im(x[18]) = tmp55i;
	Re(x[24]) = tmp29r;
	Im(x[24]) = tmp29i;
	Re(x[30]) = tmp38r;
	Im(x[30]) = tmp38i;
	Re(x[36]) = tmp47r;
	Im(x[36]) = tmp47i;
	Re(x[42]) = tmp53r;
	Im(x[42]) = tmp53i;
	Re(x[48]) = tmp30r;
	Im(x[48]) = tmp30i;
	Re(x[54]) = tmp39r;
	Im(x[54]) = tmp39i;
	Re(x[60]) = tmp45r;
	Im(x[60]) = tmp45i;
	Re(x[66]) = tmp54r;
	Im(x[66]) = tmp54i;
}

/*
*	Number of additions = 116
*	Number of multiplications = 72
*	Number of sign changes = 0
*	Number of assigns = 114
*	Total number of operations = 302
*/
void	MFFTC14(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i;

	const REAL	C13 =     0.22252093395631;	/* REALCONST	*/
	const REAL	C12 =     0.43388373911756;	/* REALCONST	*/
	const REAL	C9 =     0.62348980185873;	/* REALCONST	*/
	const REAL	C10 =     0.78183148246803;	/* REALCONST	*/
	const REAL	C11 =     0.90096886790242;	/* REALCONST	*/
	const REAL	C14 =     0.97492791218182;	/* REALCONST	*/

	tmp0i = C9*Im(x[42])-C11*Im(x[28])-C13*Im(x[14])+Im(x[0]);
	tmp0r = C9*Re(x[42])-C11*Re(x[28])-C13*Re(x[14])+Re(x[0]);
	tmp1i = Im(x[42])+Im(x[28]);
	tmp1r = Re(x[42])+Re(x[28]);
	tmp2i = -C11*Im(x[42])-C13*Im(x[28])+C9*Im(x[14])+Im(x[0]);
	tmp2r = -C11*Re(x[42])-C13*Re(x[28])+C9*Re(x[14])+Re(x[0]);
	tmp3i = tmp1i+Im(x[14]);
	tmp3r = tmp1r+Re(x[14]);
	tmp4i = -C13*Im(x[42])+C9*Im(x[28])-C11*Im(x[14])+Im(x[0]);
	tmp4r = -C13*Re(x[42])+C9*Re(x[28])-C11*Re(x[14])+Re(x[0]);
	tmp5i = C10*Re(x[42])+C12*Re(x[28])-C14*Re(x[14]);
	tmp5r = -C10*Im(x[42])-C12*Im(x[28])+C14*Im(x[14]);
	tmp6i = C12*Re(x[42])+C14*Re(x[28])+C10*Re(x[14]);
	tmp6r = -C12*Im(x[42])-C14*Im(x[28])-C10*Im(x[14]);
	tmp7i = C14*Re(x[42])-C10*Re(x[28])+C12*Re(x[14]);
	tmp7r = -C14*Im(x[42])+C10*Im(x[28])-C12*Im(x[14]);
	tmp8i = tmp0i+tmp5i;
	tmp8r = tmp0r+tmp5r;
	tmp9i = tmp0i-tmp5i;
	tmp9r = tmp0r-tmp5r;
	tmp10i = tmp2i+tmp6i;
	tmp10r = tmp2r+tmp6r;
	tmp11i = tmp2i-tmp6i;
	tmp11r = tmp2r-tmp6r;
	tmp12i = tmp4i+tmp7i;
	tmp12r = tmp4r+tmp7r;
	tmp13i = tmp4i-tmp7i;
	tmp13r = tmp4r-tmp7r;
	tmp14i = Im(x[0])+tmp3i;
	tmp14r = Re(x[0])+tmp3r;
	tmp15i = C9*Im(x[7])-C11*Im(x[21])-C13*Im(x[35]);
	tmp15r = C9*Re(x[7])-C11*Re(x[21])-C13*Re(x[35]);
	tmp16i = Im(x[7])+Im(x[21]);
	tmp16r = Re(x[7])+Re(x[21]);
	tmp17i = -C11*Im(x[7])-C13*Im(x[21])+C9*Im(x[35]);
	tmp17r = -C11*Re(x[7])-C13*Re(x[21])+C9*Re(x[35]);
	tmp18i = tmp16i+Im(x[35]);
	tmp18r = tmp16r+Re(x[35]);
	tmp19i = -C13*Im(x[7])+C9*Im(x[21])-C11*Im(x[35]);
	tmp19r = -C13*Re(x[7])+C9*Re(x[21])-C11*Re(x[35]);
	tmp20i = -C10*Re(x[7])-C12*Re(x[21])+C14*Re(x[35]);
	tmp20r = C10*Im(x[7])+C12*Im(x[21])-C14*Im(x[35]);
	tmp21i = -C12*Re(x[7])-C14*Re(x[21])-C10*Re(x[35]);
	tmp21r = C12*Im(x[7])+C14*Im(x[21])+C10*Im(x[35]);
	tmp22i = -C14*Re(x[7])+C10*Re(x[21])-C12*Re(x[35]);
	tmp22r = C14*Im(x[7])-C10*Im(x[21])+C12*Im(x[35]);
	tmp23i = tmp15i+tmp20i;
	tmp23r = tmp15r+tmp20r;
	tmp24i = tmp15i-tmp20i;
	tmp24r = tmp15r-tmp20r;
	tmp25i = tmp17i+tmp21i;
	tmp25r = tmp17r+tmp21r;
	tmp26i = tmp17i-tmp21i;
	tmp26r = tmp17r-tmp21r;
	tmp27i = tmp19i+tmp22i;
	tmp27r = tmp19r+tmp22r;
	tmp28i = tmp19i-tmp22i;
	tmp28r = tmp19r-tmp22r;
	tmp29i = tmp14i+tmp18i;
	tmp29r = tmp14r+tmp18r;
	tmp30i = tmp14i-tmp18i;
	tmp30r = tmp14r-tmp18r;
	tmp31i = tmp11i+tmp26i;
	tmp31r = tmp11r+tmp26r;
	tmp32i = tmp11i-tmp26i;
	tmp32r = tmp11r-tmp26r;
	tmp33i = tmp8i+tmp23i;
	tmp33r = tmp8r+tmp23r;
	tmp34i = tmp8i-tmp23i;
	tmp34r = tmp8r-tmp23r;
	tmp35i = tmp13i+tmp28i;
	tmp35r = tmp13r+tmp28r;
	tmp36i = tmp13i-tmp28i;
	tmp36r = tmp13r-tmp28r;
	tmp37i = tmp12i+tmp27i;
	tmp37r = tmp12r+tmp27r;
	tmp38i = tmp12i-tmp27i;
	tmp38r = tmp12r-tmp27r;
	tmp39i = tmp9i+tmp24i;
	tmp39r = tmp9r+tmp24r;
	tmp40i = tmp9i-tmp24i;
	tmp40r = tmp9r-tmp24r;
	tmp41i = tmp10i+tmp25i;
	tmp41r = tmp10r+tmp25r;
	tmp42i = tmp10i-tmp25i;
	tmp42r = tmp10r-tmp25r;
	Re(x[0]) = tmp29r;
	Im(x[0]) = tmp29i;
	Re(x[7]) = tmp32r;
	Im(x[7]) = tmp32i;
	Re(x[14]) = tmp33r;
	Im(x[14]) = tmp33i;
	Re(x[21]) = tmp36r;
	Im(x[21]) = tmp36i;
	Re(x[28]) = tmp37r;
	Im(x[28]) = tmp37i;
	Re(x[35]) = tmp40r;
	Im(x[35]) = tmp40i;
	Re(x[42]) = tmp41r;
	Im(x[42]) = tmp41i;
	Re(x[49]) = tmp30r;
	Im(x[49]) = tmp30i;
	Re(x[56]) = tmp31r;
	Im(x[56]) = tmp31i;
	Re(x[63]) = tmp34r;
	Im(x[63]) = tmp34i;
	Re(x[70]) = tmp35r;
	Im(x[70]) = tmp35i;
	Re(x[77]) = tmp38r;
	Im(x[77]) = tmp38i;
	Re(x[84]) = tmp39r;
	Im(x[84]) = tmp39i;
	Re(x[91]) = tmp42r;
	Im(x[91]) = tmp42i;
}

/*
*	Number of additions = 148
*	Number of multiplications = 72
*	Number of sign changes = 0
*	Number of assigns = 140
*	Total number of operations = 360
*/
void	MIFFTC14(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i;

	const REAL	C13 =     0.22252093395631;	/* REALCONST	*/
	const REAL	C12 =     0.43388373911756;	/* REALCONST	*/
	const REAL	C9 =     0.62348980185873;	/* REALCONST	*/
	const REAL	C10 =     0.78183148246803;	/* REALCONST	*/
	const REAL	C11 =     0.90096886790242;	/* REALCONST	*/
	const REAL	C14 =     0.97492791218182;	/* REALCONST	*/

	tmp0i = Im(x[56])+Im(x[42]);
	tmp0r = Re(x[56])+Re(x[42]);
	tmp1i = Im(x[56])-Im(x[42]);
	tmp1r = Re(x[56])-Re(x[42]);
	tmp2i = Im(x[70])+Im(x[28]);
	tmp2r = Re(x[70])+Re(x[28]);
	tmp3i = Im(x[70])-Im(x[28]);
	tmp3r = Re(x[70])-Re(x[28]);
	tmp4i = Im(x[14])+Im(x[84]);
	tmp4r = Re(x[14])+Re(x[84]);
	tmp5i = Im(x[14])-Im(x[84]);
	tmp5r = Re(x[14])-Re(x[84]);
	tmp6i = C9*tmp0i-C11*tmp2i-C13*tmp4i+Im(x[0]);
	tmp6r = C9*tmp0r-C11*tmp2r-C13*tmp4r+Re(x[0]);
	tmp7i = tmp0i+tmp2i;
	tmp7r = tmp0r+tmp2r;
	tmp8i = -C11*tmp0i-C13*tmp2i+C9*tmp4i+Im(x[0]);
	tmp8r = -C11*tmp0r-C13*tmp2r+C9*tmp4r+Re(x[0]);
	tmp9i = tmp7i+tmp4i;
	tmp9r = tmp7r+tmp4r;
	tmp10i = -C13*tmp0i+C9*tmp2i-C11*tmp4i+Im(x[0]);
	tmp10r = -C13*tmp0r+C9*tmp2r-C11*tmp4r+Re(x[0]);
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
	tmp20i = Im(x[0])+tmp9i;
	tmp20r = Re(x[0])+tmp9r;
	tmp21i = Im(x[7])+Im(x[91]);
	tmp21r = Re(x[7])+Re(x[91]);
	tmp22i = Im(x[7])-Im(x[91]);
	tmp22r = Re(x[7])-Re(x[91]);
	tmp23i = Im(x[21])+Im(x[77]);
	tmp23r = Re(x[21])+Re(x[77]);
	tmp24i = Im(x[21])-Im(x[77]);
	tmp24r = Re(x[21])-Re(x[77]);
	tmp25i = Im(x[63])+Im(x[35]);
	tmp25r = Re(x[63])+Re(x[35]);
	tmp26i = Im(x[63])-Im(x[35]);
	tmp26r = Re(x[63])-Re(x[35]);
	tmp27i = C9*tmp21i-C11*tmp23i-C13*tmp25i+Im(x[49]);
	tmp27r = C9*tmp21r-C11*tmp23r-C13*tmp25r+Re(x[49]);
	tmp28i = tmp21i+tmp23i;
	tmp28r = tmp21r+tmp23r;
	tmp29i = -C11*tmp21i-C13*tmp23i+C9*tmp25i+Im(x[49]);
	tmp29r = -C11*tmp21r-C13*tmp23r+C9*tmp25r+Re(x[49]);
	tmp30i = tmp28i+tmp25i;
	tmp30r = tmp28r+tmp25r;
	tmp31i = -C13*tmp21i+C9*tmp23i-C11*tmp25i+Im(x[49]);
	tmp31r = -C13*tmp21r+C9*tmp23r-C11*tmp25r+Re(x[49]);
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
	tmp41i = Im(x[49])+tmp30i;
	tmp41r = Re(x[49])+tmp30r;
	tmp42i = tmp20i+tmp41i;
	tmp42r = tmp20r+tmp41r;
	tmp43i = tmp20i-tmp41i;
	tmp43r = tmp20r-tmp41r;
	tmp44i = tmp17i+tmp38i;
	tmp44r = tmp17r+tmp38r;
	tmp45i = tmp17i-tmp38i;
	tmp45r = tmp17r-tmp38r;
	tmp46i = tmp14i+tmp35i;
	tmp46r = tmp14r+tmp35r;
	tmp47i = tmp14i-tmp35i;
	tmp47r = tmp14r-tmp35r;
	tmp48i = tmp19i+tmp40i;
	tmp48r = tmp19r+tmp40r;
	tmp49i = tmp19i-tmp40i;
	tmp49r = tmp19r-tmp40r;
	tmp50i = tmp18i+tmp39i;
	tmp50r = tmp18r+tmp39r;
	tmp51i = tmp18i-tmp39i;
	tmp51r = tmp18r-tmp39r;
	tmp52i = tmp15i+tmp36i;
	tmp52r = tmp15r+tmp36r;
	tmp53i = tmp15i-tmp36i;
	tmp53r = tmp15r-tmp36r;
	tmp54i = tmp16i+tmp37i;
	tmp54r = tmp16r+tmp37r;
	tmp55i = tmp16i-tmp37i;
	tmp55r = tmp16r-tmp37r;
	Re(x[0]) = tmp42r;
	Im(x[0]) = tmp42i;
	Re(x[7]) = tmp45r;
	Im(x[7]) = tmp45i;
	Re(x[14]) = tmp46r;
	Im(x[14]) = tmp46i;
	Re(x[21]) = tmp49r;
	Im(x[21]) = tmp49i;
	Re(x[28]) = tmp50r;
	Im(x[28]) = tmp50i;
	Re(x[35]) = tmp53r;
	Im(x[35]) = tmp53i;
	Re(x[42]) = tmp54r;
	Im(x[42]) = tmp54i;
	Re(x[49]) = tmp43r;
	Im(x[49]) = tmp43i;
	Re(x[56]) = tmp44r;
	Im(x[56]) = tmp44i;
	Re(x[63]) = tmp47r;
	Im(x[63]) = tmp47i;
	Re(x[70]) = tmp48r;
	Im(x[70]) = tmp48i;
	Re(x[77]) = tmp51r;
	Im(x[77]) = tmp51i;
	Re(x[84]) = tmp52r;
	Im(x[84]) = tmp52i;
	Re(x[91]) = tmp55r;
	Im(x[91]) = tmp55i;
}

/*
*	Number of additions = 112
*	Number of multiplications = 24
*	Number of sign changes = 4
*	Number of assigns = 148
*	Total number of operations = 288
*/
void	MFFTC16(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i;

	const REAL	C3 =     0.38268343236509;	/* REALCONST	*/
	const REAL	C4 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C2 =     0.92387953251129;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[32]);
	tmp0r = Re(x[0])+Re(x[32]);
	tmp1i = Im(x[0])-Im(x[32]);
	tmp1r = Re(x[0])-Re(x[32]);
	tmp2i = Im(x[0])-Re(x[32]);
	tmp2r = Re(x[0])+Im(x[32]);
	tmp3i = Im(x[0])+Re(x[32]);
	tmp3r = Re(x[0])-Im(x[32]);
	tmp4i = Im(x[16])+Im(x[48]);
	tmp4r = Re(x[16])+Re(x[48]);
	tmp5i = Im(x[16])-Im(x[48]);
	tmp5r = Re(x[16])-Re(x[48]);
	tmp6i = Im(x[16])-Re(x[48]);
	tmp6r = Re(x[16])+Im(x[48]);
	tmp7i = Im(x[16])+Re(x[48]);
	tmp7r = Re(x[16])-Im(x[48]);
	tmp8i = C4*(tmp6i-tmp6r);
	tmp8r = C4*(tmp6r+tmp6i);
	tmp9i = -C4*(tmp7i+tmp7r);
	tmp9r = -C4*(tmp7r-tmp7i);
	tmp10i = tmp0i+tmp4i;
	tmp10r = tmp0r+tmp4r;
	tmp11i = tmp0i-tmp4i;
	tmp11r = tmp0r-tmp4r;
	tmp12i = tmp2i+tmp8i;
	tmp12r = tmp2r+tmp8r;
	tmp13i = tmp2i-tmp8i;
	tmp13r = tmp2r-tmp8r;
	tmp14i = tmp1i-tmp5r;
	tmp14r = tmp1r+tmp5i;
	tmp15i = tmp1i+tmp5r;
	tmp15r = tmp1r-tmp5i;
	tmp16i = tmp3i+tmp9i;
	tmp16r = tmp3r+tmp9r;
	tmp17i = tmp3i-tmp9i;
	tmp17r = tmp3r-tmp9r;
	tmp18i = Im(x[8])+Im(x[40]);
	tmp18r = Re(x[8])+Re(x[40]);
	tmp19i = Im(x[8])-Im(x[40]);
	tmp19r = Re(x[8])-Re(x[40]);
	tmp20i = Im(x[8])-Re(x[40]);
	tmp20r = Re(x[8])+Im(x[40]);
	tmp21i = Im(x[8])+Re(x[40]);
	tmp21r = Re(x[8])-Im(x[40]);
	tmp22r = tmp18r;
	tmp22i = tmp18i;
	tmp23r = C2*tmp20r+C3*tmp20i;
	tmp23i = C2*tmp20i-C3*tmp20r;
	tmp24r = C4*(tmp19r+tmp19i);
	tmp24i = C4*(tmp19i-tmp19r);
	tmp25r = C3*tmp21r+C2*tmp21i;
	tmp25i = C3*tmp21i-C2*tmp21r;
	tmp26i = Im(x[24])+Im(x[56]);
	tmp26r = Re(x[24])+Re(x[56]);
	tmp27i = Im(x[24])-Im(x[56]);
	tmp27r = Re(x[24])-Re(x[56]);
	tmp28i = Im(x[24])-Re(x[56]);
	tmp28r = Re(x[24])+Im(x[56]);
	tmp29i = Im(x[24])+Re(x[56]);
	tmp29r = Re(x[24])-Im(x[56]);
	tmp30r = tmp26r;
	tmp30i = tmp26i;
	tmp31r = C3*tmp28r+C2*tmp28i;
	tmp31i = C3*tmp28i-C2*tmp28r;
	tmp32r = -C4*(tmp27r-tmp27i);
	tmp32i = -C4*(tmp27i+tmp27r);
	tmp33r = -C2*tmp29r-C3*tmp29i;
	tmp33i = -C2*tmp29i+C3*tmp29r;
	tmp34r = tmp22r+tmp30r;
	tmp34i = tmp22i+tmp30i;
	tmp35r = tmp23r+tmp31r;
	tmp35i = tmp23i+tmp31i;
	tmp36r = tmp24r+tmp32r;
	tmp36i = tmp24i+tmp32i;
	tmp37r = tmp25r+tmp33r;
	tmp37i = tmp25i+tmp33i;
	tmp38r = -tmp22i+tmp30i;
	tmp38i = tmp22r-tmp30r;
	tmp39r = -tmp23i+tmp31i;
	tmp39i = tmp23r-tmp31r;
	tmp40r = -tmp24i+tmp32i;
	tmp40i = tmp24r-tmp32r;
	tmp41r = -tmp25i+tmp33i;
	tmp41i = tmp25r-tmp33r;
	tmp42r = tmp10r+tmp34r;
	tmp42i = tmp10i+tmp34i;
	tmp43r = tmp12r+tmp35r;
	tmp43i = tmp12i+tmp35i;
	tmp44r = tmp14r+tmp36r;
	tmp44i = tmp14i+tmp36i;
	tmp45r = tmp16r+tmp37r;
	tmp45i = tmp16i+tmp37i;
	tmp46r = tmp11r-tmp38r;
	tmp46i = tmp11i-tmp38i;
	tmp47r = tmp13r-tmp39r;
	tmp47i = tmp13i-tmp39i;
	tmp48r = tmp15r-tmp40r;
	tmp48i = tmp15i-tmp40i;
	tmp49r = tmp17r-tmp41r;
	tmp49i = tmp17i-tmp41i;
	tmp50r = tmp10r-tmp34r;
	tmp50i = tmp10i-tmp34i;
	tmp51r = tmp12r-tmp35r;
	tmp51i = tmp12i-tmp35i;
	tmp52r = tmp14r-tmp36r;
	tmp52i = tmp14i-tmp36i;
	tmp53r = tmp16r-tmp37r;
	tmp53i = tmp16i-tmp37i;
	tmp54r = tmp11r+tmp38r;
	tmp54i = tmp11i+tmp38i;
	tmp55r = tmp13r+tmp39r;
	tmp55i = tmp13i+tmp39i;
	tmp56r = tmp15r+tmp40r;
	tmp56i = tmp15i+tmp40i;
	tmp57r = tmp17r+tmp41r;
	tmp57i = tmp17i+tmp41i;
	Re(x[0]) = tmp42r;
	Im(x[0]) = tmp42i;
	Re(x[8]) = tmp43r;
	Im(x[8]) = tmp43i;
	Re(x[16]) = tmp44r;
	Im(x[16]) = tmp44i;
	Re(x[24]) = tmp45r;
	Im(x[24]) = tmp45i;
	Re(x[32]) = tmp46r;
	Im(x[32]) = tmp46i;
	Re(x[40]) = tmp47r;
	Im(x[40]) = tmp47i;
	Re(x[48]) = tmp48r;
	Im(x[48]) = tmp48i;
	Re(x[56]) = tmp49r;
	Im(x[56]) = tmp49i;
	Re(x[64]) = tmp50r;
	Im(x[64]) = tmp50i;
	Re(x[72]) = tmp51r;
	Im(x[72]) = tmp51i;
	Re(x[80]) = tmp52r;
	Im(x[80]) = tmp52i;
	Re(x[88]) = tmp53r;
	Im(x[88]) = tmp53i;
	Re(x[96]) = tmp54r;
	Im(x[96]) = tmp54i;
	Re(x[104]) = tmp55r;
	Im(x[104]) = tmp55i;
	Re(x[112]) = tmp56r;
	Im(x[112]) = tmp56i;
	Re(x[120]) = tmp57r;
	Im(x[120]) = tmp57i;
}

/*
*	Number of additions = 144
*	Number of multiplications = 24
*	Number of sign changes = 8
*	Number of assigns = 180
*	Total number of operations = 356
*/
void	MIFFTC16(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r,
		tmp71r, tmp72r, tmp73r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i,
		tmp71i, tmp72i, tmp73i;

	const REAL	C3 =     0.38268343236509;	/* REALCONST	*/
	const REAL	C4 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C2 =     0.92387953251129;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[64]);
	tmp0r = Re(x[0])+Re(x[64]);
	tmp1i = Im(x[0])-Im(x[64]);
	tmp1r = Re(x[0])-Re(x[64]);
	tmp2i = Im(x[32])+Im(x[96]);
	tmp2r = Re(x[32])+Re(x[96]);
	tmp3i = Im(x[32])-Im(x[96]);
	tmp3r = Re(x[32])-Re(x[96]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = tmp1i+tmp3r;
	tmp6r = tmp1r-tmp3i;
	tmp7i = tmp1i-tmp3r;
	tmp7r = tmp1r+tmp3i;
	tmp8i = Im(x[16])+Im(x[80]);
	tmp8r = Re(x[16])+Re(x[80]);
	tmp9i = Im(x[16])-Im(x[80]);
	tmp9r = Re(x[16])-Re(x[80]);
	tmp10i = Im(x[48])+Im(x[112]);
	tmp10r = Re(x[48])+Re(x[112]);
	tmp11i = Im(x[48])-Im(x[112]);
	tmp11r = Re(x[48])-Re(x[112]);
	tmp12i = tmp8i+tmp10i;
	tmp12r = tmp8r+tmp10r;
	tmp13i = tmp8i-tmp10i;
	tmp13r = tmp8r-tmp10r;
	tmp14i = tmp9i+tmp11r;
	tmp14r = tmp9r-tmp11i;
	tmp15i = tmp9i-tmp11r;
	tmp15r = tmp9r+tmp11i;
	tmp16i = C4*(tmp14i+tmp14r);
	tmp16r = C4*(tmp14r-tmp14i);
	tmp17i = -C4*(tmp15i-tmp15r);
	tmp17r = -C4*(tmp15r+tmp15i);
	tmp18i = tmp4i+tmp12i;
	tmp18r = tmp4r+tmp12r;
	tmp19i = tmp4i-tmp12i;
	tmp19r = tmp4r-tmp12r;
	tmp20i = tmp6i+tmp16i;
	tmp20r = tmp6r+tmp16r;
	tmp21i = tmp6i-tmp16i;
	tmp21r = tmp6r-tmp16r;
	tmp22i = tmp5i+tmp13r;
	tmp22r = tmp5r-tmp13i;
	tmp23i = tmp5i-tmp13r;
	tmp23r = tmp5r+tmp13i;
	tmp24i = tmp7i+tmp17i;
	tmp24r = tmp7r+tmp17r;
	tmp25i = tmp7i-tmp17i;
	tmp25r = tmp7r-tmp17r;
	tmp26i = Im(x[8])+Im(x[72]);
	tmp26r = Re(x[8])+Re(x[72]);
	tmp27i = Im(x[8])-Im(x[72]);
	tmp27r = Re(x[8])-Re(x[72]);
	tmp28i = Im(x[40])+Im(x[104]);
	tmp28r = Re(x[40])+Re(x[104]);
	tmp29i = Im(x[40])-Im(x[104]);
	tmp29r = Re(x[40])-Re(x[104]);
	tmp30i = tmp26i+tmp28i;
	tmp30r = tmp26r+tmp28r;
	tmp31i = tmp26i-tmp28i;
	tmp31r = tmp26r-tmp28r;
	tmp32i = tmp27i+tmp29r;
	tmp32r = tmp27r-tmp29i;
	tmp33i = tmp27i-tmp29r;
	tmp33r = tmp27r+tmp29i;
	tmp34r = tmp30r;
	tmp34i = tmp30i;
	tmp35r = C2*tmp32r-C3*tmp32i;
	tmp35i = C2*tmp32i+C3*tmp32r;
	tmp36r = C4*(tmp31r-tmp31i);
	tmp36i = C4*(tmp31i+tmp31r);
	tmp37r = C3*tmp33r-C2*tmp33i;
	tmp37i = C3*tmp33i+C2*tmp33r;
	tmp38i = Im(x[24])+Im(x[88]);
	tmp38r = Re(x[24])+Re(x[88]);
	tmp39i = Im(x[24])-Im(x[88]);
	tmp39r = Re(x[24])-Re(x[88]);
	tmp40i = Im(x[56])+Im(x[120]);
	tmp40r = Re(x[56])+Re(x[120]);
	tmp41i = Im(x[56])-Im(x[120]);
	tmp41r = Re(x[56])-Re(x[120]);
	tmp42i = tmp38i+tmp40i;
	tmp42r = tmp38r+tmp40r;
	tmp43i = tmp38i-tmp40i;
	tmp43r = tmp38r-tmp40r;
	tmp44i = tmp39i+tmp41r;
	tmp44r = tmp39r-tmp41i;
	tmp45i = tmp39i-tmp41r;
	tmp45r = tmp39r+tmp41i;
	tmp46r = tmp42r;
	tmp46i = tmp42i;
	tmp47r = C3*tmp44r-C2*tmp44i;
	tmp47i = C3*tmp44i+C2*tmp44r;
	tmp48r = -C4*(tmp43r+tmp43i);
	tmp48i = -C4*(tmp43i-tmp43r);
	tmp49r = -C2*tmp45r+C3*tmp45i;
	tmp49i = -C2*tmp45i-C3*tmp45r;
	tmp50r = tmp34r+tmp46r;
	tmp50i = tmp34i+tmp46i;
	tmp51r = tmp35r+tmp47r;
	tmp51i = tmp35i+tmp47i;
	tmp52r = tmp36r+tmp48r;
	tmp52i = tmp36i+tmp48i;
	tmp53r = tmp37r+tmp49r;
	tmp53i = tmp37i+tmp49i;
	tmp54r = tmp34i-tmp46i;
	tmp54i = -(tmp34r-tmp46r);
	tmp55r = tmp35i-tmp47i;
	tmp55i = -(tmp35r-tmp47r);
	tmp56r = tmp36i-tmp48i;
	tmp56i = -(tmp36r-tmp48r);
	tmp57r = tmp37i-tmp49i;
	tmp57i = -(tmp37r-tmp49r);
	tmp58r = tmp18r+tmp50r;
	tmp58i = tmp18i+tmp50i;
	tmp59r = tmp20r+tmp51r;
	tmp59i = tmp20i+tmp51i;
	tmp60r = tmp22r+tmp52r;
	tmp60i = tmp22i+tmp52i;
	tmp61r = tmp24r+tmp53r;
	tmp61i = tmp24i+tmp53i;
	tmp62r = tmp19r-tmp54r;
	tmp62i = tmp19i-tmp54i;
	tmp63r = tmp21r-tmp55r;
	tmp63i = tmp21i-tmp55i;
	tmp64r = tmp23r-tmp56r;
	tmp64i = tmp23i-tmp56i;
	tmp65r = tmp25r-tmp57r;
	tmp65i = tmp25i-tmp57i;
	tmp66r = tmp18r-tmp50r;
	tmp66i = tmp18i-tmp50i;
	tmp67r = tmp20r-tmp51r;
	tmp67i = tmp20i-tmp51i;
	tmp68r = tmp22r-tmp52r;
	tmp68i = tmp22i-tmp52i;
	tmp69r = tmp24r-tmp53r;
	tmp69i = tmp24i-tmp53i;
	tmp70r = tmp19r+tmp54r;
	tmp70i = tmp19i+tmp54i;
	tmp71r = tmp21r+tmp55r;
	tmp71i = tmp21i+tmp55i;
	tmp72r = tmp23r+tmp56r;
	tmp72i = tmp23i+tmp56i;
	tmp73r = tmp25r+tmp57r;
	tmp73i = tmp25i+tmp57i;
	Re(x[0]) = tmp58r;
	Im(x[0]) = tmp58i;
	Re(x[8]) = tmp59r;
	Im(x[8]) = tmp59i;
	Re(x[16]) = tmp60r;
	Im(x[16]) = tmp60i;
	Re(x[24]) = tmp61r;
	Im(x[24]) = tmp61i;
	Re(x[32]) = tmp62r;
	Im(x[32]) = tmp62i;
	Re(x[40]) = tmp63r;
	Im(x[40]) = tmp63i;
	Re(x[48]) = tmp64r;
	Im(x[48]) = tmp64i;
	Re(x[56]) = tmp65r;
	Im(x[56]) = tmp65i;
	Re(x[64]) = tmp66r;
	Im(x[64]) = tmp66i;
	Re(x[72]) = tmp67r;
	Im(x[72]) = tmp67i;
	Re(x[80]) = tmp68r;
	Im(x[80]) = tmp68i;
	Re(x[88]) = tmp69r;
	Im(x[88]) = tmp69i;
	Re(x[96]) = tmp70r;
	Im(x[96]) = tmp70i;
	Re(x[104]) = tmp71r;
	Im(x[104]) = tmp71i;
	Re(x[112]) = tmp72r;
	Im(x[112]) = tmp72i;
	Re(x[120]) = tmp73r;
	Im(x[120]) = tmp73i;
}

/*
*	Number of additions = 160
*	Number of multiplications = 80
*	Number of sign changes = 36
*	Number of assigns = 244
*	Total number of operations = 520
*/
void	MFFTC18(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp99r, tmp100r, tmp101r, tmp102r, tmp103r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp99i, tmp100i, tmp101i, tmp102i, tmp103i;

	const REAL	C8 =     0.17364817766693;	/* REALCONST	*/
	const REAL	C11 =     0.34202014332567;	/* REALCONST	*/
	const REAL	C1 =                  0.5;	/* REALCONST	*/
	const REAL	C7 =     0.64278760968654;	/* REALCONST	*/
	const REAL	C6 =     0.76604444311898;	/* REALCONST	*/
	const REAL	C5 =     0.86602540378444;	/* REALCONST	*/
	const REAL	C10 =     0.93969262078591;	/* REALCONST	*/
	const REAL	C9 =     0.98480775301221;	/* REALCONST	*/

	tmp0i = -C1*Im(x[54]);
	tmp0r = -C1*Re(x[54]);
	tmp1i = C5*Re(x[54]);
	tmp1r = -C5*Im(x[54]);
	tmp2i = tmp0i+Im(x[0]);
	tmp2r = tmp0r+Re(x[0]);
	tmp3i = tmp2i+tmp1i;
	tmp3r = tmp2r+tmp1r;
	tmp4i = tmp2i-tmp1i;
	tmp4r = tmp2r-tmp1r;
	tmp5i = Im(x[0])+Im(x[54]);
	tmp5r = Re(x[0])+Re(x[54]);
	tmp6i = -C1*Im(x[36]);
	tmp6r = -C1*Re(x[36]);
	tmp7i = -C5*Re(x[36]);
	tmp7r = C5*Im(x[36]);
	tmp8i = tmp6i+tmp7i;
	tmp8r = tmp6r+tmp7r;
	tmp9i = tmp6i-tmp7i;
	tmp9r = tmp6r-tmp7r;
	tmp10i = -C1*Im(x[72]);
	tmp10r = -C1*Re(x[72]);
	tmp11i = C5*Re(x[72]);
	tmp11r = -C5*Im(x[72]);
	tmp12i = tmp10i+Im(x[18]);
	tmp12r = tmp10r+Re(x[18]);
	tmp13i = tmp12i+tmp11i;
	tmp13r = tmp12r+tmp11r;
	tmp14i = tmp12i-tmp11i;
	tmp14r = tmp12r-tmp11r;
	tmp15i = Im(x[18])+Im(x[72]);
	tmp15r = Re(x[18])+Re(x[72]);
	tmp16i = C6*tmp8i-C7*tmp8r;
	tmp16r = C6*tmp8r+C7*tmp8i;
	tmp17i = C8*tmp9i-C9*tmp9r;
	tmp17r = C8*tmp9r+C9*tmp9i;
	tmp18i = C8*tmp13i-C9*tmp13r;
	tmp18r = C8*tmp13r+C9*tmp13i;
	tmp19i = -C10*tmp14i-C11*tmp14r;
	tmp19r = -C10*tmp14r+C11*tmp14i;
	tmp20i = Im(x[36])+tmp15i;
	tmp20r = Re(x[36])+tmp15r;
	tmp21i = Im(x[36])-tmp15i;
	tmp21r = Re(x[36])-tmp15r;
	tmp22i = -C1*tmp20i;
	tmp22r = -C1*tmp20r;
	tmp23i = -C5*tmp21r;
	tmp23r = C5*tmp21i;
	tmp24i = tmp22i+tmp5i;
	tmp24r = tmp22r+tmp5r;
	tmp25i = tmp24i+tmp23i;
	tmp25r = tmp24r+tmp23r;
	tmp26i = tmp24i-tmp23i;
	tmp26r = tmp24r-tmp23r;
	tmp27i = tmp5i+tmp20i;
	tmp27r = tmp5r+tmp20r;
	tmp28i = tmp16i+tmp18i;
	tmp28r = tmp16r+tmp18r;
	tmp29i = tmp16i-tmp18i;
	tmp29r = tmp16r-tmp18r;
	tmp30i = -C1*tmp28i;
	tmp30r = -C1*tmp28r;
	tmp31i = -C5*tmp29r;
	tmp31r = C5*tmp29i;
	tmp32i = tmp30i+tmp3i;
	tmp32r = tmp30r+tmp3r;
	tmp33i = tmp32i+tmp31i;
	tmp33r = tmp32r+tmp31r;
	tmp34i = tmp32i-tmp31i;
	tmp34r = tmp32r-tmp31r;
	tmp35i = tmp3i+tmp28i;
	tmp35r = tmp3r+tmp28r;
	tmp36i = tmp17i+tmp19i;
	tmp36r = tmp17r+tmp19r;
	tmp37i = tmp17i-tmp19i;
	tmp37r = tmp17r-tmp19r;
	tmp38i = -C1*tmp36i;
	tmp38r = -C1*tmp36r;
	tmp39i = -C5*tmp37r;
	tmp39r = C5*tmp37i;
	tmp40i = tmp38i+tmp4i;
	tmp40r = tmp38r+tmp4r;
	tmp41i = tmp40i+tmp39i;
	tmp41r = tmp40r+tmp39r;
	tmp42i = tmp40i-tmp39i;
	tmp42r = tmp40r-tmp39r;
	tmp43i = tmp4i+tmp36i;
	tmp43r = tmp4r+tmp36r;
	tmp44i = -C1*Im(x[27]);
	tmp44r = -C1*Re(x[27]);
	tmp45i = -C5*Re(x[27]);
	tmp45r = C5*Im(x[27]);
	tmp46i = tmp44i+tmp45i;
	tmp46r = tmp44r+tmp45r;
	tmp47i = tmp44i-tmp45i;
	tmp47r = tmp44r-tmp45r;
	tmp48i = -C1*Im(x[63]);
	tmp48r = -C1*Re(x[63]);
	tmp49i = C5*Re(x[63]);
	tmp49r = -C5*Im(x[63]);
	tmp50i = tmp48i+Im(x[9]);
	tmp50r = tmp48r+Re(x[9]);
	tmp51i = tmp50i+tmp49i;
	tmp51r = tmp50r+tmp49r;
	tmp52i = tmp50i-tmp49i;
	tmp52r = tmp50r-tmp49r;
	tmp53i = Im(x[9])+Im(x[63]);
	tmp53r = Re(x[9])+Re(x[63]);
	tmp54i = -C1*Im(x[45]);
	tmp54r = -C1*Re(x[45]);
	tmp55i = -C5*Re(x[45]);
	tmp55r = C5*Im(x[45]);
	tmp56i = tmp54i+tmp55i;
	tmp56r = tmp54r+tmp55r;
	tmp57i = tmp54i-tmp55i;
	tmp57r = tmp54r-tmp55r;
	tmp58i = C6*tmp51i-C7*tmp51r;
	tmp58r = C6*tmp51r+C7*tmp51i;
	tmp59i = C8*tmp52i-C9*tmp52r;
	tmp59r = C8*tmp52r+C9*tmp52i;
	tmp60i = C8*tmp56i-C9*tmp56r;
	tmp60r = C8*tmp56r+C9*tmp56i;
	tmp61i = -C10*tmp57i-C11*tmp57r;
	tmp61r = -C10*tmp57r+C11*tmp57i;
	tmp62i = tmp53i+Im(x[45]);
	tmp62r = tmp53r+Re(x[45]);
	tmp63i = tmp53i-Im(x[45]);
	tmp63r = tmp53r-Re(x[45]);
	tmp64i = -C1*tmp62i;
	tmp64r = -C1*tmp62r;
	tmp65i = -C5*tmp63r;
	tmp65r = C5*tmp63i;
	tmp66i = tmp64i+Im(x[27]);
	tmp66r = tmp64r+Re(x[27]);
	tmp67i = tmp66i+tmp65i;
	tmp67r = tmp66r+tmp65r;
	tmp68i = tmp66i-tmp65i;
	tmp68r = tmp66r-tmp65r;
	tmp69i = Im(x[27])+tmp62i;
	tmp69r = Re(x[27])+tmp62r;
	tmp70i = tmp58i+tmp60i;
	tmp70r = tmp58r+tmp60r;
	tmp71i = tmp58i-tmp60i;
	tmp71r = tmp58r-tmp60r;
	tmp72i = -C1*tmp70i;
	tmp72r = -C1*tmp70r;
	tmp73i = -C5*tmp71r;
	tmp73r = C5*tmp71i;
	tmp74i = tmp72i+tmp46i;
	tmp74r = tmp72r+tmp46r;
	tmp75i = tmp74i+tmp73i;
	tmp75r = tmp74r+tmp73r;
	tmp76i = tmp74i-tmp73i;
	tmp76r = tmp74r-tmp73r;
	tmp77i = tmp46i+tmp70i;
	tmp77r = tmp46r+tmp70r;
	tmp78i = tmp59i+tmp61i;
	tmp78r = tmp59r+tmp61r;
	tmp79i = tmp59i-tmp61i;
	tmp79r = tmp59r-tmp61r;
	tmp80i = -C1*tmp78i;
	tmp80r = -C1*tmp78r;
	tmp81i = -C5*tmp79r;
	tmp81r = C5*tmp79i;
	tmp82i = tmp80i+tmp47i;
	tmp82r = tmp80r+tmp47r;
	tmp83i = tmp82i+tmp81i;
	tmp83r = tmp82r+tmp81r;
	tmp84i = tmp82i-tmp81i;
	tmp84r = tmp82r-tmp81r;
	tmp85i = tmp47i+tmp78i;
	tmp85r = tmp47r+tmp78r;
	tmp86i = tmp27i+tmp69i;
	tmp86r = tmp27r+tmp69r;
	tmp87i = tmp27i-tmp69i;
	tmp87r = tmp27r-tmp69r;
	tmp88i = tmp41i+tmp83i;
	tmp88r = tmp41r+tmp83r;
	tmp89i = tmp41i-tmp83i;
	tmp89r = tmp41r-tmp83r;
	tmp90i = tmp35i+tmp77i;
	tmp90r = tmp35r+tmp77r;
	tmp91i = tmp35i-tmp77i;
	tmp91r = tmp35r-tmp77r;
	tmp92i = tmp26i+tmp68i;
	tmp92r = tmp26r+tmp68r;
	tmp93i = tmp26i-tmp68i;
	tmp93r = tmp26r-tmp68r;
	tmp94i = tmp43i+tmp85i;
	tmp94r = tmp43r+tmp85r;
	tmp95i = tmp43i-tmp85i;
	tmp95r = tmp43r-tmp85r;
	tmp96i = tmp34i+tmp76i;
	tmp96r = tmp34r+tmp76r;
	tmp97i = tmp34i-tmp76i;
	tmp97r = tmp34r-tmp76r;
	tmp98i = tmp25i+tmp67i;
	tmp98r = tmp25r+tmp67r;
	tmp99i = tmp25i-tmp67i;
	tmp99r = tmp25r-tmp67r;
	tmp100i = tmp42i+tmp84i;
	tmp100r = tmp42r+tmp84r;
	tmp101i = tmp42i-tmp84i;
	tmp101r = tmp42r-tmp84r;
	tmp102i = tmp33i+tmp75i;
	tmp102r = tmp33r+tmp75r;
	tmp103i = tmp33i-tmp75i;
	tmp103r = tmp33r-tmp75r;
	Re(x[0]) = tmp86r;
	Im(x[0]) = tmp86i;
	Re(x[9]) = tmp89r;
	Im(x[9]) = tmp89i;
	Re(x[18]) = tmp90r;
	Im(x[18]) = tmp90i;
	Re(x[27]) = tmp93r;
	Im(x[27]) = tmp93i;
	Re(x[36]) = tmp94r;
	Im(x[36]) = tmp94i;
	Re(x[45]) = tmp97r;
	Im(x[45]) = tmp97i;
	Re(x[54]) = tmp98r;
	Im(x[54]) = tmp98i;
	Re(x[63]) = tmp101r;
	Im(x[63]) = tmp101i;
	Re(x[72]) = tmp102r;
	Im(x[72]) = tmp102i;
	Re(x[81]) = tmp87r;
	Im(x[81]) = tmp87i;
	Re(x[90]) = tmp88r;
	Im(x[90]) = tmp88i;
	Re(x[99]) = tmp91r;
	Im(x[99]) = tmp91i;
	Re(x[108]) = tmp92r;
	Im(x[108]) = tmp92i;
	Re(x[117]) = tmp95r;
	Im(x[117]) = tmp95i;
	Re(x[126]) = tmp96r;
	Im(x[126]) = tmp96i;
	Re(x[135]) = tmp99r;
	Im(x[135]) = tmp99i;
	Re(x[144]) = tmp100r;
	Im(x[144]) = tmp100i;
	Re(x[153]) = tmp103r;
	Im(x[153]) = tmp103i;
}

/*
*	Number of additions = 196
*	Number of multiplications = 80
*	Number of sign changes = 36
*	Number of assigns = 280
*	Total number of operations = 592
*/
void	MIFFTC18(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp120r, tmp121r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp120i, tmp121i;

	const REAL	C8 =     0.17364817766693;	/* REALCONST	*/
	const REAL	C11 =     0.34202014332567;	/* REALCONST	*/
	const REAL	C1 =                  0.5;	/* REALCONST	*/
	const REAL	C7 =     0.64278760968654;	/* REALCONST	*/
	const REAL	C6 =     0.76604444311898;	/* REALCONST	*/
	const REAL	C5 =     0.86602540378444;	/* REALCONST	*/
	const REAL	C10 =     0.93969262078591;	/* REALCONST	*/
	const REAL	C9 =     0.98480775301221;	/* REALCONST	*/

	tmp0i = Im(x[108])+Im(x[54]);
	tmp0r = Re(x[108])+Re(x[54]);
	tmp1i = Im(x[108])-Im(x[54]);
	tmp1r = Re(x[108])-Re(x[54]);
	tmp2i = -C1*tmp0i;
	tmp2r = -C1*tmp0r;
	tmp3i = C5*tmp1r;
	tmp3r = -C5*tmp1i;
	tmp4i = tmp2i+Im(x[0]);
	tmp4r = tmp2r+Re(x[0]);
	tmp5i = tmp4i+tmp3i;
	tmp5r = tmp4r+tmp3r;
	tmp6i = tmp4i-tmp3i;
	tmp6r = tmp4r-tmp3r;
	tmp7i = Im(x[0])+tmp0i;
	tmp7r = Re(x[0])+tmp0r;
	tmp8i = Im(x[36])+Im(x[144]);
	tmp8r = Re(x[36])+Re(x[144]);
	tmp9i = Im(x[36])-Im(x[144]);
	tmp9r = Re(x[36])-Re(x[144]);
	tmp10i = -C1*tmp8i;
	tmp10r = -C1*tmp8r;
	tmp11i = C5*tmp9r;
	tmp11r = -C5*tmp9i;
	tmp12i = tmp10i+Im(x[90]);
	tmp12r = tmp10r+Re(x[90]);
	tmp13i = tmp12i+tmp11i;
	tmp13r = tmp12r+tmp11r;
	tmp14i = tmp12i-tmp11i;
	tmp14r = tmp12r-tmp11r;
	tmp15i = Im(x[90])+tmp8i;
	tmp15r = Re(x[90])+tmp8r;
	tmp16i = Im(x[126])+Im(x[72]);
	tmp16r = Re(x[126])+Re(x[72]);
	tmp17i = Im(x[126])-Im(x[72]);
	tmp17r = Re(x[126])-Re(x[72]);
	tmp18i = -C1*tmp16i;
	tmp18r = -C1*tmp16r;
	tmp19i = C5*tmp17r;
	tmp19r = -C5*tmp17i;
	tmp20i = tmp18i+Im(x[18]);
	tmp20r = tmp18r+Re(x[18]);
	tmp21i = tmp20i+tmp19i;
	tmp21r = tmp20r+tmp19r;
	tmp22i = tmp20i-tmp19i;
	tmp22r = tmp20r-tmp19r;
	tmp23i = Im(x[18])+tmp16i;
	tmp23r = Re(x[18])+tmp16r;
	tmp24i = C6*tmp13i+C7*tmp13r;
	tmp24r = C6*tmp13r-C7*tmp13i;
	tmp25i = C8*tmp14i+C9*tmp14r;
	tmp25r = C8*tmp14r-C9*tmp14i;
	tmp26i = C8*tmp21i+C9*tmp21r;
	tmp26r = C8*tmp21r-C9*tmp21i;
	tmp27i = -C10*tmp22i+C11*tmp22r;
	tmp27r = -C10*tmp22r-C11*tmp22i;
	tmp28i = tmp15i+tmp23i;
	tmp28r = tmp15r+tmp23r;
	tmp29i = tmp15i-tmp23i;
	tmp29r = tmp15r-tmp23r;
	tmp30i = -C1*tmp28i;
	tmp30r = -C1*tmp28r;
	tmp31i = C5*tmp29r;
	tmp31r = -C5*tmp29i;
	tmp32i = tmp30i+tmp7i;
	tmp32r = tmp30r+tmp7r;
	tmp33i = tmp32i+tmp31i;
	tmp33r = tmp32r+tmp31r;
	tmp34i = tmp32i-tmp31i;
	tmp34r = tmp32r-tmp31r;
	tmp35i = tmp7i+tmp28i;
	tmp35r = tmp7r+tmp28r;
	tmp36i = tmp24i+tmp26i;
	tmp36r = tmp24r+tmp26r;
	tmp37i = tmp24i-tmp26i;
	tmp37r = tmp24r-tmp26r;
	tmp38i = -C1*tmp36i;
	tmp38r = -C1*tmp36r;
	tmp39i = C5*tmp37r;
	tmp39r = -C5*tmp37i;
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
	tmp47i = C5*tmp45r;
	tmp47r = -C5*tmp45i;
	tmp48i = tmp46i+tmp6i;
	tmp48r = tmp46r+tmp6r;
	tmp49i = tmp48i+tmp47i;
	tmp49r = tmp48r+tmp47r;
	tmp50i = tmp48i-tmp47i;
	tmp50r = tmp48r-tmp47r;
	tmp51i = tmp6i+tmp44i;
	tmp51r = tmp6r+tmp44r;
	tmp52i = Im(x[27])+Im(x[135]);
	tmp52r = Re(x[27])+Re(x[135]);
	tmp53i = Im(x[27])-Im(x[135]);
	tmp53r = Re(x[27])-Re(x[135]);
	tmp54i = -C1*tmp52i;
	tmp54r = -C1*tmp52r;
	tmp55i = C5*tmp53r;
	tmp55r = -C5*tmp53i;
	tmp56i = tmp54i+Im(x[81]);
	tmp56r = tmp54r+Re(x[81]);
	tmp57i = tmp56i+tmp55i;
	tmp57r = tmp56r+tmp55r;
	tmp58i = tmp56i-tmp55i;
	tmp58r = tmp56r-tmp55r;
	tmp59i = Im(x[81])+tmp52i;
	tmp59r = Re(x[81])+tmp52r;
	tmp60i = Im(x[117])+Im(x[63]);
	tmp60r = Re(x[117])+Re(x[63]);
	tmp61i = Im(x[117])-Im(x[63]);
	tmp61r = Re(x[117])-Re(x[63]);
	tmp62i = -C1*tmp60i;
	tmp62r = -C1*tmp60r;
	tmp63i = C5*tmp61r;
	tmp63r = -C5*tmp61i;
	tmp64i = tmp62i+Im(x[9]);
	tmp64r = tmp62r+Re(x[9]);
	tmp65i = tmp64i+tmp63i;
	tmp65r = tmp64r+tmp63r;
	tmp66i = tmp64i-tmp63i;
	tmp66r = tmp64r-tmp63r;
	tmp67i = Im(x[9])+tmp60i;
	tmp67r = Re(x[9])+tmp60r;
	tmp68i = Im(x[45])+Im(x[153]);
	tmp68r = Re(x[45])+Re(x[153]);
	tmp69i = Im(x[45])-Im(x[153]);
	tmp69r = Re(x[45])-Re(x[153]);
	tmp70i = -C1*tmp68i;
	tmp70r = -C1*tmp68r;
	tmp71i = C5*tmp69r;
	tmp71r = -C5*tmp69i;
	tmp72i = tmp70i+Im(x[99]);
	tmp72r = tmp70r+Re(x[99]);
	tmp73i = tmp72i+tmp71i;
	tmp73r = tmp72r+tmp71r;
	tmp74i = tmp72i-tmp71i;
	tmp74r = tmp72r-tmp71r;
	tmp75i = Im(x[99])+tmp68i;
	tmp75r = Re(x[99])+tmp68r;
	tmp76i = C6*tmp65i+C7*tmp65r;
	tmp76r = C6*tmp65r-C7*tmp65i;
	tmp77i = C8*tmp66i+C9*tmp66r;
	tmp77r = C8*tmp66r-C9*tmp66i;
	tmp78i = C8*tmp73i+C9*tmp73r;
	tmp78r = C8*tmp73r-C9*tmp73i;
	tmp79i = -C10*tmp74i+C11*tmp74r;
	tmp79r = -C10*tmp74r-C11*tmp74i;
	tmp80i = tmp67i+tmp75i;
	tmp80r = tmp67r+tmp75r;
	tmp81i = tmp67i-tmp75i;
	tmp81r = tmp67r-tmp75r;
	tmp82i = -C1*tmp80i;
	tmp82r = -C1*tmp80r;
	tmp83i = C5*tmp81r;
	tmp83r = -C5*tmp81i;
	tmp84i = tmp82i+tmp59i;
	tmp84r = tmp82r+tmp59r;
	tmp85i = tmp84i+tmp83i;
	tmp85r = tmp84r+tmp83r;
	tmp86i = tmp84i-tmp83i;
	tmp86r = tmp84r-tmp83r;
	tmp87i = tmp59i+tmp80i;
	tmp87r = tmp59r+tmp80r;
	tmp88i = tmp76i+tmp78i;
	tmp88r = tmp76r+tmp78r;
	tmp89i = tmp76i-tmp78i;
	tmp89r = tmp76r-tmp78r;
	tmp90i = -C1*tmp88i;
	tmp90r = -C1*tmp88r;
	tmp91i = C5*tmp89r;
	tmp91r = -C5*tmp89i;
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
	tmp99i = C5*tmp97r;
	tmp99r = -C5*tmp97i;
	tmp100i = tmp98i+tmp58i;
	tmp100r = tmp98r+tmp58r;
	tmp101i = tmp100i+tmp99i;
	tmp101r = tmp100r+tmp99r;
	tmp102i = tmp100i-tmp99i;
	tmp102r = tmp100r-tmp99r;
	tmp103i = tmp58i+tmp96i;
	tmp103r = tmp58r+tmp96r;
	tmp104i = tmp35i+tmp87i;
	tmp104r = tmp35r+tmp87r;
	tmp105i = tmp35i-tmp87i;
	tmp105r = tmp35r-tmp87r;
	tmp106i = tmp49i+tmp101i;
	tmp106r = tmp49r+tmp101r;
	tmp107i = tmp49i-tmp101i;
	tmp107r = tmp49r-tmp101r;
	tmp108i = tmp43i+tmp95i;
	tmp108r = tmp43r+tmp95r;
	tmp109i = tmp43i-tmp95i;
	tmp109r = tmp43r-tmp95r;
	tmp110i = tmp34i+tmp86i;
	tmp110r = tmp34r+tmp86r;
	tmp111i = tmp34i-tmp86i;
	tmp111r = tmp34r-tmp86r;
	tmp112i = tmp51i+tmp103i;
	tmp112r = tmp51r+tmp103r;
	tmp113i = tmp51i-tmp103i;
	tmp113r = tmp51r-tmp103r;
	tmp114i = tmp42i+tmp94i;
	tmp114r = tmp42r+tmp94r;
	tmp115i = tmp42i-tmp94i;
	tmp115r = tmp42r-tmp94r;
	tmp116i = tmp33i+tmp85i;
	tmp116r = tmp33r+tmp85r;
	tmp117i = tmp33i-tmp85i;
	tmp117r = tmp33r-tmp85r;
	tmp118i = tmp50i+tmp102i;
	tmp118r = tmp50r+tmp102r;
	tmp119i = tmp50i-tmp102i;
	tmp119r = tmp50r-tmp102r;
	tmp120i = tmp41i+tmp93i;
	tmp120r = tmp41r+tmp93r;
	tmp121i = tmp41i-tmp93i;
	tmp121r = tmp41r-tmp93r;
	Re(x[0]) = tmp104r;
	Im(x[0]) = tmp104i;
	Re(x[9]) = tmp107r;
	Im(x[9]) = tmp107i;
	Re(x[18]) = tmp108r;
	Im(x[18]) = tmp108i;
	Re(x[27]) = tmp111r;
	Im(x[27]) = tmp111i;
	Re(x[36]) = tmp112r;
	Im(x[36]) = tmp112i;
	Re(x[45]) = tmp115r;
	Im(x[45]) = tmp115i;
	Re(x[54]) = tmp116r;
	Im(x[54]) = tmp116i;
	Re(x[63]) = tmp119r;
	Im(x[63]) = tmp119i;
	Re(x[72]) = tmp120r;
	Im(x[72]) = tmp120i;
	Re(x[81]) = tmp105r;
	Im(x[81]) = tmp105i;
	Re(x[90]) = tmp106r;
	Im(x[90]) = tmp106i;
	Re(x[99]) = tmp109r;
	Im(x[99]) = tmp109i;
	Re(x[108]) = tmp110r;
	Im(x[108]) = tmp110i;
	Re(x[117]) = tmp113r;
	Im(x[117]) = tmp113i;
	Re(x[126]) = tmp114r;
	Im(x[126]) = tmp114i;
	Re(x[135]) = tmp117r;
	Im(x[135]) = tmp117i;
	Re(x[144]) = tmp118r;
	Im(x[144]) = tmp118i;
	Re(x[153]) = tmp121r;
	Im(x[153]) = tmp121i;
}

/*
*	Number of additions = 160
*	Number of multiplications = 40
*	Number of sign changes = 14
*	Number of assigns = 224
*	Total number of operations = 438
*/
void	MFFTC20(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i;

	const REAL	C1 =                 0.25;	/* REALCONST	*/
	const REAL	C14 =     0.55901699437495;	/* REALCONST	*/
	const REAL	C10 =     0.58778525229247;	/* REALCONST	*/
	const REAL	C8 =     0.95105651629515;	/* REALCONST	*/

	tmp0i = Im(x[40])+Im(x[80]);
	tmp0r = Re(x[40])+Re(x[80]);
	tmp1i = Im(x[40])-Im(x[80]);
	tmp1r = Re(x[40])-Re(x[80]);
	tmp2i = -C1*tmp0i;
	tmp2r = -C1*tmp0r;
	tmp3i = C14*tmp1i;
	tmp3r = C14*tmp1r;
	tmp4i = tmp2i+Im(x[0]);
	tmp4r = tmp2r+Re(x[0]);
	tmp5i = tmp4i+tmp3i;
	tmp5r = tmp4r+tmp3r;
	tmp6i = tmp4i-tmp3i;
	tmp6r = tmp4r-tmp3r;
	tmp7i = C8*Re(x[40])+C10*Re(x[80]);
	tmp7r = -C8*Im(x[40])-C10*Im(x[80]);
	tmp8i = C10*Re(x[40])-C8*Re(x[80]);
	tmp8r = -C10*Im(x[40])+C8*Im(x[80]);
	tmp9i = tmp5i+tmp7i;
	tmp9r = tmp5r+tmp7r;
	tmp10i = tmp5i-tmp7i;
	tmp10r = tmp5r-tmp7r;
	tmp11i = tmp6i+tmp8i;
	tmp11r = tmp6r+tmp8r;
	tmp12i = tmp6i-tmp8i;
	tmp12r = tmp6r-tmp8r;
	tmp13i = Im(x[0])+tmp0i;
	tmp13r = Re(x[0])+tmp0r;
	tmp14i = Im(x[10])+Im(x[90]);
	tmp14r = Re(x[10])+Re(x[90]);
	tmp15i = Im(x[10])-Im(x[90]);
	tmp15r = Re(x[10])-Re(x[90]);
	tmp16i = -C1*tmp14i;
	tmp16r = -C1*tmp14r;
	tmp17i = C14*tmp14i;
	tmp17r = C14*tmp14r;
	tmp18i = tmp16i+Im(x[50]);
	tmp18r = tmp16r+Re(x[50]);
	tmp19i = tmp18i+tmp17i;
	tmp19r = tmp18r+tmp17r;
	tmp20i = tmp18i-tmp17i;
	tmp20r = tmp18r-tmp17r;
	tmp21i = -C8*tmp15r;
	tmp21r = C8*tmp15i;
	tmp22i = -C10*tmp15r;
	tmp22r = C10*tmp15i;
	tmp23i = tmp19i+tmp21i;
	tmp23r = tmp19r+tmp21r;
	tmp24i = tmp19i-tmp21i;
	tmp24r = tmp19r-tmp21r;
	tmp25i = tmp20i+tmp22i;
	tmp25r = tmp20r+tmp22r;
	tmp26i = tmp20i-tmp22i;
	tmp26r = tmp20r-tmp22r;
	tmp27i = Im(x[50])+tmp14i;
	tmp27r = Re(x[50])+tmp14r;
	tmp28i = Im(x[60])+Im(x[20]);
	tmp28r = Re(x[60])+Re(x[20]);
	tmp29i = Im(x[60])-Im(x[20]);
	tmp29r = Re(x[60])-Re(x[20]);
	tmp30i = -C1*tmp28i;
	tmp30r = -C1*tmp28r;
	tmp31i = C14*tmp29i;
	tmp31r = C14*tmp29r;
	tmp32i = tmp30i+tmp31i;
	tmp32r = tmp30r+tmp31r;
	tmp33i = tmp30i-tmp31i;
	tmp33r = tmp30r-tmp31r;
	tmp34i = -C8*Re(x[60])-C10*Re(x[20]);
	tmp34r = C8*Im(x[60])+C10*Im(x[20]);
	tmp35i = -C10*Re(x[60])+C8*Re(x[20]);
	tmp35r = C10*Im(x[60])-C8*Im(x[20]);
	tmp36i = tmp32i+tmp34i;
	tmp36r = tmp32r+tmp34r;
	tmp37i = tmp32i-tmp34i;
	tmp37r = tmp32r-tmp34r;
	tmp38i = tmp33i+tmp35i;
	tmp38r = tmp33r+tmp35r;
	tmp39i = tmp33i-tmp35i;
	tmp39r = tmp33r-tmp35r;
	tmp40i = Im(x[70])+Im(x[30]);
	tmp40r = Re(x[70])+Re(x[30]);
	tmp41i = Im(x[70])-Im(x[30]);
	tmp41r = Re(x[70])-Re(x[30]);
	tmp42i = -C1*tmp40i;
	tmp42r = -C1*tmp40r;
	tmp43i = -C14*tmp40i;
	tmp43r = -C14*tmp40r;
	tmp44i = tmp42i+tmp43i;
	tmp44r = tmp42r+tmp43r;
	tmp45i = tmp42i-tmp43i;
	tmp45r = tmp42r-tmp43r;
	tmp46i = -C10*tmp41r;
	tmp46r = C10*tmp41i;
	tmp47i = C8*tmp41r;
	tmp47r = -C8*tmp41i;
	tmp48i = tmp44i+tmp46i;
	tmp48r = tmp44r+tmp46r;
	tmp49i = tmp44i-tmp46i;
	tmp49r = tmp44r-tmp46r;
	tmp50i = tmp45i+tmp47i;
	tmp50r = tmp45r+tmp47r;
	tmp51i = tmp45i-tmp47i;
	tmp51r = tmp45r-tmp47r;
	tmp52i = tmp13i+tmp28i;
	tmp52r = tmp13r+tmp28r;
	tmp53i = tmp13i-tmp28i;
	tmp53r = tmp13r-tmp28r;
	tmp54i = tmp27i+tmp40i;
	tmp54r = tmp27r+tmp40r;
	tmp55i = tmp27i-tmp40i;
	tmp55r = tmp27r-tmp40r;
	tmp56i = tmp52i+tmp54i;
	tmp56r = tmp52r+tmp54r;
	tmp57i = tmp52i-tmp54i;
	tmp57r = tmp52r-tmp54r;
	tmp58i = tmp53i-tmp55r;
	tmp58r = tmp53r+tmp55i;
	tmp59i = tmp53i+tmp55r;
	tmp59r = tmp53r-tmp55i;
	tmp60i = tmp10i+tmp37i;
	tmp60r = tmp10r+tmp37r;
	tmp61i = tmp10i-tmp37i;
	tmp61r = tmp10r-tmp37r;
	tmp62i = tmp24i+tmp49i;
	tmp62r = tmp24r+tmp49r;
	tmp63i = tmp24i-tmp49i;
	tmp63r = tmp24r-tmp49r;
	tmp64i = tmp60i+tmp62i;
	tmp64r = tmp60r+tmp62r;
	tmp65i = tmp60i-tmp62i;
	tmp65r = tmp60r-tmp62r;
	tmp66i = tmp61i-tmp63r;
	tmp66r = tmp61r+tmp63i;
	tmp67i = tmp61i+tmp63r;
	tmp67r = tmp61r-tmp63i;
	tmp68i = tmp12i+tmp39i;
	tmp68r = tmp12r+tmp39r;
	tmp69i = tmp12i-tmp39i;
	tmp69r = tmp12r-tmp39r;
	tmp70i = tmp26i+tmp51i;
	tmp70r = tmp26r+tmp51r;
	tmp71i = tmp26i-tmp51i;
	tmp71r = tmp26r-tmp51r;
	tmp72i = tmp68i+tmp70i;
	tmp72r = tmp68r+tmp70r;
	tmp73i = tmp68i-tmp70i;
	tmp73r = tmp68r-tmp70r;
	tmp74i = tmp69i-tmp71r;
	tmp74r = tmp69r+tmp71i;
	tmp75i = tmp69i+tmp71r;
	tmp75r = tmp69r-tmp71i;
	tmp76i = tmp11i+tmp38i;
	tmp76r = tmp11r+tmp38r;
	tmp77i = tmp11i-tmp38i;
	tmp77r = tmp11r-tmp38r;
	tmp78i = tmp25i+tmp50i;
	tmp78r = tmp25r+tmp50r;
	tmp79i = tmp25i-tmp50i;
	tmp79r = tmp25r-tmp50r;
	tmp80i = tmp76i+tmp78i;
	tmp80r = tmp76r+tmp78r;
	tmp81i = tmp76i-tmp78i;
	tmp81r = tmp76r-tmp78r;
	tmp82i = tmp77i-tmp79r;
	tmp82r = tmp77r+tmp79i;
	tmp83i = tmp77i+tmp79r;
	tmp83r = tmp77r-tmp79i;
	tmp84i = tmp9i+tmp36i;
	tmp84r = tmp9r+tmp36r;
	tmp85i = tmp9i-tmp36i;
	tmp85r = tmp9r-tmp36r;
	tmp86i = tmp23i+tmp48i;
	tmp86r = tmp23r+tmp48r;
	tmp87i = tmp23i-tmp48i;
	tmp87r = tmp23r-tmp48r;
	tmp88i = tmp84i+tmp86i;
	tmp88r = tmp84r+tmp86r;
	tmp89i = tmp84i-tmp86i;
	tmp89r = tmp84r-tmp86r;
	tmp90i = tmp85i-tmp87r;
	tmp90r = tmp85r+tmp87i;
	tmp91i = tmp85i+tmp87r;
	tmp91r = tmp85r-tmp87i;
	Re(x[0]) = tmp56r;
	Im(x[0]) = tmp56i;
	Re(x[10]) = tmp66r;
	Im(x[10]) = tmp66i;
	Re(x[20]) = tmp73r;
	Im(x[20]) = tmp73i;
	Re(x[30]) = tmp83r;
	Im(x[30]) = tmp83i;
	Re(x[40]) = tmp88r;
	Im(x[40]) = tmp88i;
	Re(x[50]) = tmp58r;
	Im(x[50]) = tmp58i;
	Re(x[60]) = tmp65r;
	Im(x[60]) = tmp65i;
	Re(x[70]) = tmp75r;
	Im(x[70]) = tmp75i;
	Re(x[80]) = tmp80r;
	Im(x[80]) = tmp80i;
	Re(x[90]) = tmp90r;
	Im(x[90]) = tmp90i;
	Re(x[100]) = tmp57r;
	Im(x[100]) = tmp57i;
	Re(x[110]) = tmp67r;
	Im(x[110]) = tmp67i;
	Re(x[120]) = tmp72r;
	Im(x[120]) = tmp72i;
	Re(x[130]) = tmp82r;
	Im(x[130]) = tmp82i;
	Re(x[140]) = tmp89r;
	Im(x[140]) = tmp89i;
	Re(x[150]) = tmp59r;
	Im(x[150]) = tmp59i;
	Re(x[160]) = tmp64r;
	Im(x[160]) = tmp64i;
	Re(x[170]) = tmp74r;
	Im(x[170]) = tmp74i;
	Re(x[180]) = tmp81r;
	Im(x[180]) = tmp81i;
	Re(x[190]) = tmp91r;
	Im(x[190]) = tmp91i;
}

/*
*	Number of additions = 208
*	Number of multiplications = 48
*	Number of sign changes = 8
*	Number of assigns = 264
*	Total number of operations = 528
*/
void	MIFFTC20(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp106r, tmp107r, tmp108r, tmp109r, tmp110r, tmp111r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp106i, tmp107i, tmp108i, tmp109i, tmp110i, tmp111i;

	const REAL	C1 =                 0.25;	/* REALCONST	*/
	const REAL	C14 =     0.55901699437495;	/* REALCONST	*/
	const REAL	C10 =     0.58778525229247;	/* REALCONST	*/
	const REAL	C8 =     0.95105651629515;	/* REALCONST	*/

	tmp0i = Im(x[160])+Im(x[40]);
	tmp0r = Re(x[160])+Re(x[40]);
	tmp1i = Im(x[160])-Im(x[40]);
	tmp1r = Re(x[160])-Re(x[40]);
	tmp2i = Im(x[120])+Im(x[80]);
	tmp2r = Re(x[120])+Re(x[80]);
	tmp3i = Im(x[120])-Im(x[80]);
	tmp3r = Re(x[120])-Re(x[80]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = -C1*tmp4i;
	tmp6r = -C1*tmp4r;
	tmp7i = C14*tmp5i;
	tmp7r = C14*tmp5r;
	tmp8i = tmp6i+Im(x[0]);
	tmp8r = tmp6r+Re(x[0]);
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
	tmp17i = Im(x[0])+tmp4i;
	tmp17r = Re(x[0])+tmp4r;
	tmp18i = Im(x[10])+Im(x[90]);
	tmp18r = Re(x[10])+Re(x[90]);
	tmp19i = Im(x[10])-Im(x[90]);
	tmp19r = Re(x[10])-Re(x[90]);
	tmp20i = Im(x[170])+Im(x[130]);
	tmp20r = Re(x[170])+Re(x[130]);
	tmp21i = Im(x[170])-Im(x[130]);
	tmp21r = Re(x[170])-Re(x[130]);
	tmp22i = tmp18i+tmp20i;
	tmp22r = tmp18r+tmp20r;
	tmp23i = tmp18i-tmp20i;
	tmp23r = tmp18r-tmp20r;
	tmp24i = -C1*tmp22i;
	tmp24r = -C1*tmp22r;
	tmp25i = C14*tmp23i;
	tmp25r = C14*tmp23r;
	tmp26i = tmp24i+Im(x[50]);
	tmp26r = tmp24r+Re(x[50]);
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
	tmp35i = Im(x[50])+tmp22i;
	tmp35r = Re(x[50])+tmp22r;
	tmp36i = Im(x[60])+Im(x[140]);
	tmp36r = Re(x[60])+Re(x[140]);
	tmp37i = Im(x[60])-Im(x[140]);
	tmp37r = Re(x[60])-Re(x[140]);
	tmp38i = Im(x[20])+Im(x[180]);
	tmp38r = Re(x[20])+Re(x[180]);
	tmp39i = Im(x[20])-Im(x[180]);
	tmp39r = Re(x[20])-Re(x[180]);
	tmp40i = tmp36i+tmp38i;
	tmp40r = tmp36r+tmp38r;
	tmp41i = tmp36i-tmp38i;
	tmp41r = tmp36r-tmp38r;
	tmp42i = -C1*tmp40i;
	tmp42r = -C1*tmp40r;
	tmp43i = C14*tmp41i;
	tmp43r = C14*tmp41r;
	tmp44i = tmp42i+Im(x[100]);
	tmp44r = tmp42r+Re(x[100]);
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
	tmp53i = Im(x[100])+tmp40i;
	tmp53r = Re(x[100])+tmp40r;
	tmp54i = Im(x[110])+Im(x[190]);
	tmp54r = Re(x[110])+Re(x[190]);
	tmp55i = Im(x[110])-Im(x[190]);
	tmp55r = Re(x[110])-Re(x[190]);
	tmp56i = Im(x[70])+Im(x[30]);
	tmp56r = Re(x[70])+Re(x[30]);
	tmp57i = Im(x[70])-Im(x[30]);
	tmp57r = Re(x[70])-Re(x[30]);
	tmp58i = tmp54i+tmp56i;
	tmp58r = tmp54r+tmp56r;
	tmp59i = tmp54i-tmp56i;
	tmp59r = tmp54r-tmp56r;
	tmp60i = -C1*tmp58i;
	tmp60r = -C1*tmp58r;
	tmp61i = C14*tmp59i;
	tmp61r = C14*tmp59r;
	tmp62i = tmp60i+Im(x[150]);
	tmp62r = tmp60r+Re(x[150]);
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
	tmp71i = Im(x[150])+tmp58i;
	tmp71r = Re(x[150])+tmp58r;
	tmp72i = tmp17i+tmp53i;
	tmp72r = tmp17r+tmp53r;
	tmp73i = tmp17i-tmp53i;
	tmp73r = tmp17r-tmp53r;
	tmp74i = tmp35i+tmp71i;
	tmp74r = tmp35r+tmp71r;
	tmp75i = tmp35i-tmp71i;
	tmp75r = tmp35r-tmp71r;
	tmp76i = tmp72i+tmp74i;
	tmp76r = tmp72r+tmp74r;
	tmp77i = tmp72i-tmp74i;
	tmp77r = tmp72r-tmp74r;
	tmp78i = tmp73i+tmp75r;
	tmp78r = tmp73r-tmp75i;
	tmp79i = tmp73i-tmp75r;
	tmp79r = tmp73r+tmp75i;
	tmp80i = tmp14i+tmp50i;
	tmp80r = tmp14r+tmp50r;
	tmp81i = tmp14i-tmp50i;
	tmp81r = tmp14r-tmp50r;
	tmp82i = tmp32i+tmp68i;
	tmp82r = tmp32r+tmp68r;
	tmp83i = tmp32i-tmp68i;
	tmp83r = tmp32r-tmp68r;
	tmp84i = tmp80i+tmp82i;
	tmp84r = tmp80r+tmp82r;
	tmp85i = tmp80i-tmp82i;
	tmp85r = tmp80r-tmp82r;
	tmp86i = tmp81i+tmp83r;
	tmp86r = tmp81r-tmp83i;
	tmp87i = tmp81i-tmp83r;
	tmp87r = tmp81r+tmp83i;
	tmp88i = tmp16i+tmp52i;
	tmp88r = tmp16r+tmp52r;
	tmp89i = tmp16i-tmp52i;
	tmp89r = tmp16r-tmp52r;
	tmp90i = tmp34i+tmp70i;
	tmp90r = tmp34r+tmp70r;
	tmp91i = tmp34i-tmp70i;
	tmp91r = tmp34r-tmp70r;
	tmp92i = tmp88i+tmp90i;
	tmp92r = tmp88r+tmp90r;
	tmp93i = tmp88i-tmp90i;
	tmp93r = tmp88r-tmp90r;
	tmp94i = tmp89i+tmp91r;
	tmp94r = tmp89r-tmp91i;
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
	tmp101i = tmp96i-tmp98i;
	tmp101r = tmp96r-tmp98r;
	tmp102i = tmp97i+tmp99r;
	tmp102r = tmp97r-tmp99i;
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
	tmp109i = tmp104i-tmp106i;
	tmp109r = tmp104r-tmp106r;
	tmp110i = tmp105i+tmp107r;
	tmp110r = tmp105r-tmp107i;
	tmp111i = tmp105i-tmp107r;
	tmp111r = tmp105r+tmp107i;
	Re(x[0]) = tmp76r;
	Im(x[0]) = tmp76i;
	Re(x[10]) = tmp86r;
	Im(x[10]) = tmp86i;
	Re(x[20]) = tmp93r;
	Im(x[20]) = tmp93i;
	Re(x[30]) = tmp103r;
	Im(x[30]) = tmp103i;
	Re(x[40]) = tmp108r;
	Im(x[40]) = tmp108i;
	Re(x[50]) = tmp78r;
	Im(x[50]) = tmp78i;
	Re(x[60]) = tmp85r;
	Im(x[60]) = tmp85i;
	Re(x[70]) = tmp95r;
	Im(x[70]) = tmp95i;
	Re(x[80]) = tmp100r;
	Im(x[80]) = tmp100i;
	Re(x[90]) = tmp110r;
	Im(x[90]) = tmp110i;
	Re(x[100]) = tmp77r;
	Im(x[100]) = tmp77i;
	Re(x[110]) = tmp87r;
	Im(x[110]) = tmp87i;
	Re(x[120]) = tmp92r;
	Im(x[120]) = tmp92i;
	Re(x[130]) = tmp102r;
	Im(x[130]) = tmp102i;
	Re(x[140]) = tmp109r;
	Im(x[140]) = tmp109i;
	Re(x[150]) = tmp79r;
	Im(x[150]) = tmp79i;
	Re(x[160]) = tmp84r;
	Im(x[160]) = tmp84i;
	Re(x[170]) = tmp94r;
	Im(x[170]) = tmp94i;
	Re(x[180]) = tmp101r;
	Im(x[180]) = tmp101i;
	Re(x[190]) = tmp111r;
	Im(x[190]) = tmp111i;
}

/*
*	Number of additions = 272
*	Number of multiplications = 200
*	Number of sign changes = 0
*	Number of assigns = 186
*	Total number of operations = 658
*/
void	MFFTC22(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
		tmp8r, tmp9r, tmp10r, tmp11r, tmp12r, tmp13r, tmp14r,
		tmp15r, tmp16r, tmp17r, tmp18r, tmp19r, tmp20r, tmp21r,
		tmp22r, tmp23r, tmp24r, tmp25r, tmp26r, tmp27r, tmp28r,
		tmp29r, tmp30r, tmp31r, tmp32r, tmp33r, tmp34r, tmp35r,
		tmp36r, tmp37r, tmp38r, tmp39r, tmp40r, tmp41r, tmp42r,
		tmp43r, tmp44r, tmp45r, tmp46r, tmp47r, tmp48r, tmp49r,
		tmp50r, tmp51r, tmp52r, tmp53r, tmp54r, tmp55r, tmp56r,
		tmp57r, tmp58r, tmp59r, tmp60r, tmp61r, tmp62r, tmp63r,
		tmp64r, tmp65r, tmp66r, tmp67r, tmp68r, tmp69r, tmp70r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
		tmp8i, tmp9i, tmp10i, tmp11i, tmp12i, tmp13i, tmp14i,
		tmp15i, tmp16i, tmp17i, tmp18i, tmp19i, tmp20i, tmp21i,
		tmp22i, tmp23i, tmp24i, tmp25i, tmp26i, tmp27i, tmp28i,
		tmp29i, tmp30i, tmp31i, tmp32i, tmp33i, tmp34i, tmp35i,
		tmp36i, tmp37i, tmp38i, tmp39i, tmp40i, tmp41i, tmp42i,
		tmp43i, tmp44i, tmp45i, tmp46i, tmp47i, tmp48i, tmp49i,
		tmp50i, tmp51i, tmp52i, tmp53i, tmp54i, tmp55i, tmp56i,
		tmp57i, tmp58i, tmp59i, tmp60i, tmp61i, tmp62i, tmp63i,
		tmp64i, tmp65i, tmp66i, tmp67i, tmp68i, tmp69i, tmp70i;

	const REAL	C19 =     0.14231483827329;	/* REALCONST	*/
	const REAL	C22 =     0.28173255684143;	/* REALCONST	*/
	const REAL	C15 =     0.41541501300189;	/* REALCONST	*/
	const REAL	C14 =      0.5406408174556;	/* REALCONST	*/
	const REAL	C17 =     0.65486073394529;	/* REALCONST	*/
	const REAL	C18 =     0.75574957435426;	/* REALCONST	*/
	const REAL	C13 =     0.84125353283118;	/* REALCONST	*/
	const REAL	C16 =     0.90963199535452;	/* REALCONST	*/
	const REAL	C21 =      0.9594929736145;	/* REALCONST	*/
	const REAL	C20 =     0.98982144188093;	/* REALCONST	*/

	tmp0i = C13*Im(x[110])+C15*Im(x[22])-C17*Im(x[44])-C19*Im(x[88])-C21*Im(x[66])+Im(x[0]);
	tmp0r = C13*Re(x[110])+C15*Re(x[22])-C17*Re(x[44])-C19*Re(x[88])-C21*Re(x[66])+Re(x[0]);
	tmp1i = Im(x[110])+Im(x[22]);
	tmp1r = Re(x[110])+Re(x[22]);
	tmp2i = C15*Im(x[110])-C17*Im(x[22])-C19*Im(x[44])-C21*Im(x[88])+C13*Im(x[66])+Im(x[0]);
	tmp2r = C15*Re(x[110])-C17*Re(x[22])-C19*Re(x[44])-C21*Re(x[88])+C13*Re(x[66])+Re(x[0]);
	tmp3i = tmp1i+Im(x[44]);
	tmp3r = tmp1r+Re(x[44]);
	tmp4i = -C17*Im(x[110])-C19*Im(x[22])-C21*Im(x[44])+C13*Im(x[88])+C15*Im(x[66])+Im(x[0]);
	tmp4r = -C17*Re(x[110])-C19*Re(x[22])-C21*Re(x[44])+C13*Re(x[88])+C15*Re(x[66])+Re(x[0]);
	tmp5i = tmp3i+Im(x[88]);
	tmp5r = tmp3r+Re(x[88]);
	tmp6i = -C19*Im(x[110])-C21*Im(x[22])+C13*Im(x[44])+C15*Im(x[88])-C17*Im(x[66])+Im(x[0]);
	tmp6r = -C19*Re(x[110])-C21*Re(x[22])+C13*Re(x[44])+C15*Re(x[88])-C17*Re(x[66])+Re(x[0]);
	tmp7i = tmp5i+Im(x[66]);
	tmp7r = tmp5r+Re(x[66]);
	tmp8i = -C21*Im(x[110])+C13*Im(x[22])+C15*Im(x[44])-C17*Im(x[88])-C19*Im(x[66])+Im(x[0]);
	tmp8r = -C21*Re(x[110])+C13*Re(x[22])+C15*Re(x[44])-C17*Re(x[88])-C19*Re(x[66])+Re(x[0]);
	tmp9i = C14*Re(x[110])-C16*Re(x[22])-C18*Re(x[44])+C20*Re(x[88])+C22*Re(x[66]);
	tmp9r = -C14*Im(x[110])+C16*Im(x[22])+C18*Im(x[44])-C20*Im(x[88])-C22*Im(x[66]);
	tmp10i = C16*Re(x[110])-C18*Re(x[22])+C20*Re(x[44])-C22*Re(x[88])-C14*Re(x[66]);
	tmp10r = -C16*Im(x[110])+C18*Im(x[22])-C20*Im(x[44])+C22*Im(x[88])+C14*Im(x[66]);
	tmp11i = C18*Re(x[110])+C20*Re(x[22])-C22*Re(x[44])+C14*Re(x[88])-C16*Re(x[66]);
	tmp11r = -C18*Im(x[110])-C20*Im(x[22])+C22*Im(x[44])-C14*Im(x[88])+C16*Im(x[66]);
	tmp12i = -C20*Re(x[110])-C22*Re(x[22])+C14*Re(x[44])+C16*Re(x[88])-C18*Re(x[66]);
	tmp12r = C20*Im(x[110])+C22*Im(x[22])-C14*Im(x[44])-C16*Im(x[88])+C18*Im(x[66]);
	tmp13i = C22*Re(x[110])+C14*Re(x[22])+C16*Re(x[44])+C18*Re(x[88])+C20*Re(x[66]);
	tmp13r = -C22*Im(x[110])-C14*Im(x[22])-C16*Im(x[44])-C18*Im(x[88])-C20*Im(x[66]);
	tmp14i = tmp0i+tmp9i;
	tmp14r = tmp0r+tmp9r;
	tmp15i = tmp0i-tmp9i;
	tmp15r = tmp0r-tmp9r;
	tmp16i = tmp2i+tmp10i;
	tmp16r = tmp2r+tmp10r;
	tmp17i = tmp2i-tmp10i;
	tmp17r = tmp2r-tmp10r;
	tmp18i = tmp4i+tmp11i;
	tmp18r = tmp4r+tmp11r;
	tmp19i = tmp4i-tmp11i;
	tmp19r = tmp4r-tmp11r;
	tmp20i = tmp6i+tmp12i;
	tmp20r = tmp6r+tmp12r;
	tmp21i = tmp6i-tmp12i;
	tmp21r = tmp6r-tmp12r;
	tmp22i = tmp8i+tmp13i;
	tmp22r = tmp8r+tmp13r;
	tmp23i = tmp8i-tmp13i;
	tmp23r = tmp8r-tmp13r;
	tmp24i = Im(x[0])+tmp7i;
	tmp24r = Re(x[0])+tmp7r;
	tmp25i = C13*Im(x[11])+C15*Im(x[99])-C17*Im(x[77])-C19*Im(x[33])-C21*Im(x[55]);
	tmp25r = C13*Re(x[11])+C15*Re(x[99])-C17*Re(x[77])-C19*Re(x[33])-C21*Re(x[55]);
	tmp26i = Im(x[11])+Im(x[99]);
	tmp26r = Re(x[11])+Re(x[99]);
	tmp27i = C15*Im(x[11])-C17*Im(x[99])-C19*Im(x[77])-C21*Im(x[33])+C13*Im(x[55]);
	tmp27r = C15*Re(x[11])-C17*Re(x[99])-C19*Re(x[77])-C21*Re(x[33])+C13*Re(x[55]);
	tmp28i = tmp26i+Im(x[77]);
	tmp28r = tmp26r+Re(x[77]);
	tmp29i = -C17*Im(x[11])-C19*Im(x[99])-C21*Im(x[77])+C13*Im(x[33])+C15*Im(x[55]);
	tmp29r = -C17*Re(x[11])-C19*Re(x[99])-C21*Re(x[77])+C13*Re(x[33])+C15*Re(x[55]);
	tmp30i = tmp28i+Im(x[33]);
	tmp30r = tmp28r+Re(x[33]);
	tmp31i = -C19*Im(x[11])-C21*Im(x[99])+C13*Im(x[77])+C15*Im(x[33])-C17*Im(x[55]);
	tmp31r = -C19*Re(x[11])-C21*Re(x[99])+C13*Re(x[77])+C15*Re(x[33])-C17*Re(x[55]);
	tmp32i = tmp30i+Im(x[55]);
	tmp32r = tmp30r+Re(x[55]);
	tmp33i = -C21*Im(x[11])+C13*Im(x[99])+C15*Im(x[77])-C17*Im(x[33])-C19*Im(x[55]);
	tmp33r = -C21*Re(x[11])+C13*Re(x[99])+C15*Re(x[77])-C17*Re(x[33])-C19*Re(x[55]);
	tmp34i = -C14*Re(x[11])+C16*Re(x[99])+C18*Re(x[77])-C20*Re(x[33])-C22*Re(x[55]);
	tmp34r = C14*Im(x[11])-C16*Im(x[99])-C18*Im(x[77])+C20*Im(x[33])+C22*Im(x[55]);
	tmp35i = -C16*Re(x[11])+C18*Re(x[99])-C20*Re(x[77])+C22*Re(x[33])+C14*Re(x[55]);
	tmp35r = C16*Im(x[11])-C18*Im(x[99])+C20*Im(x[77])-C22*Im(x[33])-C14*Im(x[55]);
	tmp36i = -C18*Re(x[11])-C20*Re(x[99])+C22*Re(x[77])-C14*Re(x[33])+C16*Re(x[55]);
	tmp36r = C18*Im(x[11])+C20*Im(x[99])-C22*Im(x[77])+C14*Im(x[33])-C16*Im(x[55]);
	tmp37i = C20*Re(x[11])+C22*Re(x[99])-C14*Re(x[77])-C16*Re(x[33])+C18*Re(x[55]);
	tmp37r = -C20*Im(x[11])-C22*Im(x[99])+C14*Im(x[77])+C16*Im(x[33])-C18*Im(x[55]);
	tmp38i = -C22*Re(x[11])-C14*Re(x[99])-C16*Re(x[77])-C18*Re(x[33])-C20*Re(x[55]);
	tmp38r = C22*Im(x[11])+C14*Im(x[99])+C16*Im(x[77])+C18*Im(x[33])+C20*Im(x[55]);
	tmp39i = tmp25i+tmp34i;
	tmp39r = tmp25r+tmp34r;
	tmp40i = tmp25i-tmp34i;
	tmp40r = tmp25r-tmp34r;
	tmp41i = tmp27i+tmp35i;
	tmp41r = tmp27r+tmp35r;
	tmp42i = tmp27i-tmp35i;
	tmp42r = tmp27r-tmp35r;
	tmp43i = tmp29i+tmp36i;
	tmp43r = tmp29r+tmp36r;
	tmp44i = tmp29i-tmp36i;
	tmp44r = tmp29r-tmp36r;
	tmp45i = tmp31i+tmp37i;
	tmp45r = tmp31r+tmp37r;
	tmp46i = tmp31i-tmp37i;
	tmp46r = tmp31r-tmp37r;
	tmp47i = tmp33i+tmp38i;
	tmp47r = tmp33r+tmp38r;
	tmp48i = tmp33i-tmp38i;
	tmp48r = tmp33r-tmp38r;
	tmp49i = tmp24i+tmp32i;
	tmp49r = tmp24r+tmp32r;
	tmp50i = tmp24i-tmp32i;
	tmp50r = tmp24r-tmp32r;
	tmp51i = tmp23i+tmp48i;
	tmp51r = tmp23r+tmp48r;
	tmp52i = tmp23i-tmp48i;
	tmp52r = tmp23r-tmp48r;
	tmp53i = tmp14i+tmp39i;
	tmp53r = tmp14r+tmp39r;
	tmp54i = tmp14i-tmp39i;
	tmp54r = tmp14r-tmp39r;
	tmp55i = tmp19i+tmp44i;
	tmp55r = tmp19r+tmp44r;
	tmp56i = tmp19i-tmp44i;
	tmp56r = tmp19r-tmp44r;
	tmp57i = tmp16i+tmp41i;
	tmp57r = tmp16r+tmp41r;
	tmp58i = tmp16i-tmp41i;
	tmp58r = tmp16r-tmp41r;
	tmp59i = tmp20i+tmp45i;
	tmp59r = tmp20r+tmp45r;
	tmp60i = tmp20i-tmp45i;
	tmp60r = tmp20r-tmp45r;
	tmp61i = tmp21i+tmp46i;
	tmp61r = tmp21r+tmp46r;
	tmp62i = tmp21i-tmp46i;
	tmp62r = tmp21r-tmp46r;
	tmp63i = tmp17i+tmp42i;
	tmp63r = tmp17r+tmp42r;
	tmp64i = tmp17i-tmp42i;
	tmp64r = tmp17r-tmp42r;
	tmp65i = tmp18i+tmp43i;
	tmp65r = tmp18r+tmp43r;
	tmp66i = tmp18i-tmp43i;
	tmp66r = tmp18r-tmp43r;
	tmp67i = tmp15i+tmp40i;
	tmp67r = tmp15r+tmp40r;
	tmp68i = tmp15i-tmp40i;
	tmp68r = tmp15r-tmp40r;
	tmp69i = tmp22i+tmp47i;
	tmp69r = tmp22r+tmp47r;
	tmp70i = tmp22i-tmp47i;
	tmp70r = tmp22r-tmp47r;
	Re(x[0]) = tmp49r;
	Im(x[0]) = tmp49i;
	Re(x[11]) = tmp52r;
	Im(x[11]) = tmp52i;
	Re(x[22]) = tmp53r;
	Im(x[22]) = tmp53i;
	Re(x[33]) = tmp56r;
	Im(x[33]) = tmp56i;
	Re(x[44]) = tmp57r;
	Im(x[44]) = tmp57i;
	Re(x[55]) = tmp60r;
	Im(x[55]) = tmp60i;
	Re(x[66]) = tmp61r;
	Im(x[66]) = tmp61i;
	Re(x[77]) = tmp64r;
	Im(x[77]) = tmp64i;
	Re(x[88]) = tmp65r;
	Im(x[88]) = tmp65i;
	Re(x[99]) = tmp68r;
	Im(x[99]) = tmp68i;
	Re(x[110]) = tmp69r;
	Im(x[110]) = tmp69i;
	Re(x[121]) = tmp50r;
	Im(x[121]) = tmp50i;
	Re(x[132]) = tmp51r;
	Im(x[132]) = tmp51i;
	Re(x[143]) = tmp54r;
	Im(x[143]) = tmp54i;
	Re(x[154]) = tmp55r;
	Im(x[154]) = tmp55i;
	Re(x[165]) = tmp58r;
	Im(x[165]) = tmp58i;
	Re(x[176]) = tmp59r;
	Im(x[176]) = tmp59i;
	Re(x[187]) = tmp62r;
	Im(x[187]) = tmp62i;
	Re(x[198]) = tmp63r;
	Im(x[198]) = tmp63i;
	Re(x[209]) = tmp66r;
	Im(x[209]) = tmp66i;
	Re(x[220]) = tmp67r;
	Im(x[220]) = tmp67i;
	Re(x[231]) = tmp70r;
	Im(x[231]) = tmp70i;
}

/*
*	Number of additions = 324
*	Number of multiplications = 200
*	Number of sign changes = 0
*	Number of assigns = 228
*	Total number of operations = 752
*/
void	MIFFTC22(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp85r, tmp86r, tmp87r, tmp88r, tmp89r, tmp90r, tmp91r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp85i, tmp86i, tmp87i, tmp88i, tmp89i, tmp90i, tmp91i;

	const REAL	C19 =     0.14231483827329;	/* REALCONST	*/
	const REAL	C22 =     0.28173255684143;	/* REALCONST	*/
	const REAL	C15 =     0.41541501300189;	/* REALCONST	*/
	const REAL	C14 =      0.5406408174556;	/* REALCONST	*/
	const REAL	C17 =     0.65486073394529;	/* REALCONST	*/
	const REAL	C18 =     0.75574957435426;	/* REALCONST	*/
	const REAL	C13 =     0.84125353283118;	/* REALCONST	*/
	const REAL	C16 =     0.90963199535452;	/* REALCONST	*/
	const REAL	C21 =      0.9594929736145;	/* REALCONST	*/
	const REAL	C20 =     0.98982144188093;	/* REALCONST	*/

	tmp0i = Im(x[132])+Im(x[110]);
	tmp0r = Re(x[132])+Re(x[110]);
	tmp1i = Im(x[132])-Im(x[110]);
	tmp1r = Re(x[132])-Re(x[110]);
	tmp2i = Im(x[22])+Im(x[220]);
	tmp2r = Re(x[22])+Re(x[220]);
	tmp3i = Im(x[22])-Im(x[220]);
	tmp3r = Re(x[22])-Re(x[220]);
	tmp4i = Im(x[44])+Im(x[198]);
	tmp4r = Re(x[44])+Re(x[198]);
	tmp5i = Im(x[44])-Im(x[198]);
	tmp5r = Re(x[44])-Re(x[198]);
	tmp6i = Im(x[88])+Im(x[154]);
	tmp6r = Re(x[88])+Re(x[154]);
	tmp7i = Im(x[88])-Im(x[154]);
	tmp7r = Re(x[88])-Re(x[154]);
	tmp8i = Im(x[176])+Im(x[66]);
	tmp8r = Re(x[176])+Re(x[66]);
	tmp9i = Im(x[176])-Im(x[66]);
	tmp9r = Re(x[176])-Re(x[66]);
	tmp10i = C13*tmp0i+C15*tmp2i-C17*tmp4i-C19*tmp6i-C21*tmp8i+Im(x[0]);
	tmp10r = C13*tmp0r+C15*tmp2r-C17*tmp4r-C19*tmp6r-C21*tmp8r+Re(x[0]);
	tmp11i = tmp0i+tmp2i;
	tmp11r = tmp0r+tmp2r;
	tmp12i = C15*tmp0i-C17*tmp2i-C19*tmp4i-C21*tmp6i+C13*tmp8i+Im(x[0]);
	tmp12r = C15*tmp0r-C17*tmp2r-C19*tmp4r-C21*tmp6r+C13*tmp8r+Re(x[0]);
	tmp13i = tmp11i+tmp4i;
	tmp13r = tmp11r+tmp4r;
	tmp14i = -C17*tmp0i-C19*tmp2i-C21*tmp4i+C13*tmp6i+C15*tmp8i+Im(x[0]);
	tmp14r = -C17*tmp0r-C19*tmp2r-C21*tmp4r+C13*tmp6r+C15*tmp8r+Re(x[0]);
	tmp15i = tmp13i+tmp6i;
	tmp15r = tmp13r+tmp6r;
	tmp16i = -C19*tmp0i-C21*tmp2i+C13*tmp4i+C15*tmp6i-C17*tmp8i+Im(x[0]);
	tmp16r = -C19*tmp0r-C21*tmp2r+C13*tmp4r+C15*tmp6r-C17*tmp8r+Re(x[0]);
	tmp17i = tmp15i+tmp8i;
	tmp17r = tmp15r+tmp8r;
	tmp18i = -C21*tmp0i+C13*tmp2i+C15*tmp4i-C17*tmp6i-C19*tmp8i+Im(x[0]);
	tmp18r = -C21*tmp0r+C13*tmp2r+C15*tmp4r-C17*tmp6r-C19*tmp8r+Re(x[0]);
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
	tmp34i = Im(x[0])+tmp17i;
	tmp34r = Re(x[0])+tmp17r;
	tmp35i = Im(x[11])+Im(x[231]);
	tmp35r = Re(x[11])+Re(x[231]);
	tmp36i = Im(x[11])-Im(x[231]);
	tmp36r = Re(x[11])-Re(x[231]);
	tmp37i = Im(x[143])+Im(x[99]);
	tmp37r = Re(x[143])+Re(x[99]);
	tmp38i = Im(x[143])-Im(x[99]);
	tmp38r = Re(x[143])-Re(x[99]);
	tmp39i = Im(x[165])+Im(x[77]);
	tmp39r = Re(x[165])+Re(x[77]);
	tmp40i = Im(x[165])-Im(x[77]);
	tmp40r = Re(x[165])-Re(x[77]);
	tmp41i = Im(x[209])+Im(x[33]);
	tmp41r = Re(x[209])+Re(x[33]);
	tmp42i = Im(x[209])-Im(x[33]);
	tmp42r = Re(x[209])-Re(x[33]);
	tmp43i = Im(x[55])+Im(x[187]);
	tmp43r = Re(x[55])+Re(x[187]);
	tmp44i = Im(x[55])-Im(x[187]);
	tmp44r = Re(x[55])-Re(x[187]);
	tmp45i = C13*tmp35i+C15*tmp37i-C17*tmp39i-C19*tmp41i-C21*tmp43i+Im(x[121]);
	tmp45r = C13*tmp35r+C15*tmp37r-C17*tmp39r-C19*tmp41r-C21*tmp43r+Re(x[121]);
	tmp46i = tmp35i+tmp37i;
	tmp46r = tmp35r+tmp37r;
	tmp47i = C15*tmp35i-C17*tmp37i-C19*tmp39i-C21*tmp41i+C13*tmp43i+Im(x[121]);
	tmp47r = C15*tmp35r-C17*tmp37r-C19*tmp39r-C21*tmp41r+C13*tmp43r+Re(x[121]);
	tmp48i = tmp46i+tmp39i;
	tmp48r = tmp46r+tmp39r;
	tmp49i = -C17*tmp35i-C19*tmp37i-C21*tmp39i+C13*tmp41i+C15*tmp43i+Im(x[121]);
	tmp49r = -C17*tmp35r-C19*tmp37r-C21*tmp39r+C13*tmp41r+C15*tmp43r+Re(x[121]);
	tmp50i = tmp48i+tmp41i;
	tmp50r = tmp48r+tmp41r;
	tmp51i = -C19*tmp35i-C21*tmp37i+C13*tmp39i+C15*tmp41i-C17*tmp43i+Im(x[121]);
	tmp51r = -C19*tmp35r-C21*tmp37r+C13*tmp39r+C15*tmp41r-C17*tmp43r+Re(x[121]);
	tmp52i = tmp50i+tmp43i;
	tmp52r = tmp50r+tmp43r;
	tmp53i = -C21*tmp35i+C13*tmp37i+C15*tmp39i-C17*tmp41i-C19*tmp43i+Im(x[121]);
	tmp53r = -C21*tmp35r+C13*tmp37r+C15*tmp39r-C17*tmp41r-C19*tmp43r+Re(x[121]);
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
	tmp69i = Im(x[121])+tmp52i;
	tmp69r = Re(x[121])+tmp52r;
	tmp70i = tmp34i+tmp69i;
	tmp70r = tmp34r+tmp69r;
	tmp71i = tmp34i-tmp69i;
	tmp71r = tmp34r-tmp69r;
	tmp72i = tmp33i+tmp68i;
	tmp72r = tmp33r+tmp68r;
	tmp73i = tmp33i-tmp68i;
	tmp73r = tmp33r-tmp68r;
	tmp74i = tmp24i+tmp59i;
	tmp74r = tmp24r+tmp59r;
	tmp75i = tmp24i-tmp59i;
	tmp75r = tmp24r-tmp59r;
	tmp76i = tmp29i+tmp64i;
	tmp76r = tmp29r+tmp64r;
	tmp77i = tmp29i-tmp64i;
	tmp77r = tmp29r-tmp64r;
	tmp78i = tmp26i+tmp61i;
	tmp78r = tmp26r+tmp61r;
	tmp79i = tmp26i-tmp61i;
	tmp79r = tmp26r-tmp61r;
	tmp80i = tmp30i+tmp65i;
	tmp80r = tmp30r+tmp65r;
	tmp81i = tmp30i-tmp65i;
	tmp81r = tmp30r-tmp65r;
	tmp82i = tmp31i+tmp66i;
	tmp82r = tmp31r+tmp66r;
	tmp83i = tmp31i-tmp66i;
	tmp83r = tmp31r-tmp66r;
	tmp84i = tmp27i+tmp62i;
	tmp84r = tmp27r+tmp62r;
	tmp85i = tmp27i-tmp62i;
	tmp85r = tmp27r-tmp62r;
	tmp86i = tmp28i+tmp63i;
	tmp86r = tmp28r+tmp63r;
	tmp87i = tmp28i-tmp63i;
	tmp87r = tmp28r-tmp63r;
	tmp88i = tmp25i+tmp60i;
	tmp88r = tmp25r+tmp60r;
	tmp89i = tmp25i-tmp60i;
	tmp89r = tmp25r-tmp60r;
	tmp90i = tmp32i+tmp67i;
	tmp90r = tmp32r+tmp67r;
	tmp91i = tmp32i-tmp67i;
	tmp91r = tmp32r-tmp67r;
	Re(x[0]) = tmp70r;
	Im(x[0]) = tmp70i;
	Re(x[11]) = tmp73r;
	Im(x[11]) = tmp73i;
	Re(x[22]) = tmp74r;
	Im(x[22]) = tmp74i;
	Re(x[33]) = tmp77r;
	Im(x[33]) = tmp77i;
	Re(x[44]) = tmp78r;
	Im(x[44]) = tmp78i;
	Re(x[55]) = tmp81r;
	Im(x[55]) = tmp81i;
	Re(x[66]) = tmp82r;
	Im(x[66]) = tmp82i;
	Re(x[77]) = tmp85r;
	Im(x[77]) = tmp85i;
	Re(x[88]) = tmp86r;
	Im(x[88]) = tmp86i;
	Re(x[99]) = tmp89r;
	Im(x[99]) = tmp89i;
	Re(x[110]) = tmp90r;
	Im(x[110]) = tmp90i;
	Re(x[121]) = tmp71r;
	Im(x[121]) = tmp71i;
	Re(x[132]) = tmp72r;
	Im(x[132]) = tmp72i;
	Re(x[143]) = tmp75r;
	Im(x[143]) = tmp75i;
	Re(x[154]) = tmp76r;
	Im(x[154]) = tmp76i;
	Re(x[165]) = tmp79r;
	Im(x[165]) = tmp79i;
	Re(x[176]) = tmp80r;
	Im(x[176]) = tmp80i;
	Re(x[187]) = tmp83r;
	Im(x[187]) = tmp83i;
	Re(x[198]) = tmp84r;
	Im(x[198]) = tmp84i;
	Re(x[209]) = tmp87r;
	Im(x[209]) = tmp87i;
	Re(x[220]) = tmp88r;
	Im(x[220]) = tmp88i;
	Re(x[231]) = tmp91r;
	Im(x[231]) = tmp91i;
}

/*
*	Number of additions = 204
*	Number of multiplications = 44
*	Number of sign changes = 38
*	Number of assigns = 284
*	Total number of operations = 570
*/
void	MFFTC24(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp113r, tmp114r, tmp115r, tmp116r, tmp117r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp113i, tmp114i, tmp115i, tmp116i, tmp117i;

	const REAL	C3 =                  0.5;	/* REALCONST	*/
	const REAL	C2 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C6 =     0.86602540378444;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[72]);
	tmp0r = Re(x[0])+Re(x[72]);
	tmp1i = Im(x[0])-Im(x[72]);
	tmp1r = Re(x[0])-Re(x[72]);
	tmp2i = Im(x[0])+Re(x[72]);
	tmp2r = Re(x[0])-Im(x[72]);
	tmp3i = Im(x[0])-Re(x[72]);
	tmp3r = Re(x[0])+Im(x[72]);
	tmp4i = Im(x[108])+Im(x[36]);
	tmp4r = Re(x[108])+Re(x[36]);
	tmp5i = Im(x[108])-Im(x[36]);
	tmp5r = Re(x[108])-Re(x[36]);
	tmp6i = Im(x[108])-Re(x[36]);
	tmp6r = Re(x[108])+Im(x[36]);
	tmp7i = Im(x[108])+Re(x[36]);
	tmp7r = Re(x[108])-Im(x[36]);
	tmp8i = C2*(tmp6i-tmp6r);
	tmp8r = C2*(tmp6r+tmp6i);
	tmp9i = -C2*(tmp7i+tmp7r);
	tmp9r = -C2*(tmp7r-tmp7i);
	tmp10i = tmp0i+tmp4i;
	tmp10r = tmp0r+tmp4r;
	tmp11i = tmp0i-tmp4i;
	tmp11r = tmp0r-tmp4r;
	tmp12i = tmp2i+tmp8i;
	tmp12r = tmp2r+tmp8r;
	tmp13i = tmp2i-tmp8i;
	tmp13r = tmp2r-tmp8r;
	tmp14i = tmp1i-tmp5r;
	tmp14r = tmp1r+tmp5i;
	tmp15i = tmp1i+tmp5r;
	tmp15r = tmp1r-tmp5i;
	tmp16i = tmp3i+tmp9i;
	tmp16r = tmp3r+tmp9r;
	tmp17i = tmp3i-tmp9i;
	tmp17r = tmp3r-tmp9r;
	tmp18i = Im(x[48])+Im(x[120]);
	tmp18r = Re(x[48])+Re(x[120]);
	tmp19i = Im(x[48])-Im(x[120]);
	tmp19r = Re(x[48])-Re(x[120]);
	tmp20i = -(Im(x[48])+Re(x[120]));
	tmp20r = -(Re(x[48])-Im(x[120]));
	tmp21i = -(Im(x[48])-Re(x[120]));
	tmp21r = -(Re(x[48])+Im(x[120]));
	tmp22i = Im(x[12])+Im(x[84]);
	tmp22r = Re(x[12])+Re(x[84]);
	tmp23i = Im(x[12])-Im(x[84]);
	tmp23r = Re(x[12])-Re(x[84]);
	tmp24i = Im(x[12])+Re(x[84]);
	tmp24r = Re(x[12])-Im(x[84]);
	tmp25i = Im(x[12])-Re(x[84]);
	tmp25r = Re(x[12])+Im(x[84]);
	tmp26i = C2*(tmp24i-tmp24r);
	tmp26r = C2*(tmp24r+tmp24i);
	tmp27i = -C2*(tmp25i+tmp25r);
	tmp27r = -C2*(tmp25r-tmp25i);
	tmp28i = tmp18i+tmp22i;
	tmp28r = tmp18r+tmp22r;
	tmp29i = tmp18i-tmp22i;
	tmp29r = tmp18r-tmp22r;
	tmp30i = tmp20i+tmp26i;
	tmp30r = tmp20r+tmp26r;
	tmp31i = tmp20i-tmp26i;
	tmp31r = tmp20r-tmp26r;
	tmp32i = tmp19i-tmp23r;
	tmp32r = tmp19r+tmp23i;
	tmp33i = tmp19i+tmp23r;
	tmp33r = tmp19r-tmp23i;
	tmp34i = tmp21i+tmp27i;
	tmp34r = tmp21r+tmp27r;
	tmp35i = tmp21i-tmp27i;
	tmp35r = tmp21r-tmp27r;
	tmp36i = Im(x[96])+Im(x[24]);
	tmp36r = Re(x[96])+Re(x[24]);
	tmp37i = Im(x[96])-Im(x[24]);
	tmp37r = Re(x[96])-Re(x[24]);
	tmp38i = Im(x[96])-Re(x[24]);
	tmp38r = Re(x[96])+Im(x[24]);
	tmp39i = Im(x[96])+Re(x[24]);
	tmp39r = Re(x[96])-Im(x[24]);
	tmp40i = Im(x[60])+Im(x[132]);
	tmp40r = Re(x[60])+Re(x[132]);
	tmp41i = Im(x[60])-Im(x[132]);
	tmp41r = Re(x[60])-Re(x[132]);
	tmp42i = -(Im(x[60])+Re(x[132]));
	tmp42r = -(Re(x[60])-Im(x[132]));
	tmp43i = -(Im(x[60])-Re(x[132]));
	tmp43r = -(Re(x[60])+Im(x[132]));
	tmp44i = C2*(tmp42i-tmp42r);
	tmp44r = C2*(tmp42r+tmp42i);
	tmp45i = -C2*(tmp43i+tmp43r);
	tmp45r = -C2*(tmp43r-tmp43i);
	tmp46i = tmp36i+tmp40i;
	tmp46r = tmp36r+tmp40r;
	tmp47i = tmp36i-tmp40i;
	tmp47r = tmp36r-tmp40r;
	tmp48i = tmp38i+tmp44i;
	tmp48r = tmp38r+tmp44r;
	tmp49i = tmp38i-tmp44i;
	tmp49r = tmp38r-tmp44r;
	tmp50i = tmp37i-tmp41r;
	tmp50r = tmp37r+tmp41i;
	tmp51i = tmp37i+tmp41r;
	tmp51r = tmp37r-tmp41i;
	tmp52i = tmp39i+tmp45i;
	tmp52r = tmp39r+tmp45r;
	tmp53i = tmp39i-tmp45i;
	tmp53r = tmp39r-tmp45r;
	tmp54i = tmp28i+tmp46i;
	tmp54r = tmp28r+tmp46r;
	tmp55i = tmp28i-tmp46i;
	tmp55r = tmp28r-tmp46r;
	tmp56i = -C3*tmp54i;
	tmp56r = -C3*tmp54r;
	tmp57i = -C6*tmp55r;
	tmp57r = C6*tmp55i;
	tmp58i = tmp56i+tmp10i;
	tmp58r = tmp56r+tmp10r;
	tmp59i = tmp58i+tmp57i;
	tmp59r = tmp58r+tmp57r;
	tmp60i = tmp58i-tmp57i;
	tmp60r = tmp58r-tmp57r;
	tmp61i = tmp10i+tmp54i;
	tmp61r = tmp10r+tmp54r;
	tmp62i = tmp34i+tmp52i;
	tmp62r = tmp34r+tmp52r;
	tmp63i = tmp34i-tmp52i;
	tmp63r = tmp34r-tmp52r;
	tmp64i = -C3*tmp62i;
	tmp64r = -C3*tmp62r;
	tmp65i = -C6*tmp63r;
	tmp65r = C6*tmp63i;
	tmp66i = tmp64i+tmp16i;
	tmp66r = tmp64r+tmp16r;
	tmp67i = tmp66i+tmp65i;
	tmp67r = tmp66r+tmp65r;
	tmp68i = tmp66i-tmp65i;
	tmp68r = tmp66r-tmp65r;
	tmp69i = tmp16i+tmp62i;
	tmp69r = tmp16r+tmp62r;
	tmp70i = tmp33i+tmp51i;
	tmp70r = tmp33r+tmp51r;
	tmp71i = tmp33i-tmp51i;
	tmp71r = tmp33r-tmp51r;
	tmp72i = -C3*tmp70i;
	tmp72r = -C3*tmp70r;
	tmp73i = -C6*tmp71r;
	tmp73r = C6*tmp71i;
	tmp74i = tmp72i+tmp15i;
	tmp74r = tmp72r+tmp15r;
	tmp75i = tmp74i+tmp73i;
	tmp75r = tmp74r+tmp73r;
	tmp76i = tmp74i-tmp73i;
	tmp76r = tmp74r-tmp73r;
	tmp77i = tmp15i+tmp70i;
	tmp77r = tmp15r+tmp70r;
	tmp78i = tmp30i+tmp48i;
	tmp78r = tmp30r+tmp48r;
	tmp79i = tmp30i-tmp48i;
	tmp79r = tmp30r-tmp48r;
	tmp80i = -C3*tmp78i;
	tmp80r = -C3*tmp78r;
	tmp81i = -C6*tmp79r;
	tmp81r = C6*tmp79i;
	tmp82i = tmp80i+tmp12i;
	tmp82r = tmp80r+tmp12r;
	tmp83i = tmp82i+tmp81i;
	tmp83r = tmp82r+tmp81r;
	tmp84i = tmp82i-tmp81i;
	tmp84r = tmp82r-tmp81r;
	tmp85i = tmp12i+tmp78i;
	tmp85r = tmp12r+tmp78r;
	tmp86i = tmp29i+tmp47i;
	tmp86r = tmp29r+tmp47r;
	tmp87i = tmp29i-tmp47i;
	tmp87r = tmp29r-tmp47r;
	tmp88i = -C3*tmp86i;
	tmp88r = -C3*tmp86r;
	tmp89i = -C6*tmp87r;
	tmp89r = C6*tmp87i;
	tmp90i = tmp88i+tmp11i;
	tmp90r = tmp88r+tmp11r;
	tmp91i = tmp90i+tmp89i;
	tmp91r = tmp90r+tmp89r;
	tmp92i = tmp90i-tmp89i;
	tmp92r = tmp90r-tmp89r;
	tmp93i = tmp11i+tmp86i;
	tmp93r = tmp11r+tmp86r;
	tmp94i = tmp35i+tmp53i;
	tmp94r = tmp35r+tmp53r;
	tmp95i = tmp35i-tmp53i;
	tmp95r = tmp35r-tmp53r;
	tmp96i = -C3*tmp94i;
	tmp96r = -C3*tmp94r;
	tmp97i = -C6*tmp95r;
	tmp97r = C6*tmp95i;
	tmp98i = tmp96i+tmp17i;
	tmp98r = tmp96r+tmp17r;
	tmp99i = tmp98i+tmp97i;
	tmp99r = tmp98r+tmp97r;
	tmp100i = tmp98i-tmp97i;
	tmp100r = tmp98r-tmp97r;
	tmp101i = tmp17i+tmp94i;
	tmp101r = tmp17r+tmp94r;
	tmp102i = tmp32i+tmp50i;
	tmp102r = tmp32r+tmp50r;
	tmp103i = tmp32i-tmp50i;
	tmp103r = tmp32r-tmp50r;
	tmp104i = -C3*tmp102i;
	tmp104r = -C3*tmp102r;
	tmp105i = -C6*tmp103r;
	tmp105r = C6*tmp103i;
	tmp106i = tmp104i+tmp14i;
	tmp106r = tmp104r+tmp14r;
	tmp107i = tmp106i+tmp105i;
	tmp107r = tmp106r+tmp105r;
	tmp108i = tmp106i-tmp105i;
	tmp108r = tmp106r-tmp105r;
	tmp109i = tmp14i+tmp102i;
	tmp109r = tmp14r+tmp102r;
	tmp110i = tmp31i+tmp49i;
	tmp110r = tmp31r+tmp49r;
	tmp111i = tmp31i-tmp49i;
	tmp111r = tmp31r-tmp49r;
	tmp112i = -C3*tmp110i;
	tmp112r = -C3*tmp110r;
	tmp113i = -C6*tmp111r;
	tmp113r = C6*tmp111i;
	tmp114i = tmp112i+tmp13i;
	tmp114r = tmp112r+tmp13r;
	tmp115i = tmp114i+tmp113i;
	tmp115r = tmp114r+tmp113r;
	tmp116i = tmp114i-tmp113i;
	tmp116r = tmp114r-tmp113r;
	tmp117i = tmp13i+tmp110i;
	tmp117r = tmp13r+tmp110r;
	Re(x[0]) = tmp61r;
	Im(x[0]) = tmp61i;
	Re(x[12]) = tmp68r;
	Im(x[12]) = tmp68i;
	Re(x[24]) = tmp75r;
	Im(x[24]) = tmp75i;
	Re(x[36]) = tmp85r;
	Im(x[36]) = tmp85i;
	Re(x[48]) = tmp92r;
	Im(x[48]) = tmp92i;
	Re(x[60]) = tmp99r;
	Im(x[60]) = tmp99i;
	Re(x[72]) = tmp109r;
	Im(x[72]) = tmp109i;
	Re(x[84]) = tmp116r;
	Im(x[84]) = tmp116i;
	Re(x[96]) = tmp59r;
	Im(x[96]) = tmp59i;
	Re(x[108]) = tmp69r;
	Im(x[108]) = tmp69i;
	Re(x[120]) = tmp76r;
	Im(x[120]) = tmp76i;
	Re(x[132]) = tmp83r;
	Im(x[132]) = tmp83i;
	Re(x[144]) = tmp93r;
	Im(x[144]) = tmp93i;
	Re(x[156]) = tmp100r;
	Im(x[156]) = tmp100i;
	Re(x[168]) = tmp107r;
	Im(x[168]) = tmp107i;
	Re(x[180]) = tmp117r;
	Im(x[180]) = tmp117i;
	Re(x[192]) = tmp60r;
	Im(x[192]) = tmp60i;
	Re(x[204]) = tmp67r;
	Im(x[204]) = tmp67i;
	Re(x[216]) = tmp77r;
	Im(x[216]) = tmp77i;
	Re(x[228]) = tmp84r;
	Im(x[228]) = tmp84i;
	Re(x[240]) = tmp91r;
	Im(x[240]) = tmp91i;
	Re(x[252]) = tmp101r;
	Im(x[252]) = tmp101i;
	Re(x[264]) = tmp108r;
	Im(x[264]) = tmp108i;
	Re(x[276]) = tmp115r;
	Im(x[276]) = tmp115i;
}

/*
*	Number of additions = 252
*	Number of multiplications = 44
*	Number of sign changes = 30
*	Number of assigns = 332
*	Total number of operations = 658
*/
void	MIFFTC24(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp141r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp141i;

	const REAL	C3 =                  0.5;	/* REALCONST	*/
	const REAL	C2 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C6 =     0.86602540378444;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[144]);
	tmp0r = Re(x[0])+Re(x[144]);
	tmp1i = Im(x[0])-Im(x[144]);
	tmp1r = Re(x[0])-Re(x[144]);
	tmp2i = Im(x[216])+Im(x[72]);
	tmp2r = Re(x[216])+Re(x[72]);
	tmp3i = Im(x[216])-Im(x[72]);
	tmp3r = Re(x[216])-Re(x[72]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = tmp1i+tmp3r;
	tmp6r = tmp1r-tmp3i;
	tmp7i = tmp1i-tmp3r;
	tmp7r = tmp1r+tmp3i;
	tmp8i = Im(x[108])+Im(x[252]);
	tmp8r = Re(x[108])+Re(x[252]);
	tmp9i = Im(x[108])-Im(x[252]);
	tmp9r = Re(x[108])-Re(x[252]);
	tmp10i = Im(x[36])+Im(x[180]);
	tmp10r = Re(x[36])+Re(x[180]);
	tmp11i = Im(x[36])-Im(x[180]);
	tmp11r = Re(x[36])-Re(x[180]);
	tmp12i = tmp8i+tmp10i;
	tmp12r = tmp8r+tmp10r;
	tmp13i = tmp8i-tmp10i;
	tmp13r = tmp8r-tmp10r;
	tmp14i = tmp9i+tmp11r;
	tmp14r = tmp9r-tmp11i;
	tmp15i = tmp9i-tmp11r;
	tmp15r = tmp9r+tmp11i;
	tmp16i = C2*(tmp14i+tmp14r);
	tmp16r = C2*(tmp14r-tmp14i);
	tmp17i = -C2*(tmp15i-tmp15r);
	tmp17r = -C2*(tmp15r+tmp15i);
	tmp18i = tmp4i+tmp12i;
	tmp18r = tmp4r+tmp12r;
	tmp19i = tmp4i-tmp12i;
	tmp19r = tmp4r-tmp12r;
	tmp20i = tmp6i+tmp16i;
	tmp20r = tmp6r+tmp16r;
	tmp21i = tmp6i-tmp16i;
	tmp21r = tmp6r-tmp16r;
	tmp22i = tmp5i+tmp13r;
	tmp22r = tmp5r-tmp13i;
	tmp23i = tmp5i-tmp13r;
	tmp23r = tmp5r+tmp13i;
	tmp24i = tmp7i+tmp17i;
	tmp24r = tmp7r+tmp17r;
	tmp25i = tmp7i-tmp17i;
	tmp25r = tmp7r-tmp17r;
	tmp26i = Im(x[192])+Im(x[48]);
	tmp26r = Re(x[192])+Re(x[48]);
	tmp27i = Im(x[192])-Im(x[48]);
	tmp27r = Re(x[192])-Re(x[48]);
	tmp28i = Im(x[120])+Im(x[264]);
	tmp28r = Re(x[120])+Re(x[264]);
	tmp29i = Im(x[120])-Im(x[264]);
	tmp29r = Re(x[120])-Re(x[264]);
	tmp30i = tmp26i+tmp28i;
	tmp30r = tmp26r+tmp28r;
	tmp31i = tmp26i-tmp28i;
	tmp31r = tmp26r-tmp28r;
	tmp32i = tmp27i+tmp29r;
	tmp32r = tmp27r-tmp29i;
	tmp33i = tmp27i-tmp29r;
	tmp33r = tmp27r+tmp29i;
	tmp34i = Im(x[12])+Im(x[156]);
	tmp34r = Re(x[12])+Re(x[156]);
	tmp35i = Im(x[12])-Im(x[156]);
	tmp35r = Re(x[12])-Re(x[156]);
	tmp36i = Im(x[228])+Im(x[84]);
	tmp36r = Re(x[228])+Re(x[84]);
	tmp37i = Im(x[228])-Im(x[84]);
	tmp37r = Re(x[228])-Re(x[84]);
	tmp38i = tmp34i+tmp36i;
	tmp38r = tmp34r+tmp36r;
	tmp39i = tmp34i-tmp36i;
	tmp39r = tmp34r-tmp36r;
	tmp40i = tmp35i+tmp37r;
	tmp40r = tmp35r-tmp37i;
	tmp41i = tmp35i-tmp37r;
	tmp41r = tmp35r+tmp37i;
	tmp42i = C2*(tmp40i+tmp40r);
	tmp42r = C2*(tmp40r-tmp40i);
	tmp43i = -C2*(tmp41i-tmp41r);
	tmp43r = -C2*(tmp41r+tmp41i);
	tmp44i = tmp30i+tmp38i;
	tmp44r = tmp30r+tmp38r;
	tmp45i = tmp30i-tmp38i;
	tmp45r = tmp30r-tmp38r;
	tmp46i = tmp32i+tmp42i;
	tmp46r = tmp32r+tmp42r;
	tmp47i = tmp32i-tmp42i;
	tmp47r = tmp32r-tmp42r;
	tmp48i = tmp31i+tmp39r;
	tmp48r = tmp31r-tmp39i;
	tmp49i = tmp31i-tmp39r;
	tmp49r = tmp31r+tmp39i;
	tmp50i = tmp33i+tmp43i;
	tmp50r = tmp33r+tmp43r;
	tmp51i = tmp33i-tmp43i;
	tmp51r = tmp33r-tmp43r;
	tmp52i = Im(x[96])+Im(x[240]);
	tmp52r = Re(x[96])+Re(x[240]);
	tmp53i = Im(x[96])-Im(x[240]);
	tmp53r = Re(x[96])-Re(x[240]);
	tmp54i = Im(x[24])+Im(x[168]);
	tmp54r = Re(x[24])+Re(x[168]);
	tmp55i = Im(x[24])-Im(x[168]);
	tmp55r = Re(x[24])-Re(x[168]);
	tmp56i = tmp52i+tmp54i;
	tmp56r = tmp52r+tmp54r;
	tmp57i = tmp52i-tmp54i;
	tmp57r = tmp52r-tmp54r;
	tmp58i = tmp53i+tmp55r;
	tmp58r = tmp53r-tmp55i;
	tmp59i = tmp53i-tmp55r;
	tmp59r = tmp53r+tmp55i;
	tmp60i = Im(x[204])+Im(x[60]);
	tmp60r = Re(x[204])+Re(x[60]);
	tmp61i = Im(x[204])-Im(x[60]);
	tmp61r = Re(x[204])-Re(x[60]);
	tmp62i = Im(x[132])+Im(x[276]);
	tmp62r = Re(x[132])+Re(x[276]);
	tmp63i = Im(x[132])-Im(x[276]);
	tmp63r = Re(x[132])-Re(x[276]);
	tmp64i = tmp60i+tmp62i;
	tmp64r = tmp60r+tmp62r;
	tmp65i = tmp60i-tmp62i;
	tmp65r = tmp60r-tmp62r;
	tmp66i = tmp61i+tmp63r;
	tmp66r = tmp61r-tmp63i;
	tmp67i = tmp61i-tmp63r;
	tmp67r = tmp61r+tmp63i;
	tmp68i = C2*(tmp66i+tmp66r);
	tmp68r = C2*(tmp66r-tmp66i);
	tmp69i = -C2*(tmp67i-tmp67r);
	tmp69r = -C2*(tmp67r+tmp67i);
	tmp70i = tmp56i+tmp64i;
	tmp70r = tmp56r+tmp64r;
	tmp71i = tmp56i-tmp64i;
	tmp71r = tmp56r-tmp64r;
	tmp72i = tmp58i+tmp68i;
	tmp72r = tmp58r+tmp68r;
	tmp73i = tmp58i-tmp68i;
	tmp73r = tmp58r-tmp68r;
	tmp74i = tmp57i+tmp65r;
	tmp74r = tmp57r-tmp65i;
	tmp75i = tmp57i-tmp65r;
	tmp75r = tmp57r+tmp65i;
	tmp76i = tmp59i+tmp69i;
	tmp76r = tmp59r+tmp69r;
	tmp77i = tmp59i-tmp69i;
	tmp77r = tmp59r-tmp69r;
	tmp78i = tmp44i+tmp70i;
	tmp78r = tmp44r+tmp70r;
	tmp79i = tmp44i-tmp70i;
	tmp79r = tmp44r-tmp70r;
	tmp80i = -C3*tmp78i;
	tmp80r = -C3*tmp78r;
	tmp81i = C6*tmp79r;
	tmp81r = -C6*tmp79i;
	tmp82i = tmp80i+tmp18i;
	tmp82r = tmp80r+tmp18r;
	tmp83i = tmp82i+tmp81i;
	tmp83r = tmp82r+tmp81r;
	tmp84i = tmp82i-tmp81i;
	tmp84r = tmp82r-tmp81r;
	tmp85i = tmp18i+tmp78i;
	tmp85r = tmp18r+tmp78r;
	tmp86i = tmp50i+tmp76i;
	tmp86r = tmp50r+tmp76r;
	tmp87i = tmp50i-tmp76i;
	tmp87r = tmp50r-tmp76r;
	tmp88i = -C3*tmp86i;
	tmp88r = -C3*tmp86r;
	tmp89i = C6*tmp87r;
	tmp89r = -C6*tmp87i;
	tmp90i = tmp88i+tmp24i;
	tmp90r = tmp88r+tmp24r;
	tmp91i = tmp90i+tmp89i;
	tmp91r = tmp90r+tmp89r;
	tmp92i = tmp90i-tmp89i;
	tmp92r = tmp90r-tmp89r;
	tmp93i = tmp24i+tmp86i;
	tmp93r = tmp24r+tmp86r;
	tmp94i = tmp49i+tmp75i;
	tmp94r = tmp49r+tmp75r;
	tmp95i = tmp49i-tmp75i;
	tmp95r = tmp49r-tmp75r;
	tmp96i = -C3*tmp94i;
	tmp96r = -C3*tmp94r;
	tmp97i = C6*tmp95r;
	tmp97r = -C6*tmp95i;
	tmp98i = tmp96i+tmp23i;
	tmp98r = tmp96r+tmp23r;
	tmp99i = tmp98i+tmp97i;
	tmp99r = tmp98r+tmp97r;
	tmp100i = tmp98i-tmp97i;
	tmp100r = tmp98r-tmp97r;
	tmp101i = tmp23i+tmp94i;
	tmp101r = tmp23r+tmp94r;
	tmp102i = tmp46i+tmp72i;
	tmp102r = tmp46r+tmp72r;
	tmp103i = tmp46i-tmp72i;
	tmp103r = tmp46r-tmp72r;
	tmp104i = -C3*tmp102i;
	tmp104r = -C3*tmp102r;
	tmp105i = C6*tmp103r;
	tmp105r = -C6*tmp103i;
	tmp106i = tmp104i+tmp20i;
	tmp106r = tmp104r+tmp20r;
	tmp107i = tmp106i+tmp105i;
	tmp107r = tmp106r+tmp105r;
	tmp108i = tmp106i-tmp105i;
	tmp108r = tmp106r-tmp105r;
	tmp109i = tmp20i+tmp102i;
	tmp109r = tmp20r+tmp102r;
	tmp110i = tmp45i+tmp71i;
	tmp110r = tmp45r+tmp71r;
	tmp111i = tmp45i-tmp71i;
	tmp111r = tmp45r-tmp71r;
	tmp112i = -C3*tmp110i;
	tmp112r = -C3*tmp110r;
	tmp113i = C6*tmp111r;
	tmp113r = -C6*tmp111i;
	tmp114i = tmp112i+tmp19i;
	tmp114r = tmp112r+tmp19r;
	tmp115i = tmp114i+tmp113i;
	tmp115r = tmp114r+tmp113r;
	tmp116i = tmp114i-tmp113i;
	tmp116r = tmp114r-tmp113r;
	tmp117i = tmp19i+tmp110i;
	tmp117r = tmp19r+tmp110r;
	tmp118i = tmp51i+tmp77i;
	tmp118r = tmp51r+tmp77r;
	tmp119i = tmp51i-tmp77i;
	tmp119r = tmp51r-tmp77r;
	tmp120i = -C3*tmp118i;
	tmp120r = -C3*tmp118r;
	tmp121i = C6*tmp119r;
	tmp121r = -C6*tmp119i;
	tmp122i = tmp120i+tmp25i;
	tmp122r = tmp120r+tmp25r;
	tmp123i = tmp122i+tmp121i;
	tmp123r = tmp122r+tmp121r;
	tmp124i = tmp122i-tmp121i;
	tmp124r = tmp122r-tmp121r;
	tmp125i = tmp25i+tmp118i;
	tmp125r = tmp25r+tmp118r;
	tmp126i = tmp48i+tmp74i;
	tmp126r = tmp48r+tmp74r;
	tmp127i = tmp48i-tmp74i;
	tmp127r = tmp48r-tmp74r;
	tmp128i = -C3*tmp126i;
	tmp128r = -C3*tmp126r;
	tmp129i = C6*tmp127r;
	tmp129r = -C6*tmp127i;
	tmp130i = tmp128i+tmp22i;
	tmp130r = tmp128r+tmp22r;
	tmp131i = tmp130i+tmp129i;
	tmp131r = tmp130r+tmp129r;
	tmp132i = tmp130i-tmp129i;
	tmp132r = tmp130r-tmp129r;
	tmp133i = tmp22i+tmp126i;
	tmp133r = tmp22r+tmp126r;
	tmp134i = tmp47i+tmp73i;
	tmp134r = tmp47r+tmp73r;
	tmp135i = tmp47i-tmp73i;
	tmp135r = tmp47r-tmp73r;
	tmp136i = -C3*tmp134i;
	tmp136r = -C3*tmp134r;
	tmp137i = C6*tmp135r;
	tmp137r = -C6*tmp135i;
	tmp138i = tmp136i+tmp21i;
	tmp138r = tmp136r+tmp21r;
	tmp139i = tmp138i+tmp137i;
	tmp139r = tmp138r+tmp137r;
	tmp140i = tmp138i-tmp137i;
	tmp140r = tmp138r-tmp137r;
	tmp141i = tmp21i+tmp134i;
	tmp141r = tmp21r+tmp134r;
	Re(x[0]) = tmp85r;
	Im(x[0]) = tmp85i;
	Re(x[12]) = tmp92r;
	Im(x[12]) = tmp92i;
	Re(x[24]) = tmp99r;
	Im(x[24]) = tmp99i;
	Re(x[36]) = tmp109r;
	Im(x[36]) = tmp109i;
	Re(x[48]) = tmp116r;
	Im(x[48]) = tmp116i;
	Re(x[60]) = tmp123r;
	Im(x[60]) = tmp123i;
	Re(x[72]) = tmp133r;
	Im(x[72]) = tmp133i;
	Re(x[84]) = tmp140r;
	Im(x[84]) = tmp140i;
	Re(x[96]) = tmp83r;
	Im(x[96]) = tmp83i;
	Re(x[108]) = tmp93r;
	Im(x[108]) = tmp93i;
	Re(x[120]) = tmp100r;
	Im(x[120]) = tmp100i;
	Re(x[132]) = tmp107r;
	Im(x[132]) = tmp107i;
	Re(x[144]) = tmp117r;
	Im(x[144]) = tmp117i;
	Re(x[156]) = tmp124r;
	Im(x[156]) = tmp124i;
	Re(x[168]) = tmp131r;
	Im(x[168]) = tmp131i;
	Re(x[180]) = tmp141r;
	Im(x[180]) = tmp141i;
	Re(x[192]) = tmp84r;
	Im(x[192]) = tmp84i;
	Re(x[204]) = tmp91r;
	Im(x[204]) = tmp91i;
	Re(x[216]) = tmp101r;
	Im(x[216]) = tmp101i;
	Re(x[228]) = tmp108r;
	Im(x[228]) = tmp108i;
	Re(x[240]) = tmp115r;
	Im(x[240]) = tmp115i;
	Re(x[252]) = tmp125r;
	Im(x[252]) = tmp125i;
	Re(x[264]) = tmp132r;
	Im(x[264]) = tmp132i;
	Re(x[276]) = tmp139r;
	Im(x[276]) = tmp139i;
}

/*
*	Number of additions = 374
*	Number of multiplications = 288
*	Number of sign changes = 0
*	Number of assigns = 222
*	Total number of operations = 884
*/
void	MFFTC26(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp78r, tmp79r, tmp80r, tmp81r, tmp82r, tmp83r, tmp84r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp78i, tmp79i, tmp80i, tmp81i, tmp82i, tmp83i, tmp84i;

	const REAL	C23 =     0.12053668025532;	/* REALCONST	*/
	const REAL	C26 =     0.23931566428756;	/* REALCONST	*/
	const REAL	C19 =     0.35460488704254;	/* REALCONST	*/
	const REAL	C16 =     0.46472317204377;	/* REALCONST	*/
	const REAL	C17 =     0.56806474673116;	/* REALCONST	*/
	const REAL	C22 =      0.6631226582408;	/* REALCONST	*/
	const REAL	C21 =      0.7485107481711;	/* REALCONST	*/
	const REAL	C18 =     0.82298386589366;	/* REALCONST	*/
	const REAL	C15 =     0.88545602565321;	/* REALCONST	*/
	const REAL	C20 =     0.93501624268541;	/* REALCONST	*/
	const REAL	C25 =     0.97094181742605;	/* REALCONST	*/
	const REAL	C24 =     0.99270887409805;	/* REALCONST	*/

	tmp0i = C15*Im(x[156])+C17*Im(x[26])-C19*Im(x[52])-C21*Im(x[104])+C23*Im(x[130])-C25*Im(x[78])+Im(x[0]);
	tmp0r = C15*Re(x[156])+C17*Re(x[26])-C19*Re(x[52])-C21*Re(x[104])+C23*Re(x[130])-C25*Re(x[78])+Re(x[0]);
	tmp1i = Im(x[156])+Im(x[26]);
	tmp1r = Re(x[156])+Re(x[26]);
	tmp2i = C17*Im(x[156])-C19*Im(x[26])-C21*Im(x[52])+C23*Im(x[104])-C25*Im(x[130])+C15*Im(x[78])+Im(x[0]);
	tmp2r = C17*Re(x[156])-C19*Re(x[26])-C21*Re(x[52])+C23*Re(x[104])-C25*Re(x[130])+C15*Re(x[78])+Re(x[0]);
	tmp3i = tmp1i+Im(x[52]);
	tmp3r = tmp1r+Re(x[52]);
	tmp4i = -C19*Im(x[156])-C21*Im(x[26])+C23*Im(x[52])-C25*Im(x[104])+C15*Im(x[130])+C17*Im(x[78])+Im(x[0]);
	tmp4r = -C19*Re(x[156])-C21*Re(x[26])+C23*Re(x[52])-C25*Re(x[104])+C15*Re(x[130])+C17*Re(x[78])+Re(x[0]);
	tmp5i = tmp3i+Im(x[104]);
	tmp5r = tmp3r+Re(x[104]);
	tmp6i = -C21*Im(x[156])+C23*Im(x[26])-C25*Im(x[52])+C15*Im(x[104])+C17*Im(x[130])-C19*Im(x[78])+Im(x[0]);
	tmp6r = -C21*Re(x[156])+C23*Re(x[26])-C25*Re(x[52])+C15*Re(x[104])+C17*Re(x[130])-C19*Re(x[78])+Re(x[0]);
	tmp7i = tmp5i+Im(x[130]);
	tmp7r = tmp5r+Re(x[130]);
	tmp8i = C23*Im(x[156])-C25*Im(x[26])+C15*Im(x[52])+C17*Im(x[104])-C19*Im(x[130])-C21*Im(x[78])+Im(x[0]);
	tmp8r = C23*Re(x[156])-C25*Re(x[26])+C15*Re(x[52])+C17*Re(x[104])-C19*Re(x[130])-C21*Re(x[78])+Re(x[0]);
	tmp9i = tmp7i+Im(x[78]);
	tmp9r = tmp7r+Re(x[78]);
	tmp10i = -C25*Im(x[156])+C15*Im(x[26])+C17*Im(x[52])-C19*Im(x[104])-C21*Im(x[130])+C23*Im(x[78])+Im(x[0]);
	tmp10r = -C25*Re(x[156])+C15*Re(x[26])+C17*Re(x[52])-C19*Re(x[104])-C21*Re(x[130])+C23*Re(x[78])+Re(x[0]);
	tmp11i = C16*Re(x[156])-C18*Re(x[26])-C20*Re(x[52])+C22*Re(x[104])+C24*Re(x[130])-C26*Re(x[78]);
	tmp11r = -C16*Im(x[156])+C18*Im(x[26])+C20*Im(x[52])-C22*Im(x[104])-C24*Im(x[130])+C26*Im(x[78]);
	tmp12i = C18*Re(x[156])-C20*Re(x[26])+C22*Re(x[52])-C24*Re(x[104])+C26*Re(x[130])+C16*Re(x[78]);
	tmp12r = -C18*Im(x[156])+C20*Im(x[26])-C22*Im(x[52])+C24*Im(x[104])-C26*Im(x[130])-C16*Im(x[78]);
	tmp13i = C20*Re(x[156])+C22*Re(x[26])-C24*Re(x[52])-C26*Re(x[104])-C16*Re(x[130])+C18*Re(x[78]);
	tmp13r = -C20*Im(x[156])-C22*Im(x[26])+C24*Im(x[52])+C26*Im(x[104])+C16*Im(x[130])-C18*Im(x[78]);
	tmp14i = -C22*Re(x[156])-C24*Re(x[26])-C26*Re(x[52])+C16*Re(x[104])-C18*Re(x[130])+C20*Re(x[78]);
	tmp14r = C22*Im(x[156])+C24*Im(x[26])+C26*Im(x[52])-C16*Im(x[104])+C18*Im(x[130])-C20*Im(x[78]);
	tmp15i = C24*Re(x[156])-C26*Re(x[26])+C16*Re(x[52])+C18*Re(x[104])-C20*Re(x[130])-C22*Re(x[78]);
	tmp15r = -C24*Im(x[156])+C26*Im(x[26])-C16*Im(x[52])-C18*Im(x[104])+C20*Im(x[130])+C22*Im(x[78]);
	tmp16i = C26*Re(x[156])+C16*Re(x[26])+C18*Re(x[52])+C20*Re(x[104])+C22*Re(x[130])+C24*Re(x[78]);
	tmp16r = -C26*Im(x[156])-C16*Im(x[26])-C18*Im(x[52])-C20*Im(x[104])-C22*Im(x[130])-C24*Im(x[78]);
	tmp17i = tmp0i+tmp11i;
	tmp17r = tmp0r+tmp11r;
	tmp18i = tmp0i-tmp11i;
	tmp18r = tmp0r-tmp11r;
	tmp19i = tmp2i+tmp12i;
	tmp19r = tmp2r+tmp12r;
	tmp20i = tmp2i-tmp12i;
	tmp20r = tmp2r-tmp12r;
	tmp21i = tmp4i+tmp13i;
	tmp21r = tmp4r+tmp13r;
	tmp22i = tmp4i-tmp13i;
	tmp22r = tmp4r-tmp13r;
	tmp23i = tmp6i+tmp14i;
	tmp23r = tmp6r+tmp14r;
	tmp24i = tmp6i-tmp14i;
	tmp24r = tmp6r-tmp14r;
	tmp25i = tmp8i+tmp15i;
	tmp25r = tmp8r+tmp15r;
	tmp26i = tmp8i-tmp15i;
	tmp26r = tmp8r-tmp15r;
	tmp27i = tmp10i+tmp16i;
	tmp27r = tmp10r+tmp16r;
	tmp28i = tmp10i-tmp16i;
	tmp28r = tmp10r-tmp16r;
	tmp29i = Im(x[0])+tmp9i;
	tmp29r = Re(x[0])+tmp9r;
	tmp30i = C15*Im(x[13])+C17*Im(x[143])-C19*Im(x[117])-C21*Im(x[65])+C23*Im(x[39])-C25*Im(x[91]);
	tmp30r = C15*Re(x[13])+C17*Re(x[143])-C19*Re(x[117])-C21*Re(x[65])+C23*Re(x[39])-C25*Re(x[91]);
	tmp31i = Im(x[13])+Im(x[143]);
	tmp31r = Re(x[13])+Re(x[143]);
	tmp32i = C17*Im(x[13])-C19*Im(x[143])-C21*Im(x[117])+C23*Im(x[65])-C25*Im(x[39])+C15*Im(x[91]);
	tmp32r = C17*Re(x[13])-C19*Re(x[143])-C21*Re(x[117])+C23*Re(x[65])-C25*Re(x[39])+C15*Re(x[91]);
	tmp33i = tmp31i+Im(x[117]);
	tmp33r = tmp31r+Re(x[117]);
	tmp34i = -C19*Im(x[13])-C21*Im(x[143])+C23*Im(x[117])-C25*Im(x[65])+C15*Im(x[39])+C17*Im(x[91]);
	tmp34r = -C19*Re(x[13])-C21*Re(x[143])+C23*Re(x[117])-C25*Re(x[65])+C15*Re(x[39])+C17*Re(x[91]);
	tmp35i = tmp33i+Im(x[65]);
	tmp35r = tmp33r+Re(x[65]);
	tmp36i = -C21*Im(x[13])+C23*Im(x[143])-C25*Im(x[117])+C15*Im(x[65])+C17*Im(x[39])-C19*Im(x[91]);
	tmp36r = -C21*Re(x[13])+C23*Re(x[143])-C25*Re(x[117])+C15*Re(x[65])+C17*Re(x[39])-C19*Re(x[91]);
	tmp37i = tmp35i+Im(x[39]);
	tmp37r = tmp35r+Re(x[39]);
	tmp38i = C23*Im(x[13])-C25*Im(x[143])+C15*Im(x[117])+C17*Im(x[65])-C19*Im(x[39])-C21*Im(x[91]);
	tmp38r = C23*Re(x[13])-C25*Re(x[143])+C15*Re(x[117])+C17*Re(x[65])-C19*Re(x[39])-C21*Re(x[91]);
	tmp39i = tmp37i+Im(x[91]);
	tmp39r = tmp37r+Re(x[91]);
	tmp40i = -C25*Im(x[13])+C15*Im(x[143])+C17*Im(x[117])-C19*Im(x[65])-C21*Im(x[39])+C23*Im(x[91]);
	tmp40r = -C25*Re(x[13])+C15*Re(x[143])+C17*Re(x[117])-C19*Re(x[65])-C21*Re(x[39])+C23*Re(x[91]);
	tmp41i = -C16*Re(x[13])+C18*Re(x[143])+C20*Re(x[117])-C22*Re(x[65])-C24*Re(x[39])+C26*Re(x[91]);
	tmp41r = C16*Im(x[13])-C18*Im(x[143])-C20*Im(x[117])+C22*Im(x[65])+C24*Im(x[39])-C26*Im(x[91]);
	tmp42i = -C18*Re(x[13])+C20*Re(x[143])-C22*Re(x[117])+C24*Re(x[65])-C26*Re(x[39])-C16*Re(x[91]);
	tmp42r = C18*Im(x[13])-C20*Im(x[143])+C22*Im(x[117])-C24*Im(x[65])+C26*Im(x[39])+C16*Im(x[91]);
	tmp43i = -C20*Re(x[13])-C22*Re(x[143])+C24*Re(x[117])+C26*Re(x[65])+C16*Re(x[39])-C18*Re(x[91]);
	tmp43r = C20*Im(x[13])+C22*Im(x[143])-C24*Im(x[117])-C26*Im(x[65])-C16*Im(x[39])+C18*Im(x[91]);
	tmp44i = C22*Re(x[13])+C24*Re(x[143])+C26*Re(x[117])-C16*Re(x[65])+C18*Re(x[39])-C20*Re(x[91]);
	tmp44r = -C22*Im(x[13])-C24*Im(x[143])-C26*Im(x[117])+C16*Im(x[65])-C18*Im(x[39])+C20*Im(x[91]);
	tmp45i = -C24*Re(x[13])+C26*Re(x[143])-C16*Re(x[117])-C18*Re(x[65])+C20*Re(x[39])+C22*Re(x[91]);
	tmp45r = C24*Im(x[13])-C26*Im(x[143])+C16*Im(x[117])+C18*Im(x[65])-C20*Im(x[39])-C22*Im(x[91]);
	tmp46i = -C26*Re(x[13])-C16*Re(x[143])-C18*Re(x[117])-C20*Re(x[65])-C22*Re(x[39])-C24*Re(x[91]);
	tmp46r = C26*Im(x[13])+C16*Im(x[143])+C18*Im(x[117])+C20*Im(x[65])+C22*Im(x[39])+C24*Im(x[91]);
	tmp47i = tmp30i+tmp41i;
	tmp47r = tmp30r+tmp41r;
	tmp48i = tmp30i-tmp41i;
	tmp48r = tmp30r-tmp41r;
	tmp49i = tmp32i+tmp42i;
	tmp49r = tmp32r+tmp42r;
	tmp50i = tmp32i-tmp42i;
	tmp50r = tmp32r-tmp42r;
	tmp51i = tmp34i+tmp43i;
	tmp51r = tmp34r+tmp43r;
	tmp52i = tmp34i-tmp43i;
	tmp52r = tmp34r-tmp43r;
	tmp53i = tmp36i+tmp44i;
	tmp53r = tmp36r+tmp44r;
	tmp54i = tmp36i-tmp44i;
	tmp54r = tmp36r-tmp44r;
	tmp55i = tmp38i+tmp45i;
	tmp55r = tmp38r+tmp45r;
	tmp56i = tmp38i-tmp45i;
	tmp56r = tmp38r-tmp45r;
	tmp57i = tmp40i+tmp46i;
	tmp57r = tmp40r+tmp46r;
	tmp58i = tmp40i-tmp46i;
	tmp58r = tmp40r-tmp46r;
	tmp59i = tmp29i+tmp39i;
	tmp59r = tmp29r+tmp39r;
	tmp60i = tmp29i-tmp39i;
	tmp60r = tmp29r-tmp39r;
	tmp61i = tmp28i+tmp58i;
	tmp61r = tmp28r+tmp58r;
	tmp62i = tmp28i-tmp58i;
	tmp62r = tmp28r-tmp58r;
	tmp63i = tmp17i+tmp47i;
	tmp63r = tmp17r+tmp47r;
	tmp64i = tmp17i-tmp47i;
	tmp64r = tmp17r-tmp47r;
	tmp65i = tmp23i+tmp53i;
	tmp65r = tmp23r+tmp53r;
	tmp66i = tmp23i-tmp53i;
	tmp66r = tmp23r-tmp53r;
	tmp67i = tmp19i+tmp49i;
	tmp67r = tmp19r+tmp49r;
	tmp68i = tmp19i-tmp49i;
	tmp68r = tmp19r-tmp49r;
	tmp69i = tmp22i+tmp52i;
	tmp69r = tmp22r+tmp52r;
	tmp70i = tmp22i-tmp52i;
	tmp70r = tmp22r-tmp52r;
	tmp71i = tmp25i+tmp55i;
	tmp71r = tmp25r+tmp55r;
	tmp72i = tmp25i-tmp55i;
	tmp72r = tmp25r-tmp55r;
	tmp73i = tmp26i+tmp56i;
	tmp73r = tmp26r+tmp56r;
	tmp74i = tmp26i-tmp56i;
	tmp74r = tmp26r-tmp56r;
	tmp75i = tmp21i+tmp51i;
	tmp75r = tmp21r+tmp51r;
	tmp76i = tmp21i-tmp51i;
	tmp76r = tmp21r-tmp51r;
	tmp77i = tmp20i+tmp50i;
	tmp77r = tmp20r+tmp50r;
	tmp78i = tmp20i-tmp50i;
	tmp78r = tmp20r-tmp50r;
	tmp79i = tmp24i+tmp54i;
	tmp79r = tmp24r+tmp54r;
	tmp80i = tmp24i-tmp54i;
	tmp80r = tmp24r-tmp54r;
	tmp81i = tmp18i+tmp48i;
	tmp81r = tmp18r+tmp48r;
	tmp82i = tmp18i-tmp48i;
	tmp82r = tmp18r-tmp48r;
	tmp83i = tmp27i+tmp57i;
	tmp83r = tmp27r+tmp57r;
	tmp84i = tmp27i-tmp57i;
	tmp84r = tmp27r-tmp57r;
	Re(x[0]) = tmp59r;
	Im(x[0]) = tmp59i;
	Re(x[13]) = tmp62r;
	Im(x[13]) = tmp62i;
	Re(x[26]) = tmp63r;
	Im(x[26]) = tmp63i;
	Re(x[39]) = tmp66r;
	Im(x[39]) = tmp66i;
	Re(x[52]) = tmp67r;
	Im(x[52]) = tmp67i;
	Re(x[65]) = tmp70r;
	Im(x[65]) = tmp70i;
	Re(x[78]) = tmp71r;
	Im(x[78]) = tmp71i;
	Re(x[91]) = tmp74r;
	Im(x[91]) = tmp74i;
	Re(x[104]) = tmp75r;
	Im(x[104]) = tmp75i;
	Re(x[117]) = tmp78r;
	Im(x[117]) = tmp78i;
	Re(x[130]) = tmp79r;
	Im(x[130]) = tmp79i;
	Re(x[143]) = tmp82r;
	Im(x[143]) = tmp82i;
	Re(x[156]) = tmp83r;
	Im(x[156]) = tmp83i;
	Re(x[169]) = tmp60r;
	Im(x[169]) = tmp60i;
	Re(x[182]) = tmp61r;
	Im(x[182]) = tmp61i;
	Re(x[195]) = tmp64r;
	Im(x[195]) = tmp64i;
	Re(x[208]) = tmp65r;
	Im(x[208]) = tmp65i;
	Re(x[221]) = tmp68r;
	Im(x[221]) = tmp68i;
	Re(x[234]) = tmp69r;
	Im(x[234]) = tmp69i;
	Re(x[247]) = tmp72r;
	Im(x[247]) = tmp72i;
	Re(x[260]) = tmp73r;
	Im(x[260]) = tmp73i;
	Re(x[273]) = tmp76r;
	Im(x[273]) = tmp76i;
	Re(x[286]) = tmp77r;
	Im(x[286]) = tmp77i;
	Re(x[299]) = tmp80r;
	Im(x[299]) = tmp80i;
	Re(x[312]) = tmp81r;
	Im(x[312]) = tmp81i;
	Re(x[325]) = tmp84r;
	Im(x[325]) = tmp84i;
}

/*
*	Number of additions = 436
*	Number of multiplications = 288
*	Number of sign changes = 0
*	Number of assigns = 272
*	Total number of operations = 996
*/
void	MIFFTC26(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp106r, tmp107r, tmp108r, tmp109r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp106i, tmp107i, tmp108i, tmp109i;

	const REAL	C23 =     0.12053668025532;	/* REALCONST	*/
	const REAL	C26 =     0.23931566428756;	/* REALCONST	*/
	const REAL	C19 =     0.35460488704254;	/* REALCONST	*/
	const REAL	C16 =     0.46472317204377;	/* REALCONST	*/
	const REAL	C17 =     0.56806474673116;	/* REALCONST	*/
	const REAL	C22 =      0.6631226582408;	/* REALCONST	*/
	const REAL	C21 =      0.7485107481711;	/* REALCONST	*/
	const REAL	C18 =     0.82298386589366;	/* REALCONST	*/
	const REAL	C15 =     0.88545602565321;	/* REALCONST	*/
	const REAL	C20 =     0.93501624268541;	/* REALCONST	*/
	const REAL	C25 =     0.97094181742605;	/* REALCONST	*/
	const REAL	C24 =     0.99270887409805;	/* REALCONST	*/

	tmp0i = Im(x[182])+Im(x[156]);
	tmp0r = Re(x[182])+Re(x[156]);
	tmp1i = Im(x[182])-Im(x[156]);
	tmp1r = Re(x[182])-Re(x[156]);
	tmp2i = Im(x[26])+Im(x[312]);
	tmp2r = Re(x[26])+Re(x[312]);
	tmp3i = Im(x[26])-Im(x[312]);
	tmp3r = Re(x[26])-Re(x[312]);
	tmp4i = Im(x[52])+Im(x[286]);
	tmp4r = Re(x[52])+Re(x[286]);
	tmp5i = Im(x[52])-Im(x[286]);
	tmp5r = Re(x[52])-Re(x[286]);
	tmp6i = Im(x[104])+Im(x[234]);
	tmp6r = Re(x[104])+Re(x[234]);
	tmp7i = Im(x[104])-Im(x[234]);
	tmp7r = Re(x[104])-Re(x[234]);
	tmp8i = Im(x[208])+Im(x[130]);
	tmp8r = Re(x[208])+Re(x[130]);
	tmp9i = Im(x[208])-Im(x[130]);
	tmp9r = Re(x[208])-Re(x[130]);
	tmp10i = Im(x[78])+Im(x[260]);
	tmp10r = Re(x[78])+Re(x[260]);
	tmp11i = Im(x[78])-Im(x[260]);
	tmp11r = Re(x[78])-Re(x[260]);
	tmp12i = C15*tmp0i+C17*tmp2i-C19*tmp4i-C21*tmp6i+C23*tmp8i-C25*tmp10i+Im(x[0]);
	tmp12r = C15*tmp0r+C17*tmp2r-C19*tmp4r-C21*tmp6r+C23*tmp8r-C25*tmp10r+Re(x[0]);
	tmp13i = tmp0i+tmp2i;
	tmp13r = tmp0r+tmp2r;
	tmp14i = C17*tmp0i-C19*tmp2i-C21*tmp4i+C23*tmp6i-C25*tmp8i+C15*tmp10i+Im(x[0]);
	tmp14r = C17*tmp0r-C19*tmp2r-C21*tmp4r+C23*tmp6r-C25*tmp8r+C15*tmp10r+Re(x[0]);
	tmp15i = tmp13i+tmp4i;
	tmp15r = tmp13r+tmp4r;
	tmp16i = -C19*tmp0i-C21*tmp2i+C23*tmp4i-C25*tmp6i+C15*tmp8i+C17*tmp10i+Im(x[0]);
	tmp16r = -C19*tmp0r-C21*tmp2r+C23*tmp4r-C25*tmp6r+C15*tmp8r+C17*tmp10r+Re(x[0]);
	tmp17i = tmp15i+tmp6i;
	tmp17r = tmp15r+tmp6r;
	tmp18i = -C21*tmp0i+C23*tmp2i-C25*tmp4i+C15*tmp6i+C17*tmp8i-C19*tmp10i+Im(x[0]);
	tmp18r = -C21*tmp0r+C23*tmp2r-C25*tmp4r+C15*tmp6r+C17*tmp8r-C19*tmp10r+Re(x[0]);
	tmp19i = tmp17i+tmp8i;
	tmp19r = tmp17r+tmp8r;
	tmp20i = C23*tmp0i-C25*tmp2i+C15*tmp4i+C17*tmp6i-C19*tmp8i-C21*tmp10i+Im(x[0]);
	tmp20r = C23*tmp0r-C25*tmp2r+C15*tmp4r+C17*tmp6r-C19*tmp8r-C21*tmp10r+Re(x[0]);
	tmp21i = tmp19i+tmp10i;
	tmp21r = tmp19r+tmp10r;
	tmp22i = -C25*tmp0i+C15*tmp2i+C17*tmp4i-C19*tmp6i-C21*tmp8i+C23*tmp10i+Im(x[0]);
	tmp22r = -C25*tmp0r+C15*tmp2r+C17*tmp4r-C19*tmp6r-C21*tmp8r+C23*tmp10r+Re(x[0]);
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
	tmp41i = Im(x[0])+tmp21i;
	tmp41r = Re(x[0])+tmp21r;
	tmp42i = Im(x[13])+Im(x[325]);
	tmp42r = Re(x[13])+Re(x[325]);
	tmp43i = Im(x[13])-Im(x[325]);
	tmp43r = Re(x[13])-Re(x[325]);
	tmp44i = Im(x[195])+Im(x[143]);
	tmp44r = Re(x[195])+Re(x[143]);
	tmp45i = Im(x[195])-Im(x[143]);
	tmp45r = Re(x[195])-Re(x[143]);
	tmp46i = Im(x[221])+Im(x[117]);
	tmp46r = Re(x[221])+Re(x[117]);
	tmp47i = Im(x[221])-Im(x[117]);
	tmp47r = Re(x[221])-Re(x[117]);
	tmp48i = Im(x[273])+Im(x[65]);
	tmp48r = Re(x[273])+Re(x[65]);
	tmp49i = Im(x[273])-Im(x[65]);
	tmp49r = Re(x[273])-Re(x[65]);
	tmp50i = Im(x[39])+Im(x[299]);
	tmp50r = Re(x[39])+Re(x[299]);
	tmp51i = Im(x[39])-Im(x[299]);
	tmp51r = Re(x[39])-Re(x[299]);
	tmp52i = Im(x[247])+Im(x[91]);
	tmp52r = Re(x[247])+Re(x[91]);
	tmp53i = Im(x[247])-Im(x[91]);
	tmp53r = Re(x[247])-Re(x[91]);
	tmp54i = C15*tmp42i+C17*tmp44i-C19*tmp46i-C21*tmp48i+C23*tmp50i-C25*tmp52i+Im(x[169]);
	tmp54r = C15*tmp42r+C17*tmp44r-C19*tmp46r-C21*tmp48r+C23*tmp50r-C25*tmp52r+Re(x[169]);
	tmp55i = tmp42i+tmp44i;
	tmp55r = tmp42r+tmp44r;
	tmp56i = C17*tmp42i-C19*tmp44i-C21*tmp46i+C23*tmp48i-C25*tmp50i+C15*tmp52i+Im(x[169]);
	tmp56r = C17*tmp42r-C19*tmp44r-C21*tmp46r+C23*tmp48r-C25*tmp50r+C15*tmp52r+Re(x[169]);
	tmp57i = tmp55i+tmp46i;
	tmp57r = tmp55r+tmp46r;
	tmp58i = -C19*tmp42i-C21*tmp44i+C23*tmp46i-C25*tmp48i+C15*tmp50i+C17*tmp52i+Im(x[169]);
	tmp58r = -C19*tmp42r-C21*tmp44r+C23*tmp46r-C25*tmp48r+C15*tmp50r+C17*tmp52r+Re(x[169]);
	tmp59i = tmp57i+tmp48i;
	tmp59r = tmp57r+tmp48r;
	tmp60i = -C21*tmp42i+C23*tmp44i-C25*tmp46i+C15*tmp48i+C17*tmp50i-C19*tmp52i+Im(x[169]);
	tmp60r = -C21*tmp42r+C23*tmp44r-C25*tmp46r+C15*tmp48r+C17*tmp50r-C19*tmp52r+Re(x[169]);
	tmp61i = tmp59i+tmp50i;
	tmp61r = tmp59r+tmp50r;
	tmp62i = C23*tmp42i-C25*tmp44i+C15*tmp46i+C17*tmp48i-C19*tmp50i-C21*tmp52i+Im(x[169]);
	tmp62r = C23*tmp42r-C25*tmp44r+C15*tmp46r+C17*tmp48r-C19*tmp50r-C21*tmp52r+Re(x[169]);
	tmp63i = tmp61i+tmp52i;
	tmp63r = tmp61r+tmp52r;
	tmp64i = -C25*tmp42i+C15*tmp44i+C17*tmp46i-C19*tmp48i-C21*tmp50i+C23*tmp52i+Im(x[169]);
	tmp64r = -C25*tmp42r+C15*tmp44r+C17*tmp46r-C19*tmp48r-C21*tmp50r+C23*tmp52r+Re(x[169]);
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
	tmp83i = Im(x[169])+tmp63i;
	tmp83r = Re(x[169])+tmp63r;
	tmp84i = tmp41i+tmp83i;
	tmp84r = tmp41r+tmp83r;
	tmp85i = tmp41i-tmp83i;
	tmp85r = tmp41r-tmp83r;
	tmp86i = tmp40i+tmp82i;
	tmp86r = tmp40r+tmp82r;
	tmp87i = tmp40i-tmp82i;
	tmp87r = tmp40r-tmp82r;
	tmp88i = tmp29i+tmp71i;
	tmp88r = tmp29r+tmp71r;
	tmp89i = tmp29i-tmp71i;
	tmp89r = tmp29r-tmp71r;
	tmp90i = tmp35i+tmp77i;
	tmp90r = tmp35r+tmp77r;
	tmp91i = tmp35i-tmp77i;
	tmp91r = tmp35r-tmp77r;
	tmp92i = tmp31i+tmp73i;
	tmp92r = tmp31r+tmp73r;
	tmp93i = tmp31i-tmp73i;
	tmp93r = tmp31r-tmp73r;
	tmp94i = tmp34i+tmp76i;
	tmp94r = tmp34r+tmp76r;
	tmp95i = tmp34i-tmp76i;
	tmp95r = tmp34r-tmp76r;
	tmp96i = tmp37i+tmp79i;
	tmp96r = tmp37r+tmp79r;
	tmp97i = tmp37i-tmp79i;
	tmp97r = tmp37r-tmp79r;
	tmp98i = tmp38i+tmp80i;
	tmp98r = tmp38r+tmp80r;
	tmp99i = tmp38i-tmp80i;
	tmp99r = tmp38r-tmp80r;
	tmp100i = tmp33i+tmp75i;
	tmp100r = tmp33r+tmp75r;
	tmp101i = tmp33i-tmp75i;
	tmp101r = tmp33r-tmp75r;
	tmp102i = tmp32i+tmp74i;
	tmp102r = tmp32r+tmp74r;
	tmp103i = tmp32i-tmp74i;
	tmp103r = tmp32r-tmp74r;
	tmp104i = tmp36i+tmp78i;
	tmp104r = tmp36r+tmp78r;
	tmp105i = tmp36i-tmp78i;
	tmp105r = tmp36r-tmp78r;
	tmp106i = tmp30i+tmp72i;
	tmp106r = tmp30r+tmp72r;
	tmp107i = tmp30i-tmp72i;
	tmp107r = tmp30r-tmp72r;
	tmp108i = tmp39i+tmp81i;
	tmp108r = tmp39r+tmp81r;
	tmp109i = tmp39i-tmp81i;
	tmp109r = tmp39r-tmp81r;
	Re(x[0]) = tmp84r;
	Im(x[0]) = tmp84i;
	Re(x[13]) = tmp87r;
	Im(x[13]) = tmp87i;
	Re(x[26]) = tmp88r;
	Im(x[26]) = tmp88i;
	Re(x[39]) = tmp91r;
	Im(x[39]) = tmp91i;
	Re(x[52]) = tmp92r;
	Im(x[52]) = tmp92i;
	Re(x[65]) = tmp95r;
	Im(x[65]) = tmp95i;
	Re(x[78]) = tmp96r;
	Im(x[78]) = tmp96i;
	Re(x[91]) = tmp99r;
	Im(x[91]) = tmp99i;
	Re(x[104]) = tmp100r;
	Im(x[104]) = tmp100i;
	Re(x[117]) = tmp103r;
	Im(x[117]) = tmp103i;
	Re(x[130]) = tmp104r;
	Im(x[130]) = tmp104i;
	Re(x[143]) = tmp107r;
	Im(x[143]) = tmp107i;
	Re(x[156]) = tmp108r;
	Im(x[156]) = tmp108i;
	Re(x[169]) = tmp85r;
	Im(x[169]) = tmp85i;
	Re(x[182]) = tmp86r;
	Im(x[182]) = tmp86i;
	Re(x[195]) = tmp89r;
	Im(x[195]) = tmp89i;
	Re(x[208]) = tmp90r;
	Im(x[208]) = tmp90i;
	Re(x[221]) = tmp93r;
	Im(x[221]) = tmp93i;
	Re(x[234]) = tmp94r;
	Im(x[234]) = tmp94i;
	Re(x[247]) = tmp97r;
	Im(x[247]) = tmp97i;
	Re(x[260]) = tmp98r;
	Im(x[260]) = tmp98i;
	Re(x[273]) = tmp101r;
	Im(x[273]) = tmp101i;
	Re(x[286]) = tmp102r;
	Im(x[286]) = tmp102i;
	Re(x[299]) = tmp105r;
	Im(x[299]) = tmp105i;
	Re(x[312]) = tmp106r;
	Im(x[312]) = tmp106i;
	Re(x[325]) = tmp109r;
	Im(x[325]) = tmp109i;
}

/*
*	Number of additions = 258
*	Number of multiplications = 108
*	Number of sign changes = 3
*	Number of assigns = 290
*	Total number of operations = 659
*/
void	MFFTC28(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp113r, tmp114r, tmp115r, tmp116r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp113i, tmp114i, tmp115i, tmp116i;

	const REAL	C13 =     0.22252093395631;	/* REALCONST	*/
	const REAL	C12 =     0.43388373911756;	/* REALCONST	*/
	const REAL	C9 =     0.62348980185873;	/* REALCONST	*/
	const REAL	C10 =     0.78183148246803;	/* REALCONST	*/
	const REAL	C11 =     0.90096886790242;	/* REALCONST	*/
	const REAL	C14 =     0.97492791218182;	/* REALCONST	*/

	tmp0i = C9*Im(x[112])-C11*Im(x[56])-C13*Im(x[168])+Im(x[0]);
	tmp0r = C9*Re(x[112])-C11*Re(x[56])-C13*Re(x[168])+Re(x[0]);
	tmp1i = Im(x[112])+Im(x[56]);
	tmp1r = Re(x[112])+Re(x[56]);
	tmp2i = -C11*Im(x[112])-C13*Im(x[56])+C9*Im(x[168])+Im(x[0]);
	tmp2r = -C11*Re(x[112])-C13*Re(x[56])+C9*Re(x[168])+Re(x[0]);
	tmp3i = tmp1i+Im(x[168]);
	tmp3r = tmp1r+Re(x[168]);
	tmp4i = -C13*Im(x[112])+C9*Im(x[56])-C11*Im(x[168])+Im(x[0]);
	tmp4r = -C13*Re(x[112])+C9*Re(x[56])-C11*Re(x[168])+Re(x[0]);
	tmp5i = -C10*Re(x[112])+C12*Re(x[56])+C14*Re(x[168]);
	tmp5r = C10*Im(x[112])-C12*Im(x[56])-C14*Im(x[168]);
	tmp6i = -C12*Re(x[112])+C14*Re(x[56])-C10*Re(x[168]);
	tmp6r = C12*Im(x[112])-C14*Im(x[56])+C10*Im(x[168]);
	tmp7i = -C14*Re(x[112])-C10*Re(x[56])-C12*Re(x[168]);
	tmp7r = C14*Im(x[112])+C10*Im(x[56])+C12*Im(x[168]);
	tmp8i = tmp0i+tmp5i;
	tmp8r = tmp0r+tmp5r;
	tmp9i = tmp0i-tmp5i;
	tmp9r = tmp0r-tmp5r;
	tmp10i = tmp2i+tmp6i;
	tmp10r = tmp2r+tmp6r;
	tmp11i = tmp2i-tmp6i;
	tmp11r = tmp2r-tmp6r;
	tmp12i = tmp4i+tmp7i;
	tmp12r = tmp4r+tmp7r;
	tmp13i = tmp4i-tmp7i;
	tmp13r = tmp4r-tmp7r;
	tmp14i = Im(x[0])+tmp3i;
	tmp14r = Re(x[0])+tmp3r;
	tmp15i = Im(x[14])+Im(x[182]);
	tmp15r = Re(x[14])+Re(x[182]);
	tmp16i = Im(x[14])-Im(x[182]);
	tmp16r = Re(x[14])-Re(x[182]);
	tmp17i = Im(x[126])+Im(x[70]);
	tmp17r = Re(x[126])+Re(x[70]);
	tmp18i = Im(x[126])-Im(x[70]);
	tmp18r = Re(x[126])-Re(x[70]);
	tmp19i = C9*tmp15i-C13*tmp17i;
	tmp19r = C9*tmp15r-C13*tmp17r;
	tmp20i = -C11*tmp15i+C9*tmp17i;
	tmp20r = -C11*tmp15r+C9*tmp17r;
	tmp21i = tmp15i+tmp17i;
	tmp21r = tmp15r+tmp17r;
	tmp22i = -C13*tmp15i-C11*tmp17i;
	tmp22r = -C13*tmp15r-C11*tmp17r;
	tmp23i = -C10*tmp16r-C14*tmp18r;
	tmp23r = C10*tmp16i+C14*tmp18i;
	tmp24i = -C12*tmp16r+C10*tmp18r;
	tmp24r = C12*tmp16i-C10*tmp18i;
	tmp25i = -C14*tmp16r+C12*tmp18r;
	tmp25r = C14*tmp16i-C12*tmp18i;
	tmp26i = tmp19i+tmp23i;
	tmp26r = tmp19r+tmp23r;
	tmp27i = tmp19i-tmp23i;
	tmp27r = tmp19r-tmp23r;
	tmp28i = tmp20i+tmp24i;
	tmp28r = tmp20r+tmp24r;
	tmp29i = tmp20i-tmp24i;
	tmp29r = tmp20r-tmp24r;
	tmp30i = tmp22i+tmp25i;
	tmp30r = tmp22r+tmp25r;
	tmp31i = tmp22i-tmp25i;
	tmp31r = tmp22r-tmp25r;
	tmp32i = C9*Im(x[84])-C11*Im(x[140])-C13*Im(x[28]);
	tmp32r = C9*Re(x[84])-C11*Re(x[140])-C13*Re(x[28]);
	tmp33i = Im(x[84])+Im(x[140]);
	tmp33r = Re(x[84])+Re(x[140]);
	tmp34i = -C11*Im(x[84])-C13*Im(x[140])+C9*Im(x[28]);
	tmp34r = -C11*Re(x[84])-C13*Re(x[140])+C9*Re(x[28]);
	tmp35i = tmp33i+Im(x[28]);
	tmp35r = tmp33r+Re(x[28]);
	tmp36i = -C13*Im(x[84])+C9*Im(x[140])-C11*Im(x[28]);
	tmp36r = -C13*Re(x[84])+C9*Re(x[140])-C11*Re(x[28]);
	tmp37i = C10*Re(x[84])-C12*Re(x[140])-C14*Re(x[28]);
	tmp37r = -C10*Im(x[84])+C12*Im(x[140])+C14*Im(x[28]);
	tmp38i = C12*Re(x[84])-C14*Re(x[140])+C10*Re(x[28]);
	tmp38r = -C12*Im(x[84])+C14*Im(x[140])-C10*Im(x[28]);
	tmp39i = C14*Re(x[84])+C10*Re(x[140])+C12*Re(x[28]);
	tmp39r = -C14*Im(x[84])-C10*Im(x[140])-C12*Im(x[28]);
	tmp40i = tmp32i+tmp37i;
	tmp40r = tmp32r+tmp37r;
	tmp41i = tmp32i-tmp37i;
	tmp41r = tmp32r-tmp37r;
	tmp42i = tmp34i+tmp38i;
	tmp42r = tmp34r+tmp38r;
	tmp43i = tmp34i-tmp38i;
	tmp43r = tmp34r-tmp38r;
	tmp44i = tmp36i+tmp39i;
	tmp44r = tmp36r+tmp39r;
	tmp45i = tmp36i-tmp39i;
	tmp45r = tmp36r-tmp39r;
	tmp46i = Im(x[42])+Im(x[154]);
	tmp46r = Re(x[42])+Re(x[154]);
	tmp47i = Im(x[42])-Im(x[154]);
	tmp47r = Re(x[42])-Re(x[154]);
	tmp48i = -C11*tmp46i+Im(x[98]);
	tmp48r = -C11*tmp46r+Re(x[98]);
	tmp49i = -C13*tmp46i+Im(x[98]);
	tmp49r = -C13*tmp46r+Re(x[98]);
	tmp50i = C9*tmp46i+Im(x[98]);
	tmp50r = C9*tmp46r+Re(x[98]);
	tmp51i = -C12*tmp47r;
	tmp51r = C12*tmp47i;
	tmp52i = -C14*tmp47r;
	tmp52r = C14*tmp47i;
	tmp53i = C10*tmp47r;
	tmp53r = -C10*tmp47i;
	tmp54i = tmp48i+tmp51i;
	tmp54r = tmp48r+tmp51r;
	tmp55i = tmp48i-tmp51i;
	tmp55r = tmp48r-tmp51r;
	tmp56i = tmp49i+tmp52i;
	tmp56r = tmp49r+tmp52r;
	tmp57i = tmp49i-tmp52i;
	tmp57r = tmp49r-tmp52r;
	tmp58i = tmp50i+tmp53i;
	tmp58r = tmp50r+tmp53r;
	tmp59i = tmp50i-tmp53i;
	tmp59r = tmp50r-tmp53r;
	tmp60i = Im(x[98])+tmp46i;
	tmp60r = Re(x[98])+tmp46r;
	tmp61i = tmp14i+tmp35i;
	tmp61r = tmp14r+tmp35r;
	tmp62i = tmp14i-tmp35i;
	tmp62r = tmp14r-tmp35r;
	tmp63i = tmp21i+tmp60i;
	tmp63r = tmp21r+tmp60r;
	tmp64i = tmp21i-tmp60i;
	tmp64r = tmp21r-tmp60r;
	tmp65i = tmp61i+tmp63i;
	tmp65r = tmp61r+tmp63r;
	tmp66i = tmp61i-tmp63i;
	tmp66r = tmp61r-tmp63r;
	tmp67i = tmp62i-tmp64r;
	tmp67r = tmp62r+tmp64i;
	tmp68i = tmp62i+tmp64r;
	tmp68r = tmp62r-tmp64i;
	tmp69i = tmp12i+tmp44i;
	tmp69r = tmp12r+tmp44r;
	tmp70i = tmp12i-tmp44i;
	tmp70r = tmp12r-tmp44r;
	tmp71i = tmp30i+tmp58i;
	tmp71r = tmp30r+tmp58r;
	tmp72i = tmp30i-tmp58i;
	tmp72r = tmp30r-tmp58r;
	tmp73i = tmp69i+tmp71i;
	tmp73r = tmp69r+tmp71r;
	tmp74i = tmp69i-tmp71i;
	tmp74r = tmp69r-tmp71r;
	tmp75i = tmp70i-tmp72r;
	tmp75r = tmp70r+tmp72i;
	tmp76i = tmp70i+tmp72r;
	tmp76r = tmp70r-tmp72i;
	tmp77i = tmp11i+tmp43i;
	tmp77r = tmp11r+tmp43r;
	tmp78i = tmp11i-tmp43i;
	tmp78r = tmp11r-tmp43r;
	tmp79i = tmp29i+tmp57i;
	tmp79r = tmp29r+tmp57r;
	tmp80i = tmp29i-tmp57i;
	tmp80r = tmp29r-tmp57r;
	tmp81i = tmp77i+tmp79i;
	tmp81r = tmp77r+tmp79r;
	tmp82i = tmp77i-tmp79i;
	tmp82r = tmp77r-tmp79r;
	tmp83i = tmp78i-tmp80r;
	tmp83r = tmp78r+tmp80i;
	tmp84i = tmp78i+tmp80r;
	tmp84r = tmp78r-tmp80i;
	tmp85i = tmp9i+tmp41i;
	tmp85r = tmp9r+tmp41r;
	tmp86i = tmp9i-tmp41i;
	tmp86r = tmp9r-tmp41r;
	tmp87i = tmp27i+tmp55i;
	tmp87r = tmp27r+tmp55r;
	tmp88i = tmp27i-tmp55i;
	tmp88r = tmp27r-tmp55r;
	tmp89i = tmp85i+tmp87i;
	tmp89r = tmp85r+tmp87r;
	tmp90i = tmp85i-tmp87i;
	tmp90r = tmp85r-tmp87r;
	tmp91i = tmp86i-tmp88r;
	tmp91r = tmp86r+tmp88i;
	tmp92i = tmp86i+tmp88r;
	tmp92r = tmp86r-tmp88i;
	tmp93i = tmp8i+tmp40i;
	tmp93r = tmp8r+tmp40r;
	tmp94i = tmp8i-tmp40i;
	tmp94r = tmp8r-tmp40r;
	tmp95i = tmp26i+tmp54i;
	tmp95r = tmp26r+tmp54r;
	tmp96i = tmp26i-tmp54i;
	tmp96r = tmp26r-tmp54r;
	tmp97i = tmp93i+tmp95i;
	tmp97r = tmp93r+tmp95r;
	tmp98i = tmp93i-tmp95i;
	tmp98r = tmp93r-tmp95r;
	tmp99i = tmp94i-tmp96r;
	tmp99r = tmp94r+tmp96i;
	tmp100i = tmp94i+tmp96r;
	tmp100r = tmp94r-tmp96i;
	tmp101i = tmp10i+tmp42i;
	tmp101r = tmp10r+tmp42r;
	tmp102i = tmp10i-tmp42i;
	tmp102r = tmp10r-tmp42r;
	tmp103i = tmp28i+tmp56i;
	tmp103r = tmp28r+tmp56r;
	tmp104i = tmp28i-tmp56i;
	tmp104r = tmp28r-tmp56r;
	tmp105i = tmp101i+tmp103i;
	tmp105r = tmp101r+tmp103r;
	tmp106i = tmp101i-tmp103i;
	tmp106r = tmp101r-tmp103r;
	tmp107i = tmp102i-tmp104r;
	tmp107r = tmp102r+tmp104i;
	tmp108i = tmp102i+tmp104r;
	tmp108r = tmp102r-tmp104i;
	tmp109i = tmp13i+tmp45i;
	tmp109r = tmp13r+tmp45r;
	tmp110i = tmp13i-tmp45i;
	tmp110r = tmp13r-tmp45r;
	tmp111i = tmp31i+tmp59i;
	tmp111r = tmp31r+tmp59r;
	tmp112i = tmp31i-tmp59i;
	tmp112r = tmp31r-tmp59r;
	tmp113i = tmp109i+tmp111i;
	tmp113r = tmp109r+tmp111r;
	tmp114i = tmp109i-tmp111i;
	tmp114r = tmp109r-tmp111r;
	tmp115i = tmp110i-tmp112r;
	tmp115r = tmp110r+tmp112i;
	tmp116i = tmp110i+tmp112r;
	tmp116r = tmp110r-tmp112i;
	Re(x[0]) = tmp65r;
	Im(x[0]) = tmp65i;
	Re(x[14]) = tmp76r;
	Im(x[14]) = tmp76i;
	Re(x[28]) = tmp82r;
	Im(x[28]) = tmp82i;
	Re(x[42]) = tmp91r;
	Im(x[42]) = tmp91i;
	Re(x[56]) = tmp97r;
	Im(x[56]) = tmp97i;
	Re(x[70]) = tmp108r;
	Im(x[70]) = tmp108i;
	Re(x[84]) = tmp114r;
	Im(x[84]) = tmp114i;
	Re(x[98]) = tmp67r;
	Im(x[98]) = tmp67i;
	Re(x[112]) = tmp73r;
	Im(x[112]) = tmp73i;
	Re(x[126]) = tmp84r;
	Im(x[126]) = tmp84i;
	Re(x[140]) = tmp90r;
	Im(x[140]) = tmp90i;
	Re(x[154]) = tmp99r;
	Im(x[154]) = tmp99i;
	Re(x[168]) = tmp105r;
	Im(x[168]) = tmp105i;
	Re(x[182]) = tmp116r;
	Im(x[182]) = tmp116i;
	Re(x[196]) = tmp66r;
	Im(x[196]) = tmp66i;
	Re(x[210]) = tmp75r;
	Im(x[210]) = tmp75i;
	Re(x[224]) = tmp81r;
	Im(x[224]) = tmp81i;
	Re(x[238]) = tmp92r;
	Im(x[238]) = tmp92i;
	Re(x[252]) = tmp98r;
	Im(x[252]) = tmp98i;
	Re(x[266]) = tmp107r;
	Im(x[266]) = tmp107i;
	Re(x[280]) = tmp113r;
	Im(x[280]) = tmp113i;
	Re(x[294]) = tmp68r;
	Im(x[294]) = tmp68i;
	Re(x[308]) = tmp74r;
	Im(x[308]) = tmp74i;
	Re(x[322]) = tmp83r;
	Im(x[322]) = tmp83i;
	Re(x[336]) = tmp89r;
	Im(x[336]) = tmp89i;
	Re(x[350]) = tmp100r;
	Im(x[350]) = tmp100i;
	Re(x[364]) = tmp106r;
	Im(x[364]) = tmp106i;
	Re(x[378]) = tmp115r;
	Im(x[378]) = tmp115i;
}

/*
*	Number of additions = 352
*	Number of multiplications = 144
*	Number of sign changes = 0
*	Number of assigns = 336
*	Total number of operations = 832
*/
void	MIFFTC28(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp134r, tmp135r, tmp136r, tmp137r, tmp138r, tmp139r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp134i, tmp135i, tmp136i, tmp137i, tmp138i, tmp139i;

	const REAL	C13 =     0.22252093395631;	/* REALCONST	*/
	const REAL	C12 =     0.43388373911756;	/* REALCONST	*/
	const REAL	C9 =     0.62348980185873;	/* REALCONST	*/
	const REAL	C10 =     0.78183148246803;	/* REALCONST	*/
	const REAL	C11 =     0.90096886790242;	/* REALCONST	*/
	const REAL	C14 =     0.97492791218182;	/* REALCONST	*/

	tmp0i = Im(x[112])+Im(x[280]);
	tmp0r = Re(x[112])+Re(x[280]);
	tmp1i = Im(x[112])-Im(x[280]);
	tmp1r = Re(x[112])-Re(x[280]);
	tmp2i = Im(x[336])+Im(x[56]);
	tmp2r = Re(x[336])+Re(x[56]);
	tmp3i = Im(x[336])-Im(x[56]);
	tmp3r = Re(x[336])-Re(x[56]);
	tmp4i = Im(x[224])+Im(x[168]);
	tmp4r = Re(x[224])+Re(x[168]);
	tmp5i = Im(x[224])-Im(x[168]);
	tmp5r = Re(x[224])-Re(x[168]);
	tmp6i = C9*tmp0i-C11*tmp2i-C13*tmp4i+Im(x[0]);
	tmp6r = C9*tmp0r-C11*tmp2r-C13*tmp4r+Re(x[0]);
	tmp7i = tmp0i+tmp2i;
	tmp7r = tmp0r+tmp2r;
	tmp8i = -C11*tmp0i-C13*tmp2i+C9*tmp4i+Im(x[0]);
	tmp8r = -C11*tmp0r-C13*tmp2r+C9*tmp4r+Re(x[0]);
	tmp9i = tmp7i+tmp4i;
	tmp9r = tmp7r+tmp4r;
	tmp10i = -C13*tmp0i+C9*tmp2i-C11*tmp4i+Im(x[0]);
	tmp10r = -C13*tmp0r+C9*tmp2r-C11*tmp4r+Re(x[0]);
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
	tmp20i = Im(x[0])+tmp9i;
	tmp20r = Re(x[0])+tmp9r;
	tmp21i = Im(x[14])+Im(x[182]);
	tmp21r = Re(x[14])+Re(x[182]);
	tmp22i = Im(x[14])-Im(x[182]);
	tmp22r = Re(x[14])-Re(x[182]);
	tmp23i = Im(x[238])+Im(x[350]);
	tmp23r = Re(x[238])+Re(x[350]);
	tmp24i = Im(x[238])-Im(x[350]);
	tmp24r = Re(x[238])-Re(x[350]);
	tmp25i = Im(x[126])+Im(x[70]);
	tmp25r = Re(x[126])+Re(x[70]);
	tmp26i = Im(x[126])-Im(x[70]);
	tmp26r = Re(x[126])-Re(x[70]);
	tmp27i = C9*tmp21i-C11*tmp23i-C13*tmp25i+Im(x[294]);
	tmp27r = C9*tmp21r-C11*tmp23r-C13*tmp25r+Re(x[294]);
	tmp28i = tmp21i+tmp23i;
	tmp28r = tmp21r+tmp23r;
	tmp29i = -C11*tmp21i-C13*tmp23i+C9*tmp25i+Im(x[294]);
	tmp29r = -C11*tmp21r-C13*tmp23r+C9*tmp25r+Re(x[294]);
	tmp30i = tmp28i+tmp25i;
	tmp30r = tmp28r+tmp25r;
	tmp31i = -C13*tmp21i+C9*tmp23i-C11*tmp25i+Im(x[294]);
	tmp31r = -C13*tmp21r+C9*tmp23r-C11*tmp25r+Re(x[294]);
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
	tmp41i = Im(x[294])+tmp30i;
	tmp41r = Re(x[294])+tmp30r;
	tmp42i = Im(x[308])+Im(x[84]);
	tmp42r = Re(x[308])+Re(x[84]);
	tmp43i = Im(x[308])-Im(x[84]);
	tmp43r = Re(x[308])-Re(x[84]);
	tmp44i = Im(x[140])+Im(x[252]);
	tmp44r = Re(x[140])+Re(x[252]);
	tmp45i = Im(x[140])-Im(x[252]);
	tmp45r = Re(x[140])-Re(x[252]);
	tmp46i = Im(x[28])+Im(x[364]);
	tmp46r = Re(x[28])+Re(x[364]);
	tmp47i = Im(x[28])-Im(x[364]);
	tmp47r = Re(x[28])-Re(x[364]);
	tmp48i = C9*tmp42i-C11*tmp44i-C13*tmp46i+Im(x[196]);
	tmp48r = C9*tmp42r-C11*tmp44r-C13*tmp46r+Re(x[196]);
	tmp49i = tmp42i+tmp44i;
	tmp49r = tmp42r+tmp44r;
	tmp50i = -C11*tmp42i-C13*tmp44i+C9*tmp46i+Im(x[196]);
	tmp50r = -C11*tmp42r-C13*tmp44r+C9*tmp46r+Re(x[196]);
	tmp51i = tmp49i+tmp46i;
	tmp51r = tmp49r+tmp46r;
	tmp52i = -C13*tmp42i+C9*tmp44i-C11*tmp46i+Im(x[196]);
	tmp52r = -C13*tmp42r+C9*tmp44r-C11*tmp46r+Re(x[196]);
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
	tmp62i = Im(x[196])+tmp51i;
	tmp62r = Re(x[196])+tmp51r;
	tmp63i = Im(x[210])+Im(x[378]);
	tmp63r = Re(x[210])+Re(x[378]);
	tmp64i = Im(x[210])-Im(x[378]);
	tmp64r = Re(x[210])-Re(x[378]);
	tmp65i = Im(x[42])+Im(x[154]);
	tmp65r = Re(x[42])+Re(x[154]);
	tmp66i = Im(x[42])-Im(x[154]);
	tmp66r = Re(x[42])-Re(x[154]);
	tmp67i = Im(x[322])+Im(x[266]);
	tmp67r = Re(x[322])+Re(x[266]);
	tmp68i = Im(x[322])-Im(x[266]);
	tmp68r = Re(x[322])-Re(x[266]);
	tmp69i = C9*tmp63i-C11*tmp65i-C13*tmp67i+Im(x[98]);
	tmp69r = C9*tmp63r-C11*tmp65r-C13*tmp67r+Re(x[98]);
	tmp70i = tmp63i+tmp65i;
	tmp70r = tmp63r+tmp65r;
	tmp71i = -C11*tmp63i-C13*tmp65i+C9*tmp67i+Im(x[98]);
	tmp71r = -C11*tmp63r-C13*tmp65r+C9*tmp67r+Re(x[98]);
	tmp72i = tmp70i+tmp67i;
	tmp72r = tmp70r+tmp67r;
	tmp73i = -C13*tmp63i+C9*tmp65i-C11*tmp67i+Im(x[98]);
	tmp73r = -C13*tmp63r+C9*tmp65r-C11*tmp67r+Re(x[98]);
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
	tmp83i = Im(x[98])+tmp72i;
	tmp83r = Re(x[98])+tmp72r;
	tmp84i = tmp20i+tmp62i;
	tmp84r = tmp20r+tmp62r;
	tmp85i = tmp20i-tmp62i;
	tmp85r = tmp20r-tmp62r;
	tmp86i = tmp41i+tmp83i;
	tmp86r = tmp41r+tmp83r;
	tmp87i = tmp41i-tmp83i;
	tmp87r = tmp41r-tmp83r;
	tmp88i = tmp84i+tmp86i;
	tmp88r = tmp84r+tmp86r;
	tmp89i = tmp84i-tmp86i;
	tmp89r = tmp84r-tmp86r;
	tmp90i = tmp85i+tmp87r;
	tmp90r = tmp85r-tmp87i;
	tmp91i = tmp85i-tmp87r;
	tmp91r = tmp85r+tmp87i;
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
	tmp97i = tmp92i-tmp94i;
	tmp97r = tmp92r-tmp94r;
	tmp98i = tmp93i+tmp95r;
	tmp98r = tmp93r-tmp95i;
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
	tmp104i = tmp100i+tmp102i;
	tmp104r = tmp100r+tmp102r;
	tmp105i = tmp100i-tmp102i;
	tmp105r = tmp100r-tmp102r;
	tmp106i = tmp101i+tmp103r;
	tmp106r = tmp101r-tmp103i;
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
	tmp112i = tmp108i+tmp110i;
	tmp112r = tmp108r+tmp110r;
	tmp113i = tmp108i-tmp110i;
	tmp113r = tmp108r-tmp110r;
	tmp114i = tmp109i+tmp111r;
	tmp114r = tmp109r-tmp111i;
	tmp115i = tmp109i-tmp111r;
	tmp115r = tmp109r+tmp111i;
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
	tmp121i = tmp116i-tmp118i;
	tmp121r = tmp116r-tmp118r;
	tmp122i = tmp117i+tmp119r;
	tmp122r = tmp117r-tmp119i;
	tmp123i = tmp117i-tmp119r;
	tmp123r = tmp117r+tmp119i;
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
	tmp129i = tmp124i-tmp126i;
	tmp129r = tmp124r-tmp126r;
	tmp130i = tmp125i+tmp127r;
	tmp130r = tmp125r-tmp127i;
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
	tmp136i = tmp132i+tmp134i;
	tmp136r = tmp132r+tmp134r;
	tmp137i = tmp132i-tmp134i;
	tmp137r = tmp132r-tmp134r;
	tmp138i = tmp133i+tmp135r;
	tmp138r = tmp133r-tmp135i;
	tmp139i = tmp133i-tmp135r;
	tmp139r = tmp133r+tmp135i;
	Re(x[0]) = tmp88r;
	Im(x[0]) = tmp88i;
	Re(x[14]) = tmp99r;
	Im(x[14]) = tmp99i;
	Re(x[28]) = tmp105r;
	Im(x[28]) = tmp105i;
	Re(x[42]) = tmp114r;
	Im(x[42]) = tmp114i;
	Re(x[56]) = tmp120r;
	Im(x[56]) = tmp120i;
	Re(x[70]) = tmp131r;
	Im(x[70]) = tmp131i;
	Re(x[84]) = tmp137r;
	Im(x[84]) = tmp137i;
	Re(x[98]) = tmp90r;
	Im(x[98]) = tmp90i;
	Re(x[112]) = tmp96r;
	Im(x[112]) = tmp96i;
	Re(x[126]) = tmp107r;
	Im(x[126]) = tmp107i;
	Re(x[140]) = tmp113r;
	Im(x[140]) = tmp113i;
	Re(x[154]) = tmp122r;
	Im(x[154]) = tmp122i;
	Re(x[168]) = tmp128r;
	Im(x[168]) = tmp128i;
	Re(x[182]) = tmp139r;
	Im(x[182]) = tmp139i;
	Re(x[196]) = tmp89r;
	Im(x[196]) = tmp89i;
	Re(x[210]) = tmp98r;
	Im(x[210]) = tmp98i;
	Re(x[224]) = tmp104r;
	Im(x[224]) = tmp104i;
	Re(x[238]) = tmp115r;
	Im(x[238]) = tmp115i;
	Re(x[252]) = tmp121r;
	Im(x[252]) = tmp121i;
	Re(x[266]) = tmp130r;
	Im(x[266]) = tmp130i;
	Re(x[280]) = tmp136r;
	Im(x[280]) = tmp136i;
	Re(x[294]) = tmp91r;
	Im(x[294]) = tmp91i;
	Re(x[308]) = tmp97r;
	Im(x[308]) = tmp97i;
	Re(x[322]) = tmp106r;
	Im(x[322]) = tmp106i;
	Re(x[336]) = tmp112r;
	Im(x[336]) = tmp112i;
	Re(x[350]) = tmp123r;
	Im(x[350]) = tmp123i;
	Re(x[364]) = tmp129r;
	Im(x[364]) = tmp129i;
	Re(x[378]) = tmp138r;
	Im(x[378]) = tmp138i;
}

/*
*	Number of additions = 304
*	Number of multiplications = 104
*	Number of sign changes = 46
*	Number of assigns = 436
*	Total number of operations = 890
*/
void	MFFTC30(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp169r, tmp170r, tmp171r, tmp172r, tmp173r, tmp174r, tmp175r,
		tmp176r, tmp177r, tmp178r, tmp179r, tmp180r, tmp181r, tmp182r,
		tmp183r, tmp184r, tmp185r, tmp186r, tmp187r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp155i, tmp156i, tmp157i, tmp158i, tmp159i, tmp160i, tmp161i,
		tmp162i, tmp163i, tmp164i, tmp165i, tmp166i, tmp167i, tmp168i,
		tmp169i, tmp170i, tmp171i, tmp172i, tmp173i, tmp174i, tmp175i,
		tmp176i, tmp177i, tmp178i, tmp179i, tmp180i, tmp181i, tmp182i,
		tmp183i, tmp184i, tmp185i, tmp186i, tmp187i;

	const REAL	C1 =                 0.25;	/* REALCONST	*/
	const REAL	C12 =                  0.5;	/* REALCONST	*/
	const REAL	C14 =     0.55901699437495;	/* REALCONST	*/
	const REAL	C10 =     0.58778525229247;	/* REALCONST	*/
	const REAL	C16 =     0.86602540378444;	/* REALCONST	*/
	const REAL	C8 =     0.95105651629515;	/* REALCONST	*/

	tmp0i = Im(x[90])+Im(x[180]);
	tmp0r = Re(x[90])+Re(x[180]);
	tmp1i = Im(x[90])-Im(x[180]);
	tmp1r = Re(x[90])-Re(x[180]);
	tmp2i = -C1*tmp0i;
	tmp2r = -C1*tmp0r;
	tmp3i = C14*tmp1i;
	tmp3r = C14*tmp1r;
	tmp4i = tmp2i+Im(x[0]);
	tmp4r = tmp2r+Re(x[0]);
	tmp5i = tmp4i+tmp3i;
	tmp5r = tmp4r+tmp3r;
	tmp6i = tmp4i-tmp3i;
	tmp6r = tmp4r-tmp3r;
	tmp7i = -C8*Re(x[90])-C10*Re(x[180]);
	tmp7r = C8*Im(x[90])+C10*Im(x[180]);
	tmp8i = -C10*Re(x[90])+C8*Re(x[180]);
	tmp8r = C10*Im(x[90])-C8*Im(x[180]);
	tmp9i = tmp5i+tmp7i;
	tmp9r = tmp5r+tmp7r;
	tmp10i = tmp5i-tmp7i;
	tmp10r = tmp5r-tmp7r;
	tmp11i = tmp6i+tmp8i;
	tmp11r = tmp6r+tmp8r;
	tmp12i = tmp6i-tmp8i;
	tmp12r = tmp6r-tmp8r;
	tmp13i = Im(x[0])+tmp0i;
	tmp13r = Re(x[0])+tmp0r;
	tmp14i = Im(x[105])+Im(x[195]);
	tmp14r = Re(x[105])+Re(x[195]);
	tmp15i = Im(x[105])-Im(x[195]);
	tmp15r = Re(x[105])-Re(x[195]);
	tmp16i = Im(x[15])+tmp14i;
	tmp16r = Re(x[15])+tmp14r;
	tmp17i = Im(x[15])-tmp14i;
	tmp17r = Re(x[15])-tmp14r;
	tmp18i = -C1*tmp16i;
	tmp18r = -C1*tmp16r;
	tmp19i = C14*tmp17i;
	tmp19r = C14*tmp17r;
	tmp20i = tmp18i+tmp19i;
	tmp20r = tmp18r+tmp19r;
	tmp21i = tmp18i-tmp19i;
	tmp21r = tmp18r-tmp19r;
	tmp22i = -C8*Re(x[15])-C10*tmp15r;
	tmp22r = C8*Im(x[15])+C10*tmp15i;
	tmp23i = -C10*Re(x[15])+C8*tmp15r;
	tmp23r = C10*Im(x[15])-C8*tmp15i;
	tmp24i = tmp20i+tmp22i;
	tmp24r = tmp20r+tmp22r;
	tmp25i = tmp20i-tmp22i;
	tmp25r = tmp20r-tmp22r;
	tmp26i = tmp21i+tmp23i;
	tmp26r = tmp21r+tmp23r;
	tmp27i = tmp21i-tmp23i;
	tmp27r = tmp21r-tmp23r;
	tmp28i = Im(x[30])+Im(x[120]);
	tmp28r = Re(x[30])+Re(x[120]);
	tmp29i = Im(x[30])-Im(x[120]);
	tmp29r = Re(x[30])-Re(x[120]);
	tmp30i = Im(x[210])+tmp28i;
	tmp30r = Re(x[210])+tmp28r;
	tmp31i = Im(x[210])-tmp28i;
	tmp31r = Re(x[210])-tmp28r;
	tmp32i = -C1*tmp30i;
	tmp32r = -C1*tmp30r;
	tmp33i = C14*tmp31i;
	tmp33r = C14*tmp31r;
	tmp34i = tmp32i+tmp33i;
	tmp34r = tmp32r+tmp33r;
	tmp35i = tmp32i-tmp33i;
	tmp35r = tmp32r-tmp33r;
	tmp36i = C8*Re(x[210])-C10*tmp29r;
	tmp36r = -C8*Im(x[210])+C10*tmp29i;
	tmp37i = C10*Re(x[210])+C8*tmp29r;
	tmp37r = -C10*Im(x[210])-C8*tmp29i;
	tmp38i = tmp34i+tmp36i;
	tmp38r = tmp34r+tmp36r;
	tmp39i = tmp34i-tmp36i;
	tmp39r = tmp34r-tmp36r;
	tmp40i = tmp35i+tmp37i;
	tmp40r = tmp35r+tmp37r;
	tmp41i = tmp35i-tmp37i;
	tmp41r = tmp35r-tmp37r;
	tmp42i = Im(x[135])+Im(x[45]);
	tmp42r = Re(x[135])+Re(x[45]);
	tmp43i = Im(x[135])-Im(x[45]);
	tmp43r = Re(x[135])-Re(x[45]);
	tmp44i = -C1*tmp42i;
	tmp44r = -C1*tmp42r;
	tmp45i = C14*tmp43i;
	tmp45r = C14*tmp43r;
	tmp46i = tmp44i+tmp45i;
	tmp46r = tmp44r+tmp45r;
	tmp47i = tmp44i-tmp45i;
	tmp47r = tmp44r-tmp45r;
	tmp48i = C8*Re(x[135])+C10*Re(x[45]);
	tmp48r = -C8*Im(x[135])-C10*Im(x[45]);
	tmp49i = C10*Re(x[135])-C8*Re(x[45]);
	tmp49r = -C10*Im(x[135])+C8*Im(x[45]);
	tmp50i = tmp46i+tmp48i;
	tmp50r = tmp46r+tmp48r;
	tmp51i = tmp46i-tmp48i;
	tmp51r = tmp46r-tmp48r;
	tmp52i = tmp47i+tmp49i;
	tmp52r = tmp47r+tmp49r;
	tmp53i = tmp47i-tmp49i;
	tmp53r = tmp47r-tmp49r;
	tmp54i = -C1*Im(x[60]);
	tmp54r = -C1*Re(x[60]);
	tmp55i = C14*Im(x[60]);
	tmp55r = C14*Re(x[60]);
	tmp56i = tmp54i+Im(x[150]);
	tmp56r = tmp54r+Re(x[150]);
	tmp57i = tmp56i+tmp55i;
	tmp57r = tmp56r+tmp55r;
	tmp58i = tmp56i-tmp55i;
	tmp58r = tmp56r-tmp55r;
	tmp59i = C8*Re(x[60]);
	tmp59r = -C8*Im(x[60]);
	tmp60i = C10*Re(x[60]);
	tmp60r = -C10*Im(x[60]);
	tmp61i = tmp57i+tmp59i;
	tmp61r = tmp57r+tmp59r;
	tmp62i = tmp57i-tmp59i;
	tmp62r = tmp57r-tmp59r;
	tmp63i = tmp58i+tmp60i;
	tmp63r = tmp58r+tmp60r;
	tmp64i = tmp58i-tmp60i;
	tmp64r = tmp58r-tmp60r;
	tmp65i = Im(x[150])+Im(x[60]);
	tmp65r = Re(x[150])+Re(x[60]);
	tmp66i = -C1*Im(x[165]);
	tmp66r = -C1*Re(x[165]);
	tmp67i = C14*Im(x[165]);
	tmp67r = C14*Re(x[165]);
	tmp68i = tmp66i+Im(x[75]);
	tmp68r = tmp66r+Re(x[75]);
	tmp69i = tmp68i+tmp67i;
	tmp69r = tmp68r+tmp67r;
	tmp70i = tmp68i-tmp67i;
	tmp70r = tmp68r-tmp67r;
	tmp71i = -C8*Re(x[165]);
	tmp71r = C8*Im(x[165]);
	tmp72i = -C10*Re(x[165]);
	tmp72r = C10*Im(x[165]);
	tmp73i = tmp69i+tmp71i;
	tmp73r = tmp69r+tmp71r;
	tmp74i = tmp69i-tmp71i;
	tmp74r = tmp69r-tmp71r;
	tmp75i = tmp70i+tmp72i;
	tmp75r = tmp70r+tmp72r;
	tmp76i = tmp70i-tmp72i;
	tmp76r = tmp70r-tmp72r;
	tmp77i = Im(x[75])+Im(x[165]);
	tmp77r = Re(x[75])+Re(x[165]);
	tmp78i = tmp65i+tmp30i;
	tmp78r = tmp65r+tmp30r;
	tmp79i = tmp65i-tmp30i;
	tmp79r = tmp65r-tmp30r;
	tmp80i = -C12*tmp78i;
	tmp80r = -C12*tmp78r;
	tmp81i = -C16*tmp79r;
	tmp81r = C16*tmp79i;
	tmp82i = tmp80i+tmp13i;
	tmp82r = tmp80r+tmp13r;
	tmp83i = tmp82i+tmp81i;
	tmp83r = tmp82r+tmp81r;
	tmp84i = tmp82i-tmp81i;
	tmp84r = tmp82r-tmp81r;
	tmp85i = tmp13i+tmp78i;
	tmp85r = tmp13r+tmp78r;
	tmp86i = tmp16i+tmp77i;
	tmp86r = tmp16r+tmp77r;
	tmp87i = tmp16i-tmp77i;
	tmp87r = tmp16r-tmp77r;
	tmp88i = -C12*tmp86i;
	tmp88r = -C12*tmp86r;
	tmp89i = -C16*tmp87r;
	tmp89r = C16*tmp87i;
	tmp90i = tmp88i+tmp42i;
	tmp90r = tmp88r+tmp42r;
	tmp91i = tmp90i+tmp89i;
	tmp91r = tmp90r+tmp89r;
	tmp92i = tmp90i-tmp89i;
	tmp92r = tmp90r-tmp89r;
	tmp93i = tmp42i+tmp86i;
	tmp93r = tmp42r+tmp86r;
	tmp94i = tmp85i+tmp93i;
	tmp94r = tmp85r+tmp93r;
	tmp95i = tmp85i-tmp93i;
	tmp95r = tmp85r-tmp93r;
	tmp96i = tmp84i+tmp92i;
	tmp96r = tmp84r+tmp92r;
	tmp97i = tmp84i-tmp92i;
	tmp97r = tmp84r-tmp92r;
	tmp98i = tmp83i+tmp91i;
	tmp98r = tmp83r+tmp91r;
	tmp99i = tmp83i-tmp91i;
	tmp99r = tmp83r-tmp91r;
	tmp100i = tmp61i+tmp38i;
	tmp100r = tmp61r+tmp38r;
	tmp101i = tmp61i-tmp38i;
	tmp101r = tmp61r-tmp38r;
	tmp102i = -C12*tmp100i;
	tmp102r = -C12*tmp100r;
	tmp103i = -C16*tmp101r;
	tmp103r = C16*tmp101i;
	tmp104i = tmp102i+tmp9i;
	tmp104r = tmp102r+tmp9r;
	tmp105i = tmp104i+tmp103i;
	tmp105r = tmp104r+tmp103r;
	tmp106i = tmp104i-tmp103i;
	tmp106r = tmp104r-tmp103r;
	tmp107i = tmp9i+tmp100i;
	tmp107r = tmp9r+tmp100r;
	tmp108i = tmp24i+tmp73i;
	tmp108r = tmp24r+tmp73r;
	tmp109i = tmp24i-tmp73i;
	tmp109r = tmp24r-tmp73r;
	tmp110i = -C12*tmp108i;
	tmp110r = -C12*tmp108r;
	tmp111i = -C16*tmp109r;
	tmp111r = C16*tmp109i;
	tmp112i = tmp110i+tmp50i;
	tmp112r = tmp110r+tmp50r;
	tmp113i = tmp112i+tmp111i;
	tmp113r = tmp112r+tmp111r;
	tmp114i = tmp112i-tmp111i;
	tmp114r = tmp112r-tmp111r;
	tmp115i = tmp50i+tmp108i;
	tmp115r = tmp50r+tmp108r;
	tmp116i = tmp107i+tmp115i;
	tmp116r = tmp107r+tmp115r;
	tmp117i = tmp107i-tmp115i;
	tmp117r = tmp107r-tmp115r;
	tmp118i = tmp106i+tmp114i;
	tmp118r = tmp106r+tmp114r;
	tmp119i = tmp106i-tmp114i;
	tmp119r = tmp106r-tmp114r;
	tmp120i = tmp105i+tmp113i;
	tmp120r = tmp105r+tmp113r;
	tmp121i = tmp105i-tmp113i;
	tmp121r = tmp105r-tmp113r;
	tmp122i = tmp63i+tmp40i;
	tmp122r = tmp63r+tmp40r;
	tmp123i = tmp63i-tmp40i;
	tmp123r = tmp63r-tmp40r;
	tmp124i = -C12*tmp122i;
	tmp124r = -C12*tmp122r;
	tmp125i = -C16*tmp123r;
	tmp125r = C16*tmp123i;
	tmp126i = tmp124i+tmp11i;
	tmp126r = tmp124r+tmp11r;
	tmp127i = tmp126i+tmp125i;
	tmp127r = tmp126r+tmp125r;
	tmp128i = tmp126i-tmp125i;
	tmp128r = tmp126r-tmp125r;
	tmp129i = tmp11i+tmp122i;
	tmp129r = tmp11r+tmp122r;
	tmp130i = tmp26i+tmp75i;
	tmp130r = tmp26r+tmp75r;
	tmp131i = tmp26i-tmp75i;
	tmp131r = tmp26r-tmp75r;
	tmp132i = -C12*tmp130i;
	tmp132r = -C12*tmp130r;
	tmp133i = -C16*tmp131r;
	tmp133r = C16*tmp131i;
	tmp134i = tmp132i+tmp52i;
	tmp134r = tmp132r+tmp52r;
	tmp135i = tmp134i+tmp133i;
	tmp135r = tmp134r+tmp133r;
	tmp136i = tmp134i-tmp133i;
	tmp136r = tmp134r-tmp133r;
	tmp137i = tmp52i+tmp130i;
	tmp137r = tmp52r+tmp130r;
	tmp138i = tmp129i+tmp137i;
	tmp138r = tmp129r+tmp137r;
	tmp139i = tmp129i-tmp137i;
	tmp139r = tmp129r-tmp137r;
	tmp140i = tmp128i+tmp136i;
	tmp140r = tmp128r+tmp136r;
	tmp141i = tmp128i-tmp136i;
	tmp141r = tmp128r-tmp136r;
	tmp142i = tmp127i+tmp135i;
	tmp142r = tmp127r+tmp135r;
	tmp143i = tmp127i-tmp135i;
	tmp143r = tmp127r-tmp135r;
	tmp144i = tmp64i+tmp41i;
	tmp144r = tmp64r+tmp41r;
	tmp145i = tmp64i-tmp41i;
	tmp145r = tmp64r-tmp41r;
	tmp146i = -C12*tmp144i;
	tmp146r = -C12*tmp144r;
	tmp147i = -C16*tmp145r;
	tmp147r = C16*tmp145i;
	tmp148i = tmp146i+tmp12i;
	tmp148r = tmp146r+tmp12r;
	tmp149i = tmp148i+tmp147i;
	tmp149r = tmp148r+tmp147r;
	tmp150i = tmp148i-tmp147i;
	tmp150r = tmp148r-tmp147r;
	tmp151i = tmp12i+tmp144i;
	tmp151r = tmp12r+tmp144r;
	tmp152i = tmp27i+tmp76i;
	tmp152r = tmp27r+tmp76r;
	tmp153i = tmp27i-tmp76i;
	tmp153r = tmp27r-tmp76r;
	tmp154i = -C12*tmp152i;
	tmp154r = -C12*tmp152r;
	tmp155i = -C16*tmp153r;
	tmp155r = C16*tmp153i;
	tmp156i = tmp154i+tmp53i;
	tmp156r = tmp154r+tmp53r;
	tmp157i = tmp156i+tmp155i;
	tmp157r = tmp156r+tmp155r;
	tmp158i = tmp156i-tmp155i;
	tmp158r = tmp156r-tmp155r;
	tmp159i = tmp53i+tmp152i;
	tmp159r = tmp53r+tmp152r;
	tmp160i = tmp151i+tmp159i;
	tmp160r = tmp151r+tmp159r;
	tmp161i = tmp151i-tmp159i;
	tmp161r = tmp151r-tmp159r;
	tmp162i = tmp150i+tmp158i;
	tmp162r = tmp150r+tmp158r;
	tmp163i = tmp150i-tmp158i;
	tmp163r = tmp150r-tmp158r;
	tmp164i = tmp149i+tmp157i;
	tmp164r = tmp149r+tmp157r;
	tmp165i = tmp149i-tmp157i;
	tmp165r = tmp149r-tmp157r;
	tmp166i = tmp62i+tmp39i;
	tmp166r = tmp62r+tmp39r;
	tmp167i = tmp62i-tmp39i;
	tmp167r = tmp62r-tmp39r;
	tmp168i = -C12*tmp166i;
	tmp168r = -C12*tmp166r;
	tmp169i = -C16*tmp167r;
	tmp169r = C16*tmp167i;
	tmp170i = tmp168i+tmp10i;
	tmp170r = tmp168r+tmp10r;
	tmp171i = tmp170i+tmp169i;
	tmp171r = tmp170r+tmp169r;
	tmp172i = tmp170i-tmp169i;
	tmp172r = tmp170r-tmp169r;
	tmp173i = tmp10i+tmp166i;
	tmp173r = tmp10r+tmp166r;
	tmp174i = tmp25i+tmp74i;
	tmp174r = tmp25r+tmp74r;
	tmp175i = tmp25i-tmp74i;
	tmp175r = tmp25r-tmp74r;
	tmp176i = -C12*tmp174i;
	tmp176r = -C12*tmp174r;
	tmp177i = -C16*tmp175r;
	tmp177r = C16*tmp175i;
	tmp178i = tmp176i+tmp51i;
	tmp178r = tmp176r+tmp51r;
	tmp179i = tmp178i+tmp177i;
	tmp179r = tmp178r+tmp177r;
	tmp180i = tmp178i-tmp177i;
	tmp180r = tmp178r-tmp177r;
	tmp181i = tmp51i+tmp174i;
	tmp181r = tmp51r+tmp174r;
	tmp182i = tmp173i+tmp181i;
	tmp182r = tmp173r+tmp181r;
	tmp183i = tmp173i-tmp181i;
	tmp183r = tmp173r-tmp181r;
	tmp184i = tmp172i+tmp180i;
	tmp184r = tmp172r+tmp180r;
	tmp185i = tmp172i-tmp180i;
	tmp185r = tmp172r-tmp180r;
	tmp186i = tmp171i+tmp179i;
	tmp186r = tmp171r+tmp179r;
	tmp187i = tmp171i-tmp179i;
	tmp187r = tmp171r-tmp179r;
	Re(x[0]) = tmp94r;
	Im(x[0]) = tmp94i;
	Re(x[15]) = tmp121r;
	Im(x[15]) = tmp121i;
	Re(x[30]) = tmp140r;
	Im(x[30]) = tmp140i;
	Re(x[45]) = tmp161r;
	Im(x[45]) = tmp161i;
	Re(x[60]) = tmp186r;
	Im(x[60]) = tmp186i;
	Re(x[75]) = tmp97r;
	Im(x[75]) = tmp97i;
	Re(x[90]) = tmp116r;
	Im(x[90]) = tmp116i;
	Re(x[105]) = tmp143r;
	Im(x[105]) = tmp143i;
	Re(x[120]) = tmp162r;
	Im(x[120]) = tmp162i;
	Re(x[135]) = tmp183r;
	Im(x[135]) = tmp183i;
	Re(x[150]) = tmp98r;
	Im(x[150]) = tmp98i;
	Re(x[165]) = tmp119r;
	Im(x[165]) = tmp119i;
	Re(x[180]) = tmp138r;
	Im(x[180]) = tmp138i;
	Re(x[195]) = tmp165r;
	Im(x[195]) = tmp165i;
	Re(x[210]) = tmp184r;
	Im(x[210]) = tmp184i;
	Re(x[225]) = tmp95r;
	Im(x[225]) = tmp95i;
	Re(x[240]) = tmp120r;
	Im(x[240]) = tmp120i;
	Re(x[255]) = tmp141r;
	Im(x[255]) = tmp141i;
	Re(x[270]) = tmp160r;
	Im(x[270]) = tmp160i;
	Re(x[285]) = tmp187r;
	Im(x[285]) = tmp187i;
	Re(x[300]) = tmp96r;
	Im(x[300]) = tmp96i;
	Re(x[315]) = tmp117r;
	Im(x[315]) = tmp117i;
	Re(x[330]) = tmp142r;
	Im(x[330]) = tmp142i;
	Re(x[345]) = tmp163r;
	Im(x[345]) = tmp163i;
	Re(x[360]) = tmp182r;
	Im(x[360]) = tmp182i;
	Re(x[375]) = tmp99r;
	Im(x[375]) = tmp99i;
	Re(x[390]) = tmp118r;
	Im(x[390]) = tmp118i;
	Re(x[405]) = tmp139r;
	Im(x[405]) = tmp139i;
	Re(x[420]) = tmp164r;
	Im(x[420]) = tmp164i;
	Re(x[435]) = tmp185r;
	Im(x[435]) = tmp185i;
}

/*
*	Number of additions = 372
*	Number of multiplications = 112
*	Number of sign changes = 42
*	Number of assigns = 496
*	Total number of operations = 1022
*/
void	MIFFTC30(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp169r, tmp170r, tmp171r, tmp172r, tmp173r, tmp174r, tmp175r,
		tmp176r, tmp177r, tmp178r, tmp179r, tmp180r, tmp181r, tmp182r,
		tmp183r, tmp184r, tmp185r, tmp186r, tmp187r, tmp188r, tmp189r,
		tmp190r, tmp191r, tmp192r, tmp193r, tmp194r, tmp195r, tmp196r,
		tmp197r, tmp198r, tmp199r, tmp200r, tmp201r, tmp202r, tmp203r,
		tmp204r, tmp205r, tmp206r, tmp207r, tmp208r, tmp209r, tmp210r,
		tmp211r, tmp212r, tmp213r, tmp214r, tmp215r, tmp216r, tmp217r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp155i, tmp156i, tmp157i, tmp158i, tmp159i, tmp160i, tmp161i,
		tmp162i, tmp163i, tmp164i, tmp165i, tmp166i, tmp167i, tmp168i,
		tmp169i, tmp170i, tmp171i, tmp172i, tmp173i, tmp174i, tmp175i,
		tmp176i, tmp177i, tmp178i, tmp179i, tmp180i, tmp181i, tmp182i,
		tmp183i, tmp184i, tmp185i, tmp186i, tmp187i, tmp188i, tmp189i,
		tmp190i, tmp191i, tmp192i, tmp193i, tmp194i, tmp195i, tmp196i,
		tmp197i, tmp198i, tmp199i, tmp200i, tmp201i, tmp202i, tmp203i,
		tmp204i, tmp205i, tmp206i, tmp207i, tmp208i, tmp209i, tmp210i,
		tmp211i, tmp212i, tmp213i, tmp214i, tmp215i, tmp216i, tmp217i;

	const REAL	C1 =                 0.25;	/* REALCONST	*/
	const REAL	C12 =                  0.5;	/* REALCONST	*/
	const REAL	C14 =     0.55901699437495;	/* REALCONST	*/
	const REAL	C10 =     0.58778525229247;	/* REALCONST	*/
	const REAL	C16 =     0.86602540378444;	/* REALCONST	*/
	const REAL	C8 =     0.95105651629515;	/* REALCONST	*/

	tmp0i = Im(x[90])+Im(x[360]);
	tmp0r = Re(x[90])+Re(x[360]);
	tmp1i = Im(x[90])-Im(x[360]);
	tmp1r = Re(x[90])-Re(x[360]);
	tmp2i = Im(x[180])+Im(x[270]);
	tmp2r = Re(x[180])+Re(x[270]);
	tmp3i = Im(x[180])-Im(x[270]);
	tmp3r = Re(x[180])-Re(x[270]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = -C1*tmp4i;
	tmp6r = -C1*tmp4r;
	tmp7i = C14*tmp5i;
	tmp7r = C14*tmp5r;
	tmp8i = tmp6i+Im(x[0]);
	tmp8r = tmp6r+Re(x[0]);
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
	tmp17i = Im(x[0])+tmp4i;
	tmp17r = Re(x[0])+tmp4r;
	tmp18i = Im(x[15])+Im(x[285]);
	tmp18r = Re(x[15])+Re(x[285]);
	tmp19i = Im(x[15])-Im(x[285]);
	tmp19r = Re(x[15])-Re(x[285]);
	tmp20i = Im(x[105])+Im(x[195]);
	tmp20r = Re(x[105])+Re(x[195]);
	tmp21i = Im(x[105])-Im(x[195]);
	tmp21r = Re(x[105])-Re(x[195]);
	tmp22i = tmp18i+tmp20i;
	tmp22r = tmp18r+tmp20r;
	tmp23i = tmp18i-tmp20i;
	tmp23r = tmp18r-tmp20r;
	tmp24i = -C1*tmp22i;
	tmp24r = -C1*tmp22r;
	tmp25i = C14*tmp23i;
	tmp25r = C14*tmp23r;
	tmp26i = tmp24i+Im(x[375]);
	tmp26r = tmp24r+Re(x[375]);
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
	tmp35i = Im(x[375])+tmp22i;
	tmp35r = Re(x[375])+tmp22r;
	tmp36i = Im(x[390])+Im(x[210]);
	tmp36r = Re(x[390])+Re(x[210]);
	tmp37i = Im(x[390])-Im(x[210]);
	tmp37r = Re(x[390])-Re(x[210]);
	tmp38i = Im(x[30])+Im(x[120]);
	tmp38r = Re(x[30])+Re(x[120]);
	tmp39i = Im(x[30])-Im(x[120]);
	tmp39r = Re(x[30])-Re(x[120]);
	tmp40i = tmp36i+tmp38i;
	tmp40r = tmp36r+tmp38r;
	tmp41i = tmp36i-tmp38i;
	tmp41r = tmp36r-tmp38r;
	tmp42i = -C1*tmp40i;
	tmp42r = -C1*tmp40r;
	tmp43i = C14*tmp41i;
	tmp43r = C14*tmp41r;
	tmp44i = tmp42i+Im(x[300]);
	tmp44r = tmp42r+Re(x[300]);
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
	tmp53i = Im(x[300])+tmp40i;
	tmp53r = Re(x[300])+tmp40r;
	tmp54i = Im(x[315])+Im(x[135]);
	tmp54r = Re(x[315])+Re(x[135]);
	tmp55i = Im(x[315])-Im(x[135]);
	tmp55r = Re(x[315])-Re(x[135]);
	tmp56i = Im(x[405])+Im(x[45]);
	tmp56r = Re(x[405])+Re(x[45]);
	tmp57i = Im(x[405])-Im(x[45]);
	tmp57r = Re(x[405])-Re(x[45]);
	tmp58i = tmp54i+tmp56i;
	tmp58r = tmp54r+tmp56r;
	tmp59i = tmp54i-tmp56i;
	tmp59r = tmp54r-tmp56r;
	tmp60i = -C1*tmp58i;
	tmp60r = -C1*tmp58r;
	tmp61i = C14*tmp59i;
	tmp61r = C14*tmp59r;
	tmp62i = tmp60i+Im(x[225]);
	tmp62r = tmp60r+Re(x[225]);
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
	tmp71i = Im(x[225])+tmp58i;
	tmp71r = Re(x[225])+tmp58r;
	tmp72i = Im(x[240])+Im(x[60]);
	tmp72r = Re(x[240])+Re(x[60]);
	tmp73i = Im(x[240])-Im(x[60]);
	tmp73r = Re(x[240])-Re(x[60]);
	tmp74i = Im(x[330])+Im(x[420]);
	tmp74r = Re(x[330])+Re(x[420]);
	tmp75i = Im(x[330])-Im(x[420]);
	tmp75r = Re(x[330])-Re(x[420]);
	tmp76i = tmp72i+tmp74i;
	tmp76r = tmp72r+tmp74r;
	tmp77i = tmp72i-tmp74i;
	tmp77r = tmp72r-tmp74r;
	tmp78i = -C1*tmp76i;
	tmp78r = -C1*tmp76r;
	tmp79i = C14*tmp77i;
	tmp79r = C14*tmp77r;
	tmp80i = tmp78i+Im(x[150]);
	tmp80r = tmp78r+Re(x[150]);
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
	tmp89i = Im(x[150])+tmp76i;
	tmp89r = Re(x[150])+tmp76r;
	tmp90i = Im(x[165])+Im(x[435]);
	tmp90r = Re(x[165])+Re(x[435]);
	tmp91i = Im(x[165])-Im(x[435]);
	tmp91r = Re(x[165])-Re(x[435]);
	tmp92i = Im(x[255])+Im(x[345]);
	tmp92r = Re(x[255])+Re(x[345]);
	tmp93i = Im(x[255])-Im(x[345]);
	tmp93r = Re(x[255])-Re(x[345]);
	tmp94i = tmp90i+tmp92i;
	tmp94r = tmp90r+tmp92r;
	tmp95i = tmp90i-tmp92i;
	tmp95r = tmp90r-tmp92r;
	tmp96i = -C1*tmp94i;
	tmp96r = -C1*tmp94r;
	tmp97i = C14*tmp95i;
	tmp97r = C14*tmp95r;
	tmp98i = tmp96i+Im(x[75]);
	tmp98r = tmp96r+Re(x[75]);
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
	tmp107i = Im(x[75])+tmp94i;
	tmp107r = Re(x[75])+tmp94r;
	tmp108i = tmp89i+tmp53i;
	tmp108r = tmp89r+tmp53r;
	tmp109i = tmp89i-tmp53i;
	tmp109r = tmp89r-tmp53r;
	tmp110i = -C12*tmp108i;
	tmp110r = -C12*tmp108r;
	tmp111i = C16*tmp109r;
	tmp111r = -C16*tmp109i;
	tmp112i = tmp110i+tmp17i;
	tmp112r = tmp110r+tmp17r;
	tmp113i = tmp112i+tmp111i;
	tmp113r = tmp112r+tmp111r;
	tmp114i = tmp112i-tmp111i;
	tmp114r = tmp112r-tmp111r;
	tmp115i = tmp17i+tmp108i;
	tmp115r = tmp17r+tmp108r;
	tmp116i = tmp35i+tmp107i;
	tmp116r = tmp35r+tmp107r;
	tmp117i = tmp35i-tmp107i;
	tmp117r = tmp35r-tmp107r;
	tmp118i = -C12*tmp116i;
	tmp118r = -C12*tmp116r;
	tmp119i = C16*tmp117r;
	tmp119r = -C16*tmp117i;
	tmp120i = tmp118i+tmp71i;
	tmp120r = tmp118r+tmp71r;
	tmp121i = tmp120i+tmp119i;
	tmp121r = tmp120r+tmp119r;
	tmp122i = tmp120i-tmp119i;
	tmp122r = tmp120r-tmp119r;
	tmp123i = tmp71i+tmp116i;
	tmp123r = tmp71r+tmp116r;
	tmp124i = tmp115i+tmp123i;
	tmp124r = tmp115r+tmp123r;
	tmp125i = tmp115i-tmp123i;
	tmp125r = tmp115r-tmp123r;
	tmp126i = tmp114i+tmp122i;
	tmp126r = tmp114r+tmp122r;
	tmp127i = tmp114i-tmp122i;
	tmp127r = tmp114r-tmp122r;
	tmp128i = tmp113i+tmp121i;
	tmp128r = tmp113r+tmp121r;
	tmp129i = tmp113i-tmp121i;
	tmp129r = tmp113r-tmp121r;
	tmp130i = tmp85i+tmp49i;
	tmp130r = tmp85r+tmp49r;
	tmp131i = tmp85i-tmp49i;
	tmp131r = tmp85r-tmp49r;
	tmp132i = -C12*tmp130i;
	tmp132r = -C12*tmp130r;
	tmp133i = C16*tmp131r;
	tmp133r = -C16*tmp131i;
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
	tmp140i = -C12*tmp138i;
	tmp140r = -C12*tmp138r;
	tmp141i = C16*tmp139r;
	tmp141r = -C16*tmp139i;
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
	tmp147i = tmp137i-tmp145i;
	tmp147r = tmp137r-tmp145r;
	tmp148i = tmp136i+tmp144i;
	tmp148r = tmp136r+tmp144r;
	tmp149i = tmp136i-tmp144i;
	tmp149r = tmp136r-tmp144r;
	tmp150i = tmp135i+tmp143i;
	tmp150r = tmp135r+tmp143r;
	tmp151i = tmp135i-tmp143i;
	tmp151r = tmp135r-tmp143r;
	tmp152i = tmp87i+tmp51i;
	tmp152r = tmp87r+tmp51r;
	tmp153i = tmp87i-tmp51i;
	tmp153r = tmp87r-tmp51r;
	tmp154i = -C12*tmp152i;
	tmp154r = -C12*tmp152r;
	tmp155i = C16*tmp153r;
	tmp155r = -C16*tmp153i;
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
	tmp162i = -C12*tmp160i;
	tmp162r = -C12*tmp160r;
	tmp163i = C16*tmp161r;
	tmp163r = -C16*tmp161i;
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
	tmp169i = tmp159i-tmp167i;
	tmp169r = tmp159r-tmp167r;
	tmp170i = tmp158i+tmp166i;
	tmp170r = tmp158r+tmp166r;
	tmp171i = tmp158i-tmp166i;
	tmp171r = tmp158r-tmp166r;
	tmp172i = tmp157i+tmp165i;
	tmp172r = tmp157r+tmp165r;
	tmp173i = tmp157i-tmp165i;
	tmp173r = tmp157r-tmp165r;
	tmp174i = tmp88i+tmp52i;
	tmp174r = tmp88r+tmp52r;
	tmp175i = tmp88i-tmp52i;
	tmp175r = tmp88r-tmp52r;
	tmp176i = -C12*tmp174i;
	tmp176r = -C12*tmp174r;
	tmp177i = C16*tmp175r;
	tmp177r = -C16*tmp175i;
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
	tmp184i = -C12*tmp182i;
	tmp184r = -C12*tmp182r;
	tmp185i = C16*tmp183r;
	tmp185r = -C16*tmp183i;
	tmp186i = tmp184i+tmp70i;
	tmp186r = tmp184r+tmp70r;
	tmp187i = tmp186i+tmp185i;
	tmp187r = tmp186r+tmp185r;
	tmp188i = tmp186i-tmp185i;
	tmp188r = tmp186r-tmp185r;
	tmp189i = tmp70i+tmp182i;
	tmp189r = tmp70r+tmp182r;
	tmp190i = tmp181i+tmp189i;
	tmp190r = tmp181r+tmp189r;
	tmp191i = tmp181i-tmp189i;
	tmp191r = tmp181r-tmp189r;
	tmp192i = tmp180i+tmp188i;
	tmp192r = tmp180r+tmp188r;
	tmp193i = tmp180i-tmp188i;
	tmp193r = tmp180r-tmp188r;
	tmp194i = tmp179i+tmp187i;
	tmp194r = tmp179r+tmp187r;
	tmp195i = tmp179i-tmp187i;
	tmp195r = tmp179r-tmp187r;
	tmp196i = tmp86i+tmp50i;
	tmp196r = tmp86r+tmp50r;
	tmp197i = tmp86i-tmp50i;
	tmp197r = tmp86r-tmp50r;
	tmp198i = -C12*tmp196i;
	tmp198r = -C12*tmp196r;
	tmp199i = C16*tmp197r;
	tmp199r = -C16*tmp197i;
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
	tmp206i = -C12*tmp204i;
	tmp206r = -C12*tmp204r;
	tmp207i = C16*tmp205r;
	tmp207r = -C16*tmp205i;
	tmp208i = tmp206i+tmp68i;
	tmp208r = tmp206r+tmp68r;
	tmp209i = tmp208i+tmp207i;
	tmp209r = tmp208r+tmp207r;
	tmp210i = tmp208i-tmp207i;
	tmp210r = tmp208r-tmp207r;
	tmp211i = tmp68i+tmp204i;
	tmp211r = tmp68r+tmp204r;
	tmp212i = tmp203i+tmp211i;
	tmp212r = tmp203r+tmp211r;
	tmp213i = tmp203i-tmp211i;
	tmp213r = tmp203r-tmp211r;
	tmp214i = tmp202i+tmp210i;
	tmp214r = tmp202r+tmp210r;
	tmp215i = tmp202i-tmp210i;
	tmp215r = tmp202r-tmp210r;
	tmp216i = tmp201i+tmp209i;
	tmp216r = tmp201r+tmp209r;
	tmp217i = tmp201i-tmp209i;
	tmp217r = tmp201r-tmp209r;
	Re(x[0]) = tmp124r;
	Im(x[0]) = tmp124i;
	Re(x[15]) = tmp151r;
	Im(x[15]) = tmp151i;
	Re(x[30]) = tmp170r;
	Im(x[30]) = tmp170i;
	Re(x[45]) = tmp191r;
	Im(x[45]) = tmp191i;
	Re(x[60]) = tmp216r;
	Im(x[60]) = tmp216i;
	Re(x[75]) = tmp127r;
	Im(x[75]) = tmp127i;
	Re(x[90]) = tmp146r;
	Im(x[90]) = tmp146i;
	Re(x[105]) = tmp173r;
	Im(x[105]) = tmp173i;
	Re(x[120]) = tmp192r;
	Im(x[120]) = tmp192i;
	Re(x[135]) = tmp213r;
	Im(x[135]) = tmp213i;
	Re(x[150]) = tmp128r;
	Im(x[150]) = tmp128i;
	Re(x[165]) = tmp149r;
	Im(x[165]) = tmp149i;
	Re(x[180]) = tmp168r;
	Im(x[180]) = tmp168i;
	Re(x[195]) = tmp195r;
	Im(x[195]) = tmp195i;
	Re(x[210]) = tmp214r;
	Im(x[210]) = tmp214i;
	Re(x[225]) = tmp125r;
	Im(x[225]) = tmp125i;
	Re(x[240]) = tmp150r;
	Im(x[240]) = tmp150i;
	Re(x[255]) = tmp171r;
	Im(x[255]) = tmp171i;
	Re(x[270]) = tmp190r;
	Im(x[270]) = tmp190i;
	Re(x[285]) = tmp217r;
	Im(x[285]) = tmp217i;
	Re(x[300]) = tmp126r;
	Im(x[300]) = tmp126i;
	Re(x[315]) = tmp147r;
	Im(x[315]) = tmp147i;
	Re(x[330]) = tmp172r;
	Im(x[330]) = tmp172i;
	Re(x[345]) = tmp193r;
	Im(x[345]) = tmp193i;
	Re(x[360]) = tmp212r;
	Im(x[360]) = tmp212i;
	Re(x[375]) = tmp129r;
	Im(x[375]) = tmp129i;
	Re(x[390]) = tmp148r;
	Im(x[390]) = tmp148i;
	Re(x[405]) = tmp169r;
	Im(x[405]) = tmp169i;
	Re(x[420]) = tmp194r;
	Im(x[420]) = tmp194i;
	Re(x[435]) = tmp215r;
	Im(x[435]) = tmp215i;
}

/*
*	Number of additions = 308
*	Number of multiplications = 84
*	Number of sign changes = 10
*	Number of assigns = 380
*	Total number of operations = 782
*/
void	MFFTC32(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp155r, tmp156r, tmp157r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp155i, tmp156i, tmp157i;

	const REAL	C3 =     0.19509032201613;	/* REALCONST	*/
	const REAL	C5 =     0.38268343236509;	/* REALCONST	*/
	const REAL	C7 =      0.5555702330196;	/* REALCONST	*/
	const REAL	C8 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C6 =     0.83146961230255;	/* REALCONST	*/
	const REAL	C4 =     0.92387953251129;	/* REALCONST	*/
	const REAL	C2 =     0.98078528040323;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[128]);
	tmp0r = Re(x[0])+Re(x[128]);
	tmp1i = Im(x[0])-Im(x[128]);
	tmp1r = Re(x[0])-Re(x[128]);
	tmp2i = Im(x[0])-Re(x[128]);
	tmp2r = Re(x[0])+Im(x[128]);
	tmp3i = Im(x[0])+Re(x[128]);
	tmp3r = Re(x[0])-Im(x[128]);
	tmp4i = Im(x[64])+Im(x[192]);
	tmp4r = Re(x[64])+Re(x[192]);
	tmp5i = Im(x[64])-Im(x[192]);
	tmp5r = Re(x[64])-Re(x[192]);
	tmp6i = Im(x[64])-Re(x[192]);
	tmp6r = Re(x[64])+Im(x[192]);
	tmp7i = Im(x[64])+Re(x[192]);
	tmp7r = Re(x[64])-Im(x[192]);
	tmp8i = C8*(tmp6i-tmp6r);
	tmp8r = C8*(tmp6r+tmp6i);
	tmp9i = -C8*(tmp7i+tmp7r);
	tmp9r = -C8*(tmp7r-tmp7i);
	tmp10i = tmp0i+tmp4i;
	tmp10r = tmp0r+tmp4r;
	tmp11i = tmp0i-tmp4i;
	tmp11r = tmp0r-tmp4r;
	tmp12i = tmp2i+tmp8i;
	tmp12r = tmp2r+tmp8r;
	tmp13i = tmp2i-tmp8i;
	tmp13r = tmp2r-tmp8r;
	tmp14i = tmp1i-tmp5r;
	tmp14r = tmp1r+tmp5i;
	tmp15i = tmp1i+tmp5r;
	tmp15r = tmp1r-tmp5i;
	tmp16i = tmp3i+tmp9i;
	tmp16r = tmp3r+tmp9r;
	tmp17i = tmp3i-tmp9i;
	tmp17r = tmp3r-tmp9r;
	tmp18i = Im(x[32])+Im(x[160]);
	tmp18r = Re(x[32])+Re(x[160]);
	tmp19i = Im(x[32])-Im(x[160]);
	tmp19r = Re(x[32])-Re(x[160]);
	tmp20i = Im(x[32])-Re(x[160]);
	tmp20r = Re(x[32])+Im(x[160]);
	tmp21i = Im(x[32])+Re(x[160]);
	tmp21r = Re(x[32])-Im(x[160]);
	tmp22r = tmp18r;
	tmp22i = tmp18i;
	tmp23r = C4*tmp20r+C5*tmp20i;
	tmp23i = C4*tmp20i-C5*tmp20r;
	tmp24r = C8*(tmp19r+tmp19i);
	tmp24i = C8*(tmp19i-tmp19r);
	tmp25r = C5*tmp21r+C4*tmp21i;
	tmp25i = C5*tmp21i-C4*tmp21r;
	tmp26i = Im(x[96])+Im(x[224]);
	tmp26r = Re(x[96])+Re(x[224]);
	tmp27i = Im(x[96])-Im(x[224]);
	tmp27r = Re(x[96])-Re(x[224]);
	tmp28i = Im(x[96])-Re(x[224]);
	tmp28r = Re(x[96])+Im(x[224]);
	tmp29i = Im(x[96])+Re(x[224]);
	tmp29r = Re(x[96])-Im(x[224]);
	tmp30r = tmp26r;
	tmp30i = tmp26i;
	tmp31r = C5*tmp28r+C4*tmp28i;
	tmp31i = C5*tmp28i-C4*tmp28r;
	tmp32r = -C8*(tmp27r-tmp27i);
	tmp32i = -C8*(tmp27i+tmp27r);
	tmp33r = -C4*tmp29r-C5*tmp29i;
	tmp33i = -C4*tmp29i+C5*tmp29r;
	tmp34r = tmp22r+tmp30r;
	tmp34i = tmp22i+tmp30i;
	tmp35r = tmp23r+tmp31r;
	tmp35i = tmp23i+tmp31i;
	tmp36r = tmp24r+tmp32r;
	tmp36i = tmp24i+tmp32i;
	tmp37r = tmp25r+tmp33r;
	tmp37i = tmp25i+tmp33i;
	tmp38r = -tmp22i+tmp30i;
	tmp38i = tmp22r-tmp30r;
	tmp39r = -tmp23i+tmp31i;
	tmp39i = tmp23r-tmp31r;
	tmp40r = -tmp24i+tmp32i;
	tmp40i = tmp24r-tmp32r;
	tmp41r = -tmp25i+tmp33i;
	tmp41i = tmp25r-tmp33r;
	tmp42r = tmp10r+tmp34r;
	tmp42i = tmp10i+tmp34i;
	tmp43r = tmp12r+tmp35r;
	tmp43i = tmp12i+tmp35i;
	tmp44r = tmp14r+tmp36r;
	tmp44i = tmp14i+tmp36i;
	tmp45r = tmp16r+tmp37r;
	tmp45i = tmp16i+tmp37i;
	tmp46r = tmp11r-tmp38r;
	tmp46i = tmp11i-tmp38i;
	tmp47r = tmp13r-tmp39r;
	tmp47i = tmp13i-tmp39i;
	tmp48r = tmp15r-tmp40r;
	tmp48i = tmp15i-tmp40i;
	tmp49r = tmp17r-tmp41r;
	tmp49i = tmp17i-tmp41i;
	tmp50r = tmp10r-tmp34r;
	tmp50i = tmp10i-tmp34i;
	tmp51r = tmp12r-tmp35r;
	tmp51i = tmp12i-tmp35i;
	tmp52r = tmp14r-tmp36r;
	tmp52i = tmp14i-tmp36i;
	tmp53r = tmp16r-tmp37r;
	tmp53i = tmp16i-tmp37i;
	tmp54r = tmp11r+tmp38r;
	tmp54i = tmp11i+tmp38i;
	tmp55r = tmp13r+tmp39r;
	tmp55i = tmp13i+tmp39i;
	tmp56r = tmp15r+tmp40r;
	tmp56i = tmp15i+tmp40i;
	tmp57r = tmp17r+tmp41r;
	tmp57i = tmp17i+tmp41i;
	tmp58i = Im(x[16])+Im(x[144]);
	tmp58r = Re(x[16])+Re(x[144]);
	tmp59i = Im(x[16])-Im(x[144]);
	tmp59r = Re(x[16])-Re(x[144]);
	tmp60i = Im(x[16])-Re(x[144]);
	tmp60r = Re(x[16])+Im(x[144]);
	tmp61i = Im(x[16])+Re(x[144]);
	tmp61r = Re(x[16])-Im(x[144]);
	tmp62i = Im(x[80])+Im(x[208]);
	tmp62r = Re(x[80])+Re(x[208]);
	tmp63i = Im(x[80])-Im(x[208]);
	tmp63r = Re(x[80])-Re(x[208]);
	tmp64i = Im(x[80])-Re(x[208]);
	tmp64r = Re(x[80])+Im(x[208]);
	tmp65i = Im(x[80])+Re(x[208]);
	tmp65r = Re(x[80])-Im(x[208]);
	tmp66i = C8*(tmp64i-tmp64r);
	tmp66r = C8*(tmp64r+tmp64i);
	tmp67i = -C8*(tmp65i+tmp65r);
	tmp67r = -C8*(tmp65r-tmp65i);
	tmp68i = tmp58i+tmp62i;
	tmp68r = tmp58r+tmp62r;
	tmp69i = tmp58i-tmp62i;
	tmp69r = tmp58r-tmp62r;
	tmp70i = tmp60i+tmp66i;
	tmp70r = tmp60r+tmp66r;
	tmp71i = tmp60i-tmp66i;
	tmp71r = tmp60r-tmp66r;
	tmp72i = tmp59i-tmp63r;
	tmp72r = tmp59r+tmp63i;
	tmp73i = tmp59i+tmp63r;
	tmp73r = tmp59r-tmp63i;
	tmp74i = tmp61i+tmp67i;
	tmp74r = tmp61r+tmp67r;
	tmp75i = tmp61i-tmp67i;
	tmp75r = tmp61r-tmp67r;
	tmp76r = tmp68r;
	tmp76i = tmp68i;
	tmp77r = C2*tmp70r+C3*tmp70i;
	tmp77i = C2*tmp70i-C3*tmp70r;
	tmp78r = C4*tmp72r+C5*tmp72i;
	tmp78i = C4*tmp72i-C5*tmp72r;
	tmp79r = C6*tmp74r+C7*tmp74i;
	tmp79i = C6*tmp74i-C7*tmp74r;
	tmp80r = C8*(tmp69r+tmp69i);
	tmp80i = C8*(tmp69i-tmp69r);
	tmp81r = C7*tmp71r+C6*tmp71i;
	tmp81i = C7*tmp71i-C6*tmp71r;
	tmp82r = C5*tmp73r+C4*tmp73i;
	tmp82i = C5*tmp73i-C4*tmp73r;
	tmp83r = C3*tmp75r+C2*tmp75i;
	tmp83i = C3*tmp75i-C2*tmp75r;
	tmp84i = Im(x[48])+Im(x[176]);
	tmp84r = Re(x[48])+Re(x[176]);
	tmp85i = Im(x[48])-Im(x[176]);
	tmp85r = Re(x[48])-Re(x[176]);
	tmp86i = Im(x[48])-Re(x[176]);
	tmp86r = Re(x[48])+Im(x[176]);
	tmp87i = Im(x[48])+Re(x[176]);
	tmp87r = Re(x[48])-Im(x[176]);
	tmp88i = Im(x[112])+Im(x[240]);
	tmp88r = Re(x[112])+Re(x[240]);
	tmp89i = Im(x[112])-Im(x[240]);
	tmp89r = Re(x[112])-Re(x[240]);
	tmp90i = Im(x[112])-Re(x[240]);
	tmp90r = Re(x[112])+Im(x[240]);
	tmp91i = Im(x[112])+Re(x[240]);
	tmp91r = Re(x[112])-Im(x[240]);
	tmp92i = C8*(tmp90i-tmp90r);
	tmp92r = C8*(tmp90r+tmp90i);
	tmp93i = -C8*(tmp91i+tmp91r);
	tmp93r = -C8*(tmp91r-tmp91i);
	tmp94i = tmp84i+tmp88i;
	tmp94r = tmp84r+tmp88r;
	tmp95i = tmp84i-tmp88i;
	tmp95r = tmp84r-tmp88r;
	tmp96i = tmp86i+tmp92i;
	tmp96r = tmp86r+tmp92r;
	tmp97i = tmp86i-tmp92i;
	tmp97r = tmp86r-tmp92r;
	tmp98i = tmp85i-tmp89r;
	tmp98r = tmp85r+tmp89i;
	tmp99i = tmp85i+tmp89r;
	tmp99r = tmp85r-tmp89i;
	tmp100i = tmp87i+tmp93i;
	tmp100r = tmp87r+tmp93r;
	tmp101i = tmp87i-tmp93i;
	tmp101r = tmp87r-tmp93r;
	tmp102r = tmp94r;
	tmp102i = tmp94i;
	tmp103r = C6*tmp96r+C7*tmp96i;
	tmp103i = C6*tmp96i-C7*tmp96r;
	tmp104r = C5*tmp98r+C4*tmp98i;
	tmp104i = C5*tmp98i-C4*tmp98r;
	tmp105r = -C3*tmp100r+C2*tmp100i;
	tmp105i = -C3*tmp100i-C2*tmp100r;
	tmp106r = -C8*(tmp95r-tmp95i);
	tmp106i = -C8*(tmp95i+tmp95r);
	tmp107r = -C2*tmp97r+C3*tmp97i;
	tmp107i = -C2*tmp97i-C3*tmp97r;
	tmp108r = -C4*tmp99r-C5*tmp99i;
	tmp108i = -C4*tmp99i+C5*tmp99r;
	tmp109r = -C7*tmp101r-C6*tmp101i;
	tmp109i = -C7*tmp101i+C6*tmp101r;
	tmp110r = tmp76r+tmp102r;
	tmp110i = tmp76i+tmp102i;
	tmp111r = tmp77r+tmp103r;
	tmp111i = tmp77i+tmp103i;
	tmp112r = tmp78r+tmp104r;
	tmp112i = tmp78i+tmp104i;
	tmp113r = tmp79r+tmp105r;
	tmp113i = tmp79i+tmp105i;
	tmp114r = tmp80r+tmp106r;
	tmp114i = tmp80i+tmp106i;
	tmp115r = tmp81r+tmp107r;
	tmp115i = tmp81i+tmp107i;
	tmp116r = tmp82r+tmp108r;
	tmp116i = tmp82i+tmp108i;
	tmp117r = tmp83r+tmp109r;
	tmp117i = tmp83i+tmp109i;
	tmp118r = -tmp76i+tmp102i;
	tmp118i = tmp76r-tmp102r;
	tmp119r = -tmp77i+tmp103i;
	tmp119i = tmp77r-tmp103r;
	tmp120r = -tmp78i+tmp104i;
	tmp120i = tmp78r-tmp104r;
	tmp121r = -tmp79i+tmp105i;
	tmp121i = tmp79r-tmp105r;
	tmp122r = -tmp80i+tmp106i;
	tmp122i = tmp80r-tmp106r;
	tmp123r = -tmp81i+tmp107i;
	tmp123i = tmp81r-tmp107r;
	tmp124r = -tmp82i+tmp108i;
	tmp124i = tmp82r-tmp108r;
	tmp125r = -tmp83i+tmp109i;
	tmp125i = tmp83r-tmp109r;
	tmp126r = tmp42r+tmp110r;
	tmp126i = tmp42i+tmp110i;
	tmp127r = tmp43r+tmp111r;
	tmp127i = tmp43i+tmp111i;
	tmp128r = tmp44r+tmp112r;
	tmp128i = tmp44i+tmp112i;
	tmp129r = tmp45r+tmp113r;
	tmp129i = tmp45i+tmp113i;
	tmp130r = tmp46r+tmp114r;
	tmp130i = tmp46i+tmp114i;
	tmp131r = tmp47r+tmp115r;
	tmp131i = tmp47i+tmp115i;
	tmp132r = tmp48r+tmp116r;
	tmp132i = tmp48i+tmp116i;
	tmp133r = tmp49r+tmp117r;
	tmp133i = tmp49i+tmp117i;
	tmp134r = tmp50r-tmp118r;
	tmp134i = tmp50i-tmp118i;
	tmp135r = tmp51r-tmp119r;
	tmp135i = tmp51i-tmp119i;
	tmp136r = tmp52r-tmp120r;
	tmp136i = tmp52i-tmp120i;
	tmp137r = tmp53r-tmp121r;
	tmp137i = tmp53i-tmp121i;
	tmp138r = tmp54r-tmp122r;
	tmp138i = tmp54i-tmp122i;
	tmp139r = tmp55r-tmp123r;
	tmp139i = tmp55i-tmp123i;
	tmp140r = tmp56r-tmp124r;
	tmp140i = tmp56i-tmp124i;
	tmp141r = tmp57r-tmp125r;
	tmp141i = tmp57i-tmp125i;
	tmp142r = tmp42r-tmp110r;
	tmp142i = tmp42i-tmp110i;
	tmp143r = tmp43r-tmp111r;
	tmp143i = tmp43i-tmp111i;
	tmp144r = tmp44r-tmp112r;
	tmp144i = tmp44i-tmp112i;
	tmp145r = tmp45r-tmp113r;
	tmp145i = tmp45i-tmp113i;
	tmp146r = tmp46r-tmp114r;
	tmp146i = tmp46i-tmp114i;
	tmp147r = tmp47r-tmp115r;
	tmp147i = tmp47i-tmp115i;
	tmp148r = tmp48r-tmp116r;
	tmp148i = tmp48i-tmp116i;
	tmp149r = tmp49r-tmp117r;
	tmp149i = tmp49i-tmp117i;
	tmp150r = tmp50r+tmp118r;
	tmp150i = tmp50i+tmp118i;
	tmp151r = tmp51r+tmp119r;
	tmp151i = tmp51i+tmp119i;
	tmp152r = tmp52r+tmp120r;
	tmp152i = tmp52i+tmp120i;
	tmp153r = tmp53r+tmp121r;
	tmp153i = tmp53i+tmp121i;
	tmp154r = tmp54r+tmp122r;
	tmp154i = tmp54i+tmp122i;
	tmp155r = tmp55r+tmp123r;
	tmp155i = tmp55i+tmp123i;
	tmp156r = tmp56r+tmp124r;
	tmp156i = tmp56i+tmp124i;
	tmp157r = tmp57r+tmp125r;
	tmp157i = tmp57i+tmp125i;
	Re(x[0]) = tmp126r;
	Im(x[0]) = tmp126i;
	Re(x[16]) = tmp127r;
	Im(x[16]) = tmp127i;
	Re(x[32]) = tmp128r;
	Im(x[32]) = tmp128i;
	Re(x[48]) = tmp129r;
	Im(x[48]) = tmp129i;
	Re(x[64]) = tmp130r;
	Im(x[64]) = tmp130i;
	Re(x[80]) = tmp131r;
	Im(x[80]) = tmp131i;
	Re(x[96]) = tmp132r;
	Im(x[96]) = tmp132i;
	Re(x[112]) = tmp133r;
	Im(x[112]) = tmp133i;
	Re(x[128]) = tmp134r;
	Im(x[128]) = tmp134i;
	Re(x[144]) = tmp135r;
	Im(x[144]) = tmp135i;
	Re(x[160]) = tmp136r;
	Im(x[160]) = tmp136i;
	Re(x[176]) = tmp137r;
	Im(x[176]) = tmp137i;
	Re(x[192]) = tmp138r;
	Im(x[192]) = tmp138i;
	Re(x[208]) = tmp139r;
	Im(x[208]) = tmp139i;
	Re(x[224]) = tmp140r;
	Im(x[224]) = tmp140i;
	Re(x[240]) = tmp141r;
	Im(x[240]) = tmp141i;
	Re(x[256]) = tmp142r;
	Im(x[256]) = tmp142i;
	Re(x[272]) = tmp143r;
	Im(x[272]) = tmp143i;
	Re(x[288]) = tmp144r;
	Im(x[288]) = tmp144i;
	Re(x[304]) = tmp145r;
	Im(x[304]) = tmp145i;
	Re(x[320]) = tmp146r;
	Im(x[320]) = tmp146i;
	Re(x[336]) = tmp147r;
	Im(x[336]) = tmp147i;
	Re(x[352]) = tmp148r;
	Im(x[352]) = tmp148i;
	Re(x[368]) = tmp149r;
	Im(x[368]) = tmp149i;
	Re(x[384]) = tmp150r;
	Im(x[384]) = tmp150i;
	Re(x[400]) = tmp151r;
	Im(x[400]) = tmp151i;
	Re(x[416]) = tmp152r;
	Im(x[416]) = tmp152i;
	Re(x[432]) = tmp153r;
	Im(x[432]) = tmp153i;
	Re(x[448]) = tmp154r;
	Im(x[448]) = tmp154i;
	Re(x[464]) = tmp155r;
	Im(x[464]) = tmp155i;
	Re(x[480]) = tmp156r;
	Im(x[480]) = tmp156i;
	Re(x[496]) = tmp157r;
	Im(x[496]) = tmp157i;
}

/*
*	Number of additions = 372
*	Number of multiplications = 84
*	Number of sign changes = 22
*	Number of assigns = 444
*	Total number of operations = 922
*/
void	MIFFTC32(Complex *x)
{
	register REAL	tmp0r, tmp1r, tmp2r, tmp3r, tmp4r, tmp5r, tmp6r, tmp7r,
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
		tmp169r, tmp170r, tmp171r, tmp172r, tmp173r, tmp174r, tmp175r,
		tmp176r, tmp177r, tmp178r, tmp179r, tmp180r, tmp181r, tmp182r,
		tmp183r, tmp184r, tmp185r, tmp186r, tmp187r, tmp188r, tmp189r;
	register REAL	tmp0i, tmp1i, tmp2i, tmp3i, tmp4i, tmp5i, tmp6i, tmp7i,
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
		tmp155i, tmp156i, tmp157i, tmp158i, tmp159i, tmp160i, tmp161i,
		tmp162i, tmp163i, tmp164i, tmp165i, tmp166i, tmp167i, tmp168i,
		tmp169i, tmp170i, tmp171i, tmp172i, tmp173i, tmp174i, tmp175i,
		tmp176i, tmp177i, tmp178i, tmp179i, tmp180i, tmp181i, tmp182i,
		tmp183i, tmp184i, tmp185i, tmp186i, tmp187i, tmp188i, tmp189i;

	const REAL	C3 =     0.19509032201613;	/* REALCONST	*/
	const REAL	C5 =     0.38268343236509;	/* REALCONST	*/
	const REAL	C7 =      0.5555702330196;	/* REALCONST	*/
	const REAL	C8 =     0.70710678118655;	/* REALCONST	*/
	const REAL	C6 =     0.83146961230255;	/* REALCONST	*/
	const REAL	C4 =     0.92387953251129;	/* REALCONST	*/
	const REAL	C2 =     0.98078528040323;	/* REALCONST	*/

	tmp0i = Im(x[0])+Im(x[256]);
	tmp0r = Re(x[0])+Re(x[256]);
	tmp1i = Im(x[0])-Im(x[256]);
	tmp1r = Re(x[0])-Re(x[256]);
	tmp2i = Im(x[128])+Im(x[384]);
	tmp2r = Re(x[128])+Re(x[384]);
	tmp3i = Im(x[128])-Im(x[384]);
	tmp3r = Re(x[128])-Re(x[384]);
	tmp4i = tmp0i+tmp2i;
	tmp4r = tmp0r+tmp2r;
	tmp5i = tmp0i-tmp2i;
	tmp5r = tmp0r-tmp2r;
	tmp6i = tmp1i+tmp3r;
	tmp6r = tmp1r-tmp3i;
	tmp7i = tmp1i-tmp3r;
	tmp7r = tmp1r+tmp3i;
	tmp8i = Im(x[64])+Im(x[320]);
	tmp8r = Re(x[64])+Re(x[320]);
	tmp9i = Im(x[64])-Im(x[320]);
	tmp9r = Re(x[64])-Re(x[320]);
	tmp10i = Im(x[192])+Im(x[448]);
	tmp10r = Re(x[192])+Re(x[448]);
	tmp11i = Im(x[192])-Im(x[448]);
	tmp11r = Re(x[192])-Re(x[448]);
	tmp12i = tmp8i+tmp10i;
	tmp12r = tmp8r+tmp10r;
	tmp13i = tmp8i-tmp10i;
	tmp13r = tmp8r-tmp10r;
	tmp14i = tmp9i+tmp11r;
	tmp14r = tmp9r-tmp11i;
	tmp15i = tmp9i-tmp11r;
	tmp15r = tmp9r+tmp11i;
	tmp16i = C8*(tmp14i+tmp14r);
	tmp16r = C8*(tmp14r-tmp14i);
	tmp17i = -C8*(tmp15i-tmp15r);
	tmp17r = -C8*(tmp15r+tmp15i);
	tmp18i = tmp4i+tmp12i;
	tmp18r = tmp4r+tmp12r;
	tmp19i = tmp4i-tmp12i;
	tmp19r = tmp4r-tmp12r;
	tmp20i = tmp6i+tmp16i;
	tmp20r = tmp6r+tmp16r;
	tmp21i = tmp6i-tmp16i;
	tmp21r = tmp6r-tmp16r;
	tmp22i = tmp5i+tmp13r;
	tmp22r = tmp5r-tmp13i;
	tmp23i = tmp5i-tmp13r;
	tmp23r = tmp5r+tmp13i;
	tmp24i = tmp7i+tmp17i;
	tmp24r = tmp7r+tmp17r;
	tmp25i = tmp7i-tmp17i;
	tmp25r = tmp7r-tmp17r;
	tmp26i = Im(x[32])+Im(x[288]);
	tmp26r = Re(x[32])+Re(x[288]);
	tmp27i = Im(x[32])-Im(x[288]);
	tmp27r = Re(x[32])-Re(x[288]);
	tmp28i = Im(x[160])+Im(x[416]);
	tmp28r = Re(x[160])+Re(x[416]);
	tmp29i = Im(x[160])-Im(x[416]);
	tmp29r = Re(x[160])-Re(x[416]);
	tmp30i = tmp26i+tmp28i;
	tmp30r = tmp26r+tmp28r;
	tmp31i = tmp26i-tmp28i;
	tmp31r = tmp26r-tmp28r;
	tmp32i = tmp27i+tmp29r;
	tmp32r = tmp27r-tmp29i;
	tmp33i = tmp27i-tmp29r;
	tmp33r = tmp27r+tmp29i;
	tmp34r = tmp30r;
	tmp34i = tmp30i;
	tmp35r = C4*tmp32r-C5*tmp32i;
	tmp35i = C4*tmp32i+C5*tmp32r;
	tmp36r = C8*(tmp31r-tmp31i);
	tmp36i = C8*(tmp31i+tmp31r);
	tmp37r = C5*tmp33r-C4*tmp33i;
	tmp37i = C5*tmp33i+C4*tmp33r;
	tmp38i = Im(x[96])+Im(x[352]);
	tmp38r = Re(x[96])+Re(x[352]);
	tmp39i = Im(x[96])-Im(x[352]);
	tmp39r = Re(x[96])-Re(x[352]);
	tmp40i = Im(x[224])+Im(x[480]);
	tmp40r = Re(x[224])+Re(x[480]);
	tmp41i = Im(x[224])-Im(x[480]);
	tmp41r = Re(x[224])-Re(x[480]);
	tmp42i = tmp38i+tmp40i;
	tmp42r = tmp38r+tmp40r;
	tmp43i = tmp38i-tmp40i;
	tmp43r = tmp38r-tmp40r;
	tmp44i = tmp39i+tmp41r;
	tmp44r = tmp39r-tmp41i;
	tmp45i = tmp39i-tmp41r;
	tmp45r = tmp39r+tmp41i;
	tmp46r = tmp42r;
	tmp46i = tmp42i;
	tmp47r = C5*tmp44r-C4*tmp44i;
	tmp47i = C5*tmp44i+C4*tmp44r;
	tmp48r = -C8*(tmp43r+tmp43i);
	tmp48i = -C8*(tmp43i-tmp43r);
	tmp49r = -C4*tmp45r+C5*tmp45i;
	tmp49i = -C4*tmp45i-C5*tmp45r;
	tmp50r = tmp34r+tmp46r;
	tmp50i = tmp34i+tmp46i;
	tmp51r = tmp35r+tmp47r;
	tmp51i = tmp35i+tmp47i;
	tmp52r = tmp36r+tmp48r;
	tmp52i = tmp36i+tmp48i;
	tmp53r = tmp37r+tmp49r;
	tmp53i = tmp37i+tmp49i;
	tmp54r = tmp34i-tmp46i;
	tmp54i = -(tmp34r-tmp46r);
	tmp55r = tmp35i-tmp47i;
	tmp55i = -(tmp35r-tmp47r);
	tmp56r = tmp36i-tmp48i;
	tmp56i = -(tmp36r-tmp48r);
	tmp57r = tmp37i-tmp49i;
	tmp57i = -(tmp37r-tmp49r);
	tmp58r = tmp18r+tmp50r;
	tmp58i = tmp18i+tmp50i;
	tmp59r = tmp20r+tmp51r;
	tmp59i = tmp20i+tmp51i;
	tmp60r = tmp22r+tmp52r;
	tmp60i = tmp22i+tmp52i;
	tmp61r = tmp24r+tmp53r;
	tmp61i = tmp24i+tmp53i;
	tmp62r = tmp19r-tmp54r;
	tmp62i = tmp19i-tmp54i;
	tmp63r = tmp21r-tmp55r;
	tmp63i = tmp21i-tmp55i;
	tmp64r = tmp23r-tmp56r;
	tmp64i = tmp23i-tmp56i;
	tmp65r = tmp25r-tmp57r;
	tmp65i = tmp25i-tmp57i;
	tmp66r = tmp18r-tmp50r;
	tmp66i = tmp18i-tmp50i;
	tmp67r = tmp20r-tmp51r;
	tmp67i = tmp20i-tmp51i;
	tmp68r = tmp22r-tmp52r;
	tmp68i = tmp22i-tmp52i;
	tmp69r = tmp24r-tmp53r;
	tmp69i = tmp24i-tmp53i;
	tmp70r = tmp19r+tmp54r;
	tmp70i = tmp19i+tmp54i;
	tmp71r = tmp21r+tmp55r;
	tmp71i = tmp21i+tmp55i;
	tmp72r = tmp23r+tmp56r;
	tmp72i = tmp23i+tmp56i;
	tmp73r = tmp25r+tmp57r;
	tmp73i = tmp25i+tmp57i;
	tmp74i = Im(x[16])+Im(x[272]);
	tmp74r = Re(x[16])+Re(x[272]);
	tmp75i = Im(x[16])-Im(x[272]);
	tmp75r = Re(x[16])-Re(x[272]);
	tmp76i = Im(x[144])+Im(x[400]);
	tmp76r = Re(x[144])+Re(x[400]);
	tmp77i = Im(x[144])-Im(x[400]);
	tmp77r = Re(x[144])-Re(x[400]);
	tmp78i = tmp74i+tmp76i;
	tmp78r = tmp74r+tmp76r;
	tmp79i = tmp74i-tmp76i;
	tmp79r = tmp74r-tmp76r;
	tmp80i = tmp75i+tmp77r;
	tmp80r = tmp75r-tmp77i;
	tmp81i = tmp75i-tmp77r;
	tmp81r = tmp75r+tmp77i;
	tmp82i = Im(x[80])+Im(x[336]);
	tmp82r = Re(x[80])+Re(x[336]);
	tmp83i = Im(x[80])-Im(x[336]);
	tmp83r = Re(x[80])-Re(x[336]);
	tmp84i = Im(x[208])+Im(x[464]);
	tmp84r = Re(x[208])+Re(x[464]);
	tmp85i = Im(x[208])-Im(x[464]);
	tmp85r = Re(x[208])-Re(x[464]);
	tmp86i = tmp82i+tmp84i;
	tmp86r = tmp82r+tmp84r;
	tmp87i = tmp82i-tmp84i;
	tmp87r = tmp82r-tmp84r;
	tmp88i = tmp83i+tmp85r;
	tmp88r = tmp83r-tmp85i;
	tmp89i = tmp83i-tmp85r;
	tmp89r = tmp83r+tmp85i;
	tmp90i = C8*(tmp88i+tmp88r);
	tmp90r = C8*(tmp88r-tmp88i);
	tmp91i = -C8*(tmp89i-tmp89r);
	tmp91r = -C8*(tmp89r+tmp89i);
	tmp92i = tmp78i+tmp86i;
	tmp92r = tmp78r+tmp86r;
	tmp93i = tmp78i-tmp86i;
	tmp93r = tmp78r-tmp86r;
	tmp94i = tmp80i+tmp90i;
	tmp94r = tmp80r+tmp90r;
	tmp95i = tmp80i-tmp90i;
	tmp95r = tmp80r-tmp90r;
	tmp96i = tmp79i+tmp87r;
	tmp96r = tmp79r-tmp87i;
	tmp97i = tmp79i-tmp87r;
	tmp97r = tmp79r+tmp87i;
	tmp98i = tmp81i+tmp91i;
	tmp98r = tmp81r+tmp91r;
	tmp99i = tmp81i-tmp91i;
	tmp99r = tmp81r-tmp91r;
	tmp100r = tmp92r;
	tmp100i = tmp92i;
	tmp101r = C2*tmp94r-C3*tmp94i;
	tmp101i = C2*tmp94i+C3*tmp94r;
	tmp102r = C4*tmp96r-C5*tmp96i;
	tmp102i = C4*tmp96i+C5*tmp96r;
	tmp103r = C6*tmp98r-C7*tmp98i;
	tmp103i = C6*tmp98i+C7*tmp98r;
	tmp104r = C8*(tmp93r-tmp93i);
	tmp104i = C8*(tmp93i+tmp93r);
	tmp105r = C7*tmp95r-C6*tmp95i;
	tmp105i = C7*tmp95i+C6*tmp95r;
	tmp106r = C5*tmp97r-C4*tmp97i;
	tmp106i = C5*tmp97i+C4*tmp97r;
	tmp107r = C3*tmp99r-C2*tmp99i;
	tmp107i = C3*tmp99i+C2*tmp99r;
	tmp108i = Im(x[48])+Im(x[304]);
	tmp108r = Re(x[48])+Re(x[304]);
	tmp109i = Im(x[48])-Im(x[304]);
	tmp109r = Re(x[48])-Re(x[304]);
	tmp110i = Im(x[176])+Im(x[432]);
	tmp110r = Re(x[176])+Re(x[432]);
	tmp111i = Im(x[176])-Im(x[432]);
	tmp111r = Re(x[176])-Re(x[432]);
	tmp112i = tmp108i+tmp110i;
	tmp112r = tmp108r+tmp110r;
	tmp113i = tmp108i-tmp110i;
	tmp113r = tmp108r-tmp110r;
	tmp114i = tmp109i+tmp111r;
	tmp114r = tmp109r-tmp111i;
	tmp115i = tmp109i-tmp111r;
	tmp115r = tmp109r+tmp111i;
	tmp116i = Im(x[112])+Im(x[368]);
	tmp116r = Re(x[112])+Re(x[368]);
	tmp117i = Im(x[112])-Im(x[368]);
	tmp117r = Re(x[112])-Re(x[368]);
	tmp118i = Im(x[240])+Im(x[496]);
	tmp118r = Re(x[240])+Re(x[496]);
	tmp119i = Im(x[240])-Im(x[496]);
	tmp119r = Re(x[240])-Re(x[496]);
	tmp120i = tmp116i+tmp118i;
	tmp120r = tmp116r+tmp118r;
	tmp121i = tmp116i-tmp118i;
	tmp121r = tmp116r-tmp118r;
	tmp122i = tmp117i+tmp119r;
	tmp122r = tmp117r-tmp119i;
	tmp123i = tmp117i-tmp119r;
	tmp123r = tmp117r+tmp119i;
	tmp124i = C8*(tmp122i+tmp122r);
	tmp124r = C8*(tmp122r-tmp122i);
	tmp125i = -C8*(tmp123i-tmp123r);
	tmp125r = -C8*(tmp123r+tmp123i);
	tmp126i = tmp112i+tmp120i;
	tmp126r = tmp112r+tmp120r;
	tmp127i = tmp112i-tmp120i;
	tmp127r = tmp112r-tmp120r;
	tmp128i = tmp114i+tmp124i;
	tmp128r = tmp114r+tmp124r;
	tmp129i = tmp114i-tmp124i;
	tmp129r = tmp114r-tmp124r;
	tmp130i = tmp113i+tmp121r;
	tmp130r = tmp113r-tmp121i;
	tmp131i = tmp113i-tmp121r;
	tmp131r = tmp113r+tmp121i;
	tmp132i = tmp115i+tmp125i;
	tmp132r = tmp115r+tmp125r;
	tmp133i = tmp115i-tmp125i;
	tmp133r = tmp115r-tmp125r;
	tmp134r = tmp126r;
	tmp134i = tmp126i;
	tmp135r = C6*tmp128r-C7*tmp128i;
	tmp135i = C6*tmp128i+C7*tmp128r;
	tmp136r = C5*tmp130r-C4*tmp130i;
	tmp136i = C5*tmp130i+C4*tmp130r;
	tmp137r = -C3*tmp132r-C2*tmp132i;
	tmp137i = -C3*tmp132i+C2*tmp132r;
	tmp138r = -C8*(tmp127r+tmp127i);
	tmp138i = -C8*(tmp127i-tmp127r);
	tmp139r = -C2*tmp129r-C3*tmp129i;
	tmp139i = -C2*tmp129i+C3*tmp129r;
	tmp140r = -C4*tmp131r+C5*tmp131i;
	tmp140i = -C4*tmp131i-C5*tmp131r;
	tmp141r = -C7*tmp133r+C6*tmp133i;
	tmp141i = -C7*tmp133i-C6*tmp133r;
	tmp142r = tmp100r+tmp134r;
	tmp142i = tmp100i+tmp134i;
	tmp143r = tmp101r+tmp135r;
	tmp143i = tmp101i+tmp135i;
	tmp144r = tmp102r+tmp136r;
	tmp144i = tmp102i+tmp136i;
	tmp145r = tmp103r+tmp137r;
	tmp145i = tmp103i+tmp137i;
	tmp146r = tmp104r+tmp138r;
	tmp146i = tmp104i+tmp138i;
	tmp147r = tmp105r+tmp139r;
	tmp147i = tmp105i+tmp139i;
	tmp148r = tmp106r+tmp140r;
	tmp148i = tmp106i+tmp140i;
	tmp149r = tmp107r+tmp141r;
	tmp149i = tmp107i+tmp141i;
	tmp150r = tmp100i-tmp134i;
	tmp150i = -(tmp100r-tmp134r);
	tmp151r = tmp101i-tmp135i;
	tmp151i = -(tmp101r-tmp135r);
	tmp152r = tmp102i-tmp136i;
	tmp152i = -(tmp102r-tmp136r);
	tmp153r = tmp103i-tmp137i;
	tmp153i = -(tmp103r-tmp137r);
	tmp154r = tmp104i-tmp138i;
	tmp154i = -(tmp104r-tmp138r);
	tmp155r = tmp105i-tmp139i;
	tmp155i = -(tmp105r-tmp139r);
	tmp156r = tmp106i-tmp140i;
	tmp156i = -(tmp106r-tmp140r);
	tmp157r = tmp107i-tmp141i;
	tmp157i = -(tmp107r-tmp141r);
	tmp158r = tmp58r+tmp142r;
	tmp158i = tmp58i+tmp142i;
	tmp159r = tmp59r+tmp143r;
	tmp159i = tmp59i+tmp143i;
	tmp160r = tmp60r+tmp144r;
	tmp160i = tmp60i+tmp144i;
	tmp161r = tmp61r+tmp145r;
	tmp161i = tmp61i+tmp145i;
	tmp162r = tmp62r+tmp146r;
	tmp162i = tmp62i+tmp146i;
	tmp163r = tmp63r+tmp147r;
	tmp163i = tmp63i+tmp147i;
	tmp164r = tmp64r+tmp148r;
	tmp164i = tmp64i+tmp148i;
	tmp165r = tmp65r+tmp149r;
	tmp165i = tmp65i+tmp149i;
	tmp166r = tmp66r-tmp150r;
	tmp166i = tmp66i-tmp150i;
	tmp167r = tmp67r-tmp151r;
	tmp167i = tmp67i-tmp151i;
	tmp168r = tmp68r-tmp152r;
	tmp168i = tmp68i-tmp152i;
	tmp169r = tmp69r-tmp153r;
	tmp169i = tmp69i-tmp153i;
	tmp170r = tmp70r-tmp154r;
	tmp170i = tmp70i-tmp154i;
	tmp171r = tmp71r-tmp155r;
	tmp171i = tmp71i-tmp155i;
	tmp172r = tmp72r-tmp156r;
	tmp172i = tmp72i-tmp156i;
	tmp173r = tmp73r-tmp157r;
	tmp173i = tmp73i-tmp157i;
	tmp174r = tmp58r-tmp142r;
	tmp174i = tmp58i-tmp142i;
	tmp175r = tmp59r-tmp143r;
	tmp175i = tmp59i-tmp143i;
	tmp176r = tmp60r-tmp144r;
	tmp176i = tmp60i-tmp144i;
	tmp177r = tmp61r-tmp145r;
	tmp177i = tmp61i-tmp145i;
	tmp178r = tmp62r-tmp146r;
	tmp178i = tmp62i-tmp146i;
	tmp179r = tmp63r-tmp147r;
	tmp179i = tmp63i-tmp147i;
	tmp180r = tmp64r-tmp148r;
	tmp180i = tmp64i-tmp148i;
	tmp181r = tmp65r-tmp149r;
	tmp181i = tmp65i-tmp149i;
	tmp182r = tmp66r+tmp150r;
	tmp182i = tmp66i+tmp150i;
	tmp183r = tmp67r+tmp151r;
	tmp183i = tmp67i+tmp151i;
	tmp184r = tmp68r+tmp152r;
	tmp184i = tmp68i+tmp152i;
	tmp185r = tmp69r+tmp153r;
	tmp185i = tmp69i+tmp153i;
	tmp186r = tmp70r+tmp154r;
	tmp186i = tmp70i+tmp154i;
	tmp187r = tmp71r+tmp155r;
	tmp187i = tmp71i+tmp155i;
	tmp188r = tmp72r+tmp156r;
	tmp188i = tmp72i+tmp156i;
	tmp189r = tmp73r+tmp157r;
	tmp189i = tmp73i+tmp157i;
	Re(x[0]) = tmp158r;
	Im(x[0]) = tmp158i;
	Re(x[16]) = tmp159r;
	Im(x[16]) = tmp159i;
	Re(x[32]) = tmp160r;
	Im(x[32]) = tmp160i;
	Re(x[48]) = tmp161r;
	Im(x[48]) = tmp161i;
	Re(x[64]) = tmp162r;
	Im(x[64]) = tmp162i;
	Re(x[80]) = tmp163r;
	Im(x[80]) = tmp163i;
	Re(x[96]) = tmp164r;
	Im(x[96]) = tmp164i;
	Re(x[112]) = tmp165r;
	Im(x[112]) = tmp165i;
	Re(x[128]) = tmp166r;
	Im(x[128]) = tmp166i;
	Re(x[144]) = tmp167r;
	Im(x[144]) = tmp167i;
	Re(x[160]) = tmp168r;
	Im(x[160]) = tmp168i;
	Re(x[176]) = tmp169r;
	Im(x[176]) = tmp169i;
	Re(x[192]) = tmp170r;
	Im(x[192]) = tmp170i;
	Re(x[208]) = tmp171r;
	Im(x[208]) = tmp171i;
	Re(x[224]) = tmp172r;
	Im(x[224]) = tmp172i;
	Re(x[240]) = tmp173r;
	Im(x[240]) = tmp173i;
	Re(x[256]) = tmp174r;
	Im(x[256]) = tmp174i;
	Re(x[272]) = tmp175r;
	Im(x[272]) = tmp175i;
	Re(x[288]) = tmp176r;
	Im(x[288]) = tmp176i;
	Re(x[304]) = tmp177r;
	Im(x[304]) = tmp177i;
	Re(x[320]) = tmp178r;
	Im(x[320]) = tmp178i;
	Re(x[336]) = tmp179r;
	Im(x[336]) = tmp179i;
	Re(x[352]) = tmp180r;
	Im(x[352]) = tmp180i;
	Re(x[368]) = tmp181r;
	Im(x[368]) = tmp181i;
	Re(x[384]) = tmp182r;
	Im(x[384]) = tmp182i;
	Re(x[400]) = tmp183r;
	Im(x[400]) = tmp183i;
	Re(x[416]) = tmp184r;
	Im(x[416]) = tmp184i;
	Re(x[432]) = tmp185r;
	Im(x[432]) = tmp185i;
	Re(x[448]) = tmp186r;
	Im(x[448]) = tmp186i;
	Re(x[464]) = tmp187r;
	Im(x[464]) = tmp187i;
	Re(x[480]) = tmp188r;
	Im(x[480]) = tmp188i;
	Re(x[496]) = tmp189r;
	Im(x[496]) = tmp189i;
}
