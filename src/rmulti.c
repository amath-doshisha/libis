#include<stdio.h>
#include<stdlib.h>
#include<mpfr.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_irmulti.h"
#include"is_cmulti.h"
#include"is_icmulti.h"
#include"is_strings.h"
#include"mt19937ar.h"

/**
 @file  rmulti.c
 @brief 多倍長精度実数型rmultiに関する関数の定義
 */

////////////////////////////////////////////////////////////

/*
 マクロ
 */

#define RAp(X,Y) ((X)=rallocate_prec(rget_prec(Y)))
#define RF(X) ((X)=rmfree(X))

////////////////////////////////////////////////////////////

/*
 大域変数
 */

/**
 @brief rmulti型の丸めモードの状態を表す変数(ユーザが直接アクセスしてはならない).
 */
mpfr_rnd_t __default_round_mode=MPFR_RNDN;


////////////////////////////////////////////////////////////

/** @name rmulti型の精度の設定に関する関数 */
/** @{ */

/**
 @brief rmulti型のメモリ割り当ての際に使用される精度のデフォルト値の設定.
 @details 設定された浮動小数点数の仮数の桁数(ビット数)のデフォルト値は@link rallocate @endlinkでrmulti型を初期化する際に使用される.
 @param[in] prec 浮動小数点数の仮数の桁数(ビット数).
 */
void set_default_prec(int prec)
{
  mpfr_set_default_prec(prec);
}

/**
 @brief rmulti型のメモリ割り当ての際に使用される精度のデフォルト値の取得.
 @return prec 浮動小数点数の仮数の桁数(ビット数).
 */
int get_default_prec()
{
  return (int)mpfr_get_default_prec();
}

/** @} */

/////////////////////////////////////////////////////

/** @name rmulti型の丸めに関する関数 */
/** @{ */

/**
 @brief     rmulti型の丸めモードの状態の設定.
 @details   丸めモードの状態 mode を整数の定数で設定する.
            近接値への丸めはMPFR_RNDN,
            上への丸めはMPFR_RNDU,
            下への丸めはMPFR_RNDD,
            零への丸めはMPFR_RNDZを指定する.
            これらの定数はMPFR(http://www.mpfr.org/)で定義されたマクロである.
 @param[in] mode rmulti型の丸めモードの状態を指定する.
 */
void set_round_mode(mpfr_rnd_t mode)
{
  __default_round_mode=mode;
}

/**
 @brief   rmulti型の丸めモードの状態の取得.
 @details 丸めモードの状態を表す整数の定数が返される.
          近接値への丸めはMPFR_RNDN,
          上への丸めはMPFR_RNDU,
          下への丸めはMPFR_RNDD,
          零への丸めはMPFR_RNDZが返される.
          これらの定数はMPFR(http://www.mpfr.org/)で定義されたマクロである.
 @return  rmulti型の丸めモードの状態を返す.
 */
mpfr_rnd_t get_round_mode()
{
  return __default_round_mode;
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name rmulti型の初期化に関する関数 */
/** @{ */

/**
 @brief  rmulti型の新規生成.
 @details 
   rmulti型の構造体のメモリを割り当てとそのメンバのメモリ割当などの初期化作業が行われる.
   初期化後にはrmulti型へのポインタが返される.
   浮動小数点数の値はNaNに設定される.
   浮動小数点数の仮数の精度(ビット数)は関数@link set_default_prec@endlinkにより設定された値が使用される.
   初期化後の仮数の精度(ビット数)の再初期化は関数@link rround@endlinkにより行う.
   初期化されたrmulti型の終了処理は関数@link rfree@endlinkにより行う.
 @return 初期化されたrmulti型へのポインタを返す．
 */
rmulti *rallocate(void)
{
  rmulti *x=NULL;
  x=(rmulti*)malloc(sizeof(rmulti));
  mpfr_init(x);
  return x;
}

/**
 @brief  rmulti型の精度を指定しての新規生成.
 @details 
   rmulti型の構造体のメモリを割り当てとそのメンバのメモリ割当などの初期化作業が行われる.
   初期化後にはrmulti型へのポインタが返される.
   浮動小数点数の値はNaNに設定される.
   浮動小数点数の仮数の精度(ビット数)は指定した prec の値が使用される.
   初期化後の仮数の精度(ビット数)の再初期化は関数@link rround@endlinkにより行う.
   初期化されたrmulti型の終了処理は関数@link rfree@endlinkにより行う.
 @param[in]  prec 浮動小数点数の仮数の桁数(ビット数)を指定する.
 @return 初期化されたrmulti型へのポインタを返す．
 */
rmulti *rallocate_prec(int prec)
{
  rmulti *x=NULL;
  x=(rmulti*)malloc(sizeof(rmulti));
  mpfr_init2(x,prec);
  return x;
}

/**
 @brief  rmulti型の複製の生成.
 @details 
   rmulti型の構造体のメモリを割り当てとそのメンバのメモリ割当などの初期化作業が行われる.
   初期化後にはrmulti型へのポインタが返される.
   浮動小数点数の値および浮動小数点数の仮数の精度(ビット数)は，変数yと同一のものが設定される.
   初期化後の仮数の精度(ビット数)の再初期化は関数@link rround@endlinkにより行う.
   初期化されたrmulti型の終了処理は関数@link rfree@endlinkにより行う.
 @param[in] y 複製されるrmulti型の指定.
 @return 初期化されたrmulti型へのポインタを返す．
 */
rmulti *rallocate_clone(rmulti *y)
{
  rmulti *x=NULL;
  RAp(x,y);
  mpfr_set(x,y,get_round_mode());
  return x;
}

/**
 @brief   rmulti型の終了処理.
 @details 初期化されたrmulti型の構造体およびそのメンバなどのメモリの開放を行う.
          rmulti型の構造体のメモリ割り当ておよび初期化は関数@link rallocate@endlink, @link rallocate_prec@endlink, @link rallocate_clone@endlink, @link rbin_load@endlinkにより行われる.
 @param[in]  x 初期化されたrmulti型ポインタ.
 @param[out] x メモリが開放される。これ以降はアクセス不可.
 @return       NULLポインタを常に返す.
 */
rmulti *rmfree(rmulti *x)
{
  if(x==NULL){ return NULL; }
  mpfr_clear(x);
  free(x);
  x=NULL;
  return x;
}


/**
 @brief   初期化済みのrmulti型の浮動小数点数の精度(ビット数)を変更し再初期化.
 @details 初期化済みのrmulti型の変数 x の浮動小数点数の精度(ビット数)を prec に変更し再初期化を行う.
          再初期化の際には浮動小数点数の値は保持される.
          ただし，変更後のビット数が元のビット数より小さい場合は関数@link set_round_mode@endlink で指定された丸めモードで丸めが行われる.
          変更後のビット数の方が大きい場合は不足する仮数の桁には零が埋められる.
          返値として丸めが発生したどうかの真偽値を返す.
 @param[in]  x    初期化済みのrmulti型.
 @param[in]  prec 望みの浮動小数点数の精度(ビット数).
 @param[out] x    浮動小数点数の精度(ビット数)が変更され再初期化されたrmulti型.
 @return          なし．
 */
void rround(rmulti *x, int prec)
{
  mpfr_prec_round(x,prec,get_round_mode());
}

/**
 @brief   rmulti型の値を複成.
 @details 書き込み先のrmulti型の変数 y の浮動小数点数の精度(ビット数)をrmulti型の変数 x の精度(ビット数)と同じに変更する.
          その後，変数 x の値を変数 y に上書きする.
          そのため，丸めは発生せず完全に同じ値がコピーされる.
          返値として丸めが発生したどうかの真偽値を返すが，この関数の場合は常に偽値を返す.
 @param[in]  x 初期化済みのrmulti型.
 @param[in,out]  y    [in]初期化済みのrmulti型.[out]変数 x と同じ精度(ビット数)が変更され，x と同じ値が上書きされる.
 */
void rclone_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  if(y==x){ return; }
  rround(y,rget_prec(x));
  mpfr_set(y,x,get_round_mode());
}

/**
 @brief   rmulti型の値の交換.
 @details rmulti型の変数 x とrmulti型の変数 y の値を交換する.
          その際に，丸めは発生しない.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が変数 y の値に変更される.
 @param[in,out]  y   [in]初期化済みのrmulti型.[out]値が変数 x の値に変更される.
 */
void rswap(rmulti *x, rmulti *y)
{
  NULL_EXC2(x,y);
  mpfr_swap(x,y);
}

/**
 @brief rmulti型のマシンイプシロンの生成.
*/ 
rmulti *rmepsilon(int prec)
{
  rmulti *eps=NULL;
  eps=rallocate_prec(prec);  // allocate
  rset_one(eps);             // eps=1
  rdiv_2exp(eps,eps,prec);   // eps=1/2^n
  return eps;
}

/** @} */

/////////////////////////////////////////////////////

/** @name rmulti型のメンバ変数に関する関数 */
/** @{ */

/**
 @brief       rmulti型の浮動小数点数の精度(ビット数)を取得.
 @param[in] x 初期化済みのrmulti型.
 @return      取得されたビット数の値.
 */
int rget_prec(rmulti *x)
{
  NULL_EXC1(x);
  return (int)mpfr_get_prec(x);
}

/**
 @brief rmulti型の浮動小数点数の指数部の取得.
 @return      取得された指数部の値.
 */
int rget_exp(rmulti *x)
{
  NULL_EXC1(x);
  return (int)mpfr_get_exp(x);
}

/**
 @brief       rmulti型の浮動小数点数の符号部の取得.
 @param[in] x 初期化済みのrmulti型.
 @return      取得された符号部の値.
 */
int rget_sgn(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_sgn(x);
}

/**
 @brief       rmulti型の浮動小数点数xの符号部の取得.x>=0のとき1，x<0のとき-1を返す．
 @param[in] x 初期化済みのrmulti型.
 @return      取得された符号部の値.
 */
int rget_sgn_not0(rmulti *x)
{
  NULL_EXC1(x);
  return (rget_sgn(x)>=0)?1:-1;
}

/**
 @brief       rmulti型の浮動小数点数の仮数部が格納されている配列のサイズの取得.
 @param[in] x 初期化済みのrmulti型.
 @return      計算されたサイズ.
 */
int rget_size(rmulti *x)
{
  int k=0;
  NULL_EXC1(x);
  k=(int)(x->_mpfr_prec/GMP_NUMB_BITS);
  if(x->_mpfr_prec % GMP_NUMB_BITS){ k++; }
  return k;
}

/**
 @brief       rmulti型の浮動小数点数の仮数部の0でない桁数の長さ[bits]の取得
 @param[in] x 初期化済みのrmulti型.
 @return      桁数の長さ[bits]
 */
int rget_length(rmulti *x)
{
  mp_limb_t mask;
  int i,j,n,done,bits;
  if(!ris_number(x) || ris_zero(x)){ return 0; }
  n=rget_size(x);
  bits=n*GMP_NUMB_BITS;
  done=0;
  for(i=0; done==0 && i<n; i++){
    mask=1;
    for(j=0; done==0 && j<GMP_NUMB_BITS; j++){
      if((x->_mpfr_d[i]) & mask){ done=1; }
      else                      { bits--; }
      mask=mask<<1;
    }
  }
  return bits;
}

/** @} */



/////////////////////////////////////////////////////

/** @name rmulti型の数種の判定に関する関数 */
/** @{ */

/**
 @brief       rmulti型がNaNであるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_nan(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_nan_p(x);
}

/**
 @brief       rmulti型がInfであるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_inf(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_inf_p(x);
}

/**
 @brief       rmulti型が数であるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_number(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_number_p(x);
}

/**
 @brief       rmulti型が零であるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_zero(rmulti *x)
{
  return mpfr_zero_p(x);
}

/**
 @brief       rmulti型が正であるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_positive(rmulti *x)
{
  NULL_EXC1(x);
  return ris_number(x) && rget_sgn(x)>0;
}

/**
 @brief       rmulti型が負であるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
 */
int ris_negative(rmulti *x)
{
  NULL_EXC1(x);
  return ris_number(x) && rget_sgn(x)<0;
}

/**
 @brief       rmulti型が整数であるかの判定.
 @param[in] x 初期化済みのrmulti型.
 @return      判定結果の真偽値.
*/
int ris_integer(rmulti *x)
{
  int value;
  rmulti *y=NULL;
  NULL_EXC1(x);
  RAp(y,x);
  rfloor_r(y,x);
  value=eq_rr(y,x);
  RF(y);
  return value;
}

/** @} */
/////////////////////////////////////////////////////////////

/** @name rmulti型の入出力に関する関数 */
/** @{ */

/**
 @brief            rmulti型の値の表示
 @param[in] x      初期化済みのrmulti型.
 @param[in] name   出力される変数名の文字列.
 @param[in] f      出力される浮動小数点数のフォーマット."e", "f", "g"のどれかを指定.
 @param[in] digits 出力される浮動小数点数の桁数を指定.
 */
void rprint(rmulti *x, const char *name, const char *f, int digits)
{
  char format[128];
  NULL_EXC2(x,f);
  sprintf(format,"%%.%dR%s\n",digits,f);
  if(name!=NULL){ printf("%s",name); }
  if(x==NULL){ printf("NULL\n"); return; }
  mpfr_printf(format,x);
}

/**
 @brief       rmulti型のメンバ変数の値の表示
 @param[in] x 初期化済みのrmulti型.
 */
void rput(rmulti *x)
{
  int i,k;
  NULL_EXC1(x);
  k=(int)(mpfr_custom_get_size(rget_prec(x))*8/mp_bits_per_limb);
  printf("prec=%4d, sgn=%+d, exp=%4d, digits=",rget_prec(x),rget_sgn(x),rget_exp(x));
  for(i=k-1; i>=0; i--){ printf("%016lx ",x->_mpfr_d[i]); }
  printf("\n");
}

/**
 @brief       rmulti型のメンバ変数の値の表示
 @param[in] x 初期化済みのrmulti型.
 @param[in] name 出力される変数名の文字列.
 */
void rput_info(rmulti *x, char *name)
{
  int i,k;
  NULL_EXC2(x,name);
  printf("GMP_NUMB_BITS=%d bits\n",GMP_NUMB_BITS);
  mpfr_printf("%s=%.50Re\n",name,x);
  printf("%s->prec=%ld bits (%ld bytes)\n",name,x->_mpfr_prec,sizeof(mpfr_prec_t));
  printf("%s->sign=%d (%ld bytes)\n",name,x->_mpfr_sign,sizeof(mpfr_sign_t));
  printf("%s->exp=%ld (%ld bytes)\n",name,x->_mpfr_exp,sizeof(mpfr_exp_t));
  printf("%s->size=%d\n",name,(k=rget_size(x)));
  for(i=0; i<k; i++){
    printf("%s->d[%d]=%16lx (%ld bytes)\n",name,i,x->_mpfr_d[i],sizeof(mp_limb_t));
  }
}

/**
 @brief         rmulti型のファイルへの保存
 @param[in] x   初期化済みのrmulti型.
 @param[in] fid 書き込み先のファイル識別子.ただし，オープン済みであること.
 */
void rbin_save(rmulti *x, FILE *fid)
{
  NULL_EXC2(x,fid);
  fwrite(&(x->_mpfr_prec),sizeof(mpfr_prec_t),1,fid);
  fwrite(&(x->_mpfr_sign),sizeof(mpfr_sign_t),1,fid);
  fwrite(&(x->_mpfr_exp),sizeof(mpfr_exp_t),1,fid);
  fwrite(x->_mpfr_d,sizeof(mp_limb_t),rget_size(x),fid);
}

/**
 @brief         rmulti型のファイルからの読み込み
 @param[in] fid 読み込むファイル識別子.ただし，オープン済みであること.
 @return        初期化され値が設定されたrmulti型.
 */
rmulti *rbin_load(FILE *fid)
{
  mpfr_prec_t prec=0;
  int n;
  size_t k;
  rmulti *x=NULL;
  NULL_EXC1(fid);
  // read prec
  k=fread(&prec,sizeof(mpfr_prec_t),1,fid);
  if(k!=1){ ERROR_AT; printf("Can't read x->_mpfr_prec."); exit(0); }
  // allocate
  x=rallocate_prec((int)prec);
  // read sign
  k=fread(&(x->_mpfr_sign),sizeof(mpfr_sign_t),1,fid);
  if(k!=1){ ERROR_AT; printf("Can't read x->_mpfr_sign."); exit(0); }
  // read exp
  k=fread(&(x->_mpfr_exp),sizeof(mpfr_exp_t),1,fid);
  if(k!=1){ ERROR_AT; printf("Can't read x->_mpfr_exp."); exit(0); }
  // read data
  n=rget_size(x);
  k=fread(x->_mpfr_d,sizeof(mp_limb_t),n,fid);
  if(k!=(size_t)n){ ERROR_AT; printf("Can't read x->_mpfr_d."); exit(0); }
  return x;
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name rmulti型の値の設定に関する関数 */
/** @{ */

/**
 @brief rmulti型の浮動小数点数を文字列から設定.
 @param[in,out]  x     [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      s     浮動小数点数に変換される文字列.
 */
void rset_s(rmulti *x, char *s)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=rget_prec(x);
  a0=callocate_prec(prec);
  a1=callocate_prec(prec);  
  icset_s(a0,a1,s);  
  irmid(x,C_R(a0),C_R(a1));
}

// y=x
void rset_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  if(y==x){ return; }
  mpfr_set(y,x,get_round_mode());
}


/**
 @brief rmulti型の浮動小数点数を倍精度浮動小数点数から設定.
 @param[in,out]  x     [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      value 浮動小数点数に変換される値.
 */
void rset_d(rmulti *x, double value)
{
  NULL_EXC1(x);
  mpfr_set_d(x,value,get_round_mode());
}

/**
 @brief rmulti型の浮動小数点数を符号あり整数から設定.
 @param[in,out]  x     [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      value 浮動小数点数に変換される値.
 */
void rset_i(rmulti *x, int value)
{
  NULL_EXC1(x);
  mpfr_set_si(x,value,get_round_mode());
}

/**
 @brief rmulti型の値をInfに設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      sgn +Infの場合はsgn=1,-Infの場合はsgn=-1,Infの場合はsgn=0を指定.
 */
void rset_inf(rmulti *x, int sgn)
{
  NULL_EXC1(x);
  mpfr_set_inf(x,sgn);
}

/**
 @brief rmulti型の値をNaNに設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_nan(rmulti *x)
{
  NULL_EXC1(x);
  mpfr_set_nan(x);
}

/**
 @brief rmulti型の値を零に設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_zero(rmulti *x)
{
  NULL_EXC1(x);
  rset_i(x,0);
}

/**
 @brief rmulti型の値を1に設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_one(rmulti *x)
{
  NULL_EXC1(x);
  rset_i(x,1);
}

/**
 @brief rmulti型の値を-1に設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_one_neg(rmulti *x)
{
  NULL_EXC1(x);
  rset_i(x,-1);
}

/**
 @brief rmulti型の値を区間(0,1)の疑似乱数値を設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_rand(rmulti *x)
{
  NULL_EXC1(x);
  /* generates a random number on [0,0xffffffff]-interval */
  rset_d(x,((double)genrand_int32())/0xffffffff);
}

/** @} */

///////////////////////////////////////////////////////////

/** @name rmulti型の型変換に関する関数 */
/** @{ */

/**
 @brief rmulti型を倍精度型に変換.
 */
double rget_d(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_get_d(x,get_round_mode());
}

/**
 @brief rmulti型を符号あり整数型に変換.
 */
int rget_i(rmulti *x)
{
  NULL_EXC1(x);
  return (int)mpfr_get_si(x,get_round_mode());
}


/** @} */

/////////////////////////////////////////////////////////////////

/** @name rmulti型の関数 y=f(x) */
/** @{ */

/**
 @brief y=-x
 */
void rneg_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  mpfr_neg(y,x,get_round_mode());
}

/**
 @brief y=-x
 */
void rneg_d(rmulti *y, double x)
{
  NULL_EXC1(y);
  mpfr_set_d(y,-x,get_round_mode());
}

/**
 @brief y=1/x
 */
void rinv_r(rmulti *y, rmulti *x){
  NULL_EXC2(y,x);
  rdiv_dr(y,1,x);
}

/**
 @brief y=1/x
 */
void rinv_d(rmulti *y, double x)
{
  NULL_EXC1(y);
  rset_one(y);
  rdiv_rd(y,y,x);
}

/** @brief rmulti型の整数への丸め y=floor(x) */
void rfloor_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_floor(y,x);
}

/** @brief rmulti型の整数への丸め y=ceil(x) */
void rceil_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_ceil(y,x);
}

/** @brief rmulti型の整数への丸め y=trunc(x) */
void rtrunc_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_trunc(y,x);
}

/**
 @brief rmulti型の絶対値 y=abs(x)
*/
void rabs_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  mpfr_abs(y,x,get_round_mode());
}

/**
 @brief rmulti型の絶対値の加算 y+=abs(x)
*/
void radd_abs_r(rmulti *y, rmulti *x)
{
  rmulti *a=NULL;
  NULL_EXC2(y,x);
  RAp(a,y);
  rabs_r(a,x);   // a=abs(x)
  radd_rr(y,y,a); // y=y+a
  RF(a);
}

/** @brief rmulti型の計算 y=sqrt(x) */
void rsqrt_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sqrt(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=exp(x) */
void rexp_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=2^x */
void rexp2_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp2(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=10^x */
void rexp10_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp10(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log(x) */
void rlog_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log_2(x) */
void rlog2_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log2(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log_10(x) */
void rlog10_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log10(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sin(x) */
void rsin_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sin(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cos(x) */
void rcos_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cos(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=tan(x) */
void rtan_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_tan(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sec(x) */
void rsec_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sec(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cosec(x) */
void rcsc_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_csc(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cot(x) */
void rcot_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cot(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arcsin(x) */
void rasin_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_asin(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arccos(x) */
void racos_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_acos(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arctan(x) */
void ratan_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_atan(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sinh(x) */
void rsinh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sinh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cosh(x) */
void rcosh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cosh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=tanh(x) */
void rtanh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_tanh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arcsinh(x) */
void rasinh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_asinh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arccosh(x) */
void racosh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_acosh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arctanh(x) */
void ratanh_r(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_atanh(y,x,get_round_mode());
}

/** @} */

/////////////////////////////////////////////////////////////////

/** @name rmulti型の関数 y=f(x,y) */
/** @{ */

/**
 @brief z=x+y
 */
void radd_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_add(z,x,y,get_round_mode());
}

/**
 @brief z=x+y
 */
int radd_rr_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int nx0,nx1,ny0,ny1,nz0,nz1,lz,e;
  NULL_EXC3(z,x,y);
  nx0=rget_exp(x);
  ny0=rget_exp(y);
  nx1=rget_exp(x)-rget_length(x)+1;
  ny1=rget_exp(y)-rget_length(y)+1;
  nz0=MAX2(nx0,ny0);
  nz1=MIN2(nx1,ny1);
  lz=nz0-nz1+2;
  rround(z,((lz/GMP_NUMB_BITS)+1)*GMP_NUMB_BITS);  
  if((e=abs(mpfr_add(z,x,y,get_round_mode())))!=0){ ERROR_AT; exit(-1); }  
  return e;
}

/**
 @brief z=x+y
 */
void radd_rd(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_add_d(z,x,y,get_round_mode());
}

/**
 @brief z=x+y
 */
void radd_dr(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_add_d(z,y,x,get_round_mode());
}

/**
 @brief z=x-y
 */
void rsub_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_sub(z,x,y,get_round_mode());
}

/**
 @brief z=x-y
 */
int rsub_rr_exact(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  int nx0,nx1,ny0,ny1,nz0,nz1,lz,e;
  NULL_EXC3(z,x,y);
  nx0=rget_exp(x);
  ny0=rget_exp(y);
  nx1=rget_exp(x)-rget_length(x)+1;
  ny1=rget_exp(y)-rget_length(y)+1;
  nz0=MAX2(nx0,ny0);
  nz1=MIN2(nx1,ny1);
  lz=nz0-nz1+2;
  rround(z,((lz/GMP_NUMB_BITS)+1)*GMP_NUMB_BITS);  
  if((e=abs(mpfr_sub(z,x,y,get_round_mode())))!=0){ ERROR_AT; exit(-1); }  
  return e;
}

/**
 @brief z=x-y
 */
void rsub_dr(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_d_sub(z,x,y,get_round_mode());
}

/**
 @brief z=x-y
 */
void rsub_rd(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_sub_d(z,x,y,get_round_mode());
}

/**
 @brief z=x*y
 */
void rmul_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_mul(z,x,y,get_round_mode());
}

/**
 @brief z=x*y
 */
int rmul_rr_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int lx,ly,lz,e=0;
  NULL_EXC3(z,x,y);
  lx=rget_length(x);
  ly=rget_length(y);
  lz=lx+ly;
  rround(z,((lz/GMP_NUMB_BITS)+1)*GMP_NUMB_BITS);
  if((e=abs(mpfr_mul(z,x,y,get_round_mode())))!=0){ ERROR_AT; exit(-1); }
  return e;
}

/**
 @brief z=x*y
 */
void rmul_rd(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_mul_d(z,x,y,get_round_mode());
}

/**
 @brief z=x*y
 */
void rmul_dr(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_mul_d(z,y,x,get_round_mode());
}

/**
 @brief z=x/y
 */
void rdiv_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_div(z,x,y,get_round_mode());
}

/**
 @brief z=x/y
 */
int rdiv_rr_rouding_check(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  return abs(mpfr_div(z,x,y,get_round_mode()));
}

/**
 @brief z=x/y
 */
void rdiv_rd(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_div_d(z,x,y,get_round_mode());
}

/**
 @brief z=x/y
 */
void rdiv_dr(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_d_div(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の2値の絶対値の和 z=abs(x)+abs(y)
 */
void radd_abs_abs_rr(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rset_zero(a);
  radd_abs_r(a,x);
  radd_abs_r(a,y);
  rset_r(z,a);
  RF(a);
}

/**
 @brief rmulti型の差の絶対値 z=abs(x-y)
 */
void rabs_sub_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  rsub_rr(z,x,y);
  rabs_r(z,z);
}

/**
 @brief rmulti型の差の絶対値 z=abs(x-y)
 */
void rabs_sub_rd(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  rsub_rd(z,x,y);
  rabs_r(z,z);
}

/** @brief rmulti型のべき乗 z=x^y */
void rpow_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_pow(z,x,y,get_round_mode());
}

/**
 @brief rmulti型のべき乗 z=x^y
 */
void rpow_rd(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rset_d(a,y);
  rpow_rr(z,x,a); // z=x^y
  RF(a);
}

/** @brief rmulti型のべき乗 z=x^y */
void rpow_r(rmulti *z, rmulti *x, int y)
{
  NULL_EXC2(z,x);
  mpfr_pow_si(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の指数部の足し算 y=x*2^n
 */
void rmul_2exp(rmulti *y, rmulti *x, int n)
{
  NULL_EXC2(x,y);
  if(n>=0){ mpfr_mul_2exp(y,x,n,get_round_mode()); }
  else    { mpfr_div_2exp(y,x,-n,get_round_mode()); }
}

/**
 @brief rmulti型の指数部の引き算 y=x/2^n
 */
void rdiv_2exp(rmulti *y, rmulti *x, int n)
{
  NULL_EXC2(x,y);
  if(n>=0){ mpfr_div_2exp(y,x,n,get_round_mode()); }
  else    { mpfr_mul_2exp(y,x,-n,get_round_mode()); }
}

/** @brief rmulti型の計算 y=arctan(x/y) */
void ratan2_rr(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_atan2(z,x,y,get_round_mode());
}

/**
 @brief 2つのrmulti型の大きい方 上丸め z=max2(x,y)
 */
void rmax2_up(rmulti *z, rmulti *x, rmulti *y)
{
  if(gt_rr(x,y)){ mpfr_set(z,x,MPFR_RNDU); }
  else        { mpfr_set(z,y,MPFR_RNDU); }
}

/**
 @brief 2つのrmulti型の小さい方 下丸め z=min2(x,y)
 */
void rmin2_down(rmulti *z, rmulti *x, rmulti *y)
{
  if(lt_rr(x,y)){ mpfr_set(z,x,MPFR_RNDD); }
  else        { mpfr_set(z,y,MPFR_RNDD); }
}

/**
 @brief rmulti型の指数部で評価 z=10^(floor(log10(abs(x))-y))
 */
void rexp10_floor_log10_abs_sub(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rabs_r(a,x);      // a=abs(x)
  rlog10_r(a,a);    // a=log10(abs(x))
  rsub_rd(a,a,y); // a=log10(abs(x))-y
  rfloor_r(a,a);    // a=floor(log10(abs(x))-y)
  rexp10_r(z,a);    // z=10^floor(log10(abs(x))-y)
  RF(a);
}

/**
 @brief rmulti型の指数部で評価 z=2^(floor(log2(abs(x))-y))
 */
void rexp2_floor_log2_abs_sub(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rabs_r(a,x);      // a=abs(x)
  rlog2_r(a,a);     // a=log2(abs(x))
  rsub_rd(a,a,y); // a=log2(abs(x))-y
  rfloor_r(a,a);    // a=floor(log2(abs(x))-y)
  rexp2_r(z,a);     // z=10^floor(log2(abs(x))-y)
  RF(a);
}


/** @} */

//////////////////////////////////////////////////////////////////////

/** @name rmulti型の関数 z=z+f(x,y) */
/** @{ */

/**
 @brief z+=x*y
*/
void radd_mul_rr(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rmul_rr(a,x,y); // a=x*y
  radd_rr(z,z,a); // z=z+x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y
*/
int radd_mul_rr_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  e+=rmul_rr_exact(a,x,y); // a=x*y
  e+=radd_rr_exact(z,z,a); // z=z+x*y
  RF(a);
  return e;
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y（作業領域モード）
*/
void radd_mul_rr_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws)
{
  NULL_EXC3(z,x,y);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in radd_mul_ws().\n",n_rws); }
  rmul_rr(rws[0],x,y); // rws[0]=x*y
  radd_rr(z,z,rws[0]); // z=z+x*y
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y
*/
void radd_mul_rd(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rmul_rd(a,x,y); // a=x*y
  radd_rr(z,z,a);   // z=z+x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul_rr(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rmul_rr(a,x,y); // a=x*y
  rsub_rr(z,z,a); // z=z-x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
int rsub_mul_rr_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  e+=rmul_rr_exact(a,x,y); // a=x*y
  e+=rsub_rr_exact(z,z,a); // z=z-x*y
  RF(a);
  return e;
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul_rr_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws)
{
  NULL_EXC3(z,x,y);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in rsub_mul_ws().\n",n_rws); }
  rmul_rr(rws[0],x,y); // a=x*y
  rsub_rr(z,z,rws[0]); // z=z-x*y
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul_rd(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rmul_rd(a,x,y); // a=x*y
  rsub_rr(z,z,a);   // z=z-x*y
  RF(a);
}


/** @} */

/////////////////////////////////////////////////////////////////

/** @name rmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief rmulti型の3値の絶対値の和 z=abs(x0)+abs(x1)+abs(x2)
*/
void radd_abs_abs_abs(rmulti *z, rmulti *x0, rmulti *x1, rmulti *x2)
{
  rmulti *a=NULL;
  NULL_EXC4(z,x0,x1,x2);
  RAp(a,z);
  rset_zero(a);
  radd_abs_r(a,x0);
  radd_abs_r(a,x1);
  radd_abs_r(a,x2);
  rset_r(z,a);
  RF(a);
}

/** @} */

//////////////////////////////////////////////////////////////////////

/** @name rmulti型の値の比較に関する関数 */
/** @{ */

/** @brief rmulti型の値の比較 x<=>y */
int cmp_rr(rmulti *x, rmulti *y)  { NULL_EXC2(x,y); return mpfr_cmp(x,y);}
/** @brief rmulti型の値の比較 x<=>y */
int cmp_rd(rmulti *x, double y) { NULL_EXC1(x); return mpfr_cmp_d(x,y);}
/** @brief rmulti型の値の比較 x==y */
int eq_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_rr(x,y)==0; }
/** @brief rmulti型の値の比較 x==y */
int eq_rd(rmulti *x, double y)  { NULL_EXC1(x); return ris_number(x) && cmp_rd(x,y)==0; }
/** @brief rmulti型の値の比較 x!=y */
int ne_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return !(ris_number(x) && ris_number(y) && cmp_rr(x,y)==0); }
/** @brief rmulti型の値の比較 x!=y */
int ne_rd(rmulti *x, double y)  { NULL_EXC1(x); return !(ris_number(x) && cmp_rd(x,y)==0); }
/** @brief rmulti型の値の比較 x>y */
int gt_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_rr(x,y)>0; }
/** @brief rmulti型の値の比較 x>y */
int gt_dr(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && cmp_rd(y,x)<0;  }
/** @brief rmulti型の値の比較 x>y */
int gt_rd(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && cmp_rd(x,y)>0;  }
/** @brief rmulti型の値の比較 x>=y */
int ge_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_rr(x,y)>=0;   }
/** @brief rmulti型の値の比較 x>=y */
int ge_dr(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && cmp_rd(y,x)<=0; }
/** @brief rmulti型の値の比較 x>=y */
int ge_rd(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && cmp_rd(x,y)>=0; }
/** @brief rmulti型の値の比較 x<y */
int lt_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_rr(x,y)<0; }
/** @brief rmulti型の値の比較 x<y */
int lt_dr(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && cmp_rd(y,x)>0;  }
/** @brief rmulti型の値の比較 x<y */
int lt_rd(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && cmp_rd(x,y)<0;  }
/** @brief rmulti型の値の比較 x<=y */
int le_rr(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_rr(x,y)<=0;   }
/** @brief rmulti型の値の比較 x<=y */
int le_dr(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && cmp_rd(y,x)>=0; }
/** @brief rmulti型の値の比較 x<=y */
int le_rd(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && cmp_rd(x,y)<=0; }

/** @brief rmulti型の値の比較 abs(x)<=>abs(y) */
int cmp_abs_rr(rmulti *x, rmulti *y)
{
  int value;  
  rmulti *ax=NULL,*ay=NULL;
  NULL_EXC2(x,y);
  RAp(ax,x); RAp(ay,y);
  rabs_r(ax,x); rabs_r(ay,y);
  value=cmp_rr(ax,ay);
  RF(ax); RF(ay);
  return value;
}

/** @brief rmulti型の値の比較 abs(x)==abs(y) */
int eq_abs_rr(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_abs_rr(x,y)==0; }
/** @brief rmulti型の値の比較 abs(x)>abs(y) */
int gt_abs_rr(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_abs_rr(x,y)>0;  }
/** @brief rmulti型の値の比較 abs(x)>=abs(y) */
int ge_abs_rr(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_abs_rr(x,y)>=0; }
/** @brief rmulti型の値の比較 abs(x)<abs(y) */
int lt_abs_rr(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_abs_rr(x,y)<0;  }
/** @brief rmulti型の値の比較 abs(x)<=abs(y) */
int le_abs_rr(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && cmp_abs_rr(x,y)<=0; }

/** @brief rmulti型の値の比較 abs(x-y)<=>z */
int cmp_dist_rrr(rmulti *x, rmulti *y, rmulti *z)
{
  rmulti *a=NULL;
  int value;
  NULL_EXC3(x,y,z);
  a=rallocate_prec(MAX2(rget_prec(x),rget_prec(y)));
  rabs_sub_rr(a,x,y);
  value=cmp_rr(a,z);
  RF(a);
  return value;
}

/** @brief rmulti型の値の比較 abs(x-y)==z */
int eq_dist_rrr(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && cmp_dist_rrr(x,y,z)==0; }
/** @brief rmulti型の値の比較 abs(x-y)>=z */
int ge_dist_rrr(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && cmp_dist_rrr(x,y,z)>=0; }
/** @brief rmulti型の値の比較 abs(x-y)>z */
int gt_dist_rrr(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && cmp_dist_rrr(x,y,z)> 0; }
/** @brief rmulti型の値の比較 abs(x-y)<=z */
int le_dist_rrr(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && cmp_dist_rrr(x,y,z)<=0; }
/** @brief rmulti型の値の比較 abs(x-y)<z */
int lt_dist_rrr(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && cmp_dist_rrr(x,y,z)< 0; }

/** @brief rmulti型の値の比較 abs(x-y)<=>z */
int cmp_dist_rdd(rmulti *x, double y, double z)
{
  rmulti *a=NULL;
  int value;
  NULL_EXC1(x);
  RAp(a,x);
  rabs_sub_rd(a,x,y);
  value=cmp_rd(a,z);
  RF(a);
  return value;
}

/** @brief rmulti型の値の比較 abs(x-y)==z */
int eq_dist_rdd(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && cmp_dist_rdd(x,y,z)==0; }
/** @brief rmulti型の値の比較 abs(x-y)>=z */
int ge_dist_rdd(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && cmp_dist_rdd(x,y,z)>=0; }
/** @brief rmulti型の値の比較 abs(x-y)>z */
int gt_dist_rdd(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && cmp_dist_rdd(x,y,z)> 0; }
/** @brief rmulti型の値の比較 abs(x-y)<=z */
int le_dist_rdd(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && cmp_dist_rdd(x,y,z)<=0; }
/** @brief rmulti型の値の比較 abs(x-y)<z */
int lt_dist_rdd(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && cmp_dist_rdd(x,y,z)< 0; }

/** @} */

//////////////////////////////////////////////////////////////////

//EOF
