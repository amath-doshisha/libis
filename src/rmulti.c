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
 @details rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
	  rmulti型の3次元配列の関数に関する定義は@link rmat3.c@endlinkを参照のこと.
 */

////////////////////////////////////////////////////////////

/*
 マクロ
 */

#define RAp(X,Y) ((X)=rallocate_prec(rget_prec(Y)))
#define RF(X) ((X)=rfree(X))

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
rmulti *rfree(rmulti *x)
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
void rclone(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  if(y==x){ return; }
  rround(y,rget_prec(x));
  mpfr_set(y,x,get_round_mode());
}

/**
 @brief   rmulti型の値のコピー y=x.
 @details rmulti型の変数 x の値をrmulti型の変数 y に上書きする.
          ただし，y の浮動小数点数の精度(ビット数)が x の精度(ビット数)より小さい場合は関数@link set_round_mode@endlink で指定された丸めモードで丸めが行われる.
          y の精度(ビット数)の方が大きい場合は不足する仮数の桁には零が埋められる.
 @param[in]  x 初期化済みのrmulti型.
 @param[in]  y 初期化済みのrmulti型.
 @param[out] y 値が上書きされたrmulti型.場合によっては精度(ビット数)が変更される.
 @return       丸めが発生した場合は真,発生しなければ偽の値が返される.
 */
void rcopy(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  if(y==x){ return; }
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
  rfloor(y,x);
  value=req(y,x);
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
  /*
  strings *list=NULL;
  int ret;
  NULL_EXC2(x,value);
  list=strings_split_number(value);
  if(list!=NULL && strings_size(list)>=1 && strings_at(list,0)!=NULL){
    ret=mpfr_set_str(x,strings_at(list,0),10,get_round_mode());
    if(ret){ mpfr_set_nan(x); }
  }else{ mpfr_set_nan(x); }  
  list=strings_del(list);
  */
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
 @brief rmulti型の浮動小数点数を符号なし整数から設定.
 @param[in,out]  x     [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      value 浮動小数点数に変換される値.
 */
void rset_ui(rmulti *x, ulong value)
{
  NULL_EXC1(x);
  mpfr_set_ui(x,value,get_round_mode());
}

/**
 @brief rmulti型の浮動小数点数を符号あり整数から設定.
 @param[in,out]  x     [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 @param[in]      value 浮動小数点数に変換される値.
 */
void rset_si(rmulti *x, long value)
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
  rset_si(x,0);
}

/**
 @brief rmulti型の値を1に設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_one(rmulti *x)
{
  NULL_EXC1(x);
  rset_si(x,1);
}

/**
 @brief rmulti型の値を-1に設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_one_neg(rmulti *x)
{
  NULL_EXC1(x);
  rset_si(x,-1);
}

/**
 @brief rmulti型の値を区間(0,1)の疑似乱数値を設定.
 @param[in,out]  x   [in]初期化済みのrmulti型.[out]値が設定されたrmulti型.
 */
void rset_rand(rmulti *x)
{
  NULL_EXC1(x);
  /* generates a random number on [0,0xffffffff]-interval */
  rset_ui(x,genrand_int32());
  rdiv_ui2(x,x,0xffffffff);
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
 @brief rmulti型を符号なし整数型に変換.
 */
ulong rget_ui(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_get_ui(x,get_round_mode());
}

/**
 @brief rmulti型を符号あり整数型に変換.
 */
long rget_si(rmulti *x)
{
  NULL_EXC1(x);
  return mpfr_get_si(x,get_round_mode());
}


/** @} */

/////////////////////////////////////////////////////////////////

/** @name rmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief rmulti型の符号反転 y=-x
*/
void rneg(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  mpfr_neg(y,x,get_round_mode());
}

/** @brief rmulti型の逆数 y=1/x */
void rinv(rmulti *y, rmulti *x){
  NULL_EXC2(y,x);
  rdiv_d1(y,1,x);
}

/** @brief rmulti型の逆数 y=1/x */
void rinv_d(rmulti *y, double x)
{
  NULL_EXC1(y);
  rset_one(y);
  rdiv_d2(y,y,x);
}

/** @brief rmulti型の整数への丸め y=floor(x) */
void rfloor(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_floor(y,x);
}

/** @brief rmulti型の整数への丸め y=ceil(x) */
void rceil(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_ceil(y,x);
}

/** @brief rmulti型の整数への丸め y=trunc(x) */
void rtrunc(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_trunc(y,x);
}

/**
 @brief rmulti型の絶対値 y=abs(x)
*/
void rabs(rmulti *y, rmulti *x)
{
  NULL_EXC2(x,y);
  mpfr_abs(y,x,get_round_mode());
}

/**
 @brief rmulti型の絶対値の加算 y+=abs(x)
*/
void radd_abs(rmulti *y, rmulti *x)
{
  rmulti *a=NULL;
  NULL_EXC2(y,x);
  RAp(a,y);
  rabs(a,x);   // a=abs(x)
  radd(y,y,a); // y=y+a
  RF(a);
}

/** @brief rmulti型の計算 y=sqrt(x) */
void rsqrt(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sqrt(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=exp(x) */
void rexp(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=2^x */
void rexp2(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp2(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=10^x */
void rexp10(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_exp10(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log(x) */
void rlog(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log_2(x) */
void rlog2(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log2(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=log_10(x) */
void rlog10(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_log10(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sin(x) */
void rsin(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sin(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cos(x) */
void rcos(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cos(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=tan(x) */
void rtan(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_tan(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sec(x) */
void rsec(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sec(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cosec(x) */
void rcsc(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_csc(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cot(x) */
void rcot(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cot(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arcsin(x) */
void rasin(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_asin(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arccos(x) */
void racos(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_acos(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arctan(x) */
void ratan(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_atan(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=sinh(x) */
void rsinh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_sinh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=cosh(x) */
void rcosh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_cosh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=tanh(x) */
void rtanh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_tanh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arcsinh(x) */
void rasinh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_asinh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arccosh(x) */
void racosh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_acosh(y,x,get_round_mode());
}

/** @brief rmulti型の計算 y=arctanh(x) */
void ratanh(rmulti *y, rmulti *x)
{
  NULL_EXC2(y,x);
  mpfr_atanh(y,x,get_round_mode());
}

/** @} */

/////////////////////////////////////////////////////////////////

/** @name rmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief rmulti型の足し算 z=x+y
*/
void radd(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_add(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の足し算 z=x+y
*/
int radd_exact(rmulti *z, rmulti *x, rmulti *y)
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
 @brief rmulti型の足し算 z=x+y
*/
void radd_d(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_add_d(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の足し算 z=x+y
*/
void radd_ui(rmulti *z, rmulti *x, ulong y)
{
  NULL_EXC2(z,x);
  mpfr_add_ui(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の足し算 z=x+y
*/
void radd_si(rmulti *z, rmulti *x, long y)
{
  NULL_EXC2(z,x);
  mpfr_add_si(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_sub(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
int rsub_exact(rmulti *z, rmulti *x, rmulti *y)
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
 @brief rmulti型の引き算 z=x-y
*/
void rsub_d1(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_d_sub(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub_d2(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_sub_d(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub_ui1(rmulti *z, ulong x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_ui_sub(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub_ui2(rmulti *z, rmulti *x, ulong y)
{
  NULL_EXC2(z,x);
  mpfr_sub_ui(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub_si1(rmulti *z, long x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_si_sub(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の引き算 z=x-y
*/
void rsub_si2(rmulti *z, rmulti *x, long y)
{
  NULL_EXC2(z,x);
  mpfr_sub_si(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の掛け算 z=x*y
*/
void rmul(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_mul(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の掛け算 z=x*y
*/
int rmul_exact(rmulti *z, rmulti *x, rmulti *y)
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
 @brief rmulti型の掛け算 z=x*y
*/
void rmul_d(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_mul_d(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の掛け算 z=x*y
*/
void rmul_ui(rmulti *z, rmulti *x, ulong y)
{
  NULL_EXC2(z,x);
  mpfr_mul_ui(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の掛け算 z=x*y
*/
void rmul_si(rmulti *z, rmulti *x, long y)
{
  NULL_EXC2(z,x);
  mpfr_mul_si(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y
*/
void radd_mul(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rmul(a,x,y); // a=x*y
  radd(z,z,a); // z=z+x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y
*/
int radd_mul_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  e+=rmul_exact(a,x,y); // a=x*y
  e+=radd_exact(z,z,a); // z=z+x*y
  RF(a);
  return e;
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y（作業領域モード）
*/
void radd_mul_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws)
{
  NULL_EXC3(z,x,y);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in radd_mul_ws().\n",n_rws); }
  rmul(rws[0],x,y); // rws[0]=x*y
  radd(z,z,rws[0]); // z=z+x*y
}

/**
 @brief rmulti型の掛け算の加算 z+=x*y
*/
void radd_mul_d(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rmul_d(a,x,y); // a=x*y
  radd(z,z,a);   // z=z+x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rmul(a,x,y); // a=x*y
  rsub(z,z,a); // z=z-x*y
  RF(a);
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
int rsub_mul_exact(rmulti *z, rmulti *x, rmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  e+=rmul_exact(a,x,y); // a=x*y
  e+=rsub_exact(z,z,a); // z=z-x*y
  RF(a);
  return e;
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws)
{
  NULL_EXC3(z,x,y);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in rsub_mul_ws().\n",n_rws); }
  rmul(rws[0],x,y); // a=x*y
  rsub(z,z,rws[0]); // z=z-x*y
}

/**
 @brief rmulti型の掛け算の減算 z-=x*y
*/
void rsub_mul_d(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rmul_d(a,x,y); // a=x*y
  rsub(z,z,a);   // z=z-x*y
  RF(a);
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_div(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
int rdiv_rouding_check(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  return abs(mpfr_div(z,x,y,get_round_mode()));
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_d1(rmulti *z, double x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_d_div(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_d2(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  mpfr_div_d(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_ui1(rmulti *z, ulong x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_ui_div(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_ui2(rmulti *z, rmulti *x, ulong y)
{
  NULL_EXC2(z,x);
  mpfr_div_ui(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_si1(rmulti *z, long x, rmulti *y)
{
  NULL_EXC2(z,y);
  mpfr_si_div(z,x,y,get_round_mode());
}

/** @brief rmulti型の割り算 z=x/y */
void rdiv_si2(rmulti *z, rmulti *x, long y)
{
  NULL_EXC2(z,x);
  mpfr_div_si(z,x,y,get_round_mode());
}

/**
 @brief rmulti型の2値の絶対値の和 z=abs(x)+abs(y)
*/
void radd_abs_abs(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rset_zero(a);
  radd_abs(a,x);
  radd_abs(a,y);
  rcopy(z,a);
  RF(a);
}

/**
 @brief rmulti型の差の絶対値 z=abs(x-y)
*/
void rabs_sub(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  rsub(z,x,y);
  rabs(z,z);
}

/**
 @brief rmulti型の差の絶対値 z=abs(x-y)
*/
void rabs_sub_d(rmulti *z, rmulti *x, double y)
{
  NULL_EXC2(z,x);
  rsub_d2(z,x,y);
  rabs(z,z);
}

/**
 @brief rmulti型の商の分母が絶対値 z=x/abs(y)
*/
void rdiv_abs(rmulti *z, rmulti *x, rmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC3(z,x,y);
  RAp(a,z);
  rabs(a,y);   // a=abs(y)
  rdiv(z,x,a); // z=x/abs(y)
  RF(a);
}

/** @brief rmulti型のべき乗 z=x^y */
void rpow(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_pow(z,x,y,get_round_mode());
}

/**
 @brief rmulti型のべき乗 z=x^y
*/
void rpow_d2(rmulti *z, rmulti *x, double y)
{
  rmulti *a=NULL;
  NULL_EXC2(z,x);
  RAp(a,z);
  rset_d(a,y);
  rpow(z,x,a); // z=x^y
  RF(a);
}

/**
 @brief rmulti型のべき乗 z=x^y
*/
void rpow_ui(rmulti *z, rmulti *x, ulong y)
{
  NULL_EXC2(z,x);
  mpfr_pow_ui(z,x,y,get_round_mode());
}

/** @brief rmulti型のべき乗 z=x^y */
void rpow_si(rmulti *z, rmulti *x, long y)
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
void ratan2(rmulti *z, rmulti *x, rmulti *y)
{
  NULL_EXC3(z,x,y);
  mpfr_atan2(z,x,y,get_round_mode());
}

/**
 @brief 2つのrmulti型の大きい方 上丸め z=max2(x,y)
*/
void rmax2_up(rmulti *z, rmulti *x, rmulti *y)
{
  if(rgt(x,y)){ mpfr_set(z,x,MPFR_RNDU); }
  else        { mpfr_set(z,y,MPFR_RNDU); }
}

/**
 @brief 2つのrmulti型の小さい方 下丸め z=min2(x,y)
*/
void rmin2_down(rmulti *z, rmulti *x, rmulti *y)
{
  if(rlt(x,y)){ mpfr_set(z,x,MPFR_RNDD); }
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
  rabs(a,x);      // a=abs(x)
  rlog10(a,a);    // a=log10(abs(x))
  rsub_d2(a,a,y); // a=log10(abs(x))-y
  rfloor(a,a);    // a=floor(log10(abs(x))-y)
  rexp10(z,a);    // z=10^floor(log10(abs(x))-y)
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
  rabs(a,x);      // a=abs(x)
  rlog2(a,a);     // a=log2(abs(x))
  rsub_d2(a,a,y); // a=log2(abs(x))-y
  rfloor(a,a);    // a=floor(log2(abs(x))-y)
  rexp2(z,a);     // z=10^floor(log2(abs(x))-y)
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
  radd_abs(a,x0);
  radd_abs(a,x1);
  radd_abs(a,x2);
  rcopy(z,a);
  RF(a);
}

/** @} */

//////////////////////////////////////////////////////////////////////

/** @name rmulti型の値の比較に関する関数 */
/** @{ */

/** @brief rmulti型の値の比較 x<=>y */
int rcmp(rmulti *x, rmulti *y)  { NULL_EXC2(x,y); return mpfr_cmp(x,y);}
/** @brief rmulti型の値の比較 x<=>y */
int rcmp_ui(rmulti *x, ulong y) { NULL_EXC1(x); return mpfr_cmp_ui(x,y);}
/** @brief rmulti型の値の比較 x<=>y */
int rcmp_si(rmulti *x, long y)  { NULL_EXC1(x); return mpfr_cmp_si(x,y);}
/** @brief rmulti型の値の比較 x<=>y */
int rcmp_d(rmulti *x, double y) { NULL_EXC1(x); return mpfr_cmp_d(x,y);}
/** @brief rmulti型の値の比較 x==y */
int req(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rcmp(x,y)==0; }
/** @brief rmulti型の値の比較 x==y */
int req_d(rmulti *x, double y)  { NULL_EXC1(x); return ris_number(x) && rcmp_d(x,y)==0; }
/** @brief rmulti型の値の比較 x==y */
int req_ui(rmulti *x, ulong y)  { NULL_EXC1(x); return ris_number(x) && rcmp_ui(x,y)==0; }
/** @brief rmulti型の値の比較 x==y */
int req_si(rmulti *x, long y)   { NULL_EXC1(x); return ris_number(x) && rcmp_si(x,y)==0; }
/** @brief rmulti型の値の比較 x!=y */
int rne(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return !(ris_number(x) && ris_number(y) && rcmp(x,y)==0); }
/** @brief rmulti型の値の比較 x!=y */
int rne_d(rmulti *x, double y)  { NULL_EXC1(x); return !(ris_number(x) && rcmp_d(x,y)==0); }
/** @brief rmulti型の値の比較 x!=y */
int rne_ui(rmulti *x, ulong y)  { NULL_EXC1(x); return !(ris_number(x) && rcmp_ui(x,y)==0); }
/** @brief rmulti型の値の比較 x!=y */
int rne_si(rmulti *x, long y)   { NULL_EXC1(x); return !(ris_number(x) && rcmp_si(x,y)==0); }
/** @brief rmulti型の値の比較 x>y */
int rgt(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rcmp(x,y)>0; }
/** @brief rmulti型の値の比較 x>y */
int rgt_d1(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_d(y,x)<0;  }
/** @brief rmulti型の値の比較 x>y */
int rgt_d2(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && rcmp_d(x,y)>0;  }
/** @brief rmulti型の値の比較 x>y */
int rgt_ui1(ulong x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_ui(y,x)<0;  }
/** @brief rmulti型の値の比較 x>y */
int rgt_ui2(rmulti *x, ulong y) { NULL_EXC1(x); return ris_number(x) && rcmp_ui(x,y)>0;  }
/** @brief rmulti型の値の比較 x>y */
int rgt_si1(long x, rmulti *y)  { NULL_EXC1(y); return ris_number(y) && rcmp_si(y,x)<0;  }
/** @brief rmulti型の値の比較 x>y */
int rgt_si2(rmulti *x, long y)  { NULL_EXC1(x); return ris_number(x) && rcmp_si(x,y)>0;  }
/** @brief rmulti型の値の比較 x>=y */
int rge(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rcmp(x,y)>=0;   }
/** @brief rmulti型の値の比較 x>=y */
int rge_d1(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_d(y,x)<=0; }
/** @brief rmulti型の値の比較 x>=y */
int rge_d2(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && rcmp_d(x,y)>=0; }
/** @brief rmulti型の値の比較 x>=y */
int rge_ui1(ulong x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_ui(y,x)<=0; }
/** @brief rmulti型の値の比較 x>=y */
int rge_ui2(rmulti *x, ulong y) { NULL_EXC1(x); return ris_number(x) && rcmp_ui(x,y)>=0; }
/** @brief rmulti型の値の比較 x>=y */
int rge_si1(long x, rmulti *y)  { NULL_EXC1(y); return ris_number(y) && rcmp_si(y,x)<=0; }
/** @brief rmulti型の値の比較 x>=y */
int rge_si2(rmulti *x, long y)  { NULL_EXC1(x); return ris_number(x) && rcmp_si(x,y)>=0; }
/** @brief rmulti型の値の比較 x<y */
int rlt(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rcmp(x,y)<0; }
/** @brief rmulti型の値の比較 x<y */
int rlt_d1(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_d(y,x)>0;  }
/** @brief rmulti型の値の比較 x<y */
int rlt_d2(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && rcmp_d(x,y)<0;  }
/** @brief rmulti型の値の比較 x<y */
int rlt_ui1(ulong x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_ui(y,x)>0; }
/** @brief rmulti型の値の比較 x<y */
int rlt_ui2(rmulti *x, ulong y) { NULL_EXC1(x); return ris_number(x) && rcmp_ui(x,y)<0; }
/** @brief rmulti型の値の比較 x<y */
int rlt_si1(long x, rmulti *y)  { NULL_EXC1(y); return ris_number(y) && rcmp_si(y,x)>0; }
/** @brief rmulti型の値の比較 x<y */
int rlt_si2(rmulti *x, long y)  { NULL_EXC1(x); return ris_number(x) && rcmp_si(x,y)<0; }
/** @brief rmulti型の値の比較 x<=y */
int rle(rmulti *x, rmulti *y)   { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rcmp(x,y)<=0;   }
/** @brief rmulti型の値の比較 x<=y */
int rle_d1(double x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_d(y,x)>=0; }
/** @brief rmulti型の値の比較 x<=y */
int rle_d2(rmulti *x, double y) { NULL_EXC1(x); return ris_number(x) && rcmp_d(x,y)<=0; }
/** @brief rmulti型の値の比較 x<=y */
int rle_ui1(ulong x, rmulti *y) { NULL_EXC1(y); return ris_number(y) && rcmp_ui(y,x)>=0; }
/** @brief rmulti型の値の比較 x<=y */
int rle_ui2(rmulti *x, ulong y) { NULL_EXC1(x); return ris_number(x) && rcmp_ui(x,y)<=0; }
/** @brief rmulti型の値の比較 x<=y */
int rle_si1(long x, rmulti *y)  { NULL_EXC1(y); return ris_number(y) && rcmp_si(y,x)>=0; }
/** @brief rmulti型の値の比較 x<=y */
int rle_si2(rmulti *x, long y)  { NULL_EXC1(x); return ris_number(x) && rcmp_si(x,y)<=0; }

/** @brief rmulti型の値の比較 abs(x)<=>abs(y) */
int rabs_cmp(rmulti *x, rmulti *y)
{
  int value;  
  rmulti *ax=NULL,*ay=NULL;
  NULL_EXC2(x,y);
  RAp(ax,x); RAp(ay,y);
  rabs(ax,x); rabs(ay,y);
  value=rcmp(ax,ay);
  RF(ax); RF(ay);
  return value;
}

/** @brief rmulti型の値の比較 abs(x)==abs(y) */
int rabs_eq(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rabs_cmp(x,y)==0; }
/** @brief rmulti型の値の比較 abs(x)>abs(y) */
int rabs_gt(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rabs_cmp(x,y)>0;  }
/** @brief rmulti型の値の比較 abs(x)>=abs(y) */
int rabs_ge(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rabs_cmp(x,y)>=0; }
/** @brief rmulti型の値の比較 abs(x)<abs(y) */
int rabs_lt(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rabs_cmp(x,y)<0;  }
/** @brief rmulti型の値の比較 abs(x)<=abs(y) */
int rabs_le(rmulti *x, rmulti *y) { NULL_EXC2(x,y); return ris_number(x) && ris_number(y) && rabs_cmp(x,y)<=0; }

/** @brief rmulti型の値の比較 abs(x-y)<=>z */
int rdist_cmp(rmulti *x, rmulti *y, rmulti *z)
{
  rmulti *a=NULL;
  int value;
  NULL_EXC3(x,y,z);
  a=rallocate_prec(MAX2(rget_prec(x),rget_prec(y)));
  rabs_sub(a,x,y);
  value=rcmp(a,z);
  RF(a);
  return value;
}

/** @brief rmulti型の値の比較 abs(x-y)==z */
int rdist_eq(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && rdist_cmp(x,y,z)==0; }
/** @brief rmulti型の値の比較 abs(x-y)>=z */
int rdist_ge(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && rdist_cmp(x,y,z)>=0; }
/** @brief rmulti型の値の比較 abs(x-y)>z */
int rdist_gt(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && rdist_cmp(x,y,z)> 0; }
/** @brief rmulti型の値の比較 abs(x-y)<=z */
int rdist_le(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && rdist_cmp(x,y,z)<=0; }
/** @brief rmulti型の値の比較 abs(x-y)<z */
int rdist_lt(rmulti *x, rmulti *y, rmulti *z){ NULL_EXC3(x,y,z); return ris_number(x) && ris_number(y) && ris_number(z) && rdist_cmp(x,y,z)< 0; }

/** @brief rmulti型の値の比較 abs(x-y)<=>z */
int rdist_cmp_d(rmulti *x, double y, double z)
{
  rmulti *a=NULL;
  int value;
  NULL_EXC1(x);
  RAp(a,x);
  rabs_sub_d(a,x,y);
  value=rcmp_d(a,z);
  RF(a);
  return value;
}

/** @brief rmulti型の値の比較 abs(x-y)==z */
int rdist_eq_d(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && rdist_cmp_d(x,y,z)==0; }
/** @brief rmulti型の値の比較 abs(x-y)>=z */
int rdist_ge_d(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && rdist_cmp_d(x,y,z)>=0; }
/** @brief rmulti型の値の比較 abs(x-y)>z */
int rdist_gt_d(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && rdist_cmp_d(x,y,z)> 0; }
/** @brief rmulti型の値の比較 abs(x-y)<=z */
int rdist_le_d(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && rdist_cmp_d(x,y,z)<=0; }
/** @brief rmulti型の値の比較 abs(x-y)<z */
int rdist_lt_d(rmulti *x, double y, double z){ NULL_EXC1(x); return ris_number(x) && rdist_cmp_d(x,y,z)< 0; }

/** @} */

//////////////////////////////////////////////////////////////////

//EOF
