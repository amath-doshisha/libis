#ifndef IS_ROUNDING_H
#define IS_ROUNDING_H

/* 丸め制御のマクロ */
#ifdef __LINUX__
#include <fpu_control.h>
#define SET_ROUND_NEAR _FPU_SETCW(cwnear) /* 最近点への丸め */
#define SET_ROUND_UP   _FPU_SETCW(cwup)   /* 上への丸め */
#define SET_ROUND_DOWN _FPU_SETCW(cwdown) /* 下への丸め */
extern fpu_control_t cwnear;
extern fpu_control_t cwup;
extern fpu_control_t cwdown;
#elif __MACOSX__
#include <fenv.h>
#define SET_ROUND_NEAR (fesetround(FE_TONEAREST)) /* 最近点への丸め */
#define SET_ROUND_UP   (fesetround(FE_UPWARD))   /* 上への丸め */
#define SET_ROUND_DOWN (fesetround(FE_DOWNWARD)) /* 下への丸め */
#else
// do nothing
#define SET_ROUND_NEAR (0) /* 最近点への丸め */
#define SET_ROUND_UP   (0)   /* 上への丸め */
#define SET_ROUND_DOWN (0) /* 下への丸め */
#endif

#endif
//EOF
