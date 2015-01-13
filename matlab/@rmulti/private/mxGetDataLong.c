#include "mex.h"

long mxGetDataLong(const mxArray *value)
{
    if(mxIsDouble(value)){ return (long)(*(double*)mxGetData(value)); }
    else if(mxIsInt8(value)){ return (long)(*(int8_t*)mxGetData(value)); }
    else if(mxIsInt16(value)){ return (long)(*(int16_t*)mxGetData(value)); }
    else if(mxIsInt32(value)){ return (long)(*(int32_t*)mxGetData(value)); }
    else if(mxIsInt64(value)){ return (long)(*(int64_t*)mxGetData(value)); }
    else{ return 0; }
}
