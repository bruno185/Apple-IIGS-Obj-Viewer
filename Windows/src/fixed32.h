#ifndef FIXED32_H
#define FIXED32_H

#include <stdint.h>
#include <stdlib.h>

typedef int32_t Fixed32;
typedef int64_t Fixed64;

#define FIXED_SHIFT 16
#define FIXED_SCALE (1<<FIXED_SHIFT)
#define FIXED_ONE FIXED_SCALE

static inline Fixed32 float_to_fixed(float f) { return (Fixed32)(f * (float)FIXED_SCALE); }
static inline float fixed_to_float(Fixed32 x) { return ((float)x) / (float)FIXED_SCALE; }

static inline Fixed32 FIXED_ADD(Fixed32 a, Fixed32 b) { return (Fixed32)(a + b); }
static inline Fixed32 FIXED_SUB(Fixed32 a, Fixed32 b) { return (Fixed32)(a - b); }
static inline Fixed32 FIXED_MUL_64(Fixed32 a, Fixed32 b) { return (Fixed32)(( (Fixed64)a * (Fixed64)b ) >> FIXED_SHIFT); }
static inline Fixed32 FIXED_DIV_64(Fixed32 a, Fixed32 b) { return (Fixed32)(( (Fixed64)a << FIXED_SHIFT ) / (Fixed64)b); }

#endif // FIXED32_H
