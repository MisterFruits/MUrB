/*
 * myIntrinsics.h
 *
 *  Created on: 25 apr. 2014
 *      Author: Adrien CASSAGNE
 */

#ifndef MY_INTRINSICS_H
#define MY_INTRINSICS_H

#ifdef __ARM_NEON__
#include "arm_neon.h"
#else
#include "immintrin.h"
#endif

#ifdef __ARM_NEON__
#ifdef NBODY_DOUBLE /* NEON-128 double */
	/* ces instructions n'existent pas ! */
	#define VECTOR_SIZE 2

#elif defined(NBODY_FLOAT) /* NEON-128 float */
	/* intrinsics NEON headers
	float32x4_t vld1q_f32   (const float32_t *);        // load
	float32x4_t vdupq_n_f32 (float32_t)                 // set 1
	float32x4_t vaddq_f32   (float32x4_t, float32x4_t); // add
	float32x4_t vsubq_f32   (float32x4_t, float32x4_t); // sub
	float32x4_t vmulq_f32   (float32x4_t, float32x4_t); // mul
	void        vst1q_f32   (float32_t *, float32x4_t); // store
	*/
	
	#define vec                    float32x4_t
	#define vec_load(mem_addr)     vld1q_f32   (mem_addr)
	#define vec_store(mem_addr, a) vst1q_f32   (mem_addr, a)
	#define vec_set1(a)            vdupq_n_f32 (a)
	#define vec_add(a, b)          vaddq_f32   (a, b)
	#define vec_sub(a, b)          vsubq_f32   (a, b)
	#define vec_mul(a, b)          vmulq_f32   (a, b)

	#define VECTOR_SIZE 4

#endif
#elif defined(__AVX__)
#ifdef NBODY_DOUBLE /* AVX-256 double */
	/* intrinsics AVX headers
	__m256d _mm256_load_pd  (double const *mem_addr);
	__m256d _mm256_set1_pd  (double a);
	__m256d _mm256_add_pd   (__m256d a, __m256d b);
	__m256d _mm256_sub_pd   (__m256d a, __m256d b);
	__m256d _mm256_mul_pd   (__m256d a, __m256d b);
	__m256d _mm256_div_pd   (__m256d a, __m256d b);
        __m256d _mm256_min_pd (__m256d a, __m256d b)
        __m256d _mm256_sqrt_pd  (__m256d a)
	__m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c);
        __m256d _mm256_permute4x64_pd (__m256d a, const int imm)
	void    _mm256_store_pd (double * mem_addr, __m256d a);
	*/

	#define vec                    __m256d
	#define vec_load(mem_addr)     _mm256_load_pd (mem_addr)
	#define vec_store(mem_addr, a) _mm256_store_pd(mem_addr, a)
	#define vec_set1(a)            _mm256_set1_pd (a)
	#define vec_add(a, b)          _mm256_add_pd  (a, b)
	#define vec_sub(a, b)          _mm256_sub_pd  (a, b)
	#define vec_mul(a, b)          _mm256_mul_pd  (a, b)
	#define vec_div(a, b)          _mm256_div_pd  (a, b)

	#define vec_min(a, b)          _mm256_min_pd  (a, b)
	#define vec_sqrt(a)            _mm256_sqrt_pd (a)

#ifdef __AVX2__
	#define vec_fmadd(a, b, c)     _mm256_fmadd_pd(a, b, c)
	#define vec_permute(a, b)      _mm256_permute4x64_pd(a, b)
#endif

	#define VECTOR_SIZE 4

#elif defined(NBODY_FLOAT) /* AVX-256 float */
	/* intrinsics AVX headers
	__m256 _mm256_load_ps  (float const *mem_addr);
	__m256 _mm256_set1_ps  (float a);
	__m256 _mm256_add_ps   (__m256 a, __m256 b);
	__m256 _mm256_sub_ps   (__m256 a, __m256 b);
	__m256 _mm256_mul_ps   (__m256 a, __m256 b);
	__m256 _mm256_div_ps   (__m256 a, __m256 b);
	__m256 _mm256_min_ps   (__m256 a, __m256 b);
        __m256 _mm256_sqrt_ps  (__m256 a);
	__m256 _mm256_fmadd_ps (__m256 a, __m256 b, __m256 c);
	void   _mm256_store_ps (float * mem_addr, __m256 a);
	*/

	#define vec                    __m256
	#define vec_load(mem_addr)     _mm256_load_ps (mem_addr)
	#define vec_store(mem_addr, a) _mm256_store_ps(mem_addr, a)
	#define vec_set1(a)            _mm256_set1_ps (a)
	#define vec_add(a, b)          _mm256_add_ps  (a, b)
	#define vec_sub(a, b)          _mm256_sub_ps  (a, b)
	#define vec_mul(a, b)          _mm256_mul_ps  (a, b)
	#define vec_div(a, b)          _mm256_div_ps  (a, b)

	#define vec_min(a, b)          _mm256_min_ps  (a, b)
        #define vec_sqrt(a)            _mm256_sqrt_ps (a)

#ifdef __AVX2__
	#define vec_fmadd(a, b, c)     _mm256_fmadd_ps(a, b, c)
#endif

	#define VECTOR_SIZE 8

#endif
#endif
#elif defined(__SSE2__)
#ifdef NBODY_DOUBLE /* SSE2-128 double */
	/* intrinsics SSE2 headers
	__m128d _mm_load_pd  (double const *mem_addr);
	__m128d _mm_set1_pd  (double a);
	__m128d _mm_add_pd   (__m128d a, __m128d b);
	__m128d _mm_sub_pd   (__m128d a, __m128d b);
	__m128d _mm_mul_pd   (__m128d a, __m128d b);
	void    _mm_store_pd (double * mem_addr, __m128d a);
	*/

	#define vec                    __m128d
	#define vec_load(mem_addr)     _mm_load_pd (mem_addr)
	#define vec_store(mem_addr, a) _mm_store_pd(mem_addr, a)
	#define vec_set1(a)            _mm_set1_pd (a)
	#define vec_add(a, b)          _mm_add_pd  (a, b)
	#define vec_sub(a, b)          _mm_sub_pd  (a, b)
	#define vec_mul(a, b)          _mm_mul_pd  (a, b)

	#define VECTOR_SIZE 2

#elif defined(NBODY_FLOAT) /* SSE-128 float */
	/* intrinsics SSE headers
	__m128 _mm_load_ps  (float const *mem_addr);
	__m128 _mm_set1_ps  (float a);
	__m128 _mm_add_ps   (__m128 a, __m128 b);
	__m128 _mm_sub_ps   (__m128 a, __m128 b);
	__m128 _mm_mul_ps   (__m128 a, __m128 b);
	void   _mm_store_ps (float * mem_addr, __m128 a);
	*/

	#define vec                    __m128
	#define vec_load(mem_addr)     _mm_load_ps (mem_addr)
	#define vec_store(mem_addr, a) _mm_store_ps(mem_addr, a)
	#define vec_set1(a)            _mm_set1_ps (a)
	#define vec_add(a, b)          _mm_add_ps  (a, b)
	#define vec_sub(a, b)          _mm_sub_ps  (a, b)
	#define vec_mul(a, b)          _mm_mul_ps  (a, b)

	#define VECTOR_SIZE 4

#endif

#endif /* MY_INTRINSICS_H */
