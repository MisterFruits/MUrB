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
#else // we assume that in this case we use an x86 processor with SSE or AVX instructions :-)
#include "immintrin.h"
#endif

#ifdef __ARM_NEON__
	#ifdef NBODY_DOUBLE /* NEON-128 double */
		/* those instructions does not really exist */
		#define VECTOR_SIZE 2
		#define REQUIRED_ALIGNEMENT 16

	#elif defined(NBODY_FLOAT) /* NEON-128 float */
		/* intrinsics NEON headers
		float32x4_t vld1q_f32    (const float32_t *);        // load
		float32x4_t vdupq_n_f32  (float32_t)                 // set 1
		float32x4_t vaddq_f32    (float32x4_t, float32x4_t); // add
		float32x4_t vsubq_f32    (float32x4_t, float32x4_t); // sub
		float32x4_t vmulq_f32    (float32x4_t, float32x4_t); // mul
		float32x4_t vrecpeq_f32  (float32x4_t);              // 1.0 / float32x4_t
		float32x4_t vminq_f32    (float32x4_t, float32x4_t); // min
		float32x4_t vmaxq_f32    (float32x4_t, float32x4_t); // max
		float32x4_t vrsqrteq_f32 (float32x4_t);              // 1 / sqrt(float32x4_t)
		void        vst1q_f32    (float32_t *, float32x4_t); // store
		*/

		#define vec                    float32x4_t
		#define vec_load(mem_addr)     vld1q_f32   (mem_addr)
		#define vec_store(mem_addr, a) vst1q_f32   (mem_addr, a)
		#define vec_set1(a)            vdupq_n_f32 (a)
		#define vec_add(a, b)          vaddq_f32   (a, b)
		#define vec_sub(a, b)          vsubq_f32   (a, b)
		#define vec_mul(a, b)          vmulq_f32   (a, b)
		#define vec_div(a, b)          vec_mul     (a, vrecpeq_f32(b))
		#define vec_min(a, b)          vminq_f32   (a, b)
		#define vec_max(a, b)          vmaxq_f32   (a, b)
		#define vec_rsqrt(a)           vrsqrteq_f32(a)
		#define vec_sqrt(a)            vrecpeq_f32 (vec_rsqrt(a))
		#define vec_fmadd(a, b, c)     vec_add     (c, vec_mul(a, b))
		#define vec_rot(a)             vextq_f32   (a, a, 1)

		#define VECTOR_SIZE 4
		#define REQUIRED_ALIGNEMENT 16

	#endif
#elif defined(__AVX__)
	#ifdef NBODY_DOUBLE /* AVX-256 double */
		/* intrinsics AVX headers
		__m256d _mm256_load_pd        (double const *mem_addr);
		__m256d _mm256_set1_pd        (double a);
		__m256d _mm256_add_pd         (__m256d a, __m256d b);
		__m256d _mm256_sub_pd         (__m256d a, __m256d b);
		__m256d _mm256_mul_pd         (__m256d a, __m256d b);
		__m256d _mm256_div_pd         (__m256d a, __m256d b);
		__m256d _mm256_min_pd         (__m256d a, __m256d b);
		__m256d _mm256_max_pd         (__m256d a, __m256d b);
		__m256d _mm256_sqrt_pd        (__m256d a);
		__m256d _mm256_fmadd_pd       (__m256d a, __m256d b, __m256d c);
		__m256d _mm256_permute4x64_pd (__m256d a, const int imm)
		  void  _mm256_store_pd       (double * mem_addr, __m256d a);
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
		#define vec_max(a, b)          _mm256_max_pd  (a, b)
		#define vec_sqrt(a)            _mm256_sqrt_pd (a)
		#define vec_rsqrt(a)           vec_div(vec_set1(1), vec_sqrt(a))

		#ifdef __AVX2__
			#define vec_fmadd(a, b, c) _mm256_fmadd_pd(a, b, c)
			// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
			#define vec_rot(a)         _mm256_permute4x64_pd (a, _MM_SHUFFLE(0, 3, 2, 1))
		#else
			#define vec_fmadd(a, b, c) vec_add(c, vec_mul(a, b))
			// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
			//
			//   -> _mm256_permute_pd(a, _MM_SHUFFLE(1, 1, 1, 1)) # rotation per lane of 128 bits
			//          l0      l1           l0      l1
			//      in[3, 2, | 1, 0] => out[2, 3, | 0, 1]
			//
			//   -> _mm256_permute2f128_pd(a, a, _MM_SHUFFLE(0, 0, 0, 1)) # switch lanes
			//          l0     l1             l0     l1
 			//      in[3, 2, | 1, 0] => out[1, 0, | 3, 2]
			//
			//   -> _mm256_blend_pd(a, b, _MM_SHUFFLE(0, 0, 2, 2))
			//      ina[3a, 2a, 1a, 0a] and inb[3b, 2b, 1b, 0b] => out[3b, 2a, 1b, 0a]
			#define vec_rot(a)         _mm256_blend_pd(                                                               \
			                                    _mm256_permute_pd(a, _MM_SHUFFLE(1, 1, 1, 1)),                        \
			                                    _mm256_permute2f128_pd(_mm256_permute_pd(a, _MM_SHUFFLE(1, 1, 1, 1)), \
			                                                           _mm256_permute_pd(a, _MM_SHUFFLE(1, 1, 1, 1)), \
			                                                           _MM_SHUFFLE(0, 0, 0, 1)),                      \
			                                    _MM_SHUFFLE(0, 0, 2, 2))
		#endif
	
		#define VECTOR_SIZE 4
		#define REQUIRED_ALIGNEMENT 32

	#elif defined(NBODY_FLOAT) /* AVX-256 float */
		/* intrinsics AVX headers
		__m256 _mm256_load_ps        (float const *mem_addr);
		__m256 _mm256_set1_ps        (float a);
		__m256 _mm256_add_ps         (__m256 a, __m256 b);
		__m256 _mm256_sub_ps         (__m256 a, __m256 b);
		__m256 _mm256_mul_ps         (__m256 a, __m256 b);
		__m256 _mm256_div_ps         (__m256 a, __m256 b);
		__m256 _mm256_min_ps         (__m256 a, __m256 b);
		__m256 _mm256_max_ps         (__m256 a, __m256 b);
		__m256 _mm256_sqrt_ps        (__m256 a);
		__m256 _mm256_rsqrt_ps       (__m256 a);
		__m256 _mm256_fmadd_ps       (__m256 a, __m256 b, __m256 c);
		__m256 _mm256_permute4x64_ps (__m256d a, const int imm)
		  void _mm256_store_ps       (float * mem_addr, __m256 a);
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
		#define vec_max(a, b)          _mm256_max_ps  (a, b)
		#define vec_sqrt(a)            _mm256_sqrt_ps (a)
		#define vec_rsqrt(a)           _mm256_rsqrt_ps(a)

		#ifdef __AVX2__
			#define vec_fmadd(a, b, c) _mm256_fmadd_ps(a, b, c)
			// make a rotation in:[7, 6, 5, 4, 3, 2 , 1, 0] => out:[0, 7, 6, 5, 4, 3, 2, 1]
			#define vec_rot(a)         _mm256_permute8x32_ps (a, _mm256_setr_epi32(0, 7, 6, 5, 4, 3, 2, 1))
		#else
			#define vec_fmadd(a, b, c) vec_add(c, vec_mul(a, b))
			// make a rotation in:[7, 6, 5, 4, 3, 2 , 1, 0] => out:[0, 7, 6, 5, 4, 3, 2, 1]
			//
			//   -> _mm256_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)) # rotation per lane of 128 bits
			//           lane 0        lane 1             lane 0        lane 1
			//      in[7, 6, 5, 4, | 3, 2, 1, 0] => out[4, 5, 6, 7, | 0, 3, 2, 1]
			//
			//   -> _mm256_permute2f128_ps(a, a, _MM_SHUFFLE(0, 0, 0, 1)) # switch lanes
			//           lane 0        lane 1             lane 0        lane 1
 			//      in[7, 6, 5, 4, | 3, 2, 1, 0] => out[3, 2, 1, 0, | 7, 6, 5, 4]
			//
			//   -> _mm256_blend_ps(a, b, _MM_SHUFFLE(2, 0, 2, 0))
			//      ina[7a, 6a, 5a, 4a, 3a, 2a, 1a, 0a] and inb[7b, 6b, 5b, 4b, 3b, 2b, 1b, 0b] => out[7b, 6a, 5a, 4a, 3b, 2a, 1a, 0a]
			#define vec_rot(a)         _mm256_blend_ps(                                                               \
			                                    _mm256_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)),                        \
			                                    _mm256_permute2f128_ps(_mm256_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)), \
			                                                           _mm256_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)), \
			                                                           _MM_SHUFFLE(0, 0, 0, 1)),                      \
			                                     _MM_SHUFFLE(2, 0, 2, 0))
		#endif

		#define VECTOR_SIZE 8
		#define REQUIRED_ALIGNEMENT 32

	#endif
#elif defined(__SSE2__)
	#ifdef NBODY_DOUBLE /* SSE2-128 double */
		/* intrinsics SSE2 headers
		__m128d _mm_load_pd  (double const *mem_addr);
		__m128d _mm_set1_pd  (double a);
		__m128d _mm_add_pd   (__m128d a, __m128d b);
		__m128d _mm_sub_pd   (__m128d a, __m128d b);
		__m128d _mm_mul_pd   (__m128d a, __m128d b);
		__m128d _mm_div_pd   (__m128d a, __m128d b);
		__m128d _mm_min_pd   (__m128d a, __m128d b);
		__m128d _mm_max_pd   (__m128d a, __m128d b);
		__m128  _mm_sqrt_pd  (__m128d a);
		  void  _mm_store_pd (double * mem_addr, __m128d a);
		*/

		#define vec                    __m128d
		#define vec_load(mem_addr)     _mm_load_pd (mem_addr)
		#define vec_store(mem_addr, a) _mm_store_pd(mem_addr, a)
		#define vec_set1(a)            _mm_set1_pd (a)
		#define vec_add(a, b)          _mm_add_pd  (a, b)
		#define vec_sub(a, b)          _mm_sub_pd  (a, b)
		#define vec_mul(a, b)          _mm_mul_pd  (a, b)
		#define vec_div(a, b)          _mm_div_pd  (a, b)
		#define vec_min(a, b)          _mm_min_pd  (a, b)
		#define vec_max(a, b)          _mm_max_pd  (a, b)
		#define vec_sqrt(a)            _mm_sqrt_pd (a)
		#define vec_rsqrt(a)           vec_div(vec_set1(1), vec_sqrt(a))
		#define vec_fmadd(a, b, c)     vec_add(c, vec_mul(a, b))
		// make a rotation in:[1, 0] => out:[0, 1]
		#define vec_rot(a)             (__m128d) _mm_shuffle_epi32 ((__m128i) a, _MM_SHUFFLE(1, 0, 3, 2))

		#define VECTOR_SIZE 2
		#define REQUIRED_ALIGNEMENT 16

	#elif defined(NBODY_FLOAT) /* SSE-128 float */
		/* intrinsics SSE headers
		__m128  _mm_load_ps  (float const *mem_addr);
		__m128  _mm_set1_ps  (float a);
		__m128  _mm_add_ps   (__m128 a, __m128 b);
		__m128  _mm_sub_ps   (__m128 a, __m128 b);
		__m128  _mm_mul_ps   (__m128 a, __m128 b);
		__m128d _mm_div_ps   (__m128 a, __m128 b);
		__m128d _mm_min_ps   (__m128 a, __m128 b);
		__m128d _mm_max_ps   (__m128 a, __m128 b);
		__m128  _mm_sqrt_ps  (__m128 a)
		__m128  _mm_rsqrt_ps (__m128 a)
		  void  _mm_store_ps (float * mem_addr, __m128 a);
		__m128i _mm_shuffle_epi32 (__m128i a, int imm);
		*/

		#define vec                    __m128
		#define vec_load(mem_addr)     _mm_load_ps (mem_addr)
		#define vec_store(mem_addr, a) _mm_store_ps(mem_addr, a)
		#define vec_set1(a)            _mm_set1_ps (a)
		#define vec_add(a, b)          _mm_add_ps  (a, b)
		#define vec_sub(a, b)          _mm_sub_ps  (a, b)
		#define vec_mul(a, b)          _mm_mul_ps  (a, b)
		#define vec_div(a, b)          _mm_div_ps  (a, b)
		#define vec_min(a, b)          _mm_min_ps  (a, b)
		#define vec_max(a, b)          _mm_max_ps  (a, b)
		#define vec_sqrt(a)            _mm_sqrt_ps (a)
		#define vec_rsqrt(a)           _mm_rsqrt_ps(a)
		#define vec_fmadd(a, b, c)     vec_add(c, vec_mul(a, b))
		// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
		#define vec_rot(a)             (__m128) _mm_shuffle_epi32 ((__m128i) a, _MM_SHUFFLE(0, 3, 2, 1))

		#define VECTOR_SIZE 4
		#define REQUIRED_ALIGNEMENT 16

	#endif
#endif

#endif /* MY_INTRINSICS_H */
