/*
 * myIntrinsicsPlusPlus.h
 *
 *  Created on: 5 oct. 2014
 *      Author: Adrien Cassagne
 */

#ifndef MY_INTRINSICS_PLUS_PLUS_H_
#define MY_INTRINSICS_PLUS_PLUS_H_

#ifdef __ARM_NEON__
#include "arm_neon.h"
#else // we assume that in this case we use an x86 processor with SSE or AVX instructions :-)
#include "immintrin.h"
#endif
#include <cassert>

namespace mipp //My Intrinsics Plus Plus => mipp
{
// ------------------------------------------------------------------------------------------ myIntrinsics vector sizes
// --------------------------------------------------------------------------------------------------------------------

#ifdef __ARM_NEON__ // ----------------------------------------------------------------------------------- ARM NEON-128
	#define REQUIRED_ALIGNEMENT 16
	static const unsigned RequiredAlignement = REQUIRED_ALIGNEMENT;

	#define VECTOR_SIZE_BIT 128
	static const unsigned VectorSizeBit = VECTOR_SIZE_BIT;

	using vec = float32x4_t;

#elif defined(__AVX__) // --------------------------------------------------------------------------------- X86 AVX-256
	#define REQUIRED_ALIGNEMENT 32
	static const unsigned RequiredAlignement = REQUIRED_ALIGNEMENT;

	#define VECTOR_SIZE_BIT 256
	static const unsigned VectorSizeBit = VECTOR_SIZE_BIT;

	using vec = __m256;

#elif defined(__SSE2__) // -------------------------------------------------------------------------------- X86 SSE-128
	#define REQUIRED_ALIGNEMENT 16
	static const unsigned RequiredAlignement = REQUIRED_ALIGNEMENT;

	#define VECTOR_SIZE_BIT 128
	static const unsigned VectorSizeBit = VECTOR_SIZE_BIT;

	using vec = __m128;
#endif

template <typename T>
static inline unsigned int vectorSize() {
	return VECTOR_SIZE_BIT / (8 * sizeof(T));
}

// -------------------------------------------------------------------------------------------- myIntrinsics prototypes
// --------------------------------------------------------------------------------------------------------------------

template <typename T>
inline vec load(const T *mem_addr) {
	std::cerr << "mipp::load is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline void store(T *mem_addr, const vec v) {
	std::cerr << "mipp::store is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec set1(const T val) {
	std::cerr << "mipp::set1 is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec add(const vec v1, const vec v2) {
	std::cerr << "mipp::add is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec sub(const vec v1, const vec v2) {
	std::cerr << "mipp::sub is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec mul(const vec v1, const vec v2) {
	std::cerr << "mipp::mul is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec div(const vec v1, const vec v2) {
	std::cerr << "mipp::div is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec min(const vec v1, const vec v2) {
	std::cerr << "mipp::min is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec max(const vec v1, const vec v2) {
	std::cerr << "mipp::max is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec sqrt(const vec v1) {
	std::cerr << "mipp::sqrt is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec rsqrt(const vec v1) {
	std::cerr << "mipp::rsqrt is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec fmadd(const vec v1, const vec v2, const vec v3) {
	std::cerr << "mipp::fmadd is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

template <typename T>
inline vec rot(const vec v1) {
	std::cerr << "mipp::rot is undefined! Program halting..." << std::endl;
	exit(-1);
	return nullptr;
}

// --------------------------------------------------------------------------------------- myIntrinsics implementations
// --------------------------------------------------------------------------------------------------------------------

#ifdef __ARM_NEON__ // ----------------------------------------------------------------------------------- ARM NEON-128
	// ----------------------------------------------------------------------------------------------------------------

	/* intrinsics NEON-128 headers (float)
	float32x4_t vld1q_f32    (const float32_t *);        // load
	float32x4_t vdupq_n_f32  (float32_t)                 // set 1
	float32x4_t vaddq_f32    (float32x4_t, float32x4_t); // add
	float32x4_t vsubq_f32    (float32x4_t, float32x4_t); // sub
	float32x4_t vmulq_f32    (float32x4_t, float32x4_t); // mul
	float32x4_t vrecpeq_f32  (float32x4_t);              // 1.0 / float32x4_t
	float32x4_t vminq_f32    (float32x4_t, float32x4_t); // min
	float32x4_t vmaxq_f32    (float32x4_t, float32x4_t); // max
	float32x4_t vrsqrteq_f32 (float32x4_t);              // 1.0 / sqrt(float32x4_t)
	void        vst1q_f32    (float32_t *, float32x4_t); // store
	*/

	// ----------------------------------------------------------------------------------------------------------- load
	template <>
	inline vec load<float>(const float *mem_addr) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		return vld1q_f32(mem_addr);
	}

	// ---------------------------------------------------------------------------------------------------------- store
	template <>
	inline void store<float>(float *mem_addr, const vec v) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		vst1q_f32(mem_addr, v);
	}

	// ----------------------------------------------------------------------------------------------------------- set1
	template <>
	inline vec set1<float>(const float val) {
		return vdupq_n_f32(val);
	}

	// ------------------------------------------------------------------------------------------------------------ add
	template <>
	inline vec add<float>(const vec v1, const vec v2) {
		return vaddq_f32(v1, v2);
	}

	// ------------------------------------------------------------------------------------------------------------ sub
	template <>
	inline vec sub<float>(const vec v1, const vec v2) {
		return vsubq_f32(v1, v2);
	}

	// ------------------------------------------------------------------------------------------------------------ mul
	template <>
	inline vec mul<float>(const vec v1, const vec v2) {
		return vmulq_f32(v1, v2);
	}

	// ------------------------------------------------------------------------------------------------------------ div
	template <>
	inline vec div<float>(const vec v1, const vec v2) {
		return mipp::mul<float>(v1, vrecpeq_f32(v2));
	}

	// ------------------------------------------------------------------------------------------------------------ min
	template <>
	inline vec min<float>(const vec v1, const vec v2) {
		return vminq_f32(v1, v2);
	}

	// ------------------------------------------------------------------------------------------------------------ max
	template <>
	inline vec max<float>(const vec v1, const vec v2) {
		return vmaxq_f32(v1, v2);
	}

	// ---------------------------------------------------------------------------------------------------------- rsqrt
	template <>
	inline vec rsqrt<float>(const vec v1) {
		return vrsqrteq_f32(v1);
	}

	// ----------------------------------------------------------------------------------------------------------- sqrt
	template <>
	inline vec sqrt<float>(const vec v1) {
		return vrecpeq_f32(mipp::rsqrt<float>(v1));
	}

	// ---------------------------------------------------------------------------------------------------------- fmadd
	template <>
	inline vec fmadd<float>(const vec v1, const vec v2, const vec v3) {
		return mipp::add<float>(v3, mipp::mul<float>(v1, v2));
	}

	// ------------------------------------------------------------------------------------------------------------ rot
	template <>
	inline vec rot<float>(const vec v1) {
		// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
		return vextq_f32(v1, v1, 1);
	}

#elif defined(__AVX__)// ---------------------------------------------------------------------------------- X86 AVX-256
	// ----------------------------------------------------------------------------------------------------------------

	/* intrinsics AVX headers (float)                               intrinsics AVX headers (double)
	__m256 _mm256_load_ps        (float const *mem_addr);        __m256d _mm256_load_pd        (double const *mem_addr);
	__m256 _mm256_set1_ps        (float a);                      __m256d _mm256_set1_pd        (double a);
	__m256 _mm256_add_ps         (__m256 a, __m256 b);           __m256d _mm256_add_pd         (__m256d a, __m256d b);
	__m256 _mm256_sub_ps         (__m256 a, __m256 b);           __m256d _mm256_sub_pd         (__m256d a, __m256d b);
	__m256 _mm256_mul_ps         (__m256 a, __m256 b);           __m256d _mm256_mul_pd         (__m256d a, __m256d b);
	__m256 _mm256_div_ps         (__m256 a, __m256 b);           __m256d _mm256_div_pd         (__m256d a, __m256d b);
	__m256 _mm256_min_ps         (__m256 a, __m256 b);           __m256d _mm256_min_pd         (__m256d a, __m256d b);
	__m256 _mm256_max_ps         (__m256 a, __m256 b);           __m256d _mm256_max_pd         (__m256d a, __m256d b);
	__m256 _mm256_sqrt_ps        (__m256 a);                     __m256d _mm256_sqrt_pd        (__m256d a);
	__m256 _mm256_rsqrt_ps       (__m256 a);
	__m256 _mm256_fmadd_ps       (__m256 a, __m256 b, __m256 c); __m256d _mm256_fmadd_pd       (__m256d a, __m256d b, __m256d c);
	__m256 _mm256_permute4x64_ps (__m256d a, const int imm)      __m256d _mm256_permute4x64_pd (__m256d a, const int imm)
	  void _mm256_store_ps       (float * mem_addr, __m256 a);     void  _mm256_store_pd       (double * mem_addr, __m256d a);
	*/

	// ----------------------------------------------------------------------------------------------------------- load
	template <>
	inline vec load<float>(const float *mem_addr) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		return _mm256_load_ps(mem_addr);
	}

	template <>
	inline vec load<double>(const double *mem_addr) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		return (__m256) _mm256_load_pd(mem_addr);
	}

	// ---------------------------------------------------------------------------------------------------------- store
	template <>
	inline void store<float>(float *mem_addr, const vec v) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		_mm256_store_ps(mem_addr, v);
	}

	template <>
	inline void store<double>(double *mem_addr, const vec v) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		_mm256_store_pd(mem_addr, (__m256d) v);
	}

	// ----------------------------------------------------------------------------------------------------------- set1
	template <>
	inline vec set1<float>(const float val) {
		return _mm256_set1_ps(val);
	}

	template <>
	inline vec set1<double>(const double val) {
		return (__m256) _mm256_set1_pd(val);
	}

	// ------------------------------------------------------------------------------------------------------------ add
	template <>
	inline vec add<float>(const vec v1, const vec v2) {
		return _mm256_add_ps(v1, v2);
	}

	template <>
	inline vec add<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_add_pd((__m256d) v1, (__m256d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ sub
	template <>
	inline vec sub<float>(const vec v1, const vec v2) {
		return _mm256_sub_ps(v1, v2);
	}

	template <>
	inline vec sub<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_sub_pd((__m256d) v1, (__m256d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ mul
	template <>
	inline vec mul<float>(const vec v1, const vec v2) {
		return _mm256_mul_ps(v1, v2);
	}

	template <>
	inline vec mul<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_mul_pd((__m256d) v1, (__m256d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ div
	template <>
	inline vec div<float>(const vec v1, const vec v2) {
		return _mm256_div_ps(v1, v2);
	}

	template <>
	inline vec div<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_div_pd((__m256d) v1, (__m256d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ min
	template <>
	inline vec min<float>(const vec v1, const vec v2) {
		return _mm256_min_ps(v1, v2);
	}

	template <>
	inline vec min<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_min_pd((__m256d) v1, (__m256d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ max
	template <>
	inline vec max<float>(const vec v1, const vec v2) {
		return _mm256_max_ps(v1, v2);
	}

	template <>
	inline vec max<double>(const vec v1, const vec v2) {
		return (__m256) _mm256_max_pd((__m256d) v1, (__m256d) v2);
	}

	// ----------------------------------------------------------------------------------------------------------- sqrt
	template <>
	inline vec sqrt<float>(const vec v1) {
		return _mm256_sqrt_ps(v1);
	}

	template <>
	inline vec sqrt<double>(const vec v1) {
		return (__m256) _mm256_sqrt_pd((__m256d) v1);
	}

	// ---------------------------------------------------------------------------------------------------------- rsqrt
	template <>
	inline vec rsqrt<float>(const vec v1) {
		return _mm256_rsqrt_ps(v1);
	}

	template <>
	inline vec rsqrt<double>(const vec v1) {
		return mipp::div<double>(mipp::set1<double>(1.0), mipp::sqrt<double>(v1));
	}

	// ---------------------------------------------------------------------------------------------------------- fmadd
	#ifdef __AVX2__
		template <>
		inline vec fmadd<float>(const vec v1, const vec v2, const vec v3) {
			return _mm256_fmadd_ps(v1, v2, v3);
		}

		template <>
		inline vec fmadd<double>(const vec v1, const vec v2, const vec v3) {
			return (__m256) _mm256_fmadd_pd((__m256d) v1, (__m256d) v2, (__m256d) v3);
		}
	#else
		template <>
		inline vec fmadd<float>(const vec v1, const vec v2, const vec v3) {
			return mipp::add<float>(v3, mipp::mul<float>(v1, v2));
		}

		template <>
		inline vec fmadd<double>(const vec v1, const vec v2, const vec v3) {
			return mipp::add<double>(v3, mipp::mul<double>(v1, v2));
		}
	#endif

	// ------------------------------------------------------------------------------------------------------------ rot
	#ifdef __AVX2TODO__
		template <>
		inline vec rot<float>(const vec v1) {
			// make a rotation in:[7, 6, 5, 4, 3, 2 , 1, 0] => out:[0, 7, 6, 5, 4, 3, 2, 1]
			// TODO: this intrinsic does not work well :-(
			return _mm256_permute8x32_ps (v1, _mm256_setr_epi32(0, 7, 6, 5, 4, 3, 2, 1));
		}

		template <>
		inline vec rot<double>(const vec v1) {
			// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
			return (__m256) _mm256_permute4x64_pd ((__m256d) v1, _MM_SHUFFLE(0, 3, 2, 1));
		}
	#else
		template <>
		inline vec rot<float>(const vec v1) {
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
			vec rTmp = _mm256_permute_ps(v1, _MM_SHUFFLE(0, 3, 2, 1));
			return _mm256_blend_ps(rTmp, _mm256_permute2f128_ps(rTmp, rTmp, _MM_SHUFFLE(0, 0, 0, 1)),
			                       _MM_SHUFFLE(2, 0, 2, 0));
		}

		template <>
		inline vec rot<double>(const vec v1) {
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
			vec rTmp = (__m256) _mm256_permute_pd((__m256d) v1, _MM_SHUFFLE(1, 1, 1, 1));
			return (__m256) _mm256_blend_pd((__m256d) rTmp, _mm256_permute2f128_pd((__m256d) rTmp, (__m256d) rTmp, _MM_SHUFFLE(0, 0, 0, 1)),
			                                _MM_SHUFFLE(0, 0, 2, 2));
		}
	#endif

#elif defined(__SSE2__)// --------------------------------------------------------------------------------- X86 SSE-128
	// ----------------------------------------------------------------------------------------------------------------

	/* intrinsics SSE2-128 headers (float)                intrinsics SSE2-128 headers (double)
	__m128  _mm_load_ps  (float const *mem_addr);      __m128d _mm_load_pd  (double const *mem_addr);
	__m128  _mm_set1_ps  (float a);                    __m128d _mm_set1_pd  (double a);
	__m128  _mm_add_ps   (__m128 a, __m128 b);         __m128d _mm_add_pd   (__m128d a, __m128d b);
	__m128  _mm_sub_ps   (__m128 a, __m128 b);         __m128d _mm_sub_pd   (__m128d a, __m128d b);
	__m128  _mm_mul_ps   (__m128 a, __m128 b);         __m128d _mm_mul_pd   (__m128d a, __m128d b);
	__m128d _mm_div_ps   (__m128 a, __m128 b);         __m128d _mm_div_pd   (__m128d a, __m128d b);
	__m128d _mm_min_ps   (__m128 a, __m128 b);         __m128d _mm_min_pd   (__m128d a, __m128d b);
	__m128d _mm_max_ps   (__m128 a, __m128 b);         __m128d _mm_max_pd   (__m128d a, __m128d b);
	__m128  _mm_sqrt_ps  (__m128 a)_m128               __m128d _mm_sqrt_pd  (__m128d a);
	__m128  _mm_rsqrt_ps (__m128 a)
	  void  _mm_store_ps (float * mem_addr, __m128 a);   void  _mm_store_pd (double * mem_addr, __m128d a);

	   intrinsics AVX headers (int)
	__m128i _mm_shuffle_epi32 (__m128i a, int imm);
	*/

	// ----------------------------------------------------------------------------------------------------------- load
	template <>
	inline vec load<float>(const float *mem_addr) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		return _mm_load_ps(mem_addr);
	}

	template <>
	inline vec load<double>(const double *mem_addr) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		return (__m128) _mm_load_pd(mem_addr);
	}

	// ---------------------------------------------------------------------------------------------------------- store
	template <>
	inline void store<float>(float *mem_addr, const vec v) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		_mm_store_ps(mem_addr, v);
	}

	template <>
	inline void store<double>(double *mem_addr, const vec v) {
		assert(alignof(mem_addr) == (mipp::RequiredAlignement));
		_mm_store_pd(mem_addr, (__m128d) v);
	}

	// ----------------------------------------------------------------------------------------------------------- set1
	template <>
	inline vec set1<float>(const float val) {
		return _mm_set1_ps(val);
	}

	template <>
	inline vec set1<double>(const double val) {
		return (__m128) _mm_set1_pd(val);
	}

	// ------------------------------------------------------------------------------------------------------------ add
	template <>
	inline vec add<float>(const vec v1, const vec v2) {
		return _mm_add_ps(v1, v2);
	}

	template <>
	inline vec add<double>(const vec v1, const vec v2) {
		return (__m128) _mm_add_pd((__m128d) v1, (__m128d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ sub
	template <>
	inline vec sub<float>(const vec v1, const vec v2) {
		return _mm_sub_ps(v1, v2);
	}

	template <>
	inline vec sub<double>(const vec v1, const vec v2) {
		return (__m128) _mm_sub_pd((__m128d) v1, (__m128d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ mul
	template <>
	inline vec mul<float>(const vec v1, const vec v2) {
		return _mm_mul_ps(v1, v2);
	}

	template <>
	inline vec mul<double>(const vec v1, const vec v2) {
		return (__m128) _mm_mul_pd((__m128d) v1, (__m128d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ div
	template <>
	inline vec div<float>(const vec v1, const vec v2) {
		return _mm_div_ps(v1, v2);
	}

	template <>
	inline vec div<double>(const vec v1, const vec v2) {
		return (__m128) _mm_div_pd((__m128d) v1, (__m128d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ min
	template <>
	inline vec min<float>(const vec v1, const vec v2) {
		return _mm_min_ps(v1, v2);
	}

	template <>
	inline vec min<double>(const vec v1, const vec v2) {
		return (__m128) _mm_min_pd((__m128d) v1, (__m128d) v2);
	}

	// ------------------------------------------------------------------------------------------------------------ max
	template <>
	inline vec max<float>(const vec v1, const vec v2) {
		return _mm_max_ps(v1, v2);
	}

	template <>
	inline vec max<double>(const vec v1, const vec v2) {
		return (__m128) _mm_max_pd((__m128d) v1, (__m128d) v2);
	}

	// ----------------------------------------------------------------------------------------------------------- sqrt
	template <>
	inline vec sqrt<float>(const vec v1) {
		return _mm_sqrt_ps(v1);
	}

	template <>
	inline vec sqrt<double>(const vec v1) {
		return (__m128) _mm_sqrt_pd((__m128d) v1);
	}

	// ---------------------------------------------------------------------------------------------------------- rsqrt
	template <>
	inline vec rsqrt<float>(const vec v1) {
		return _mm_rsqrt_ps(v1);
	}

	template <>
	inline vec rsqrt<double>(const vec v1) {
		return mipp::div<double>(mipp::set1<double>(1.0), mipp::sqrt<double>(v1));
	}

	// ---------------------------------------------------------------------------------------------------------- fmadd
	template <>
	inline vec fmadd<float>(const vec v1, const vec v2, const vec v3) {
		return mipp::add<float>(v3, mipp::mul<float>(v1, v2));
	}

	template <>
	inline vec fmadd<double>(const vec v1, const vec v2, const vec v3) {
		return mipp::add<double>(v3, mipp::mul<double>(v1, v2));
	}

	// ------------------------------------------------------------------------------------------------------------ rot
	template <>
	inline vec rot<float>(const vec v1) {
		// make a rotation in:[3, 2 , 1, 0] => out:[0, 3, 2, 1]
		return (__m128) _mm_shuffle_epi32 ((__m128i) v1, _MM_SHUFFLE(0, 3, 2, 1));
	}

	template <>
	inline vec rot<double>(const vec v1) {
		// make a rotation in:[1, 0] => out:[0, 1]
		return (__m128) _mm_shuffle_epi32 ((__m128i) v1, _MM_SHUFFLE(1, 0, 3, 2));
	}

#endif
}

#endif /* MY_INTRINSICS_PLUS_PLUS_H_ */
