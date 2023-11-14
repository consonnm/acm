#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
#define fr(i, l, r) for (int i = (l); i <= (r); ++i)
#define rf(i, l, r) for (int i = (l); i >= (r); --i)
#define fo(i, l, r) for (int i = (l); i < (r); ++i)
#define foredge(i, u, v) for (int i = fir[u], v; v = to[i], i; i = nxt[i])
#define filein(file) freopen(file ".in", "r", stdin)
#define fileout(file) freopen(file ".out", "w", stdout)

namespace POLY
{
/*#pragma GCC optimize("-Ofast","-funroll-all-loops","-ffast-math")
#pragma GCC optimize("-fno-math-errno")
#pragma GCC optimize("-funsafe-math-optimizations")
#pragma GCC optimize("-freciprocal-math")
#pragma GCC optimize("-fno-trapping-math")
#pragma GCC optimize("-ffinite-math-only")
#pragma GCC optimize("-fno-stack-protector")
#pragma GCC target ("avx2","sse4.2","fma")*/
#define pb push_back
#define mp make_pair
	typedef long long ll;
#define fi first
#define se second
#define SS 524288			// max length of DFT
#define PS (SS * 20 + 1000) // pool for temp arrays
#define PS2 (SS * 7 + 1000) // another pool, see ocmul
#define FS 666666			// length of fac, rfac
	const int MOD = 998244353;
	ll qp(ll a, ll b)
	{
		ll ans = 1;
		a %= MOD;
		while (b)
		{
			if (b & 1)
				ans = ans * a % MOD;
			a = a * a % MOD;
			b >>= 1;
		}
		return ans;
	}
	int getK(int n)
	{
		int s = 1;
		while (s < n)
			s <<= 1;
		return s;
	}
	typedef unsigned us;
	typedef unsigned long long ull;
	//+ ntt
	namespace RawNTT
	{
#define USE_INTRINSICS
		us pool[SS * 8 + 10000] __attribute__((aligned(64))), *ptr = pool;
		us *p0[SS], *p1[SS], *q0[SS], *q1[SS], rv[SS];
		inline void bit_flip(us *p, int t)
		{
			for (int i = 0, j = 0; i < t; ++i)
			{
				if (i > j)
					swap(p[i], p[j]);
				for (int l = t >> 1; (j ^= l) < l; l >>= 1)
					;
			}
		}
		void prep(int n)
		{
			static int t = 1;
			rv[1] = 1;
			for (; t < n; t <<= 1)
			{
				int g = qp(3, (MOD - 1) / (t * 2));
				rv[t + t] = rv[t] * (ull)((MOD + 1) / 2) % MOD;
				us *p, *q;
				p = p0[t] = ptr;
				ptr += max(t, 16);
				p[0] = 1;
				for (int m = 1; m < t; ++m)
					p[m] = p[m - 1] * (ull)g % us(MOD);
				bit_flip(p, t);
				q = q0[t] = ptr;
				ptr += max(t, 16);
				for (int i = 0; i < t; ++i)
					q[i] = (ull(p[i]) << 32) / MOD;
				g = qp(g, MOD - 2);
				p = p1[t] = ptr;
				ptr += max(t, 16);
				p[0] = 1;
				for (int m = 1; m < t; ++m)
					p[m] = p[m - 1] * (ull)g % us(MOD);
				bit_flip(p, t);
				q = q1[t] = ptr;
				ptr += max(t, 16);
				for (int i = 0; i < t; ++i)
					q[i] = (ull(p[i]) << 32) / MOD;
			}
		}
		typedef unsigned long long ull;
		inline us my_mul(us a, us b, us c) { return b * (ull)a - ((ull(a) * c) >> 32) * ull(998244353); }
#ifdef USE_INTRINSICS
#warning intrinsics is on!
#include <immintrin.h>
		inline __m128i my_mullo_epu32(const __m128i &a, const __m128i &b)
		{
#if defined(__SSE4_2__)
			return _mm_mullo_epi32(a, b);
#else
			__m128i a13 = _mm_shuffle_epi32(a, 0xF5), b13 = _mm_shuffle_epi32(b, 0xF5),
					prod02 = _mm_mul_epu32(a, b), prod13 = _mm_mul_epu32(a13, b13),
					prod = _mm_unpacklo_epi64(_mm_unpacklo_epi32(prod02, prod13),
											  _mm_unpackhi_epi32(prod02, prod13));
			return prod;
#endif
		}
		inline __m128i my_mulhi_epu32(const __m128i &a, const __m128i &b)
		{
			__m128i a13 = _mm_shuffle_epi32(a, 0xF5), b13 = _mm_shuffle_epi32(b, 0xF5),
					prod02 = _mm_mul_epu32(a, b), prod13 = _mm_mul_epu32(a13, b13),
					prod = _mm_unpackhi_epi64(_mm_unpacklo_epi32(prod02, prod13),
											  _mm_unpackhi_epi32(prod02, prod13));
			return prod;
		}
#endif
		void ntt(us *x, int n, bool f = true)
		{
			prep(n);
			int t = n;
			for (int m = 1; m < n; m <<= 1)
			{
				t >>= 1;
				us *p = p0[m], *q = q0[m];
#ifdef USE_INTRINSICS
				if (t == 1 || t == 2)
				{
#endif
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
						for (int j = 0; j < t; ++j)
						{
							us u = xa[j] - (xa[j] >= us(MOD + MOD)) * us(MOD + MOD);
							us v = my_mul(xb[j], p[i], q[i]);
							xa[j] = u + v;
							xb[j] = u - v + us(MOD + MOD);
						}
#ifdef USE_INTRINSICS
				}
				else if (t == 4)
				{
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
					{
						const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]),
									  mm = _mm_set1_epi32(MOD + MOD), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(MOD);
						for (int j = 0; j < t; j += 4)
						{
							__m128i u = _mm_loadu_si128((__m128i *)(xa + j));
							u = _mm_sub_epi32(u, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u, mm),
																			_mm_cmpgt_epi32(m0, u)),
															   mm));
							__m128i v = _mm_loadu_si128((__m128i *)(xb + j));
							v = _mm_sub_epi32(my_mullo_epu32(v, p4), my_mullo_epu32(my_mulhi_epu32(v, q4), m1));
							_mm_storeu_si128((__m128i *)(xa + j), _mm_add_epi32(u, v));
							_mm_storeu_si128((__m128i *)(xb + j), _mm_add_epi32(_mm_sub_epi32(u, v), mm));
						}
					}
				}
				else
				{
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
					{
						const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]),
									  mm = _mm_set1_epi32(MOD + MOD), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(MOD);
						for (int j = 0; j < t; j += 8)
						{
							__m128i u0 = _mm_loadu_si128((__m128i *)(xa + j)),
									u1 = _mm_loadu_si128((__m128i *)(xa + j + 4)), v0 = _mm_loadu_si128((__m128i *)(xb + j)),
									v1 = _mm_loadu_si128((__m128i *)(xb + j + 4));
							u0 = _mm_sub_epi32(u0,
											   _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u0, mm), _mm_cmpgt_epi32(m0, u0)), mm));
							u1 = _mm_sub_epi32(u1, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u1, mm),
																			  _mm_cmpgt_epi32(m0, u1)),
																 mm));
							v0 = _mm_sub_epi32(my_mullo_epu32(v0, p4),
											   my_mullo_epu32(my_mulhi_epu32(v0, q4), m1));
							v1 = _mm_sub_epi32(my_mullo_epu32(v1, p4),
											   my_mullo_epu32(my_mulhi_epu32(v1, q4), m1));
							_mm_storeu_si128((__m128i *)(xa + j),
											 _mm_add_epi32(u0, v0));
							_mm_storeu_si128((__m128i *)(xa + j + 4), _mm_add_epi32(u1, v1));
							_mm_storeu_si128((__m128i *)(xb + j), _mm_add_epi32(_mm_sub_epi32(u0, v0), mm));
							_mm_storeu_si128((__m128i *)(xb + j + 4), _mm_add_epi32(_mm_sub_epi32(u1, v1), mm));
						}
					}
				}
#endif
			}
			for (int i = 0; i < n; ++i)
				x[i] -= (x[i] >= us(MOD + MOD)) * us(MOD + MOD),
					x[i] -= (x[i] >= us(MOD)) * us(MOD);
			if (f)
				bit_flip(x, n);
		}
		void intt(us *x, int n, bool f = true)
		{
			prep(n);
			int t = 1;
			if (f)
				bit_flip(x, n);
			for (int m = (n >> 1); m; m >>= 1)
			{
				us *p = p1[m], *q = q1[m];
#ifdef USE_INTRINSICS
				if (t == 1 || t == 2)
				{
#endif
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
						for (int j = 0; j < t; ++j)
						{
							us u = xa[j], v = xb[j];
							xa[j] = u + v - (u + v >= us(MOD + MOD)) * us(MOD + MOD);
							xb[j] = my_mul(u - v + us(MOD + MOD), p[i], q[i]);
						}
#ifdef USE_INTRINSICS
				}
				else if (t == 4)
				{
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
					{
						const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]),
									  mm = _mm_set1_epi32(MOD + MOD), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(MOD);
						for (int j = 0; j < t; j += 4)
						{
							__m128i u = _mm_loadu_si128((__m128i *)(xa + j)),
									v = _mm_loadu_si128((__m128i *)(xb + j)), uv = _mm_add_epi32(u, v);
							_mm_storeu_si128((__m128i *)(xa + j), _mm_sub_epi32(uv, _mm_and_si128(_mm_or_si128(
																									  _mm_cmpgt_epi32(uv, mm), _mm_cmpgt_epi32(m0, uv)),
																								  mm)));
							uv = _mm_add_epi32(
								_mm_sub_epi32(u, v), mm);
							_mm_storeu_si128((__m128i *)(xb + j), _mm_sub_epi32(
																	  my_mullo_epu32(uv, p4), my_mullo_epu32(my_mulhi_epu32(uv, q4), m1)));
						}
					}
				}
				else
				{
					us *xa = x, *xb = x + t;
					for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
					{
						const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]), mm = _mm_set1_epi32(MOD + MOD),
									  m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(MOD);
						for (int j = 0; j < t; j += 8)
						{
							__m128i u0 = _mm_loadu_si128((__m128i *)(xa + j)), u1 = _mm_loadu_si128((__m128i *)(xa + j + 4)),
									v0 = _mm_loadu_si128((__m128i *)(xb + j)), v1 = _mm_loadu_si128((__m128i *)(xb + j + 4)),
									uv0 = _mm_add_epi32(u0, v0), uv1 = _mm_add_epi32(u1, v1);
							_mm_storeu_si128((__m128i *)(xa + j), _mm_sub_epi32(uv0, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(uv0, mm),
																												_mm_cmpgt_epi32(m0, uv0)),
																								   mm)));
							_mm_storeu_si128((__m128i *)(xa + j + 4),
											 _mm_sub_epi32(uv1, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(uv1, mm), _mm_cmpgt_epi32(m0, uv1)), mm)));
							uv0 = _mm_add_epi32(_mm_sub_epi32(u0, v0), mm);
							uv1 = _mm_add_epi32(
								_mm_sub_epi32(u1, v1), mm);
							_mm_storeu_si128((__m128i *)(xb + j), _mm_sub_epi32(
																	  my_mullo_epu32(uv0, p4), my_mullo_epu32(my_mulhi_epu32(uv0, q4), m1)));
							_mm_storeu_si128((__m128i *)(xb + j + 4), _mm_sub_epi32(my_mullo_epu32(uv1, p4), my_mullo_epu32(my_mulhi_epu32(uv1, q4), m1)));
						}
					}
				}
#endif
				t <<= 1;
			}
			us rn = rv[n];
			for (int i = 0; i < n; ++i)
				x[i] = x[i] * (ull)rn % MOD;
		}
	}
	//+ modint
	union mi // modint, treat as POD
	{
		us w;
		mi() { w = 0; }
		mi(us u) { w = u; }
		mi(int u)
		{
			u %= MOD;
			w = u + ((u < 0) ? MOD : 0);
		}
		explicit operator us() const { return w; }
		explicit operator int() const { return w; }
	};
	mi operator+(const mi &a, const mi &b)
	{
		return mi{a.w + b.w - ((a.w + b.w >= MOD) ? (MOD) : 0)};
	}
	mi operator-(const mi &a, const mi &b)
	{
		return mi{a.w - b.w + ((a.w < b.w) ? (MOD) : 0)};
	}
	mi operator*(const mi &a, const mi &b)
	{
		return mi{us((ull)a.w * b.w % MOD)};
	}
	mi operator/(const mi &a, const mi &b)
	{
		return mi{us((ull)a.w * qp(b.w, MOD - 2) % MOD)};
	}
	mi inv(const mi &a) { return mi{us(qp(a.w, MOD - 2))}; }
	bool operator==(const mi &a, const mi &b) { return a.w == b.w; }
	bool operator!=(const mi &a, const mi &b) { return a.w != b.w; }
	mi &operator+=(mi &a, const mi &b)
	{
		a.w = a.w + b.w - ((a.w + b.w >= MOD) ? MOD : 0);
		return a;
	}
	mi &operator-=(mi &a, const mi &b)
	{
		a.w = a.w - b.w + ((a.w < b.w) ? MOD : 0);
		return a;
	}
	mi operator-(const mi &a) { return mi{a.w ? (MOD - a.w) : 0}; }
	mi &operator++(mi &a)
	{
		a.w = a.w + 1 - ((a.w + 1 >= MOD) ? MOD : 0);
		return a;
	}
	mi &operator--(mi &a)
	{
		a.w = a.w - 1 + (a.w ? 0 : MOD);
		return a;
	}
	void ntt(mi *x, int n, bool f = true) { RawNTT::ntt((us *)x, n, f); } // make sure this works
	void intt(mi *x, int n, bool f = true) { RawNTT::intt((us *)x, n, f); }
	void fft(mi *x, int n, bool r, bool f = true)
	{
		if (r)
			intt(x, n, f);
		else
			ntt(x, n, f);
	}
	void cp(mi *t, const mi *s, int K)
	{
		if (s)
			memcpy(t, s, sizeof(mi) * K);
		else
			memset(t, 0, sizeof(mi) * K);
	}
	void cp(mi *t, mi s, int K)
	{
		if (s.w == 0)
			memset(t, 0, sizeof(mi) * K);
		else
			for (int i = 0; i < K; ++i)
				t[i] = s;
	}
	mi qp(mi a, ll b)
	{
		mi x = 1;
		while (b)
		{
			if (b & 1)
				x = x * a;
			a = a * a;
			b >>= 1;
		}
		return x;
	}
	string to_string(mi f) { return to_string((int)f); }
	string pretty_guess(mi x, int max_dem = 1000)
	{
		string s = to_string((int)x);
		auto upd = [&](string v)
		{
			if (v.size() < s.size())
				s = v;
		};
		upd("-" + to_string((int)(-x)));
		for (int i = 1; i <= max_dem; ++i)
		{
			mi w = x * i;
			upd(to_string((int)w) + "/" + to_string(i));
			upd("-" + to_string((int)(-w)) + "/" + to_string(i));
		}
		return s;
	}
	ostream &operator<<(ostream &os, const mi &m)
	{
		os << m.w;
		return os;
	}
	istream &operator>>(istream &is, mi &m)
	{
		int x;
		is >> x;
		m = x;
		return is;
	}
	//+ basic ops
	mi fac[FS], rfac[FS];
	struct Fac_Initer
	{
		Fac_Initer()
		{
			fac[0] = 1;
			for (int i = 1; i < FS; ++i)
				fac[i] = fac[i - 1] * i;
			rfac[FS - 1] = inv(fac[FS - 1]);
			for (int i = FS - 1; i; --i)
				rfac[i - 1] = rfac[i] * i;
		}
	} fac__initer__;
	mi mempool[PS], *pt = mempool;
	mi *alc(int t, bool c = 0)
	{
		if (c)
			cp(pt, 0, t);
		pt += t;
		assert(pt < mempool + PS);
		return pt - t;
	}
	void ginv_K(mi *x, mi *o, int K)
	{
		if (K == 1)
		{
			o[0] = inv(x[0]);
			return;
		}
		ginv_K(x, o, K >> 1);
		mi *fo = alc(K, 1), *fx = alc(K), *fw = alc(K);
		cp(fo, o, (K >> 1));
		fft(fo, K, 0);
		cp(fx, x, K);
		fft(fx, K, 0);
		for (int i = 0; i < K; ++i)
			fw[i] = fx[i] * fo[i];
		fft(fw, K, 1);
		cp(fw, fw + (K >> 1), K >> 1);
		cp(fw + (K >> 1), 0, K >> 1);
		ntt(fw, K);
		for (int i = 0; i < K; ++i)
			fw[i] = fw[i] * fo[i];
		fft(fw, K, 1);
		for (int i = 0; i < (K >> 1); ++i)
			o[i + (K >> 1)] = -fw[i];
		pt -= K + K + K;
	}
	void ginv(mi *x, mi *o, int n)
	{
		int K = getK(n);
		mi *fx = alc(K, 1), *fo = alc(K);
		cp(fx, x, n);
		ginv_K(fx, fo, K);
		cp(o, fo, n);
		pt -= K + K;
	}
	void gdiv(mi *a, mi *b, mi *d, int n, int m)
	{
		int s = getK(max(n, m));
		mi *ra = alc(s + s, 1), *rb = alc(s + s, 1);
		for (int i = 0; i < n; ++i)
			ra[i] = a[n - 1 - i];
		for (int i = 0; i < m; ++i)
			rb[i] = b[m - 1 - i];
		ginv(rb, rb, s);
		fft(ra, s + s, 0);
		fft(rb, s + s, 0);
		for (int i = 0; i < s + s; ++i)
			rb[i] = ra[i] * rb[i];
		fft(rb, s + s, 1);
		for (int i = 0; i <= n - m; ++i)
			d[i] = rb[n - m - i];
		pt -= s * 4;
	}
	void gdiv(mi *a, mi *b, mi *d, mi *r, int n, int m)
	{
		gdiv(a, b, d, n, m);
		int s = getK(n + 1);
		mi *bb = alc(s, 1), *dd = alc(s, 1);
		cp(bb, b, m);
		cp(dd, d, n - m + 1);
		fft(bb, s, 0);
		fft(dd, s, 0);
		for (int i = 0; i < s; ++i)
			bb[i] = -bb[i] * dd[i];
		fft(bb, s, 1);
		for (int i = 0; i < m - 1; ++i)
			r[i] = a[i] + bb[i];
		pt -= s * 2;
	}
	void gln(mi *a, mi *b, int n)
	{
		int s = getK(n + n);
		mi *ra = alc(s, 1);
		ginv(a, ra, n);
		mi *rb = alc(s, 1);
		for (int i = 0; i + 1 < n; ++i)
			rb[i] = a[i + 1] * (i + 1);
		fft(ra, s, 0);
		fft(rb, s, 0);
		for (int i = 0; i < s; ++i)
			ra[i] = ra[i] * rb[i];
		fft(ra, s, 1);
		b[0] = 0;
		for (int i = 1; i < n; ++i)
			b[i] = ra[i - 1] * rfac[i] * fac[i - 1];
		pt -= s * 2;
	}
	mi sqrt_f0;
	void gsqrt_K(mi *f, mi *g, mi *h, int K, bool ch = 1)
	{
		static mi gh[SS];
		if (K == 1)
		{
			assert(sqrt_f0 * sqrt_f0 - f[0] == 0);
			g[0] = sqrt_f0;
			h[0] = inv(sqrt_f0);
			gh[0] = sqrt_f0;
			return;
		}
		gsqrt_K(f, g, h, K >> 1);
		mi *fh = alc(K, 1), *gg = alc(K >> 1), *rr = alc(K, 1);
		cp(fh, h, (K >> 1));
		fft(fh, K, 0);
		cp(gg, gh, K >> 1);
		for (int i = 0; i < (K >> 1); ++i)
			gg[i] = gg[i] * gg[i];
		fft(gg, K >> 1, 1);
		for (int i = 0; i < (K >> 1); ++i)
			rr[i + (K >> 1)] = gg[i] - f[i] - f[i + (K >> 1)];
		fft(rr, K, 0);
		for (int i = 0; i < K; ++i)
			rr[i] = rr[i] * fh[i] * ((MOD + 1) / 2);
		fft(rr, K, 1);
		for (int i = (K >> 1); i < K; ++i)
			g[i] = -rr[i];
		if (ch)
		{
			mi *fg = alc(K), *fw = alc(K);
			cp(fg, g, K);
			fft(fg, K, 0);
			for (int i = 0; i < K; ++i)
				fw[i] = fg[i] * fh[i];
			fft(fw, K, 1);
			for (int i = 0; i < (K >> 1); ++i)
				fw[i] = fw[i + (K >> 1)];
			cp(fw + (K >> 1), 0, K >> 1);
			fft(fw, K, 0);
			for (int i = 0; i < K; ++i)
				fw[i] = fw[i] * fh[i];
			fft(fw, K, 1);
			for (int i = 0; i < (K >> 1); ++i)
				h[i + (K >> 1)] = -fw[i];
			cp(gh, fg, K);
			pt -= K + K;
		}
		pt -= K + K + (K >> 1);
	}
	void gsqrt(mi *f, mi *g, int n)
	{
		int s = getK(n);
		mi *mf = alc(s, 1), *mg = alc(s), *mh = alc(s);
		cp(mf, f, n);
		gsqrt_K(mf, mg, mh, s, 0);
		cp(g, mg, n);
		pt -= s + s + s;
	}
	void gexp_K(mi *f, mi *g, mi *h, int K, bool ch = 1)
	{
		if (K == 1)
		{
			g[0] = h[0] = 1;
			return;
		}
		gexp_K(f, g, h, K >> 1);
		mi *gg = alc(K >> 1), *hh = alc(K >> 1), *fh = 0,
		   *dg = alc(K >> 1), *t1 = alc(K, 1), *t2 = alc(K, 1);
		dg[(K >> 1) - 1] = 0;
		for (int i = 0; i + 1 < (K >> 1); ++i)
			dg[i] = g[i + 1] * (i + 1);
		cp(gg, g, K >> 1);
		mi c = 0;
		for (int i = 0; i < (K >> 1); ++i)
			c = c + dg[i] * h[((K >> 1) - 1) - i];
		if (!ch)
			cp(hh, h, K >> 1), fft(hh, K >> 1, 0);
		else
		{
			fh = alc(K, 1);
			cp(fh, h, (K >> 1));
			fft(fh, K, 0);
			for (int i = 0; i < K; i += 2)
				hh[i >> 1] = fh[i];
		}
		fft(gg, K >> 1, 0);
		fft(dg, K >> 1, 0);
		for (int i = 0; i < (K >> 1); ++i)
			gg[i] = gg[i] * hh[i];
		fft(gg, K >> 1, 1);
		for (int i = 0; i < (K >> 1); ++i)
			t1[i + (K >> 1)] = (i == 0) - gg[i];
		for (int i = 0; i + 1 < (K >> 1); ++i)
			t2[i] = f[i + 1] * (i + 1);
		fft(t1, K, 0);
		fft(t2, K, 0);
		for (int i = 0; i < K; ++i)
			t1[i] = t1[i] * t2[i];
		fft(t1, K, 1);
		for (int i = 0; i < (K >> 1); ++i)
			t1[i] = 0;
		for (int i = 0; i + 1 < K; ++i)
			t1[i] = t1[i] - f[i + 1] * (i + 1);
		for (int i = 0; i < (K >> 1); ++i)
			dg[i] = dg[i] * hh[i];
		fft(dg, K >> 1, 1);
		mi r;
		for (int i = 0; i < (K >> 1); ++i)
			r = (i + 1 == (K >> 1)) ? c : (f[i + 1] * (i + 1)),
			t1[i] = t1[i] + r, t1[i + (K >> 1)] = t1[i + (K >> 1)] + dg[i] - r;
		t2[0] = 0;
		for (int i = 0; i + 1 < K; ++i)
			t2[i + 1] = t1[i] * rfac[i + 1] * fac[i];
		cp(t1, g, K >> 1);
		cp(t1 + (K >> 1), 0, K >> 1);
		fft(t1, K, 0);
		fft(t2, K, 0);
		for (int i = 0; i < K; ++i)
			t1[i] = t1[i] * t2[i];
		fft(t1, K, 1);
		for (int i = (K >> 1); i < K; ++i)
			g[i] = -t1[i];
		pt -= K * 2 + (K >> 1) * 3;
		if (ch)
		{
			mi *fg = alc(K), *fw = alc(K);
			cp(fg, g, K);
			fft(fg, K, 0);
			for (int i = 0; i < K; ++i)
				fw[i] = fg[i] * fh[i];
			fft(fw, K, 1);
			for (int i = 0; i < (K >> 1); ++i)
				fw[i] = fw[i + (K >> 1)];
			cp(fw + (K >> 1), 0, K >> 1);
			fft(fw, K, 0);
			for (int i = 0; i < K; ++i)
				fw[i] = fw[i] * fh[i];
			fft(fw, K, 1);
			for (int i = 0; i < (K >> 1); ++i)
				h[i + (K >> 1)] = -fw[i];
			pt -= K + K + K;
		}
	}
	void gexp(mi *f, mi *g, int n)
	{
		int s = getK(n);
		mi *mf = alc(s, 1), *mg = alc(s), *mh = alc(s);
		cp(mf, f, n);
		gexp_K(mf, mg, mh, s, 0);
		cp(g, mg, n);
		pt -= s + s + s;
	}
	//+ mod_sqrt
	namespace QR
	{
		typedef pair<ll, ll> pll;
		ll pll_s;
		inline pll mul(pll a, pll b, ll p)
		{
			pll ans;
			ans.fi = a.fi * b.fi % p + a.se * b.se % p * pll_s % p;
			ans.se = a.fi * b.se % p + a.se * b.fi % p;
			ans.fi %= p;
			ans.se %= p;
			return ans;
		}
		pll qp(pll a, ll b, ll c)
		{
			pll ans(1, 0);
			while (b)
			{
				if (b & 1)
					ans = mul(ans, a, c);
				a = mul(a, a, c);
				b >>= 1;
			}
			return ans;
		}
		ll qp(ll a, ll b, ll c)
		{
			ll ans = 1;
			while (b)
			{
				if (b & 1)
					ans = ans * a % c;
				a = a * a % c;
				b >>= 1;
			}
			return ans;
		}
		int mod_sqrt(ll a, ll p = MOD)
		{
			if (!a)
				return 0;
			if (p == 2)
				return 1;
			ll w, q;
			while (1)
			{
				w = rand() % p;
				q = w * w - a;
				q = (q % p + p) % p;
				if (qp(q, (p - 1) / 2, p) != 1)
					break;
			}
			pll_s = q;
			pll rst = qp(pll(w, 1), (p + 1) / 2, p);
			ll ans = (rst.fi % p + p) % p;
			return min(ans, p - ans);
		}
	}
	using QR::mod_sqrt;
//+ poly
#include <functional>
	int default_shrink = -1; // mod x^n
	struct poly
	{
		vector<mi> coeff;
		int shrink_len;
		void rev()
		{
			fit_shrink();
			reverse(coeff.begin(), coeff.end());
		}
		void insert(mi x)
		{
			coeff.insert(coeff.begin(), x);
			shrink();
		}
		mi &operator[](int x)
		{
			if ((x < 0) || (shrink_len != -1 && x >= shrink_len))
				throw out_of_range("invalid offset");
			if ((int)coeff.size() < x + 1)
				coeff.resize(x + 1);
			return coeff[x];
		}
		mi operator[](int x) const
		{
			if ((x < 0) || (shrink_len != -1 && x >= shrink_len))
				throw out_of_range("invalid offset");
			if ((int)coeff.size() < x + 1)
				return mi(0);
			return coeff[x];
		}
		mi get(int x) const
		{
			if ((x < 0) || (shrink_len != -1 && x >= shrink_len))
				return 0;
			if ((int)coeff.size() < x + 1)
				return mi(0);
			return coeff[x];
		}
		explicit poly(int shrink_len_ = default_shrink) : shrink_len(shrink_len_)
		{
		}
		poly(vector<mi> coeff_, int shrink_len_ = default_shrink) : coeff(coeff_), shrink_len(shrink_len_)
		{
			this->shrink();
		}
		poly(vector<int> coeff_, int shrink_len_ = default_shrink) : shrink_len(shrink_len_)
		{
			this->coeff.resize(coeff_.size());
			for (int i = 0; i < (int)coeff.size(); ++i)
				this->coeff[i] = coeff_[i];
			this->shrink();
		}
		void clean_maybe()
		{
			if (is_poly())
				while (coeff.size() && coeff.back() == 0)
					coeff.pop_back();
		}
		void clean()
		{
			assert(is_poly());
			clean_maybe();
		}
		void fit_shrink()
		{
			assert(is_series());
			coeff.resize(shrink_len);
		}
		void set_shrink(int shrink_len_ = default_shrink)
		{
			this->shrink_len = shrink_len_;
			this->shrink();
		}
		void dump(char e = 0, bool g = 1, int l = 9) const
		{
			auto format = [&](mi num)
			{
				return g ? pretty_guess(num) : to_string(num);
			};
			int u = (int)coeff.size() - 1;
			while (u >= 0 && coeff[u] == 0)
				--u;
			if (u < 0)
			{
				printf("{}");
			}
			else
			{
				for (int j = 0; j <= u && j <= l; ++j)
					printf("%c%s", "{,"[j != 0], format(coeff[j]).c_str());
				if (u > l)
					printf("...%s(x^%d)", format(coeff[u]).c_str(), u);
				printf("}");
			}
			if (shrink_len == -1)
				printf(" (poly)");
			else
				printf(" (mod x^%d)", shrink_len);
			if (e)
				putchar(e);
		}
		const mi *coeff_ptr() const
		{
			if (!coeff.size())
				return 0;
			return coeff.data();
		}
		mi *coeff_ptr()
		{
			if (!coeff.size())
				return 0;
			return coeff.data();
		}
		int size() const
		{
			return coeff.size();
		}
		void reserve(int l)
		{
			if (shrink_len != -1)
				l = min(l, shrink_len);
			if (l > (int)coeff.size())
				coeff.resize(l);
		}
		void print_shrink(char e)
		{
			fit_shrink();
			for (int i = 0; i < shrink_len; ++i)
			{
				if (i)
					printf(" ");
				printf("%d", (int)coeff[i]);
			}
			if (e)
				printf("%c", e);
		}
		void print_len(int s, char e)
		{
			for (int i = 0; i < s; ++i)
			{
				if (i)
					printf(" ");
				printf("%d", (int)get(i));
			}
			if (e)
				printf("%c", e);
		}
		void shrink()
		{
			if (shrink_len != -1 && (int)coeff.size() > shrink_len)
				coeff.resize(shrink_len);
		}
		bool is_poly() const { return shrink_len == -1; }
		bool is_series() const { return shrink_len != -1; }
		mi eval(mi x)
		{
			assert(is_poly());
			mi w = 0;
			for (int i = size() - 1; i >= 0; --i)
				w = w * x + coeff[i];
			return w;
		}
	};
	poly polyi(mi x)
	{
		return poly(vector<mi>{x}, -1);
	}
	poly operator"" _p(unsigned long long int x)
	{
		return poly(vector<mi>{int(x % MOD)}, -1);
	}
	poly operator"" _p(const char *str, std::size_t len)
	{
		poly ans(-1);
		int sgn = 1, phase = 0, coeff = 0, touch = 0;
		ll cnum = 0;
		auto clean = [&]()
		{
			if (phase == -1)
				ans[1] += sgn * coeff;
			else if (phase == 0)
				ans[0] += sgn * (int)cnum;
			else if (phase == 1)
				ans[cnum] += sgn * coeff;
			else
				assert(0);
			phase = cnum = touch = 0;
		};
		for (int i = 0; i < (int)len; ++i)
		{
			if (str[i] == '+')
				clean(), sgn = 1;
			else if (str[i] == '-')
				clean(), sgn = -1;
			else if (isdigit(str[i]))
			{
				assert(phase == 0 || phase == 1);
				if (phase == 0)
					touch = 1, cnum = (cnum * 10LL + str[i] - 48) % MOD;
				else
					cnum = cnum * 10LL + str[i] - 48, assert(cnum < 1e8);
			}
			else if (str[i] == 'x')
			{
				assert(str[i + 1] == '^' || str[i + 1] == '+' || str[i + 1] == '-' || str[i + 1] == 0);
				phase = -1;
				coeff = touch ? cnum : 1;
				cnum = 0;
			}
			else if (str[i] == '^')
			{
				assert(phase == -1);
				phase = 1;
			}
		}
		clean();
		return ans;
	}
	//+ poly ops
	void share_shrink(poly &a, poly &b)
	{
		int l = max(a.shrink_len, b.shrink_len);
		a.set_shrink(l);
		b.set_shrink(l);
	}
	poly ginv(poly p)
	{
		p.fit_shrink();
		ginv(p.coeff_ptr(), p.coeff_ptr(), p.shrink_len);
		return p;
	}
	poly gln(poly p)
	{
		p.fit_shrink();
		gln(p.coeff_ptr(), p.coeff_ptr(), p.shrink_len);
		return p;
	}
	poly gsqrt(poly p, mi f0 = mi(1))
	{
		p.fit_shrink();
		sqrt_f0 = f0;
		gsqrt(p.coeff_ptr(), p.coeff_ptr(), p.shrink_len);
		return p;
	}
	poly gexp(poly p)
	{
		p.fit_shrink();
		gexp(p.coeff_ptr(), p.coeff_ptr(), p.shrink_len);
		return p;
	}
	int merge_shrink(int s1, int s2)
	{
		if (s1 == -1)
			return s2;
		if (s2 == -1)
			return s1;
		assert(s1 == s2); // usually s1=s2
		return min(s1, s2);
	}
	poly operator+(const poly &a, const poly &b)
	{
		poly c(merge_shrink(a.shrink_len, b.shrink_len));
		c.reserve(max(a.size(), b.size()));
		for (int i = 0; i < c.size(); ++i)
			c[i] = a[i] + b[i];
		return c;
	}
	poly operator-(const poly &a, const poly &b)
	{
		poly c(merge_shrink(a.shrink_len, b.shrink_len));
		c.reserve(max(a.size(), b.size()));
		for (int i = 0; i < c.size(); ++i)
			c[i] = a[i] - b[i];
		return c;
	}
	poly operator-(poly a)
	{
		for (auto &s : a.coeff)
			s = -s;
		return a;
	}
	poly operator*(mi v, poly a)
	{
		for (auto &t : a.coeff)
			t = t * v;
		return a;
	}
	poly operator*(poly a, mi v)
	{
		for (auto &t : a.coeff)
			t = t * v;
		return a;
	}
	poly operator+(poly a, mi b)
	{
		a.reserve(1);
		if (a.size())
			a[0] += b;
		return a;
	}
	poly operator-(poly a, mi b)
	{
		a.reserve(1);
		if (a.size())
			a[0] -= b;
		return a;
	}
	poly operator*(const poly &a, const poly &b)
	{
		if (!a.size())
			return a;
		if (!b.size())
			return b;
		poly c(merge_shrink(a.shrink_len, b.shrink_len));
		c.reserve(a.size() + b.size() - 1);
		int as = min(a.size(), c.size()),
			bs = min(b.size(), c.size());
		int K = getK(as + bs - 1);
		mi *da = alc(K, 1), *db = alc(K, 1);
		for (int i = 0; i < as; ++i)
			da[i] = a[i];
		for (int i = 0; i < bs; ++i)
			db[i] = b[i];
		fft(da, K, 0, 0);
		fft(db, K, 0, 0);
		for (int i = 0; i < K; ++i)
			da[i] = da[i] * db[i];
		fft(da, K, 1, 0);
		for (int i = 0; i < c.size(); ++i)
			c[i] = da[i];
		pt -= K * 2;
		return c;
	}
	pair<poly, poly> gdiv(poly a, poly b)
	{ //(quotient, remainder)
		assert(a.is_poly() && b.is_poly());
		int n = a.size(), m = b.size();
		assert(m > 0);
		if (n < m)
			return make_pair(poly(-1), a);
		poly d(-1), r(-1);
		d.reserve(n - m + 1);
		r.reserve(m - 1);
		gdiv(a.coeff_ptr(), b.coeff_ptr(), d.coeff_ptr(), r.coeff_ptr(), n, m);
		return make_pair(d, r);
	}
	poly operator%(const poly &a, const poly &b)
	{
		return gdiv(a, b).se;
	}
	poly gdiv_quo(poly a, poly b)
	{ // quotient
		assert(a.is_poly() && b.is_poly());
		int n = a.size(), m = b.size();
		assert(m > 0);
		if (n < m)
			return poly(-1);
		poly d(-1);
		d.reserve(n - m + 1);
		gdiv(a.coeff_ptr(), b.coeff_ptr(), d.coeff_ptr(), n, m);
		return d;
	}
	poly gint(poly a)
	{ // shrink actually changes
		a.reserve(a.size() + 1);
		for (int i = a.size() - 1; i >= 1; --i)
			a[i] = a[i - 1] * rfac[i] * fac[i - 1];
		if (a.size())
			a[0] = 0;
		return a;
	}
	poly gde(poly a)
	{ // shrink actually changes
		if (!a.size())
			return a;
		for (int i = 1; i < a.size(); ++i)
			a[i - 1] = a[i] * i;
		a[a.size() - 1] = 0;
		a.clean_maybe();
		return a;
	}
	poly gnewton( // solve G(f)=0
		function<poly(const poly &)> g,
		function<poly(const poly &)> gp,
		int f0, int len = default_shrink)
	{
		poly f(1);
		f[0] = f0;
		while (f.shrink_len != len)
		{
			int old_len = f.shrink_len;
			int new_len = min(old_len * 2, len);
			f.set_shrink(new_len);
			poly s = g(f);
			s.fit_shrink();
			poly h = f;
			h.set_shrink(new_len - old_len);
			h = ginv(gp(h));
			s.coeff.erase(s.coeff.begin(), s.coeff.begin() + old_len);
			s.set_shrink(new_len - old_len);
			s = h * s;
			s.set_shrink(new_len);
			s.coeff.insert(s.coeff.begin(), old_len, mi(0));
			f = f - s;
		}
		return f;
	}
	poly gnewton_d( // solve f'=G(f), gp=(G,G')
		function<pair<poly, poly>(const poly &)> gp,
		int f0, int len = default_shrink)
	{
		poly f(1);
		f[0] = f0;
		while (f.shrink_len != len)
		{
			int old_len = f.shrink_len;
			int new_len = min(old_len * 2, len);
			f.set_shrink(new_len);
			auto fp = gp(f);
			poly gpf = fp.second, r = gexp(mi(-1) * gint(gpf));
			f = gint((fp.first - gpf * f) * r);
			f[0] = f0;
			f = f * ginv(r);
		}
		return f;
	}
	poly gnewton_d(
		function<poly(const poly &)> g,
		function<poly(const poly &)> gp,
		int f0, int len = default_shrink)
	{
		return gnewton_d([&](const poly &s)
						 { return make_pair(g(s), gp(s)); },
						 f0, len);
	}
	poly gcompinv( // G^<-1>, gp=(G,G')
		function<pair<poly, poly>(const poly &)> gp,
		int f0, int len = default_shrink)
	{
		poly f(1);
		f[0] = f0;
		while (f.shrink_len != len)
		{
			int old_len = f.shrink_len;
			int new_len = min(old_len * 2, len);
			f.set_shrink(new_len);
			auto fp = gp(f);
			auto gf = fp.first;
			gf = mi(-1) * gf;
			gf.reserve(2);
			++gf[1];
			f = f + gf * ginv(fp.second);
		}
		return f;
	}
	poly gcompinv(
		function<poly(const poly &)> g,
		function<poly(const poly &)> gp,
		int f0, int len = default_shrink)
	{
		return gcompinv([&](const poly &s)
						{ return make_pair(g(s), gp(s)); },
						f0, len);
	}
	poly prod_recurse(poly *a, int n)
	{
		if (n == 1)
			return *a;
		int m = n >> 1;
		return prod_recurse(a, m) * prod_recurse(a + m, n - m);
	}
	poly prod(vector<poly> p)
	{
		if (!p.size())
			return poly(-1);
		sort(p.begin(), p.end(), [&](const poly &a, const poly &b)
			 { return a.size() < b.size(); });
		return prod_recurse(p.data(), p.size());
	}
	poly gcorner(poly a, vector<mi> b)
	{
		a.reserve(b.size());
		for (int i = 0; i < a.size() && i < (int)b.size(); ++i)
			a[i] = b[i];
		return a;
	}
	poly gshl(poly a, int b)
	{
		a.coeff.insert(a.coeff.begin(), b, mi(0));
		a.shrink();
		return a;
	}
	poly gshr(poly a, int b)
	{
		if (a.size() < b)
			a.coeff.clear();
		else
			a.coeff.erase(a.coeff.begin(), a.coeff.begin() + b);
		a.shrink();
		return a;
	}
	poly gamp(const poly &a, int u)
	{ // A(x)=a(x^u)
		assert(a.is_series() && u >= 1);
		poly b(a.shrink_len);
		for (int i = 0; i * u < a.shrink_len; ++i)
			b[i * u] = a[i];
		return b;
	}
	//+ less-used poly ops
	mi linear_eval(poly a, poly b, ll n)
	{ //[x^n](a(x)/b(x))
		assert(a.is_poly() && b.is_poly() && b.size() >= 1);
		while (n)
		{
			poly nb = b;
			for (int i = 1; i < nb.size(); i += 2)
				nb[i] = -nb[i];
			auto clip = [&](poly p)
			{
				p.reserve(1);
				int u = p.size() - 1;
				for (int i = 1; i <= u / 2; ++i)
					p[i] = p[i + i];
				p.coeff.resize(u / 2 + 1);
				return p;
			};
			poly s = a * nb, t = b * nb;
			if (n & 1)
				s.reserve(1), s.coeff.erase(s.coeff.begin());
			a = clip(s);
			b = clip(t);
			n >>= 1;
		}
		return a.get(0) * inv(b.get(0));
	}
	vector<mi> BM(vector<mi> x)
	{
		vector<mi> ls, cur;
		int lf = 0;
		mi ldt;
		for (int i = 0; i < int(x.size()); ++i)
		{
			mi t = -x[i];
			for (int j = 0; j < int(cur.size()); ++j)
				t = t + x[i - j - 1] * cur[j];
			if (t == 0)
				continue;
			if (!cur.size())
			{
				cur.resize(i + 1);
				lf = i;
				ldt = t;
				continue;
			}
			mi k = -t * inv(ldt);
			vector<mi> c(i - lf - 1);
			c.pb(-k);
			for (int j = 0; j < int(ls.size()); ++j)
				c.pb(ls[j] * k);
			if (c.size() < cur.size())
				c.resize(cur.size());
			for (int j = 0; j < int(cur.size()); ++j)
				c[j] = c[j] + cur[j];
			if (i - lf + (int)ls.size() >= (int)cur.size())
				ls = cur, lf = i, ldt = t;
			cur = c;
		}
		return cur;
	}
	pair<poly, poly> bm_poly(vector<mi> x)
	{
		vector<mi> f = BM(x);
		int k = f.size();
		f.insert(f.begin(), mi(-1));
		x.resize(k);
		poly r(f, -1), s(x, -1);
		poly u = r * s;
		u.coeff.resize(k);
		return make_pair(u, r);
	}
	mi linear_eval(vector<mi> x, ll n)
	{
		auto s = bm_poly(x);
		return linear_eval(s.first, s.second, n);
	}
	vector<poly> eval_tmp;
	void eval_build(int x, mi *a, int n)
	{
		if ((int)eval_tmp.size() < x + 1)
			eval_tmp.resize(x + 1);
		if (n == 1)
		{
			eval_tmp[x] = poly(vector<mi>{mi(1), -*a}, -1);
			return;
		}
		int m = (n + 1) >> 1;
		eval_build(x + x, a, m);
		eval_build(x + x + 1, a + m, n - m);
		eval_tmp[x] = eval_tmp[x + x] * eval_tmp[x + x + 1];
	}
	poly mul_transpose(const poly &a, const poly &b)
	{ // transposed mul
		assert(a.is_poly() && b.size() > 0 && b.is_poly());
		if (a.size() < b.size())
			return poly(-1);
		poly c(-1);
		c.reserve(a.size() - b.size() + 1);
		int K = getK(a.size());
		mi *da = alc(K, 1), *db = alc(K, 1);
		for (int i = 0; i < a.size(); ++i)
			da[i] = a[i];
		for (int i = 0; i < b.size(); ++i)
			db[i] = b[b.size() - 1 - i];
		fft(da, K, 0);
		fft(db, K, 0);
		for (int i = 0; i < K; ++i)
			da[i] = da[i] * db[i];
		fft(da, K, 1);
		for (int i = 0; i < c.size(); ++i)
			c[i] = da[i + b.size() - 1];
		pt -= K * 2;
		return c;
	}
	void eval_recurse(int x, poly p, mi *o, int n)
	{
		if (n == 1)
		{
			*o = p.get(0);
			return;
		}
		int m = (n + 1) >> 1;
		eval_recurse(x + x, mul_transpose(p, eval_tmp[x + x + 1]), o, m);
		eval_recurse(x + x + 1, mul_transpose(p, eval_tmp[x + x]), o + m, n - m);
	}
	vector<mi> multipoint_eval(poly p, vector<mi> q)
	{
		assert(p.is_poly());
		if (!q.size())
			return q;
		eval_build(1, q.data(), q.size());
		int d = p.size();
		p.set_shrink(d);
		p.rev();
		poly o = eval_tmp[1];
		o.set_shrink(d);
		p = p * ginv(o);
		p.rev();
		p.set_shrink(-1);
		p.coeff.resize(q.size());
		vector<mi> s(q.size());
		eval_recurse(1, p, s.data(), q.size());
		eval_tmp.clear();
		return s;
	}
	vector<poly> interp_tmp;
	vector<mi> interp_y;
	void interp_build(int x, mi *a, int n)
	{
		if ((int)interp_tmp.size() < x + 1)
			interp_tmp.resize(x + 1);
		if (n == 1)
		{
			interp_tmp[x] = poly(vector<mi>{-*a, mi(1)}, -1);
			return;
		}
		int m = (n + 1) >> 1;
		interp_build(x + x, a, m);
		interp_build(x + x + 1, a + m, n - m);
		interp_tmp[x] = interp_tmp[x + x] * interp_tmp[x + x + 1];
	}
	poly interp_calc(int x, int o, int n)
	{
		if (n == 1)
			return poly(vector<mi>{interp_y[o]}, -1);
		int m = (n + 1) >> 1;
		return interp_calc(x + x, o, m) * interp_tmp[x + x + 1] + interp_calc(x + x + 1, o + m, n - m) * interp_tmp[x + x];
	}
	poly multipoint_interp(vector<mi> x, vector<mi> y)
	{
		assert(x.size() == y.size());
		interp_build(1, x.data(), x.size());
		interp_y = multipoint_eval(gde(interp_tmp[1]), x);
		for (int i = 0; i < (int)y.size(); ++i)
			interp_y[i] = y[i] * inv(interp_y[i]);
		poly ans = interp_calc(1, 0, x.size());
		interp_tmp.clear();
		interp_y.clear();
		return ans;
	}

}
using POLY::mi;
using POLY::poly;

const int N = 1 << 17 | 1, MOD = 998244353;
ll qpow(ll a, int x)
{
	ll z = 1;
	for (; x; x >>= 1, a = a * a % MOD)
		if (x & 1)
			z = z * a % MOD;
	return z;
}
typedef pair<int, int> pii;
#define fi first
#define se second
int m, q, a[N], ans[N];
vector<pii> buc[N];
poly s1[N << 2], s2[N << 2], s3[N << 2];
poly prod1(int l, int r)
{
	if (l == r)
		return poly(vector<int>{-l, 1});
	int mid = l + r >> 1;
	return prod1(l, mid) * prod1(mid + 1, r);
}
poly prod2(int id, int l, int r)
{
	if (l == r - 1)
		return poly(vector<int>{-buc[id][l].fi, 1});
	int mid = l + r >> 1;
	return prod2(id, l, mid) * prod2(id, mid, r);
}
#define lc (k << 1)
#define rc (lc | 1)
#define left lc, l, mid
#define right rc, mid + 1, r
void build(int k, int l, int r)
{
	if (l == r)
	{
		ll tmp = 1;
		fr(i, a[l - 1] + 1, a[l]) tmp = 1ll * tmp * i % MOD;
		s1[k] = prod1(a[l - 1], a[l] - 1) * mi(int(qpow(tmp, MOD - 2)));
		if (buc[l].size())
			s2[k] = prod2(l, 0, buc[l].size());
		else
			s2[k] = poly(vector<int>{1});
		s3[k] = s1[k];
		return;
	}
	int mid = l + r >> 1;
	build(left);
	build(right);
	s1[k] = s1[lc] * s1[rc];
	s2[k] = s2[lc] * s2[rc];
	s3[k] = s3[lc] + s1[lc] * s3[rc];
}
void divide(int k, int l, int r, poly s, poly mul)
{
	if (s2[k].size() == 0)
		return;
	s = s % s2[k];
	mul = mul % s2[k];
	if (l == r)
	{
		s = (s + mul * s1[k]) % s2[k];
		vector<mi> X;
		for (auto x : buc[l])
			X.push_back(mi(x.fi));
		auto Y = multipoint_eval(s, X);
		fo(i, 0, Y.size()) ans[buc[l][i].se] = Y[i].w;
	}
	int mid = l + r >> 1;
	divide(left, s, mul);
	divide(right, s + mul * s3[lc], mul * s1[lc]);
}
int main()
{
	scanf("%d%d", &m, &q);
	fr(i, 1, m) scanf("%d", a + i);
	fr(id, 1, q)
	{
		int n, R;
		scanf("%d%d", &n, &R);
		buc[R].push_back(pii(n, id));
	}
	build(1, 1, m);
	divide(1, 1, m, poly(vector<int>{0}), poly(vector<int>{1}));
	fr(i, 1, q) printf("%d\n", ans[i]);
	return 0;
}
