#include <bits/stdc++.h>
using namespace std;
int n, L[26], R[26], cnt[26], tot, pos, t;
char s[1000002], Mn[1000002];
int main()
{
	scanf("%d", &t);
	while (t--)
	{
		scanf("%s", s + 1), n = strlen(s + 1);
		for (int i = 0; i < 26; ++i)
			L[i] = n + 1, R[i] = 0, cnt[i] = 0;
		for (int i = 1; i <= n; ++i)
			s[i] -= 'a', L[s[i]] = min(L[s[i]], i), R[s[i]] = max(R[s[i]], i), ++cnt[s[i]];
		int n1 = 0, n2 = 0;
		for (int i = 0; i < 26; ++i)
			if (L[i] <= R[i])
				n1 += (R[i] - L[i] + 1) == cnt[i], ++n2;
		if (n1 == n2 && n1 <= 2)
		{
			puts("-1");
			continue;
		}
		pos = n;
		while (s[pos - 1] == s[pos])
			--pos;
		Mn[n] = s[n];
		for (int i = n - 1; i; --i)
			Mn[i] = min(Mn[i + 1], s[i]);
		for (int i = 0; i < 26; ++i)
			if (L[i] <= R[i])
			{
				if (R[i] - L[i] + 1 != cnt[i])
				{
					int ss = 0, mx = 0;
					for (int j = 1; j <= n; ++j)
						if (s[j] == i)
							++ss, mx = max(mx, ss);
						else
							ss = 0;
					for (int j = 1; j <= mx + 1; ++j)
						putchar(i + 'a');
					puts("");
				}
				else
				{
					if (R[i] + 1 != pos && R[i] != n)
					{
						for (int j = 1; j <= cnt[i]; ++j)
							putchar(i + 'a');
						int cc = R[i] + 1;
						while (s[cc] == Mn[cc] && cc < pos - 1)
							putchar(s[cc] + 'a'), ++cc;
						char mn = 26;
						for (int j = cc + 1; j <= n; ++j)
							if (s[j] != s[cc])
								mn = min(mn, s[j]);
						putchar(mn + 'a');
						puts("");
					}
					else if (R[i] == n)
					{
						for (int j = i + 1; j < 26; ++j)
							if (L[j] <= R[j])
							{
								if (R[j] != L[i] - 1)
									putchar(j + 'a'), putchar(i + 'a'), puts("");
								else if (R[j] - L[j] + 1 != cnt[j])
								{
									int pp = R[j], cc = 1;
									while (s[pp - 1] == s[pp])
										++cc, --pp;
									int ss = 0, mx = 0;
									for (int k = 1; k <= n; ++k)
										if (s[k] == j)
											++ss, mx = max(mx, ss);
										else
											ss = 0;
									for (int k = 1; k <= cc + 1; ++k)
										putchar(j + 'a');
									if (cc ^ mx)
										putchar(i + 'a');
									puts("");
								}
								else
								{
									char mn = 26;
									for (int k = 1; k < L[j]; ++k)
										mn = min(mn, s[k]);
									putchar(mn + 'a'), putchar(i + 'a'), puts("");
								}
								break;
							}
					}
					else
					{
						char mn = 26;
						for (int k = 1; k < L[i]; ++k)
							mn = min(mn, s[k]);
						if (mn != s[L[i] - 1])
							putchar(mn + 'a'), putchar(i + 'a'), puts("");
						else if (L[i] != R[i] || mn > s[R[i] + 1])
						{
							putchar(mn + 'a');
							for (int k = 1; k < cnt[i]; ++k)
								putchar(i + 'a');
							putchar(s[R[i] + 1] + 'a'), puts("");
						}
						else
						{
							int ia = 1;
							for (int k = L[mn]; k < L[i]; ++k)
								ia &= s[k] == mn;
							if (ia)
							{
								if (mn != s[n])
								{
									for (int k = L[mn]; k < L[i]; ++k)
										putchar(s[k] + 'a');
									putchar(s[n] + 'a');
									puts("");
								}
								else
								{
									int mx = max(L[i] - L[mn], n - R[i]);
									for (int k = 1; k <= mx + 1; ++k)
										putchar(s[n] + 'a');
									puts("");
								}
							}
							else
							{
								int pp = L[i] - 1, cc = 1;
								while (s[pp] == s[pp - 1])
									--pp, ++cc;
								int ss = 0, mx = 0;
								for (int k = 1; k <= n; ++k)
									if (s[k] == s[L[i] - 1])
										++ss, mx = max(mx, ss);
									else
										ss = 0;
								for (int k = 1; k <= cc + 1; ++k)
									putchar(s[L[i] - 1] + 'a');
								if (cc != mx)
									putchar(i + 'a');
								puts("");
							}
						}
					}
				}
				break;
			}
	}
}
