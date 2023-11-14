
#include <bits/stdc++.h>
#define LL long long
using namespace std;
LL M = 998244353;

LL qpow(LL x, LL y, LL mod = M) {
    LL re = 1;
    x %= mod;
    while (y > 0) {
        if (y & 1) re = re * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return re;
}

void MAIN() {
    LL n, m, k;
    cin >> n >> m >> k;
    if (n != k)
        cout << qpow(m, n - k, M) << endl;
    else
        cout << qpow(m, n, M) << endl;
    return;
}

int main() {
    ios::sync_with_stdio(false);
    int ttt = 1;
    cin >> ttt;
    for (int i = 1; i <= ttt; ++i) {
        MAIN();
    }
    return 0;
}
/*


*/
