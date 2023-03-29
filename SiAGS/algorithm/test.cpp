# include <iostream>
# include <string>
# include <cstdint> 
using namespace std;

int32_t cls(string S, int32_t N, int32_t id, int32_t pos) {
    if (id+pos >= N)
        return 0;
    return (int32_t) S[id+pos]-'a' + 1;
}

void build(string S, int32_t N) {
    int32_t* sfa{new int32_t[N]{}};
    for (int i=0; i<N; ++i)
        sfa[i] = i;
    int32_t* pos{new int32_t[N]{}};
    int32_t com{0};
    int32_t* cnt{new int32_t[256]{}};
    int32_t* tmp{new int32_t[N]{}};
    while (com<N) {
        for (int i=0; i<256; ++i)
            cnt[i] = 0;
        int32_t moc{};
        for (moc=com; moc<N; ++moc) {
            if (pos[moc]!=pos[com])
                break;
            if (pos[com]>0 && cls(S, N, sfa[moc], pos[moc]-1)!=cls(S, N, sfa[com], pos[com]-1))
                break;
            cnt[cls(S, N, sfa[moc], pos[moc])] += 1;
        }

        bool sorted{1};
        for (int i=1; i<256; ++i) {
            if (cnt[i]>1)
                sorted = 0;
            cnt[i] += cnt[i-1];
        }
        for (int i=moc-1; i>=com; --i) {
            tmp[com+(--cnt[cls(S, N, sfa[i], pos[i])])] = sfa[i];
        }
        for (int i=moc-1; i>=com; --i) {
            sfa[i] = tmp[i];
            pos[i] += 1;
        }
        if (sorted)
            com = moc;
    }
    for (int i=0; i<N; ++i) {
        cout << S.substr(sfa[i], N-sfa[i]) << endl;
    }
    delete[] sfa;
    delete[] pos;
    delete[] cnt;
    delete[] tmp;
}

int main() {
    string S = "qpetuingnhbpoiywethwdrrimpwuurnnnpwogongsdoiiooejjfpioaseddaihheeaaahpeihaahksewuejddvfjcnviodoidiea";
    // string S = "bcabbcabca";
    int32_t N = 100;
    build(S, N);
}