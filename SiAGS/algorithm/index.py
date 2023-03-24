MAXN = pow(2, 16)   # pow(2, 32) in real implementation
N = 0
ref = [0 for i in range(MAXN)]   # 2 bit * pow(2, 32) = 1G
grp = [0 for i in range(MAXN)]   # 4 bit * pow(2, 32) = 2G
grc = [0 for i in range(17)]
pmt1 = [0 for i in range(MAXN//16)]   # 4 byte * pow(2, 28) = 1G
pmt2 = [0 for i in range(MAXN//16)]   # 4 byte * pow(2, 28) = 1G
cls1 = [0 for i in range(MAXN//16)]   # 4 byte * pow(2, 28) = 1G
cls2 = [0 for i in range(MAXN//16)]   # 4 byte * pow(2, 28) = 1G
cnt = [0 for i in range(MAXN//16)]   # 4 byte * pow(2, 28) = 1G
bwt = [0 for i in range(MAXN//16)]   # 4 bit * pow(2, 28) = 128M
occ = [[0 for i in range(MAXN//16)] for j in range(4)] # 4 byte * pow(2, 28) * 4 = 4G

def read():
    seq = "ATCGATAGTCGTAGC"
    global N
    N = 15
    c2i = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(N):
        ref[i] = c2i[seq[i]]

# Divide into 16 groups (AA ~ TT) and conquer
# O(N)
def group():
    grp[0] += ref[0] * 4
    grc[0] += 1
    for i in range(1, N):
        grp[i-1] += ref[i]
        grp[i] += ref[i] * 4
        grc[grp[i-1]+1] += 1
    grc[ref[N-1]*4] += 1
    for i in range(1, 17):
        grc[i] += grc[i-1]

# Build the suffix array
# O(Nlog(N))
def build():
    for i in range(N):
        pmt1[i] = i
        cls1[i] = ref[i] + 1
    pmt1[N] = N
    cls1[N] = 0
    M = N + 1
    classes = 5
    l = 1
    while l//2 < M:
        for i in range(M):
            pmt2[i] = (pmt1[i] + M - l//2) % M
        for i in range(M):
            cls2[i] = cls1[i]
            cls1[i] = 0
            cnt[i] = 0
        for i in range(M):
            cnt[cls2[pmt2[i]]] += 1
        for i in range(1, classes):
            cnt[i] += cnt[i-1]
        for i in range(M-1, -1, -1):
            cnt[cls2[pmt2[i]]] -= 1
            pmt1[cnt[cls2[pmt2[i]]]] = pmt2[i]
        cls1[pmt1[0]] = 0
        classes = 1
        for i in range(1, M):
            if cls2[pmt1[i]] != cls2[pmt1[i-1]]:
                classes += 1
            elif cls2[(pmt1[i]+l//2)%M] != cls2[(pmt1[i-1]+l//2)%M]:
                classes += 1
            cls1[pmt1[i]] = classes - 1
        l *= 2

def index():
    read()
    group()
    build()

if __name__ == '__main__':
    index()