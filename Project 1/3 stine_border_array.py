x = "abca"
m = "abcannabcnnabc"


def borderarray(x):
    ba = list([0] * len(x))
    for i in range(1,len(x)):
        b = ba[i-1]
        while x[b] != x[i] and b > 0: # I'm not sure that we get to the last index and extend the empty string
            b -= 1
        if x[b] == x[i]:
            ba[i] = b + 1

    return ba

def bordersearch(x,m):
    ba = borderarray(x)
    print(ba)
    bm = list([0] * len(m))
    for i in range(len(m)):
        b = ba[len(x)-1]
        while x[b] != m[i] and b > 0: # I'm not sure that we get to the last index and extend the empty string
            b -= 1
        if x[b] == m[i]:
            bm[i] = bm[i-1] + 1

    return bm


# It does not catch the second match


print(bordersearch(x,m))
