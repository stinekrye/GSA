x = "abcab"


def borderarray(x):
    ba = list([0] * len(x))
    for i in range(1,len(x)):
        b = ba[i-1]
        j = i
        while x[b] != x[j] and len(ba[0:j]) > 0: # I'm not sure that we get to the last index and extend the empty string
            j -= 1
            b = ba[j-1]
        if x[b] == x[i]:
            ba[i] = b + 1

    return ba

print(borderarray(x))
