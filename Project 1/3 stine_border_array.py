x = "abcabc"


def borderarray(x):
    ba = list([0] * len(x))
    for i in range(1,len(x)):
        b = ba[i-1]
        while x[b] != x[i] and b > 0: # I'm not sure that we get to the last index and extend the empty string
            b -= 1
        if x[b] == x[i]:
            ba[i] = b + 1

    return ba

print(borderarray(x))
