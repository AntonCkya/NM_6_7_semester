import matplotlib.pyplot as plt

def read_time(s: str):
    if 'ms' in s:
        return float(s[:-3])
    elif 'm' in s:
        mins, secs = s.split('m')
        return float(mins) * 1000 * 60 + float(secs[:-2]) * 1000
    else:
        return float(s[:-2]) * 1000

d = {
    1: [-1] * 13,
    5: [-1] * 13,
    10: [-1] * 13,
    100: [-1] * 13
}

for i in range(10, 75, 5):
    f = open(f"{i}.out", 'r')
    pos = (i // 5) - 2
    f.readline()
    t1 = read_time(f.readline()[4:])
    t5 = read_time(f.readline()[4:])
    t10 = read_time(f.readline()[5:])
    t100 = read_time(f.readline()[6:])
    d[1][pos] = t1
    d[5][pos] = t5
    d[10][pos] = t10
    d[100][pos] = t100
    print(f"{i}.out readed")

plt.figure(figsize=(8, 4))
x = range(10, 75, 5)
for n, t in d.items():
    print("n:", n, "t:", t, "len(t):", len(t))
    plt.plot(x, t, label = str(n) + " threads")
plt.legend()
plt.xlabel("n")
plt.ylabel("t, ms")
plt.show()
