import random

for n in range(10, 75, 5):
    f = open(f"{n}.txt", "w")
    arr = [[-1 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            arr[i][j] = random.randint(1, 10)
            arr[j][i] = arr[i][j]
    filee = ""
    filee += str(n) + "\n"
    for i in arr:
        for j in i:
            filee += str(j) + " "
        filee += "\n"
    filee += '0.001\n'
    filee += '4\n'
    filee += '1\n'
    filee += '5\n'
    filee += '10\n'
    filee += '100\n'
    f.write(filee)
    f.close()
