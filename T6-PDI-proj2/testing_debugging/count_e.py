import sys

filename = sys.argv[1]

c, h, n, o, S, Se = 0, 0, 0, 0, 0, 0

with open(filename, "r") as f:
    for line in f.readlines()[2:]:
        word = line.split()
        #print(word)
        if word[0] == "C":
            c += 1
        elif word[0] == "H":
            h += 1
        elif word[0] == "O":
            o += 1
        elif word[0] == "N":
            n += 1
        elif word[0] == "Se":
            Se += 1
        elif word[0] == "S":
            S += 1

print("c, h, o, n, S, Se")
print(c, h, o, n, S, Se)
print(c * 4 + h + o*6 +n*5 + S*6+Se*6)
