import matplotlib.pyplot as plt
import random

def subword_complexity(w,k):
    allwords = dict()
    F = 0
    for i in range(len(w)-k+1):
        if(not (w[i:i+k] in allwords)):
            F = F+1
            allwords[w[i:i+k]] = True
    return F

#print(subword_complexity("aaaaaaabbbbbbbabba",4))

def random_dna(length):
    w = ""
    for i in range(length):
        w = w + random.choice(['A', 'T', 'C', 'G'])
    return w

def fibonacci_word(length):
    w="a"
    v=""
    while(len(w)<length):
        for x in w:
            if(x=='a'):
                v = v + "ab"
            else:
                v = v + "a"
        w = v
        v = ""
    return w[0:length]

#print("    " + fibonacci_word(18))

with open("part_of_genome.txt", "r") as file:
    dna = file.read().replace("\n", "")

x1 = []
y1 = []
for i in range(1,30):
    x1.append(i)
    y1.append(subword_complexity(dna, i))

fig, ax = plt.subplots()

ax.plot(x1,y1,linewidth=2)

plt.show()


rand = random_dna(len(dna))

x2 = []
y2 = []
for i in range(1,30):
    x2.append(i)
    y2.append(subword_complexity(rand, i))

fig, ax = plt.subplots()

ax.plot(x2,y2,linewidth=2)

plt.show()


fib = fibonacci_word(len(dna))

x3 = []
y3 = []
for i in range(1,30):
    x3.append(i)
    y3.append(subword_complexity(fib, i))

fig, ax = plt.subplots()

ax.plot(x3,y3,linewidth=2)

plt.show()


#The two first plots look very alike with an exponential curve followed by a logarithmic one, while the third one is linear. The difference between the
#two last curves is probably due to the fact that the fibonacci word is well constructed with precise rules, while the random one is not.
#This brings the intuition that the DNA does not really follow a pattern or is not constructed by simple rules on words (or by rules that have the same properties as random
#generation on subword complexity), and can be thought of as a specific random string that represents something (for the moment and with only this first observation). 
