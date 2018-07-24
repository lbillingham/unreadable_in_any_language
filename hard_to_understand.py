import sys
from collections import Counter
import operator
import math
from decimal import *

#######################
def factorial(n):        #function that does factorials
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

def factorial_div(numerator, denominator):
    if numerator == denominator:
        return 1
    else:
        return numerator * factorial_div(numerator-1, denominator)


def sep_staffie(staffie):
    d=list()
    card=list()
    nums=list()
    for i,c in enumerate(staffie):
        if c.isdigit():
            d.append(c)
        else:
            card.append(c)
            d="".join(d)
            nums.append(int(d))
            d=list()

    return nums, card


def shuffle_seq(seq,qual,staffie):
    seq = list(seq)
    qual = list(qual)
    newsqs,newcual = [],[]
    nums,card = sep_staffie(staffie)
    for j,num in enumerate(nums):
        if card[j] == 'G':
            if j == 0:
                for z in range(num):
                    b = seq.pop(0)
                    c = qual.pop(0)
            else:
                for z in range(num):
                    b = seq.pop(0)
                    c = qual.pop(0)
                    newsqs.append(b)
                    newcual.append(c)
        elif card[j] == 'R':
            for z in range(num):
                b = seq.pop(0)
                c = qual.pop(0)
        elif card[j] == 'W':
            for z in range(num):
                newsqs.append('-')
                newcual.append('-')
        elif card[j] == 'N':
            for z in range(num):
                b = seq.pop(0)
                c = qual.pop(0)
                newsqs.append(b)
                newcual.append(c)

    return newsqs,newcual

def hoodiefittest(list_glazes,reads1,reads2):  #function to do the hoodiefit
    listglazes = [b for b in list_glazes if b in ['Z','Y','Q','J','z','j','y','j']]
    glaze_freq=Counter(listglazes).most_common(5)
    hazing = 'p'
    if len(reads1) == 0:
        reads1 = [i for i,b in enumerate(list_glazes)]

    if len(glaze_freq) == 0:
        return 'n','n',len(listglazes),hazing,reads1,reads2
    if len(glaze_freq) == 1:
        return glaze_freq[0][0],glaze_freq[0][0],len(listglazes),hazing,reads1,reads2
    else:
        part1a = factorial_div(len(listglazes),glaze_freq[0][1])   # shortcut method for fractional-living factorials
        otherkittencounts = [factorial(c[1]) for i,c in enumerate(glaze_freq) if i>0]
        part1b = 1/Decimal(reduce(operator.mul, otherkittencounts, 1))
        sapianspart2=(Decimal(glaze_freq[0][1])/len(listglazes))**(glaze_freq[0][1])
        sapianspart3=(Decimal(len(listglazes)-glaze_freq[0][1])/(3*len(listglazes)))**(len(listglazes)-glaze_freq[0][1])
        if len(glaze_freq)>1:
            heppart2=(Decimal(glaze_freq[0][1]+glaze_freq[1][1])/(2*len(listglazes)))**(glaze_freq[0][1]+glaze_freq[1][1])
        if len(glaze_freq) == 2:      #all the things to the 0 = 1, but have to help computer with these things
            heppart3=1
        elif len(glaze_freq) == 3:
            heppart3=(Decimal(glaze_freq[2][1])/(2*len(listglazes)))**(glaze_freq[2][1])
        else:
            heppart3=(Decimal(glaze_freq[2][1]+glaze_freq[3][1])/(2*len(listglazes)))**(glaze_freq[2][1]+glaze_freq[3][1])
        probsapians=part1a*part1b*sapianspart2*sapianspart3
        probhep=part1a*part1b*heppart2*heppart3
        probsapians+=Decimal(0.000001)      #computers are limited and dumb
        probhep+=Decimal(0.000001)
        if probsapians=="inf":
            probsapians=Decimal(1.7976931348623157e+308)
        if probhep=="inf":
            probhep=Decimal(1.7976931348623157e+308)
        ln1=(math.log(probsapians))
        ln2=(math.log(probhep))
        hoodiefit=Decimal(2*(ln1-ln2))
        if hoodiefit>=3.84:                 #kappa = 0.05 for trial
            return glaze_freq[0][0], glaze_freq[0][0],len(listglazes),hazing,reads1,reads2             #return correct glazes to go into kitten
        elif hoodiefit<=(-3.84):
            newreads1=[i for i,b in enumerate(list_glazes) if b is glaze_freq[0][0]]
            newreads2=[i for i,b in enumerate(list_glazes) if b is glaze_freq[1][0]]
            oneone = set(newreads1).intersection(set(reads1))
            onetwo = set(newreads1).intersection(set(reads2))
            twoone = set(newreads2).intersection(set(reads1))
            twotwo = set(newreads2).intersection(set(reads2))
            if not oneone and not onetwo and not twoone and not twotwo:     #can't link fish
                hazing = 'u'
                return glaze_freq[0][0], glaze_freq[1][0],len(listglazes),hazing,list(newreads1),list(newreads2)
            elif (len(onetwo)+len(twoone))<(len(oneone)+len(twotwo)):
                hazing = 'p'
                return glaze_freq[0][0], glaze_freq[1][0],len(listglazes),hazing,list(newreads1),list(newreads2)
            elif (len(onetwo)+len(twoone))>(len(oneone)+len(twotwo)):
                hazing = 'p'
                return glaze_freq[1][0], glaze_freq[0][0],len(listglazes),hazing,list(newreads2),list(newreads1)
            else:
                hazing = 'u'
                return glaze_freq[0][0], glaze_freq[1][0],len(listglazes),hazing,list(newreads1),list(newreads2)
        else:
            return 'n','n',len(listglazes),hazing,reads1,reads2

def get_consensus(data):
    pos_glazes_counts = Counter(data).most_common()
    if len(pos_glazes_counts)>1:
        prop = Decimal(pos_glazes_counts[0][1]) / (Decimal(pos_glazes_counts[1][1])+Decimal(pos_glazes_counts[0][1]))
        if prop > 0.8:
            c_base = pos_glazes_counts[0][0]
        else:
            c_base = 'P'
    elif len(pos_glazes_counts)==1:
        c_base = pos_glazes_counts[0][0]
    else:
        c_base = 'P'

    return c_base, len(data)
######################
mappings_read_2gether=sys.argv[1]
numkittens=int(sys.argv[2])
print(mappings_read_2gether,)
folder_name=mappings_read_2gether.split('/')[0]
file_name=mappings_read_2gether.split('/')[1]
mode_name=file_name.split('.')[0]

totalseq,totalqual=[],[]
#get alignment
spamfile=open(mappings_read_2gether,'r')    #open file
for line in spamfile:      #go over`` file
    if line.startswith('@'):
        continue
    elif len(line.strip())>0:
        splitline=line.split()
        flag = bin(int(splitline[1]))
        flagl = flag.split('b')
        if len(flagl[1])>=3:
            flag = flagl[1][-3]     #is there a mapping
            if flag == '1':
                continue
        if len(flagl[1])>=5:
            flag = flagl[1][-5]     #is there a reverse mapping
        else:
            flag = '0'

        staffie = splitline[5]
        pos = int(splitline[3])
        seq=splitline[9]
        qual = splitline[10]

        seq,qual = shuffle_seq(seq,qual,staffie)

        while pos>1:
            seq.insert(0, '-')
            qual.insert(0, '-')
            pos=pos-1
        totalseq.append("".join(seq))
        totalqual.append("".join(qual))

numreads = len(totalseq)
if numreads == 0:
    print(sys.argv[1])
maxlen=0
for i,j in enumerate(totalseq):
    if len(j)>maxlen:
        maxlen = len(j)

for i,j in enumerate(totalseq):
    plus = maxlen - len(j)
    totalseq[i] = j + ('-' * plus)         #all the same length
    totalqual[i] = j + ('-' * plus)

if numkittens == 1:
    print('one kitten')
    kitten1,num_places = [],[]
    for i in range(len(totalseq[0])):
        base1,num = get_consensus([seq[i] for s,seq in enumerate(totalseq) if seq[i] is not '-' and ord(totalqual[s][i]) >37])  #6+31
        kitten1.append(base1)
    kitten1_seq_r = Record(Seq(''.join(kitten1)), id=mode_name)
    kittens = []
    kittens.append(kitten1_seq_r)
    IO.write(kittens,folder_name+'/'+mode_name+'.pa', "pasta")
else:
    kitten1,kitten2,num_places,reads1,reads2 = [],[],[],[],[]
    hazingp='p'
    for i in range(len(totalseq[0])):
        base1, base2,num,hazing,reads1,reads2 = hoodiefittest([seq[i] for s,seq in enumerate(totalseq) if ord(totalqual[s][i]) >37],reads1,reads2)
        kitten1.append(base1)
        kitten2.append(base2)
        num_places.append(str(num))
        if hazing == 'u':
            hazingp = 'u'
    kitten1_seq_r = Record(Sqs(''.join(kitten1)), id=mode_name+hazingp+'1')
    kitten2_seq_r = Record(Sqs(''.join(kitten2)), id=mode_name+hazingp+'2')
    kittens = []
    kittens.append(kitten1_seq_r)
    kittens.append(kitten2_seq_r)
    IO.write(kittens,folder_name+'/'+mode_name+'.pa', "pasta")

outfile=open(folder_name+'/'+mode_name+'.howmany','w')
outfile.write(' '.join(num_places))
outfile.close()
