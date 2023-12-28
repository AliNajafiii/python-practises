phe1=['ttt','ttc']
leu2=['tta','ttg','ctt','ctc','cta','ctg']
ile3=['att','atc','ata']
met4=['atg']
val5=['gtt','gtc','gta','gtg']
ser6=['tct','tcc','tca','tcg','agt','agc']
pro7=['cct','ccc','cca','ccg']
thr8=['act','acc','aca','acg']
ala9=['gct','gcc','gca','gcg']
tyr10=['tat','tac']
his11=['cat','cac']
gln12=['caa','cag']
asn13=['aat','aac']
lys14=['aaa','aag']
asp15=['gat','gac']
glt16=['gaa','gag']
cys17=['tgt','tgc']
trp18=['tgg']
arg19=['cgt','cgc','cga','cgg','aga','agg']
gly20=['ggt','ggc','gga','ggg']
stop=['taa','tag','tga']

gcode=[phe1,leu2,ile3,met4,val5,ser6,pro7,thr8,ala9,tyr10,
        his11,gln12,asn13,lys14,asp15,glt16,cys17,trp18,arg19,gly20]
aacode=['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
ncode=["a","t","c","g"]
gseq=input("Enter a gene sequence: ")
bcl2_h=input("Enter a human gene sequence for tlanslocation mutation: ")
#gseq='atgactgaatataaacttgtggtagttggagctggtggcgtaggcaagagtgccttgacgatacagctaattcagaatcattttgtggacgaatatgatccaacaatagaggattcctacaggaagcaagtagtaattgatggagaaacctgtctcttggatattctcgacacagcaggtcaagaggagtacagtgcaatgagggaccagtacatgaggactggggagggctttctttgtgtatttgccataaataatactaaatcatttgaagatattcaccattatagagaacaaattaaaagagttaaggactctgaagatgtacctatggtcctagtaggaaataaatgtgatttgccttctagaacagtagacacaaaacaggctcaggacttagcaagaagttatggaattccttttattgaaacatcagcaaagacaagacagagagtggaggatgctttttatacattggtgagagagatccgacaatacagattgaaaaaaatcagcaaagaagaaaagactcctggctgtgtgaaaattaaaaaatgcattataatgtaa'
#bcl2_h='atggacgaggacgttttgcctggagaggtgttggccattgaagggatattcatggcctgtggattaaacgaacctgagtacctgtaccatcctttgctcagccctattaagctatacatcacaggcttaatgcgagacaaggagtctttattcgaggccatgttggctaatgtgagatttcacagcaccaccggtataaaccagcttgggttgagcatgctgcaggttagcggcgatggaaacatgaactgggggcgagccctggctatactgacctttggcagttttgtggcccagaagttatccaacgaacctcacctgcgagactttgctttggccgttttacctgtatatgcgtatgaagcaatcggaccccagtggtttcgcgctcgcggaggctggcgaggcctgaaggcgtattgtacacaggtgcttaccagaagaaggggacggagaatgacagcgctattgggaagcattgcattattggccactatattggcagcggtcgcgatgagcaggagataa'
gseqc=list(gseq)
gseqc2=gseqc.copy()

point=int(input('How many point mutation do you want? '))
inv=input('How many inversion do you want at each time and how many bases of each and separated by slash ? for example two times and 3 bases: 2,3/ : ')
ins=input('How many insertion do you want at each time and how many bases of each and separated by slash ? for example two times and 3 bases: 2,3/ : ')
dele=input('How many deletion do you want at each time and how many bases of each and separated by slash ? for example two times and 3 bases: 2,3/ : ')
trl=input("Do you want talanslocation? Enter y or n: ")

#point mutation
import random
rndind= random.sample(range(len(gseq)),point)
for i in range(point):
    ncodec=ncode.copy()
    print("gseq[",i,"] in native gene is:",gseqc[rndind[i]])
    ncodec.remove(gseqc[rndind[i]])
    rndn=random.sample(ncodec,1)
    gseqc[rndind[i]]=rndn[0]
    print("gseq[",i,"] after point mutation:",gseqc[rndind[i]])
    fr=0
    gseqc=''.join(gseqc)
    for f in range(0,len(gseqc)-2,3):
        if gseqc[f:f+3]==gseqc2[f:f+3]:
            fr+=1
        if gseqc[f:f+3] in stop:
            s=f+3
            break
    if gseqc2==gseqc:
        print("silent mutation")
        print()
    elif len(gseqc2)==s and gseqc!=gseqc2:
        print("missense mutation")
        print()
    elif s<len(gseqc2):
        print("nonsense mutation")
        print("There are stop codons in bases [",f,"] to [",s-1,"]")
        print()
    gseqc=list(gseqc)
    gseqc2=gseqc.copy()
print("***************************************************************")

1#inversion
invmblist=[]
invmb=""
for i in range(len(inv)):
    if inv[i]!="/":
        invmb+=inv[i]
    if inv[i]=="/" or i==(len(inv))-1:
        invmb1=""
        invmblist1=[]
        for j in range(len(invmb)):
            if invmb[j]!=",":
                invmb1+=invmb[j]
            if (invmb[j]==",") or (j==(len(invmb))-1):
                invmblist1.append(int(invmb1))
                invmb1=""
        invmblist.append(invmblist1)
        invmb=""
#print("invmblist",invmblist)
for i in range(len(invmblist)):
    rndinv=random.sample(range(len(gseqc)-invmblist[i][1]),invmblist[i][0])
    for j in rndinv:
        print(j,"to",j+invmblist[i][1]-1,"in native gene is:",gseqc[j:j+invmblist[i][1]])
        invrev=list(reversed(gseqc[j:j+invmblist[i][1]]))
        for v in range(j,j+invmblist[i][1]):
            gseqc[v]=invrev[0]
            invrev.pop(0)
        print(j,"to",j+invmblist[i][1]-1,"inversion",gseqc[j:j+invmblist[i][1]])
        fr=0
        gseqc=''.join(gseqc)
        for f in range(0,len(gseqc)-2,3):
            if gseqc[f:f+3]!=gseqc2[f:f+3]:
                fr+=1
            if gseqc[f:f+3] in stop:
                s=f+3
                break
        if gseqc2==gseqc:
            print("silent mutation")
            print()
        elif len(gseqc2)==s and gseqc!=gseqc2:
            print("missense mutation")
            print()
        elif s<len(gseqc2):
            print("nonsense mutation")
            print("There are stop codons in bases [",f,"] to [",s-1,"]")
            print()
        gseqc=list(gseqc)
        gseqc2=gseqc.copy()
print("***************************************************************")

#insertion
insblist=[]
insb=""
for i in range(len(ins)):
    if ins[i]!="/":
        insb+=ins[i]
    if ins[i]=="/" or i==(len(ins))-1:
        insb1=""
        insblist1=[]
        for j in range(len(insb)):
            if insb[j]!=",":
                insb1+=insb[j]
            if (insb[j]==",") or (j==(len(insb))-1):
                insblist1.append(int(insb1))
                insb1=""
        insblist.append(insblist1)
        insb=""
#print("insblist",insblist)
ncodec=ncode.copy()
inscode=[]
for i in range(len(insblist)):
    rndins=random.sample(range(len(gseqc)),insblist[i][0])
    for j in rndins:
        for s in range(insblist[i][1]):
            rndins1= random.sample(ncodec,1)
            inscode.append(rndins1)
        print(j,"to",j+insblist[i][1]-1,"insertion",inscode)
        for ri in range(j,j+insblist[i][1]):
            inse="".join(inscode[0])
            gseqc.insert(ri,inse)
            inscode.pop(0)
        fr=0
        gseqc=''.join(gseqc)
        for f in range(0,len(gseqc2)-2,3):
            if gseqc[f:f+3]==gseqc2[f:f+3]:
                fr+=1
            if gseqc[f:f+3] in stop:
                s=f+3
                break
        if gseqc2==gseqc:
            print("silent mutation")
            print()
        elif len(gseqc2)==s and gseqc!=gseqc2 and fr==1:
            print("missense mutation")
            print()
        elif s<len(gseqc2) and fr<=(insblist[i][1]//3)+1:
            print("nonsense mutation")
            print("There are stop codons in bases [",f,"] to [",s-1,"]")
            print()
        elif fr>1:
            print("frameshift mutation")
            print("There are stop codons in bases [",f,"] to [",s-1,"]")
            print()
        gseqc=list(gseqc)
        gseqc2=gseqc.copy()

print("***************************************************************")

#deletion
deleblist=[]
deleb=""
for i in range(len(dele)):
    if dele[i]!="/":
        deleb+=dele[i]
    if dele[i]=="/" or i==(len(dele))-1:
        deleb1=""
        deleblist1=[]
        for j in range(len(deleb)):
            if deleb[j]!=",":
                deleb1+=deleb[j]
            if (deleb[j]==",") or (j==(len(deleb))-1):
                deleblist1.append(int(deleb1))
                deleb1=""
        deleblist.append(deleblist1)
        deleb=""
#print("deleblist",deleblist)
for i in range(len(deleblist)):
    rnddele=random.sample(range(len(gseqc)-deleblist[i][1]),deleblist[i][0])
    for j in rnddele:
        delecode=[]
        if j+deleblist[i][1]-1>len(gseqc):
            de=(j+deleblist[i][1]-1)-(len(gseqc))
        else:
            de=deleblist[i][1]
        for d in range(de):
            delecode.append(gseqc.pop(j))
        print(j,"to",j+deleblist[i][1]-1,"deletion",delecode)
        fr=0
        gseqc=''.join(gseqc)
        for f in range(0,len(gseqc)-2,3):
            if gseqc[f:f+3]==gseqc2[f:f+3]:
                fr+=1
            if gseqc[f:f+3] in stop:
                s=f+3
                break
        if gseqc2==gseqc:
            print("silent mutation")
            print()
        elif len(gseqc2)==s and gseqc!=gseqc2 and fr==1:
            print("missense mutation")
            print()
        elif s<len(gseqc2) and fr<=2:
            print("nonsense mutation")
            print("There are stop codons in bases [",f,"] to [",s-1,"]")
            print()
        elif fr>1:
            print("frameshift mutation")
            print("There are stop codons in bases [",f,"] to [",s-1,"]")
            print()
        gseqc=list(gseqc)
        gseqc2=gseqc.copy()
        
print("***************************************************************")
#translocation
if trl=="y":
    rndg=random.randint(1,len(gseqc))
    rndb=random.randint(1,len(bcl2_h))
    gseqc="".join(gseqc)
    trlcode=gseqc[:rndg+1]+ bcl2_h[rndb:]
    print("translocation:",0,"to", rndg,"of gseq and",rndb,"to end of bcl2_human")
    fr=0
    gseqc=''.join(gseqc)
    lenmin=min(len(gseqc),len(gseqc2))
    for f in range(0,lenmin-2,3):
        if trlcode[f:f+3]==gseqc2[f:f+3]:
            fr+=1
        if trlcode[f:f+3] in stop:
                s=f+3
                break
    if gseqc2==trlcode:
        print("silent mutation")
    elif len(gseqc2)==s and trlcode!=gseqc2 and fr==1:
        print("missense mutation")
    elif s<len(gseqc2) and fr<=1:
        print("nonsense mutation")
        print("There are stop codons in bases [",f,"] to [",s-1,"]")
        print()
    elif fr>1:
        print("frameshift mutation")
        print("There are stop codons in bases [",f,"] to [",s-1,"]")
        print()
gseqc=gseqc[:s+1]

protein=""
for i in range(3,len(gseqc),3):
    t=gseqc[i-3:i]
    for j in gcode:
        if t in j:
            protein+=aacode[gcode.index(j)]
            break
    if t in stop:
        break
print(protein)
