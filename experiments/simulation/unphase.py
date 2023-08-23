import sys

F=open(sys.argv[1],"r")

for line in F:
    sl=line.strip().split("\t")
    if len(sl)>1:
        if (sl[6]=="1|1"):
            sl[6]="1/1"
        elif (sl[6]=="0|0"):
            sl[6]="0/0"
        else:
            sl[6]="0/1"
        newline=sl[0]
        for s in sl[1:]:
            newline+="\t"+s
        print(newline)
    else:
        print(line,end="")

F.close()