def parseLine(line):
    if line[0]=="#":
        return None
    sl=line.split()
    if sl[2][:8]!="svim.DUP":
        return None
    NewAlt="<DUP>"
    info=sl[7]
    isl=sl[7].split(";")
    first=True
    newinfo=""
    for s in isl:
        if not first:
            newinfo+=";"
        first=False
        ss=s.split("=")
        if ss[0]=="SVTYPE":
            newinfo+="SVTYPE=DUP"
        else:
            newinfo+=s
    newline=sl[0]
    for i in range(1,len(sl)):
        if (i==4):
            newline+="\t"+NewAlt
        elif (i==7):
            newline+="\t"+newinfo
        else:
            newline+="\t"+sl[i]
    return newline

import sys
for line in sys.stdin:
    newline=parseLine(line)
    if newline==None:
        print(line,end="")
    else:
        print(newline)
    