def parseLine(line):
    if line[0]=="#":
        return None
    sl=line.split()
    if sl[2][:8]!="svim.INV":
        return None
    pos=int(sl[1])
    info=sl[7]
    isl=sl[7].split(";")
    for s in isl:
        ss=s.split("=")
        if ss[0]=="END":
            end=int(ss[1])
            break
    newinfo=info+";SVLEN="+str(end-pos)
    newline=sl[0]
    for i in range(1,len(sl)):
        if (i==7):
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
    