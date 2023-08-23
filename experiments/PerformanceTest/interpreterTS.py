Lines=[]
import sys
with open(sys.argv[1]) as f:
    for line in f:
        Lines.append(line)

Results={}
Programs=["cuteSV","sniffles", "svim", "kled"]
DataSets=["CLR","ONT","CCS"]
def getDataSet(DataSets,Line):
    for d in DataSets:
        if d in Line:
            return d
    return None
def getProgram(Programs, Line):
    for p in Programs:
        if p in Line:
            return p
    return None
def getT(Program, Line):
    if Program=="sniffles":
        sl=Line.split()
        for i in range(len(sl)):
            if sl[i]=="--threads":
                return int(sl[i+1].strip('"'))
    elif Program=="cuteSV":
        sl=Line.split()
        for i in range(len(sl)):
            if sl[i]=="-t":
                return int(sl[i+1].strip('"'))
    elif Program=="kled":
        sl=Line.split()
        for i in range(len(sl)):
            if sl[i]=="-t":
                return int(sl[i+1].strip('"'))

def searchTmemc(Lines, Index, N=10):
    for i in range(max(0,Index-N), min(len(Lines),Index+22+N)):
        line=Lines[i]
        if line[:len("Run every")]=="Run every":
            return line
    return None

def getTimeInSeconds(Time):
    sl=Time.split(":")
    if len(sl)==1:
        return int(float(sl[0]))
    elif len(sl)==2:
        return int(sl[0])*60+int(float(sl[1]))
    elif len(sl)==3:
        return int(sl[0])*3600+int(sl[1])*60+int(float(sl[2]))
    return None

for i in range(len(Lines)):
    line=Lines[i]
    if "Command being timed:" in line:
        Program=getProgram(Programs,line)
        DataSet=getDataSet(DataSets,line)
        if Program==None or DataSet==None:
            continue
        T=getT(Program,line)
        TmemcLine=searchTmemc(Lines,i)
        if TmemcLine==None:
            continue
        TimeLine=Lines[i+4]
        Time=TimeLine.split()[-1]
        Time=getTimeInSeconds(Time)
        for item in TmemcLine.split():
            item=item.strip(".")
            if item[-2:]=="kb":
                Space="%.2f"%(float(item[:-2])/1024)
        if DataSet not in Results.keys():
            Results[DataSet]={}
        if Program not in Results[DataSet]:
            Results[DataSet][Program]={}
        if T not in Results[DataSet][Program].keys():
            Results[DataSet][Program][T]={}
        Results[DataSet][Program][T]["Time"]=Time
        Results[DataSet][Program][T]["Space"]=Space

Ts=[]
Ps=[]
Ds=[]
for DataSet in DataSets:
    for Program in Programs:
        for T in [1,2,4,8,16]:
            if DataSet in Results.keys() and Program in Results[DataSet].keys() and T in Results[DataSet][Program].keys():
                if T not in Ts:
                    Ts.append(T)
                if Program not in Ps:
                    Ps.append(Program)
                if DataSet not in Ds:
                    Ds.append(DataSet)
print(Ts)
for DataSet in Ds:
    print(DataSet+":")
    print("Time(s):")
    for Program in Ps:
        print(Program,end="\t")
        for T in Ts:
            if DataSet in Results.keys() and Program in Results[DataSet].keys() and T in Results[DataSet][Program].keys():
                print(Results[DataSet][Program][T]["Time"],end="\t")
        print()
    print("Space(MB):")
    for Program in Ps:
        print(Program,end="\t")
        for T in Ts:
            if DataSet in Results.keys() and Program in Results[DataSet].keys() and T in Results[DataSet][Program].keys():
                print(Results[DataSet][Program][T]["Space"],end="\t")
        print()
