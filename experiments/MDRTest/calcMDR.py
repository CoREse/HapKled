import sys
#gt: record.samples[i]['GT'] (tuple)
def parse_gt(gt):
    if gt == None or gt[0] == None:
        return 0
    if gt[0] + gt[1] == 0:
        return 0
    if gt[0] + gt[1] == 1:
        return 1
    if gt[0] + gt[1] == 2:
        return 2

# False -> inconsistency
def check(son, fa, ma):
    if fa == 0 and ma == 0:
        if son == 1 or son == 2:
            return False
    if fa == 0 and ma == 1 and son == 2:
        return False
    if fa == 1 and ma == 0 and son == 2:
        return False
    if fa == 0 and ma == 2:
        if son == 0 or son == 2:
            return False
    if fa == 2 and ma == 0:
        if son == 0 or son == 2:
            return False
    if fa == 1 and ma == 2 and son == 0:
        return False
    if fa == 2 and ma == 1 and son == 0:
        return False
    if fa == 2 and ma == 2:
        if son == 0 or son == 1:
            return False
    return True

from pysam import VariantFile
def main(vcf_file, output_file):
    vcf_reader = VariantFile(vcf_file, 'r')
    bcf_out = VariantFile(output_file, 'w', header=vcf_reader.header)
    cnt = 0
    tot = 0
    InconsistentErrorCount=0
    ConsistentSVCount=0
    for record in vcf_reader.fetch():
        tot += 1
        son = parse_gt(record.samples[0]['GT'])
        father = parse_gt(record.samples[1]['GT'])
        mother = parse_gt(record.samples[2]['GT'])
        if son !=0 and father!=0 and mother!=0:
            nonzerotot+=1
        res = check(son, father, mother)
        if res == False:
            if son==0 or father==0 or mother==0:
                InconsistentErrorCount+=1
            cnt += 1
            bcf_out.write(record)
    print("InconsistentError:",InconsistentErrorCount)
    print("ConsistentSVDiscordance:",cnt-InconsistentErrorCount)
    print("ConsistentSVCount:",ConsistentSVCount)
    print("Non zero MDR:",(cnt-InconsistentErrorCount)/ConsistentSVCount)
    print("%s\t%s\t%s"%(ConsistentSVCount,cnt-InconsistentErrorCount,(cnt-InconsistentErrorCount)/ConsistentSVCount))
    return cnt, tot

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])