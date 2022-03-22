# 建立权重字典：Dict_OR & Dict_SC （建立字典Dict存储每个SNP各风险等位基因的信息）
import pandas as pd
from math import log



def get_Weight_Dict(snp_Results_tuple, matched_rs):

    df_SNP = pd.DataFrame(list(snp_Results_tuple), columns=['Risk_allele','OR', 'MAF'])#将snp_Results转换成df
    # 1. Dict_OR
    Dict_OR = dict() 
    #风险基因&OR值存入字典Dict(不包含用户不提供的SNP位点信息)
    #Dict_OR={'rs9592461': [['A', [7.062, 1.025]]], ['G', [1.021451]]], 'rs645592': [['T', [5.571]]]}
    for index, row in df_SNP.iterrows():
        rs = str(row['Risk_allele']).split('-')[0]
        AGCT = str(row['Risk_allele']).split('-')[1][-1]
        if rs in matched_rs:
            if rs in Dict_OR : #rs在字典中存在
                if AGCT not in [i[0] for i in Dict_OR[rs]]:  #但AGCT尚未存在时
                    Dict_OR[rs].append([AGCT,[float(row['OR'])]])
                else:
                    for i in Dict_OR[rs]:
                        if i[0]==AGCT:
                            i[1].append(float(row['OR']))
            elif not(row['OR'] is None):
                Dict_OR[rs]=[[AGCT,[float(row['OR'])]]]

    #根据字典Dict计算各SNP等位风险基因对应的风险权重
    for key, value in Dict_OR.items():
        #value,eg.:[['A', [7.062, 1.025]]], ['G', [1.021451]]]
        n = len(value)
        if n > 0:
            for x in value:
                #x,eg.:['A', [7.062, 1.025]]
                AGCT = x[0]
                OR_list = x[1]
                OR_num=len(x[1])
                OR=sum(OR_list)/OR_num
                OR_log=log(OR)
                x[1]=OR_log
                #x,eg.:['A', 0.3]
            Dict_OR[key]=value
    print('Dict_OR： ', Dict_OR)  #{'rs1656369': [['A', 1.76541],['G', 1.76541549]], 'rs10811965': [['T', 1.7279318]]}


    # 2.Dict_SC
    Dict_SC = dict() 
    for index, row in df_SNP.iterrows():
        rs = str(row['Risk_allele']).split('-')[0]
        AGCT = str(row['Risk_allele']).split('-')[1][-1]
        if rs in matched_rs :
            if rs in Dict_SC and AGCT not in Dict_SC[rs]:
                Dict_SC[rs].append(AGCT)
            else:
                Dict_SC[rs]=[AGCT] 
    print('Dict_SC： ', Dict_SC) #{'rs1234:['A','G'], 'rs1235':['C']}

    #3. 返回字典
    return_value = {'Dict_OR': Dict_OR, 'Dict_SC': Dict_SC}
    return return_value