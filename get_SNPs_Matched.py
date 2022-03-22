#将用户数据与数据库中查找出的风险基因位点做匹配
import datetime
import pandas as pd



def get_SNPs_Matched(snp_Results, df_user):

    
    df_SNP = pd.DataFrame(list(snp_Results), columns=['Risk_allele','OR', 'MAF'])#将snp_Results转换成df
    unique_rs_in_DB = sorted(set([str(x).split('-')[0] for x in df_SNP['Risk_allele']]))#提取GWAS中的SNP（去重）
    user_rs_index = []  #记录匹配上SNP在df_user表中的序号,方便后续计算时查找
    matched_rs = [] #记录匹配上的SNP

    start=datetime.datetime.now()

    # for rs in unique_rs_in_DB:
    #     for index, row in df_user.iterrows():
    #         if rs==row['rs_ID']:
    #             matched_rs.append(rs)
    #             user_rs_index.append(index)

    list_user_rs = df_user['rs_ID'].tolist()
    for rs in unique_rs_in_DB:
        try:
            index_i = list_user_rs.index(rs)
            matched_rs.append(rs)
            user_rs_index.append(index_i)
        except Exception:
            continue


    end=datetime.datetime.now()
    print('匹配所花费时间: %s Seconds'%(end-start))
    num_SNP_matched = len(matched_rs)
    print("能匹配上的SNP个数：  ", num_SNP_matched)
    return_value = {'matched_rs': matched_rs, 'num_SNP_matched': num_SNP_matched, 'user_rs_index': user_rs_index}
    return return_value



