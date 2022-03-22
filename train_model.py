import SimpleCount
import OddsRatio
import EV_GRS
import pandas as pd
import numpy as np
import vcf2excel
import datetime
import db_Tools
import get_SNPs_Matched
import get_Weight_Dict





def train_model(user_id, user_file_dir, **db_args):   #输入用户名, 被测试样本文件的路径
    #0. 创建输出表格
    Report_data=pd.DataFrame(columns=['Chinese', 'English', 'OR_score','SC_score','EV_score', 'SNP_info',\
        'num_SNP_in_GWAS', 'num_SNP_matched'])
    

    #1. 导入用户数据
    start=datetime.datetime.now()
    df_user=vcf2excel.vcf_process(user_file_dir)
    user_rs=df_user['rs_ID']#按照rs_ID排序(方便后续匹配查找)
    user_rs_num = len(user_rs)
    print('用户表里的信息总数：', user_rs_num)
    end=datetime.datetime.now()
    print('导入用户数据所用时间: %s Seconds'%(end-start))

    #获取所有可检测的疾病
    all_Disease = db_Tools.get_all_Disease(**db_args)
    disease_Results = all_Disease['disease_Results']
    num_disease_in_GWAS = all_Disease['num_disease_in_GWAS']

    #2. 循环计算每种疾病的患病风险
    for disease_index, disease in enumerate(disease_Results):
        print('-'*30 ,disease_index+1, ' / ', num_disease_in_GWAS,'正在运算疾病： ', disease[0],'-'*30)
        # 2.1 数据库中查询所需的SNP位点
        all_SNPs = db_Tools.find_risk_SNPs(disease[0], **db_args)
        snp_Results_tuple = all_SNPs['snp_Results']
        num_SNP_in_GWAS = all_SNPs['num_SNP_in_GWAS']
        
        if num_SNP_in_GWAS == 0:    #将数据库中包含SNP位点的个数为零的疾病直接写入
            Report_data=Report_data.append({'Chinese': disease[1], 'English': disease[0],
            'OR_score': None, 'SC_score': None, 'EV_score': None, 'SNP_info': 'No risk SNPs in Database',\
                'num_SNP_in_GWAS': 0, 'num_SNP_matched': 0}, ignore_index=True)
            continue
        
        # 2.2 将用户数据与数据库中查找出的风险基因位点做匹配
        SNPs_matched = get_SNPs_Matched.get_SNPs_Matched(snp_Results_tuple, df_user)
        matched_rs = SNPs_matched['matched_rs']
        num_SNP_matched = SNPs_matched['num_SNP_matched']
        user_rs_index = SNPs_matched['user_rs_index']
        
        if num_SNP_matched==0:  #若匹配上的SNP位点为0，直接写入
            Report_data=Report_data.append({'Chinese': disease[1], 'English': disease[0],
                'OR_score': None, 'SC_score': None,'EV_score': None, 'SNP_info': 'No SNPs matched',\
                    'num_SNP_in_GWAS': num_SNP_in_GWAS, 'num_SNP_matched': 0}, ignore_index=True)
            continue
        print(matched_rs)  # 能配对上的SNP

        # 2.3 建立权重字典：Dict_OR & Dict_SC （建立字典Dict存储每个SNP各风险等位基因的信息）
        all_weighted_Dict = get_Weight_Dict.get_Weight_Dict(snp_Results_tuple, matched_rs)
        Dict_OR = all_weighted_Dict['Dict_OR']
        Dict_SC = all_weighted_Dict['Dict_SC']

        # 2.4 运行模型
        if len(Dict_OR)>0:
            OR_score = OddsRatio.Method1_OR(Dict_OR, df_user, user_rs_index)
        else:
            OR_score = None
        if len(Dict_SC)>0:
            SC_score= SimpleCount.Method2_SC(Dict_SC, df_user, user_rs_index)
        else:
            SC_score = None
        EV_score= None

        # 2.5 将运行结果保存至数据库中
        kw = {'OR_score': OR_score, 'SC_score': SC_score, 'EV_score': EV_score}
        db_Tools.insert_Record(disease[0], user_id, db_args,  **kw)

        # 2.6 将成功运行的疾病信息输出至报告中
        Report_data=Report_data.append({'Chinese': disease[1], 'English': disease[0],
            'OR_score': OR_score, 'SC_score': SC_score, 'EV_score': EV_score, \
                'SNP_info': Dict_OR, 'num_SNP_in_GWAS': num_SNP_in_GWAS, 'num_SNP_matched': num_SNP_matched}, ignore_index=True)

    # 3. 输出生成报告所需文件
    #Report_data.to_excel('output_data/report_data/'+ 'report_data' + user_id + '.xlsx',index=False)



























