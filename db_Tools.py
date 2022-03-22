#获取所有可检测的疾病
import pymysql
import pandas as pd




def get_all_Disease(**db_args):
    db = pymysql.connect(**db_args)
    cursor = db.cursor()
    try:   
        sql = "SELECT * FROM disease_name"
        cursor.execute(sql)  #返回的是查询到的个数，是个int类型的数字
        disease_Results = cursor.fetchall()
        num_disease_in_GWAS = len(disease_Results)
        db.close()
    except Exception as e:
        db.close()
        print(e)
        print("Error: unable to fecth data")
    
    return_value = {'disease_Results': disease_Results, 'num_disease_in_GWAS': num_disease_in_GWAS}
    return return_value


def find_risk_SNPs(disease_EN, **db_args):
    #1. 数据库中查询所需的SNP位点, 按risk_Allele进行升序排列   
    try:
        db = pymysql.connect(**db_args)
        cursor = db.cursor()
        sql = "SELECT risk_Allele, OddsRatio, MAF  FROM SNP WHERE disease_EN = '" + disease_EN.replace("'", "\\'") + "' order by risk_Allele"
        print(sql)
        cursor.execute(sql)  #返回的是查询到的个数，是个int类型的数字      
        snp_Results = cursor.fetchall()
        db.close()
    except Exception as e:
        print(e)
        db.close()
    num_SNP_in_GWAS = len(snp_Results)
    print('数据库中包含SNP位点的个数： ', num_SNP_in_GWAS)
    print(snp_Results)
    return_value = {'snp_Results': snp_Results, 'num_SNP_in_GWAS': num_SNP_in_GWAS}
    return return_value


def insert_Record(disease_EN, user_ID, db_args, **kw):
    try:
        db = pymysql.connect(**db_args)
        cursor = db.cursor()
        sql = "INSERT INTO scoreRecord(disease_EN, user_ID, OR_score, SC_score, EV_score) VALUES (%s, %s, %s, %s, %s)" 
        cursor.execute(sql, (disease_EN, user_ID, kw['OR_score'], kw['SC_score'], kw['EV_score']))
        db.commit()
        db.close()
    except Exception as e:
        print(e)
        db.rollback()
        db.close()
