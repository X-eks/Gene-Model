import pymysql
import train_model




if __name__ == '__main__':

    #数据库配置文件
    db_args = {'host': "47.107.250.206", 'user': "genemodel", 'password': "123456",\
         'database':"genemodel",  'charset' : 'utf8'}

    ## 循环运行文件
    for i in range(2,21):
        user_id = str(i+1000)
        uni_dir='D:/2_Professional Courses/132_基因项目/GWAS数据/新增_20基因数据/'
        user_file_dir=uni_dir+str(i)+'.vcf'
        print("目前正在运行文件：  "+ str(i)+'.vcf')
        train_model.train_model(user_id, user_file_dir, **db_args)

    ## 单独运行文件
    # user_id = input('请输入用户名： \n')
    # user_file_dir=input('请输入用户文件路径： \n')

    ## 单独运行01号文件
    # user_id = '1001'
    # user_file_dir= 'D:/2_Professional Courses/132_基因项目/GWAS数据/新增_20基因数据/1.vcf'
    # train_model.train_model(user_id, user_file_dir, **db_args)




# 1001
#model.SC_model('D:/2_Professional Courses/132_基因项目/GWAS数据/新增_20基因数据/1.vcf')