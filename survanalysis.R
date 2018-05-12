library(survival)

paths = c("C:/Users/zchen/Desktop/research/cancer/CESC_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/LGG_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/READ_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/PAAD_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/LIHC_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/OV_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/KIRC_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/GBM_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/SARC_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/KIRP_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/DLBC_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/PRAD_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/TGCT_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/THCA_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/KICH_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/CHOL_survdata/",
           "C:/Users/zchen/Desktop/research/cancer/UVM_survdata/")


for (path in paths){
setwd(path)
print(path)
pairs = list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)
x = vector("numeric",length = length(pairs))
i = 1
for( pair in pairs){
    tmp = read.table(pair,quote = "\"",comment.char = "")
    diff = survdiff(Surv(V2,V3)~V4,data = tmp)
    x[i] = pchisq(diff$chisq,length(diff$n)-1,lower.tail = FALSE)
    if(x[i]<0.05){
      print(pair)
    }
    i = i+1
  }
}
