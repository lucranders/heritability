library(dplyr)
library(ggplot2)
library(rpart)
library(rpart.plot)

refPath_ = yaml::read_yaml('conf/local/paths.yml')
refPathResults_ = refPath_$pathSaveFilesSimple
results_ = read.table(refPathResults_ , sep = '|' , header = TRUE)
results_ %>% head(3)

selectedCols = list(columns_ = results_  %>% colnames) %>% as.data.frame %>% 
                filter(!grepl('_final|size|randEffects',columns_))

vals_ = unique(results_[selectedCols$columns_[1]])
for ( i in 2:nrow(selectedCols)){
    aux_ = unique(results_[selectedCols$columns_[i]])
    vals_ = merge(vals_, aux_ , by = NULL)
}

vals_ = as.data.frame(vals_)
vals_ %>% head(3)

threshold_ = .97
infos_ = merge(vals_ , results_ , all.x = T , all.y = F)
infos_ = infos_ %>% 
    mutate(h2 = Genes_final/(Genes_final + Residuals_final) ) %>% 
    mutate(flgConverge_ = ifelse(is.na(h2) | h2 > threshold_ , 0 , 1 ))

infosSelectedCols = list(columns_ = infos_  %>% colnames) %>% as.data.frame %>% 
                filter(!grepl('_final|size|randEffects|method|formulaF|h2|VIF',columns_))
infosSelect_ = infos_ %>% select(infosSelectedCols$columns_)
for (col_ in infosSelectedCols$columns_){
    if(col_ == 'outliers'){
        infosSelect_[is.na(infosSelect_[,col_]),col_] = 1    
    }
    infosSelect_[,col_] = as.factor(infosSelect_[,col_])
}
model <- glm(flgConverge_ ~.,family=binomial(link='logit'),data=infosSelect_)



tree <- rpart(flgConverge_ ~., data = infosSelect_)
jpeg(file="images/tree.jpeg" , quality = 100 , pointsize = 20 , width = 1080 , height = 1080)
rpart.plot(tree)
dev.off()

infosSelect_ = infosSelect_ %>% mutate(flgPr_ = predict(tree,infosSelect_)[,2])
infosSelect_ %>% head(3)

final_ = infosSelect_ %>% 
select(!c(Genes_init , Residuals_init , flgConverge_)) %>%
filter(flgPr_ > .9)  %>% 
select(!flgPr_) %>% distinct()
