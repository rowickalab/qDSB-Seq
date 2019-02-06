#!/usr/bin/Rscript
# take cutting efficiency for each site, BLESS reads on each site, wig for BLESS reads
# then calculate BLESS reads for each DSB, and total studied DSBs
# Author: Yingjie Zhu yizhu@utmb.edu

calc_fcut <- function(data,summit){

    sum_uncut <- sum(data$Uncut_reads)
    sum_cut   <- sum(data$Left_reads) + sum(data$Right_reads)

    fcut <- round(sum_cut/(sum_cut+2*sum_uncut),6)
    if( summit == 1 )
        fcut <- round(sum_cut/(sum_cut+sum_uncut),6)

    fcut
}

calc_fcut_by_site <- function(data,summit){

    sum_uncut <- data$Uncut_reads
    sum_cut   <- data$Left_reads + data$Right_reads


    fcut <- round(sum_cut/(sum_cut+2*sum_uncut),6)

    if( summit == 1 )
        fcut <- round(sum_cut/(sum_cut+sum_uncut),6)

    fcut
}


calc_fbg <- function(data){

    fbg <- round(quantile(data$Cutting_efficiency,probs=0.5,type=1,na.rm=TRUE),6)
    print(fbg)

}

calc_fbg_sd <- function(data){

    fbg_sd <- round(sd(data$Cutting_efficiency),6)

}

calc_fcut_sd <- function(data_fcut){

    fcut_sd <- round(sd(data_fcut$Cutting_efficiency_rmbg),6)
}

calc_fcut_sd_by_fbg <- function(fcut,fbg,fbg_sd){

    fcut_min <- fcut - fbg - fbg_sd 
    fcut_max <- fcut - fbg + fbg_sd 

    fcut_sd <- round(sd(c(fcut_min,fcut_max)),6)
}

# filter out sites based on low gDNA reads
filter_fcut_low_coverage <- function(data_fcut,lowest){

    data_fcut <- data_fcut[data_fcut$Total_reads>=lowest,]

}

# filter out sites based on percentage and values
filter_fcut_outliers_q10_q90 <- function(data_fcut){

    q10 <- quantile(data_fcut$Cutting_efficiency_rmbg,probs=0.1,type=1,na.rm=TRUE)
    q90 <- quantile(data_fcut$Cutting_efficiency_rmbg,probs=0.9,type=1,na.rm=TRUE)

    data_fcut <- data_fcut[data_fcut$Cutting_efficiency_rmbg > q10, ]
    data_fcut <- data_fcut[data_fcut$Cutting_efficiency_rmbg < q90, ]

}

# filter out sites based on percentage of fcut
filter_fcut_outliers_N_20perc <- function(data_fcut){

    n=nrow(data_fcut)
    n10=round(n*0.1)
    data_fcut <- data_fcut[order(data_fcut$Cutting_efficiency_rmbg),] 
    data_fcut <- data_fcut[(n10+1):(n-n10),]

}

# filter out sites based on BLESS reads
# (BLESS_reads - median)/median
filter_fcut_outliers_BLESS_reads <- function(data_fcut){
    
    score <- (data_fcut$BLESS_reads - median(data_fcut$BLESS_reads)) / median(data_fcut$BLESS_reads)

    data_fcut <- data_fcut[which(score>=-0.5 & score<=0.5),]
}

# count reads in wig
count_wig <- function(data){
    sum(abs(data$V4))   
}

# count reads in sequence file
count_sequence <- function(data){
    nrow(data)   
}

# calculate DSBs
calc_DSBs <- function(sum_of_labeled_reads_studied,sum_of_labeled_reads_enz,fcut,nsites,index_12DSB,mixprop){

    DSBs <- (sum_of_labeled_reads_studied * fcut * nsites * index_12DSB * mixprop ) / sum_of_labeled_reads_enz

}

# calculate DSBs perl million cells
calc_DSBs_perM <- function(labeled_reads_studied,sum_of_labeled_reads_enz,fcut,nsites,index_12DSB,mixprop){

    DSBs_perM <- (labeled_reads_studied * fcut * nsites * index_12DSB * mixprop ) / sum_of_labeled_reads_enz * 1000000

}
