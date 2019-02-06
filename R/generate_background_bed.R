
######################## INPUT FILES AND PARAMETERS  #######################

Args <- commandArgs()

if( length(Args) < 8 ){
    print("Usage: Rscript generate_background_bed.R <chr length> <sampling number> <exclude list> <output>")
    quit()
}

file_chr         <- Args[6]
sampling_numbers <- as.numeric(Args[7])
output           <- Args[8]
file_exclude     <- Args[9]

print(file_exclude)

library(regioneR)
subtract_regions <- function(regiona,regionb){
  subtract <- subtractRegions(regiona,regionb)
  return(data.frame(chr=seqnames(subtract),start=start(subtract),end=end(subtract)))
}

writeTable <- function(df,output,col.names=F){
  write.table(df,output,quote=F,col.names = col.names,row.names = F,sep="\t")
}

chr_length  <- read.table(file_chr,header=F)

df_regions <- data.frame(chr=chr_length$V1,start=1,end=chr_length$V2)

if( ! is.na(file_exclude) ){
  exclude     <- read.table(file_exclude,header=F)
  df_regions <- subtract_regions(data.frame(chr=chr_length$V1,start=1,end=chr_length$V2),data.frame(chr=exclude$V1,start=exclude$V2,end=exclude$V3))
}

print(paste("Genome_size:",sum(chr_length$V2)))
region_size = sum(df_regions$end - df_regions$start + 1)
print(paste("Region_size:",region_size))

# simulate 2-DSB in fork
index_sampled_positions <- sample(c(1:region_size),size=sampling_numbers,replace=FALSE)

k = 0 
sampled_chr  <- c() 
sampled_pos  <- c() 
for( i in 1:nrow(chr_length) ){
  chr   = chr_length[i,1]
  left  = 1
  right = chr_length[i,2]
  
  k_start = k + 1
  k_end   = k + (right - left + 1)
  k = k_end
  
  for( s in index_sampled_positions){
    if( s >= k_start & s <= k_end){
      d = s - k_start
      p = left + d
      
      sampled_chr <- c(sampled_chr,as.character(chr))
      sampled_pos <- c(sampled_pos,p )
    }
  }
}


df_results <- data.frame(chr=sampled_chr,start=sampled_pos,end=sampled_pos)

writeTable(df_results,output)
