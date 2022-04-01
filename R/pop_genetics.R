calculate_Fst = function(snpdata, groups=NA, from="Country"){
    if(is.na(groups)){
        stop("Please provide a vector of groups to be compared")
    }
    if(any(!(groups %in% unique(snpdata$meta[[from]])))){
        stop("not all specified groups belong to the ", from, "column of the metadata table")
    }
    for(i in 1:(length(groups)-1)){
        idx1 = which(snpdata$meta[[from]]==groups[i])
        idx1 = paste(idx1, collapse = ",")
        for(j in (i+1):length(groups)){
            idx2 = which(snpdata$meta[[from]]==groups[j])
            idx2 = paste(idx2, collapse = ",")
            system(sprintf("wcFst --target %s --background %s --file %s --deltaaf 0.1 --type PL", idx1, id2, snpdata$vcf))
        }
    }
}

calculate_iR = function(snpdata){

}

calculate_LD = function(snpdata){

}

get_ibs = function(snpdata){

}

get_ibd = function(snpdata){

}
