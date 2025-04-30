max_row_norm <- function(networks){
  T = length(networks)
  num_countries = nrow(networks[[1]])
  rowsums = rep(0,num_countries)
  for(t in 1:T){ # FOR EVERY NETWORK
    for(i in 1:num_countries){ # FOR EVERY 
      #print(unname(colSums(networks[[t]])[i]))
      if(abs(unname(rowSums(networks[[t]])[i])) > abs(rowsums[i])){
        rowsums[i] = rowSums(networks[[t]])[i]
      }
    }
  }
  for(t in 1:T){ # FOR EVERY NETWORK
    for(i in 1:num_countries){ # FOR EVERY COLUMN
      networks[[t]][i,] = networks[[t]][i,]/rowsums[i] #DIV ROW BY 'MAX COLSUM' OVER T
    }
  }
  print(rowsums)
  # colsums = rep(0,num_countries)
  # for(t in 1:T){ # FOR EVERY NETWORK
  #   for(i in 1:num_countries){ # FOR EVERY 
  #     #print(unname(colSums(networks[[t]])[i]))
  #     if(abs(unname(colSums(networks[[t]])[i])) > abs(colsums[i])){
  #       colsums[i] = colSums(networks[[t]])[i]
  #     }
  #   }
  # }
  # for(t in 1:T){ # FOR EVERY NETWORK
  #   for(i in 1:num_countries){ # FOR EVERY COLUMN
  #     networks[[t]][,i] = networks[[t]][,i]/colsums[i] #DIV ROW BY 'MAX COLSUM' OVER T
  #   }
  # }
  # print(colsums)
  return(networks)
}
