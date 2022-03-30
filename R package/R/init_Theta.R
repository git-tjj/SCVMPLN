# init_Theta
init_Theta<-function (VMPLN_list,lambda_use,penalize_diagonal,Theta_threshold,core_num) 
{
  ######
  celltypes_num <- ncol(VMPLN_list$U_mat)
  ###
  ##
  VMPLN_list$Theta_mat_list<-Theta_fun(U_mat = VMPLN_list$U_mat,
                                     Theta_mat_list = lapply(X = 1:celltypes_num,FUN = function(g){return(as(diag(ncol(VMPLN_list$mu_mat)),"sparseMatrix"))}),
                                     m_mat_list = VMPLN_list$m_mat_list,
                                     s_sq_mat_list = VMPLN_list$s_sq_mat_list,
                                     p_GRN = length(VMPLN_list$gene_GRN_index_use),
                                     lambda_use = lambda_use,
                                     zero_GRN = VMPLN_list$zero_GRN,
                                     stop_threshold = Theta_threshold,
                                     var_index_Theta = rep(1,celltypes_num),
                                     penalize_diagonal = penalize_diagonal,
                                     zero_GRN_use = VMPLN_list$zero_GRN_use,
                                     core_num = 1)
  ######
  ##
  log_det_vec<-NULL
  if((ncol(VMPLN_list$obs_mat) - length(VMPLN_list$gene_GRN_index_use)) == 0){
    log_det_vec <- log_det(matrix_list = VMPLN_list$Theta_mat_list,celltypes_num = celltypes_num,core_num = min(core_num,celltypes_num))
  }else{
    log_det_vec <- log_det_block(matrix_list = VMPLN_list$Theta_mat_list,celltypes_num = celltypes_num,p_GRN = length(VMPLN_list$gene_GRN_index_use),core_num = min(core_num,celltypes_num))
  }
  ######################################################################################
  VMPLN_list[["log_det_vec"]]<-log_det_vec
  
  ######################################################################################
  
  return(VMPLN_list)
}

