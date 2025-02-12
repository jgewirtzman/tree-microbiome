#clr transforms data
#run after Tree_microbiome_all.R

object <- trans_norm$new(dataset = amplicon_data)
object_clr <- object$norm(method = "clr")
object_lefse <- trans_diff$new(dataset = object_clr, method = "lefse", alpha = 0.01, lefse_subgroup = NULL, rm_un=TRUE)