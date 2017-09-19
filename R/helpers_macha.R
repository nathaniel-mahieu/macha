macha.check_validity = function(macha) {
  
  names = names(macha)
  
  if ("s" %in% names) {
    cat("macha$s:")
    cat("\n Sequential scans -", all(diff(macha$s$s) == 1))
    cat("\n No NA values -", !any(is.na(macha$s)))
  } else {
    cat("macha$s not found")
  }
  cat("\n\n")
  
  if ("k_r" %in% names) {
    cat("macha$k_r:")
    cat("\n One scan per ROI -", macha$k_r[macha$k,,on="k",nomatch=0][,.(dups = sum(duplicated(s)) > 0),by="r"]$dups %>% { !any(.) })
    cat("\n No NA values -", !any(is.na(macha$k_r)))
  } else {
    cat("macha$k_r not found")
  }
  cat("\n\n")
  
  if ("r_mchan" %in% names) {
    cat("macha$r_mchan:")
    cat("\n One scan per mchan -", macha$r_mchan[macha$k_r[macha$k,,on="k",nomatch=0],,on="r",nomatch=0][,.(dups = sum(duplicated(s)) > 0),by="mchan"]$dups %>% { !any(.) })
    cat("\n No NA values -", !any(is.na(macha$r_mchan)))
  } else {
    cat("macha$r_mchan not found")
  }
  cat("\n\n")
  
  if ("roi_cache" %in% names) {
    cat("macha$roi_cache:")
    cat("\n One scan per roi -", macha$roi_cache[,.(dups = sum(duplicated(s)) > 0),by="r"]$dups %>% { !any(.) })
    cat("\n No missing scans -", macha$roi_cache[,.(diffs = diff(s) <= 1),by="r"]$diffs %>% all)
    #cat("\n No NA values -", !any(is.na(macha$roi_cache)), sep="")
  } else {
    cat("macha$roi_cache not found")
  }
  cat("\n\n")
  
}
