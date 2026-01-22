# this function creates a simple progress bar
create_progress_bar <- function(pct, width = 44) {
	filled <- round(width * pct / 100)
	bar <- paste0(rep("#", filled), collapse = "")
	empty <- paste0(rep(" ", width - filled), collapse = "")
	sprintf("[%s%s] %5.1f%%", bar, empty, pct)
}
