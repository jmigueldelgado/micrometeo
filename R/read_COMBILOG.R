#' read logger file from theodor friedrichs
#' @import dplyr
#' @import readr
#' @import lubridate
#' @param path is the absolute path to the logger file
#' @export
read_COMBILOG <- function(path)
{
    line <- readLines(path,n=6)[6]
    header <- unlist(strsplit(line,split="\t"))    
    tbl <- read_tsv(path,skip=8,col_names=header) %>%
        mutate(datetime=mdy_hms(TimeDate)) %>%
        select(-Name,-TimeDate)
    return(tbl)
}
