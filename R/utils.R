

# Function to check the data
#
verifData <- function(u1, u2)
{
  if (length(u1) != length(u2))
  {
    stop(paste0("u1 and u2 have different lengths: ",
                length(u1),
                " and ", length(u2), ".") )
  }

  which_not01 = which(u1 < 0 | u1 > 1 | u2 < 0 | u2 > 1)
  if (length(which_not01) > 0)
  {
    stop(paste0("Data in u1 and u2 must belong to [0,1]."))
  }

}


